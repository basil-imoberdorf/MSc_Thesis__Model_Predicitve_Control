# Standard libraries
import os
import sys
import glob
import time
import signal
from datetime import datetime

# Libraries that have to be installed
import do_mpc
from casadi import *
import megaind
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
# Note: this script must be executed within the virtual environment on the Raspberry Pi to import do-mpc and casadi

#######################################################################################################################
# The MPC is build with a watchdog FileHandler as parent class to dynamically read measurements and pass it to the MPC
class CustomMPC(FileSystemEventHandler):
    def __init__(self):
        self.last_modified_time = 0  # helper variable used to prevent calling on_modified() tice when new values arrive
        self.debounce_delay = 1  # Time in seconds to wait before accepting new changes, on_modified fires twice usually
        self.measurements = np.zeros(10)  # only initialization, these zeros will be overwritten
        self.k = 0  # variable to keep track of the current time stap, will increase for by 1 every iteration
        # Time step for the algorithm
        self.time_step = 2 / 60  # [h] every 2 minutes a measurement becomes available therefore time step = 2/60 h
        self.n_mpc_horizon = 20
        self.n_mhe_horizon = 10
        self.is_mpc_active = False
        self.current_measurements = 50*np.ones(13)  # initialize measurement vector with glucose value above 1 g/L
        # Set default parameters, these variables will be updated upon successful parameter estimation...
        self.mu_max_estimated = 0.69  # Maximal specific growth rate [h-1]
        self.Y_xs_estimated = 0.4     # Growth yield coefficient on glucose (during exp. feed) [g_cdw/g_glucose]
        self.K_s_estimated = 0.015    # Monod substrate saturation constant   [h-1]
        # Known model parameters. Feed glucose concentration was measured using HPLC-RID the day before start
        self.S_in = 690.81    # Glucose concentration in feed solution [g_glucose/L_feed]
        self.ms = 0.025       # Specific maintenance coefficient [g_glucose/(g_cdw*h)]
        # These variables will be updated before the mpc starts to subtract them from the current measurements
        self.added_vol_feed_pre_mpc = 0  # only initialization
        self.added_vol_base_pre_mpc = 0  # only initialization
        self.added_vol_acid_pre_mpc = 0  # only initialization
        
        # Parameters for the specific growth rate reference trajectory with shape of an exponential decay
        self.T = 16  # define the time between start and end value of the reference trajectory, starting when MPC starts
        self.mu_set_start = 0.4  # [h-1] initial mu_set
        self.mu_set_end = 0.1  # [h-1]   final mu_set
        self.mu_decay = 0.24  # shape with this parameter how sharp the decay is. Higher values == sharper decay
        self.mu_set_range = abs(self.mu_set_start - self.mu_set_end)  # [h-1]
        self.y_min = np.exp(-self.mu_decay * self.T)
        self.y_max = np.exp(-self.mu_decay * 0)  # is always 1 but hey it looks cooler

    def init_mpc_mhe(self, mu_max, Y_xs, K_s):
        # Initialize Model
        model = do_mpc.model.Model(model_type='continuous')  # model_type is either 'discrete' or 'continuous'
        # Introduce new variables to the model class. Define variable type, name and shape (optional).
        X = model.set_variable('_x', 'X')  # Cell Dry Weight [g/L]
        S = model.set_variable('_x', 'S')  # Glucose Concentration [g/L]
        V = model.set_variable('_x', 'V')  # Working Volume [g/L]
        F_feed = model.set_variable('_u', 'F_feed')    # Feed flow rate [g/L]
        F_base = model.set_variable('_tvp', 'F_base')  # Base flow rate [g/L]
        F_acid = model.set_variable('_tvp', 'F_acid')  # Acid flow rate [g/L]
        F_feed = model.set_expression('F_feed', F_feed)  # to store Feed flow rate
        F_base = model.set_expression('F_base', F_base)  # to store Base flow rate
        F_acid = model.set_expression('F_acid', F_acid)  # to store Acid flow rate

        # Variables for time-varying set-point tracking
        model.set_variable('_tvp', 'mu_set')

        # Define which variables are used by the MHE
        model.set_meas('X_meas', X, meas_noise=True)
        model.set_meas('S_meas', S, meas_noise=True)
        model.set_meas('V_meas', V, meas_noise=False)  # I assume the balances to be correct
        # Flow rate measurements, I assume no noise on flow rates (MHE has already enough burden)
        model.set_meas('F_f_meas', F_feed, meas_noise=False)
        model.set_meas('F_b_meas', F_base, meas_noise=False)
        model.set_meas('F_a_meas', F_acid, meas_noise=False)

        # Auxiliary terms
        mu = mu_max * S / (K_s + S)  # mu = mu_max*S/(K_s + S + (S**2/K_i))
        F_tot = F_feed + F_base + F_acid  # total flow rate into the bioreactor, sampling is not accounted for
        # Using set_expression makes it possible to store and plot auxiliary variables
        model.set_expression('F_tot', F_tot)
        model.set_expression('mu', mu)

        # Right-hand-side (rhs) of the ordinary differential equations
        model.set_rhs('X', mu*X - F_tot/V*X, process_noise=True)
        model.set_rhs('S', -(mu/Y_xs + self.ms)*X + F_feed/V*self.S_in - F_tot/V*S, process_noise=True)
        model.set_rhs('V', F_tot, process_noise=True)
        model.setup()  # Build the model
        # --------------------------------------------------------------------------------------------------------------
        # Initialize Moving Horizon State Estimator
        self.MHE = do_mpc.estimator.MHE(model)  # no model parameters are estimated, only states
        setup_mhe = {
            'n_horizon': self.n_mhe_horizon,  # Number of time steps for the estimation window
            't_step': self.time_step,  # Discrete sampling time in hours
            'state_discretization': 'collocation',
            'store_full_solution': True,
            'meas_from_data': True  # this means that the MHE takes the measurements from his own data storage variable
                                    # The MHE appends each incoming measurement to a "data" object. See documentation
        }
        self.MHE.set_param(**setup_mhe)

        # Define the diagonal weighting/penalty matrices
        P_v = 0.01 * np.diag(np.array([1, 1]))  # measurement noise penalty matrix
        P_w = 1 * np.diag(np.array([1, 1, 1]))  # process noise penalty matrix
        P_x = 1 * np.diag(np.array([1, 1, 1]))  # to penalize deviations from the states

        # Define objective function for MHE
        self.MHE.set_default_objective(P_x=P_x, P_v=P_v, P_w=P_w)  # see documentation for the details

        # State and parameter constraints/bounds:
        self.MHE.bounds['lower', '_x', 'X'] = 0.0
        self.MHE.bounds['lower', '_x', 'S'] = -0.01
        self.MHE.bounds['lower', '_x', 'V'] = 1.4

        self.MHE.bounds['upper', '_x', 'X'] = 120.0
        self.MHE.bounds['upper', '_x', 'S'] = 4.0
        self.MHE.bounds['upper', '_x', 'V'] = 3.0

        # Set up the functionality for time-varying parameters
        self.tvp_template_mhe = self.MHE.get_tvp_template()
        # Create dummy function simply returning the tvp_template structure
        self.tvp_fun_mhe = lambda t_now: self.tvp_template_mhe  # https://github.com/do-mpc/do-mpc/discussions/240
        self.MHE.set_tvp_fun(self.tvp_fun_mhe)
        self.MHE.settings.supress_ipopt_output()  # increase speed
        self.MHE.setup()
        # --------------------------------------------------------------------------------------------------------------
        # Build Model Predictive Controller
        self.MPC = do_mpc.controller.MPC(model)
        setup_mpc = {
            'n_horizon': self.n_mpc_horizon,
            'n_robust': 0,
            'open_loop': False,
            't_step': self.time_step,
            'state_discretization': 'collocation',
            'collocation_type': 'radau',
            'collocation_deg': 2,
            'collocation_ni': 2,
            'store_full_solution': True,
        }
        self.MPC.set_param(**setup_mpc)  # setting parameters can also be done like this: MPC.settings.n_horizon = 20

        # Define objective (cost) function
        mterm = (model.aux['mu'] - model.tvp['mu_set'])**2  # terminal cost, Meyer term
        lterm = (model.aux['mu'] - model.tvp['mu_set'])**2  # stage cost, Lagrange term
        self.MPC.set_objective(mterm=mterm, lterm=lterm)  # to omit a term use: mterm=DM(0) or lterm=DM(0)
        self.MPC.set_rterm(F_feed=50.0)  # set here the diagonal elements of R from rterm = delta_u'*R*delta_u

        # Define Constraints (lower/upper bounds) on control inputs and states
        self.MPC.bounds['lower', '_x', 'X'] = 0.0
        self.MPC.bounds['lower', '_x', 'S'] = -0.01
        self.MPC.bounds['lower', '_x', 'V'] = 1.4
        # upper bounds of the states
        self.MPC.bounds['upper', '_x', 'X'] = 120.0
        self.MPC.bounds['upper', '_x', 'S'] = 4.0
        self.MPC.bounds['upper', '_x', 'V'] = 3.0
        # upper and lower bounds of the control input
        self.MPC.bounds['lower', '_u', 'F_feed'] = 0.0
        self.MPC.bounds['upper', '_u', 'F_feed'] = 0.102  # L/h, pump capped at 16 RPM for this

        # Time-Varying-Parameters for MHE
        self.tvp_template_mpc = self.MPC.get_tvp_template()
        # Create dummy function simply returning the tvp_template structure
        tvp_fun_mpc = lambda t_now: self.tvp_template_mpc  # https://github.com/do-mpc/do-mpc/discussions/240
        self.MPC.set_tvp_fun(tvp_fun_mpc)
        self.MPC.settings.supress_ipopt_output()  # increase speed
        self.MPC.setup()
        # --------------------------------------------------------------------------------------------------------------
        # Initial Conditions
        X_0 = self.current_measurements[12]  # Biomass conc. via CO2 integral instead of via Base, because CDW via CO2 integral is considered to be more accurate at that time of the fermentation
        S_0 = self.current_measurements[2]   # Substrate conc. via PLS [g/L_culture]
        V_0 = 1.5   # Culture working volume inside the bioreactor [L]
        x0 = np.array([[X_0], [S_0], [V_0]])
        self.MPC.x0 = x0
        self.MHE.x0 = x0
        self.MPC.set_initial_guess()
        self.MHE.set_initial_guess()
        
        print('\n---------------------------------------')
        print(' MPC and MHE successfully initialized!')
        print(f"\t x_0 = {np.round(X_0, 2)} g/L")
        print(f"\t s_0 = {np.round(S_0, 2)} g/L")
        print(f"\t V_0 = {np.round(V_0, 2)} L")
        print('---------------------------------------\n')
        # --------------------------------------------------------------------------------------------------------------

    def read_measurements(self):
        print(f" {os.path.basename(self.input_filename)} has been modified.")
        try:  # to read the updated CSV file and extract the measurements to pass to the state estimator
            measurements_df = pd.read_csv(self.input_filename)  # No with-statement needed as pandas handles closing
            self.measurements = (measurements_df.iloc[:, 1:].to_numpy())  # convert to numpy array
            # Extract also the last row (most recent measurement) excluding the first column which holds the datetime
            self.current_measurements = measurements_df.iloc[-1, 1:].to_numpy()  # to pass the measurements to the MHE
            self.current_measurements = np.array(self.current_measurements, dtype=np.float64)  # ensure numpy floats
            print(f" Current measurements: {np.round(self.current_measurements, 2)}\n")
        except Exception as exe:
            # in case the file cannot be read use the last measure
            self.current_measurements = self.measurements[-1, :]
            print(f" Error reading file! Previous values used instead!\n {exe}\n")

    def estimate_model_parameters(self):
        # Extract only values form the last 3 hours to avoid fitting to lag phase data
        measured_time = self.measurements[-90:, 0]  # 90 values correspond to 3 hours since the time step is 2 min
        measured_cdw = self.measurements[-90:, -1]  # the last entry in the measurements is CDW pred. via CO2 integral
        measured_glucose = self.measurements[-90:, 2]
        measured_states = [measured_cdw, measured_glucose]

        # Define initial guesses for the parameters
        init_param = np.array([0.6871, 0.015, 0.40186])  # Initial guesses for [mu_max, Ks, Yxs]
        # Define bounds for the parameters
        bounds = [(0.6871 - 0.0332, 0.6871 + 0.0332),  # Bounds for mu_max are +- 1 std form the estimates via CDW, OD600 and offgas
                  (0.015, 0.1),  # Bounds for Ks, its very important that Ks is not too low to control the system
                  (0.40186 - 0.0043, 0.40186 + 0.0043)]  # Bounds for Yxs are +- 2 std

        # Obtain the optimal parameters using Sci Py minimize()
        result = minimize(cost_fun,
                          init_param,
                          args=(measured_states, measured_time),  # Pass additional arguments
                          bounds=bounds,
                          method='nelder-mead',  # non-gradient based solver
                          options={'disp': True})

        # Extract optimal parameters from the result
        self.mu_max_estimated, self.K_s_estimated, self.Y_xs_estimated = result.x

        print(f"\n Estimated Model Parameters:\n --> mu_max:\t {np.round(self.mu_max_estimated, 2)} h-1")
        print(f" --> Y_xs:\t {np.round(self.Y_xs_estimated, 2)} g/g")
        print(f" --> K_s:\t {np.round(self.K_s_estimated, 3)} g/L")
        # Save the estimated parameters to textfile to be able to look at them while the MPC is running
        try:
            with open("model_parameter_estimates.txt", "w") as file:
                # Write the estimated values to the file
                file.write(f"mu_max_estimated:\t {self.mu_max_estimated}\n")
                file.write(f"Y_xs_estimated:\t\t {self.Y_xs_estimated}\n")
                file.write(f"K_s_estimated:\t\t {self.K_s_estimated}\n")
                file.write(f"-------------------------------------\nTime of estimation:\t {pd.Timestamp.now()}\n")
        except Exception as e:
            print(f" Could not save estimated model parameters to .txt file!\n {e} \n")

    def on_modified(self, event):
        try:  # this try statement allows the script to be started even when the inputs file is not yet created
            self.input_filename = glob.glob(os.path.join(path, 'Measurement_Inputs*'))[0]
        except Exception as e:
            print(f"No CSV file named Measurement_Inputs_ found yet..\n {e} \n")
            return  # Return early if the file is not yet created
        if event.src_path != self.input_filename:
            return  # Return early if the event is not on the input file (i.e. when the output file is modified)
        print("\n-----------------------------------------------------------------------------------")  # indicate new measurement
        # Both writing and saving the file causes this function to trigger, but we only want one call to the function
        time.sleep(0.01)  # on_modified() detects also opening file, but the writer is not yet finished so wait a bit
        current_time = time.time()  # variable to mitigate the double call: After writing has finished it gets triggered
        if current_time - self.last_modified_time > self.debounce_delay:  # check time difference since last call
            self.last_modified_time = current_time  # update time
            self.read_measurements()  # Read the new measurements

        # Estimate model parameters and start MPC as soon as the glucose prediction is below 1 g/L
        if self.current_measurements[2] < 1 and not self.is_mpc_active:  # "and not" is there to only trigger it once
            try:
                self.estimate_model_parameters()
            except Exception as e:
                print(f"Could not fit parameters. Default values used!\n {e}\n")
                
            try:
                self.init_mpc_mhe(mu_max=self.mu_max_estimated, Y_xs=self.Y_xs_estimated, K_s=self.K_s_estimated)
                self.is_mpc_active = True
            except Exception as e:
                print(f"Could not initialize MPC!")
                raise

            # Copy the amount added pre MPC start to new variable to subtract it from the new measurements
            self.added_vol_feed_pre_mpc = self.current_measurements[6]
            self.added_vol_base_pre_mpc = self.current_measurements[7]
            self.added_vol_acid_pre_mpc = self.current_measurements[8]
        elif not self.is_mpc_active:
            print(f" Waiting for glucose to drop below 1 g/L. Current glucose conc. {np.round(self.current_measurements[2],2)} g/L")
            print("-----------------------------------------------------------------------------------\n")  # indicate step end

        # Enter MHE & MPC only as soon as glucose has dropped below a given threshold
        if self.is_mpc_active:
            self.k += 1  # increment time step index
            print(f" MPC Step Nr:\t {self.k}\n")
            try:
                # Extract latest estimates of CDW and working volume
                if self.k > 1:
                    current_cdw_estimate = self.MHE.data['_x', 'X'][-1]
                    current_vol_estimate = self.MHE.data['_x', 'V'][-1]
                    current_DOT_measurement = self.current_measurements[9]  # extract the DOT to prevent O2-limitation

                # Check conditions and update time-varying-parameters for the current time step accordingly
                if self.k > 1 and (current_cdw_estimate >= 80 or current_vol_estimate >= 2.6 or current_DOT_measurement < 5):  # the first condition is there to protect the second from failing in the first round when there is no data yet
                    self.tvp_template_mpc['_tvp', :, 'mu_set'] = 0
                    self.tvp_template_mhe['_tvp', :, 'mu_set'] = 0
                else:
                    # Create the (relative to MPC start) time vector spanning the current range of the horizon
                    horizon_time_mpc = self.time_step * np.arange(self.k, self.k + self.n_mpc_horizon + 1)
                    horizon_time_mhe = self.time_step * np.arange(self.k, self.k + self.n_mhe_horizon + 1)
                    # Get the raw shape of the mu set-point over the current horizon
                    y_mpc = np.exp(-self.mu_decay * horizon_time_mpc)
                    y_mhe = np.exp(-self.mu_decay * horizon_time_mhe)
                    # Min-Max scaling but with Min and Max of the whole reference trajectory
                    y_mpc_scaled = (y_mpc - self.y_min) / (self.y_max - self.y_min)
                    y_mhe_scaled = (y_mhe - self.y_min) / (self.y_max - self.y_min)
                    # Final mapping for the current horizon's reference trajectory
                    mu_set_mpc_current_horizon = self.mu_set_range * y_mpc_scaled + self.mu_set_end
                    mu_set_mhe_current_horizon = self.mu_set_range * y_mhe_scaled + self.mu_set_end
                    # Convert numpy array to list and paste the values to the TVP template
                    self.tvp_template_mpc['_tvp', :, 'mu_set'] = mu_set_mpc_current_horizon.tolist()
                    self.tvp_template_mhe['_tvp', :, 'mu_set'] = mu_set_mhe_current_horizon.tolist()

                # Paste the flow rates of the feed and base to the time varying parameter templates
                self.tvp_template_mpc['_tvp', :, 'F_base'] = self.measurements[-(self.n_mpc_horizon+1):, 4].tolist()
                self.tvp_template_mpc['_tvp', :, 'F_acid'] = self.measurements[-(self.n_mpc_horizon+1):, 5].tolist()
                
                # Paste the flow rates of the feed and base to the time varying parameter templates
                self.tvp_template_mhe['_tvp', :, 'F_base'] = self.measurements[-(self.n_mhe_horizon+1):, 4].tolist()
                self.tvp_template_mhe['_tvp', :, 'F_acid'] = self.measurements[-(self.n_mhe_horizon+1):, 5].tolist()

                # Estimate the current working volume based on the feed, base, and acid volumes added since MPC start
                V_feed_added = self.current_measurements[6] - self.added_vol_feed_pre_mpc
                V_base_added = self.current_measurements[7] - self.added_vol_base_pre_mpc
                V_acid_added = self.current_measurements[8] - self.added_vol_acid_pre_mpc
                V_estimate = 1.5 + V_feed_added + V_base_added + V_acid_added

                # Pass measurements to the MHE for state estimation: [CDW, Glucose, Volume, F_Feed, F_Base, F_Acid]
                measurements_for_mhe = np.hstack([self.current_measurements[1:3], V_estimate, self.current_measurements[3:6]])
                state_estimates = self.MHE.make_step(measurements_for_mhe)
                # Compute the optimal feed flow rate
                optimal_flow_rate = self.MPC.make_step(state_estimates)  # [L_feed/h]

            except Exception as worst_case:
                if self.k > 1:
                    optimal_flow_rate = self.ms / self.S_in * self.MHE.data['_x', 'X'][-1] * self.MHE.data['_x', 'V'][-1]
                    print(f" Error in do-mpc: {worst_case} \n Maintenance flow is applied: {optimal_flow_rate} L/h!\n")          
                else:  # special fail case: if the mpc fails directly in the first step
                    optimal_flow_rate = np.array([0])
                    print(f" Error in first do-mpc step: {worst_case} \n Flow rate set to: {optimal_flow_rate} L/h!\n")

            # Safety "if" statement that overrides the MPC output when the DOT falls below 1%
            if self.k > 1 and current_DOT_measurement < 1:
                optimal_flow_rate = np.array([0])  # Bypassing the MPC in case he does not react quick enough
                print(f" !!!Attention!!! \n Setting flow to 0 L/h because the current DOT is only {np.round(current_DOT_measurement, 2)} %\n !!!Bypassing MPC!!!\n")
            # Safety "if" statement that overrides the MPC output when the DOT falls below 1%
            if self.k > 1 and V_estimate > 2.75:
                optimal_flow_rate = np.array([0])  # Bypassing the MPC in case he does not react quick enough
                print(f" !!!Attention!!! \n Setting flow to 0 L/h because the estimated working volume is {np.round(V_estimate, 2)} L\n !!!Bypassing MPC!!!\n")

            # Send the command to the DAC which will apply the corresponding 0-10V signal to the pump
            send_command_to_pump(optimal_flow_rate.item())

            try:
                # Extract latest MHE and MPC data and append to CSV files
                mhe_relative_time_h = self.MHE.data['_time'][-1]
                mhe_state_estimates = self.MHE.data['_x'][-1, :]
                mhe_control_signal = self.MHE.data['_u'][-1, :]
                mhe_auxiliary_var = self.MHE.data['_aux'][-1, :]  # auxiliary variables [F_feed, F_base, F_acid, F_tot, mu]
                mhe_tvp_var = self.MHE.data['_tvp'][-1, :]  # time varying parameters [F_base, F_acid, mu_set]
                mhe_data_vector = np.hstack([mhe_relative_time_h, mhe_state_estimates, mhe_control_signal, mhe_auxiliary_var, mhe_tvp_var])

                mpc_relative_time_h = self.MPC.data['_time'][-1]
                mpc_state_estimates = self.MPC.data['_x'][-1, :]
                mpc_control_signal = self.MPC.data['_u'][-1, :]
                mpc_auxiliary_var = self.MPC.data['_aux'][-1, :]  # auxiliary variables [F_feed, F_base, F_acid, F_tot, mu]
                mpc_tvp_var = self.MPC.data['_tvp'][-1, :]  # time varying parameters [F_base, F_acid, mu_set]
                mpc_data_vector = np.hstack([mpc_relative_time_h, mpc_state_estimates, mpc_control_signal, mpc_auxiliary_var, mpc_tvp_var])
            except Exception as e:
                print(f" Could not extract data!\n{e}\n")
                mpc_data_vector = -99*np.ones(13)  # put -99 so one can easily identify missing values
                mhe_data_vector = -99*np.ones(13) 
                
            save_to_csv(mhe_data_vector, mhe_output_filename_string)
            save_to_csv(mpc_data_vector, mpc_output_filename_string)
            print("-----------------------------------------------------------------------------------\n")  # indicate step end


#######################################################################################################################
# Define the function to compute the ODEs
def fermentation_odes(state, t, mu_max, Ks, Yxs):
    X, S = state  # Biomass (X) and Substrate (S)
    mu = mu_max * S / (Ks + S)  # Monod kinetics
    dXdt = mu * X
    dSdt = -(mu / Yxs + 0.025) * X
    return [dXdt, dSdt]


#######################################################################################################################
# Simulate the ODEs for given parameters
def simulate_cdw_and_glucose(t, mu_max, Ks, Yxs, initial_state):
    state = odeint(fermentation_odes, initial_state, t, args=(mu_max, Ks, Yxs))
    return state[:, 0], state[:, 1]  # Return CDW (Biomass) and Glucose


#######################################################################################################################
# Define the cost function for the parameter estimation via Sci Py minimize()
def cost_fun(parameters, measured_states, measured_time):
    mu_max, Ks, Yxs = parameters  # Extract parameters to optimize
    cdw, glucose = measured_states  # Extract states

    # Simulate the ODE model
    initial_state = [cdw[0], glucose[0]]  # Initial [Biomass, Substrate]
    model_time = np.linspace(min(measured_time), max(measured_time), 10000)  # High-resolution time grid
    model_cdw, model_glucose = simulate_cdw_and_glucose(model_time, mu_max, Ks, Yxs, initial_state)

    # Interpolate model results to match the measured time points
    interp_cdw = np.interp(measured_time, model_time, model_cdw)
    interp_glucose = np.interp(measured_time, model_time, model_glucose)

    # Compute residuals (differences between measured and interpolated values)
    residuals_cdw = cdw - interp_cdw
    residuals_glucose = glucose - interp_glucose

    # Normalize residuals by the range of measured values
    cdw_range = np.max(cdw) - np.min(cdw)
    glucose_range = np.max(glucose) - np.min(glucose)

    normalized_residuals_cdw = residuals_cdw / cdw_range
    normalized_residuals_glucose = residuals_glucose / glucose_range

    # Define the cost function as the sum of squared normalized residuals
    cost = np.sum(normalized_residuals_cdw**2 + normalized_residuals_glucose**2)
    return cost


#######################################################################################################################
def send_command_to_pump(vol_flow_rate_command):
    """Write a single value to the file, opening and closing the file each time."""
    # Convert the volumetric flow rate calculated by the MPC into a mass flow rate, on which the pump is calibrated
    mass_flow_rate_command = vol_flow_rate_command * feed_density
    # Now we can convert the mass flow rate to the desired corresponding RPM
    rpm = feed_calibration_slope * mass_flow_rate_command  # convert mass flow rate to rpm

    print(f" Optimal Vol. Flow Rate:\t {round(vol_flow_rate_command, 3)} L/h")
    print(f" Optimal Mass Flow Rate:\t {round(mass_flow_rate_command, 3)} g/h\n")

    if rpm <= 0:  # handle the case where a negative flow rate is written to prevent a crash from the DAC
        print(f" Warning! Cannot set negative flow rate: {round(mass_flow_rate_command, 2)} g/h")
        print(f" Current8 Pump Command:\t\t 0 RPM")
        volts = 0
    elif rpm < max_rpm:
        print(f" Current Pump Command:\t\t {round(rpm, 2)} RPM")
        # volts = rpm / max_rpm * 10  # use this if you trust the internal RPM-to-Voltage calibration
        # volts = 0.39271*rpm + 0.06329 # when max RPM = 25
        volts = 0.620998 * rpm + 0.055793  # when max RPM = 16

    else:
        print(f" Warning! Desired flow rate cannot be met with current max RPM!")
        print(f" Current Pump Command:\t\t {round(rpm, 2)} RPM")
        volts = 10

    try:
        megaind.set0_10Out(0, 1, volts)  # Send the voltage to the pump via the DAC
        print(f" Current DAC Output:\t\t {round(volts, 2)} V")
        time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
        print(f" New command send! --> Task completed at: {time_date}!\n")
    except Exception as e:
        print(f"\n Error! Could not send desired flow rate command to the pump!\n {e}")


#######################################################################################################################
def save_to_csv(numpy_data_array, filename_string):
    """
    This function saves the outputs by the MPC and MHE to a .csv
    filename_string must be defined within the script
    :param numpy_data_array: The function needs a numpy array of shape (n,). One value for each measured quantity
    :param filename_string: Desired file name, if it already exists data will be appended to it otherwise it is created
    :return: None
    """
    try:  # in case something happens... better use try
        # Force data to be a row array to be able to add rows to a pandas df.
        column_data_array = numpy_data_array.reshape(1, -1)
        # If the file does not yet exist (or gets deleted) it will be (re)created here
        if not os.path.exists(filename_string):  # This is usually called only done on the first call of save_to_csv()
            # Create the CSV file and initialize with the column names and values for the first row
            column_names = ["MHE Time [h]", "CDW [g/L]", "Glucose [g/L]", "Volume [L]", "_u Feed Flow MPC", "Flow Rate [L/h]", "F_feed",
                            "F_base", "F_acid", "F_tot", "mu", "_tvp_F_base", "_tvp_F_acid", "_tvp_mu_set"]
            new_data_df = pd.DataFrame(column_data_array, columns=column_names)
            new_data_df.insert(0, 'TimeStamp', pd.Timestamp.now())  # add datetime in ISO 8601 format
            new_data_df.to_csv(filename_string, index=False)  # create the CSV
            # Signal successful creation of the CSV file
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f"\n CSV file, {filename_string}, created and initialized!\n --> Task completed at: {time_date}!")
        else:
            # Append values to the existing file
            new_data_df = pd.DataFrame(column_data_array)
            new_data_df.insert(0, 'TimeStamp', pd.Timestamp.now())
            with open(filename_string, mode='a', newline='') as file:
                new_data_df.to_csv(file, index=False, header=False)
            # Signal successful data treatment
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" Current outputs appended to {filename_string}\n --> Task completed at {time_date}!")
    except Exception as e:
        print(f" Could not store data to .csv!\n {e} \n")


#######################################################################################################################
# Safety measure to prevent accidental termination of the script (does not prevent termination if cmd window is closed!)
def signal_handler(sig, frame):
    """Handle the signal and ask for user confirmation."""
    user_input = input("\n Are you sure you want to exit? (y/n): ").strip().lower()
    if user_input == 'y':
        print("Exiting...")
        # Raise a KeyboardInterrupt to be caught by the except block
        raise KeyboardInterrupt
    else:
        print(" Continuing...\n")
        # Re-register the signal handler to continue
        signal.signal(signal.SIGINT, signal_handler)


#######################################################################################################################
def get_input_with_default(prompt, default):
    """Prompt the user for numeric input with a default value."""
    while True:
        user_input = input(f"{prompt} {default} or: ").strip()
        # Check if the input is empty (default value) when the user hits enter
        if not user_input:
            print(f" --> Default value used: {default} ")
            return default
        # Validate the input as a number
        try:
            numeric_value = float(user_input)  # Change to int(user_input) if only integers are allowed
            print(f" --> New value set: {numeric_value} ")
            return numeric_value
        except ValueError:
            print("Invalid input! Please enter a number.")


#######################################################################################################################
# Set up the watchdog:
if __name__ == "__main__":
    # Set the path to the folder to be monitored and where the CSV with the measurements will appear
    #path = "C:/Users/basil/Documents/Studium/Master_Applied_Biosciences/X03 - Master's Thesis/Scripts and Code/Python_Control_Algorithm"
    path = "/home/natto/MSc/Python_Scripts/MPC_Algorithm"
    mhe_output_filename_string = 'MHE_Outputs_' + datetime.now().strftime("%d_%m_%Y") + '.csv'
    mpc_output_filename_string = 'MPC_Outputs_' + datetime.now().strftime("%d_%m_%Y") + '.csv'
    
    print("-----------------------------------------------------------------------------------") 
    # Prompt the user for the density of the feed
    # default_density = 1185  # g_Feed/L_Feed of a 500 g_glucose/L_feed solution + 3.4 g/L MgSO4*7H2O
    default_density = 1262  # g_Feed/L_Feed of a 700 g_glucose/L_feed solution + 3.4 g/L MgSO4*7H2O

    feed_density = get_input_with_default("\n Set density of feed solution to [g/L]:", round(default_density, 2))

    # Prompt the user for the feed calibration and the max rpm for the pump
    default_feed_calibration_slope = 0.1236677  # [RPM/(g_feed/L_feed)]
    default_max_rpm = 16  # [RPM]
    feed_calibration_slope = get_input_with_default("\n Set feed calibration slope to:", default_feed_calibration_slope)
    max_rpm = get_input_with_default("\n Set maximal pump RPM to:", default_max_rpm)

    # Register the signal handler to protect the script from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Create an event handler for the watchdog which will activated when new data arrives
    event_handler = CustomMPC()
    observer = Observer()
    observer.schedule(event_handler, path=path, recursive=False)
    observer.start()
    # Print the time when the watchdog was activated
    start_time_date = datetime.now()
    print('\n\n ----\n ------')
    print(' MPC Watchdog activated! Start at:', start_time_date.strftime("%H:%M:%S on %d-%m-%Y"))
    print(' ------\n ----\n')

    try:
        while True:
            time.sleep(1)  # Define how long the watchdog "sleeps" before checking again for new files
    except KeyboardInterrupt:  # Press Ctrl + C in the command window to stop monitoring the folder
        print('\n Watchdog goes to sleep..zZz..')
        observer.stop()  # Stop monitoring the folder
        observer.join()  # Wait for the observer's thread to finish its execution before proceeding
        print(' Watchdog has fallen asleep! Target folder is no longer monitored!')






