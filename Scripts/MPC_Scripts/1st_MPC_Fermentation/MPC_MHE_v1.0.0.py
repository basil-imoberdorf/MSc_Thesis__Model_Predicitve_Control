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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
# Note: this script must be executed within the virtual environment on the Raspberry Pi to import do-mpc and casadi


#######################################################################################################################
# The MPC is build with a watchdog FileHandler as parent class to dynamically read measurements and pass it to the MPC
class CustomMPC(FileSystemEventHandler):
    def __init__(self):
        self.input_filename = glob.glob(os.path.join(path, 'Measurement_*'))[0]  # find the file that starts with 'Meas'
        self.last_modified_time = 0  # helper variable used to prevent calling on_modified() tice when new values arrive
        self.debounce_delay = 1  # Time in seconds to wait before accepting new changes, on_modified fires twice usually
        self.measurements = pd.read_csv(self.input_filename).iloc[:, 1:].to_numpy()  # load all available measurements
        self.k = 0  # variable to keep track of the current time stap, will increase by 1 every iteration
        
        # Parameters for the dynamic mu_max set point: linear decrease form mu_start to mu_final at t_final
        self.t_final = 20  # [h]
        self.mu_set_start = 0.3  # [h-1]
        self.mu_set_final = 0.1  # [h-1]
        self.mu_slope = (self.mu_set_start - self.mu_set_final)/self.t_final
        # Time step for the algorithm
        self.time_step = 2 / 60  # [h] every 2 minutes a measurement becomes available therefore time step = 2/60 h
        self.n_mpc_horizon = 20      

        # --------------------------------------------------------------------------------------------------------------
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
        # Flow rate measurements, I assume no noise on flow rates (MHE has already enough burden)
        model.set_meas('F_f_meas', F_feed, meas_noise=False)
        model.set_meas('F_b_meas', F_base, meas_noise=False)
        model.set_meas('F_a_meas', F_acid, meas_noise=False)

        # Known model parameters estimated and hardcoded manually after the batch phase in a separate file
        mu_max = 0.5024  # max specific growth rate [h-1]
        Y_xs = 0.4284   # Growth yield coefficient on substrate [g_bio/g_subst]
        K_s = 1     # Monod substrate saturation constant   [h-1]
        ms = 0.025      # Specific maintenance coefficient [g/(gh)]
        S_in = 500.0    # concentration of substrate in feed solution [g_subst/L_feed]

        # Auxiliary terms
        mu = mu_max * S / (K_s + S)  # mu = mu_max*S/(K_s + S + (S**2/K_i))
        F_tot = F_feed + F_base + F_acid  # total flow rate into the bioreactor, sampling is not accounted for
        # Using set_expression makes it possible to store and plot auxiliary variables
        mu = model.set_expression('mu', mu)
        F_tot = model.set_expression('F_tot', F_tot)

        # Right-hand-side (rhs) of the ordinary differential equations
        model.set_rhs('X', mu*X - F_tot/V*X, process_noise=True)
        model.set_rhs('S', -(mu/Y_xs + ms)*X + F_feed/V*S_in - F_tot/V*S, process_noise=True)
        model.set_rhs('V', F_tot, process_noise=True)
        model.setup()  # Build the model
        # --------------------------------------------------------------------------------------------------------------
        # Initialize Moving Horizon State Estimator
        self.MHE = do_mpc.estimator.MHE(model)  # no model parameters are estimated, only states
        setup_mhe = {
            'n_horizon': 10,  # Number of time steps for the estimation window
            't_step': self.time_step,  # Sampling time in hours
            'state_discretization': 'collocation',
            'store_full_solution': True,
            'meas_from_data': True  # this means that the MHE takes the measurements from his own data storage variable
                                    # The MHE appends each incoming measurement to a "data" object. See documentation
        }
        self.MHE.set_param(**setup_mhe)

        # Define the diagonal weighting/penalty matrices
        P_v = 0.01 * np.diag(np.array([1, 1]))  # measurement noise penalty matrix
        P_w = 1 * np.diag(np.array([1, 1, 1]))  # process noise penalty matrix
        P_x = 1 * np.diag(np.array([1, 1, 1000]))  # to penalize deviations from the states

        # Define objective function for MHE
        self.MHE.set_default_objective(P_x=P_x, P_v=P_v, P_w=P_w)  # see documentation for the details

        # State and parameter constraints/bounds:
        self.MHE.bounds['lower', '_x', 'X'] = 0
        self.MHE.bounds['lower', '_x', 'S'] = 0
        self.MHE.bounds['lower', '_x', 'V'] = 1.5

        # mhe.bounds['upper', '_x', 'X'] = 50
        self.MHE.bounds['upper', '_x', 'S'] = 25
        self.MHE.bounds['upper', '_x', 'V'] = 2.5

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
        self.MPC.set_param(**setup_mpc)  # setting parameters can also be set like this: mpc.settings.n_horizon = 20

        # Define objective (cost) function
        mterm = (model.aux['mu'] - model.tvp['mu_set'])**2  # terminal cost, Meyer term
        lterm = (model.aux['mu'] - model.tvp['mu_set'])**2  # stage cost, Lagrange term
        self.MPC.set_objective(mterm=mterm, lterm=lterm)  # to omit a term use: mterm=DM(0) or lterm=DM(0)
        self.MPC.set_rterm(F_feed=10.0)  # set here the diagonal elements of R from rterm = delta_u'*R*delta_u

        # Define Constraints (lower/upper bounds) on control inputs and states
        self.MPC.bounds['lower', '_x', 'X'] = 0.0
        self.MPC.bounds['lower', '_x', 'S'] = -0.01
        self.MPC.bounds['lower', '_x', 'V'] = 1.5
        # upper bounds of the states
        self.MPC.bounds['upper', '_x', 'S'] = 2.0
        self.MPC.bounds['upper', '_x', 'V'] = 2.5
        # upper and lower bounds of the control input
        self.MPC.bounds['lower', '_u', 'F_feed'] = 0.0
        self.MPC.bounds['upper', '_u', 'F_feed'] = 0.17

        # Time-Varying-Parameters for MHE
        self.tvp_template_mpc = self.MPC.get_tvp_template()
        # Create dummy function simply returning the tvp_template structure
        tvp_fun_mpc = lambda t_now: self.tvp_template_mpc  # GOAT --> https://github.com/do-mpc/do-mpc/discussions/240
        self.MPC.set_tvp_fun(tvp_fun_mpc)
        self.MPC.settings.supress_ipopt_output()  # increase speed
        self.MPC.setup()
        # --------------------------------------------------------------------------------------------------------------
        # Initial Conditions hardcoded after batch phase
        X_0 = 8.27  # Biomass concentration    [g/L_culture]
        S_0 = 0.1  # Substrate Concentration  [g/L_culture]
        V_0 = 1.5  # Culture working volume inside the bioreactor [L]
        x0 = np.array([[X_0], [S_0], [V_0]])
        self.MPC.x0 = x0
        self.MHE.x0 = x0
        self.MPC.set_initial_guess()
        self.MHE.set_initial_guess()
        # --------------------------------------------------------------------------------------------------------------

    def on_modified(self, event):
        if event.src_path != self.input_filename:
            return  # Return early if the event is not for the input file (when the output file is modified or else)

        # Both writing and saving the file causes this function to trigger. We only need one activation
        time.sleep(0.01)  # on_modified() detects also opening file, but the writer is not yet finished so wait a bit
        current_time = time.time()  # variable to mitigate the double call: After writing has finished it get triggered
        if current_time - self.last_modified_time > self.debounce_delay:  # check time difference since last call
            self.last_modified_time = current_time  # update time
            print("\n-------------------------------------------------------------------------")  # indicate step start
            print(f" {os.path.basename(self.input_filename)} has been modified.")

            try:  # to read the updated CSV file and extract the measurements to pass to the state estimator
                measurements_df = pd.read_csv(self.input_filename)  # No with-statement needed as pandas handles closing
                # Extract the last row (most recent measurement) excluding the first column which holds the timestamp
                self.measurements = (measurements_df.iloc[:, 1:].to_numpy())  # to pass the base and acid flow to the MHE
                current_measurements = measurements_df.iloc[-1, 1:].to_numpy()  # to pass the measurements to the MHE
                current_measurements = np.array(current_measurements, dtype=np.float64)  # auf nummer sicher gehen
                print(f" Current measurements: {np.round(current_measurements,2)}")
            except Exception as exe:
                current_measurements = self.measurements[-1, :]  # in case the file cannot be rad use the last measure
                print(f" Error reading file:\n {exe}\n")

            try:
                self.k += 1
                # Update time-varying-parameters for the current time step
                if self.k * self.time_step > 10 and self.MHE.data['_x', 'X'][-1] >= 72:  # the first condition is also there to protect the second from failing in the first round when there is not data yet
                    self.tvp_template_mpc['_tvp', :, 'mu_set'] = 0.0
                    self.tvp_template_mhe['_tvp', :, 'mu_set'] = 0.0
                else:  # create the linearly decreasing mu reference for the current horizon
                    self.tvp_template_mpc['_tvp', :, 'mu_set'] = [self.mu_slope * i * self.time_step + self.mu_set_start for i in range(self.k, self.k + self.n_mpc_horizon + 1)]
                    self.tvp_template_mpc['_tvp', :, 'mu_set'] = [self.mu_slope * i * self.time_step + self.mu_set_start for i in range(self.k, self.k + self.n_mpc_horizon + 1)]

                # Pass measurements to the MHE for state estimation
                states_prediction = self.MHE.make_step(current_measurements)

                # Get the optimal control input from the MPC
                optimal_flow_rate = self.MPC.make_step(states_prediction)

            except Exception as worst_case:  # bruh please dont
                print(f" Error in do-mpc:\n {worst_case} \n")

            # Write the optimal control input to the FlowRateCommand.txt file to send the command to the DAC via autoDAC
            write_command(optimal_flow_rate.item())  # extract the scalar value from the np array

            # Signal successful data treatment
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" Setpoint written to Command File at: {time_date}")

            # Save MHE and MPC data to CSV file
            save_to_csv(np.vstack([optimal_flow_rate, states_prediction]))
            print("------------------------------------------------------------------\n")  # indicate step end


#######################################################################################################################
def write_command(vol_flow_rate_command):
    """Write a single value to the file, opening and closing the file each time."""
    # Convert the volumetric flow rate calculated by the MPC into a mass flow rate, on which the pump is calibrated
    mass_flow_rate_command = vol_flow_rate_command*feed_density
    with open("FlowRateCommand.txt", "w") as file:
        file.write(f"{mass_flow_rate_command}")
    print(f"\n Vol. Flow Rate computed: {np.round(vol_flow_rate_command, 3)} L/h")
    print(f" Mass Flow Rate Setpoint: {np.round(mass_flow_rate_command, 3)} g/h\n")


#######################################################################################################################
def save_to_csv(numpy_data_array):
    """
    This function saves the aggregated measurements received from the BioSCADA platform and the ReactRaman to a .csv
    filename_string must be defined within the script
    :param numpy_data_array: The function needs a numpy array of shape (4,). One value for each measured quantity
    :return: None
    """
    try:  # in case something happens... better use try
        # Force data to be a row array to be able to add rows to a pandas df.
        column_data_array = numpy_data_array.reshape(1, -1)
        # If the file does not yet exist (or gets deleted) it will be (re)created here
        if not os.path.exists(output_filename_string):  # This is usually called only done on the first call of save_to_csv()
            # Create the CSV file and initialize with the column names and values for the first row
            column_names = ["Flow Rate [L/h]", "CDW [g/L]", "Glucose [g/L]", "Volume [L]"]
            new_data_df = pd.DataFrame(column_data_array, columns=column_names)
            new_data_df.insert(0, 'TimeStamp', pd.Timestamp.now())  # add datetime in ISO 8601 format
            new_data_df.to_csv(output_filename_string, index=False)  # create the CSV
            # Signal successful creation of the CSV file
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f"\n CSV file, {output_filename_string}, created and initialized!\n --> Task completed at: {time_date}!")
        else:
            # Append values to the existing file
            new_data_df = pd.DataFrame(column_data_array)
            new_data_df.insert(0, 'TimeStamp', pd.Timestamp.now())
            with open(output_filename_string, mode='a', newline='') as file:
                new_data_df.to_csv(file, index=False, header=False)
            # Signal successful data treatment
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(
                f" Current MPC and MHE outputs appended to {output_filename_string}\n --> Task completed at {time_date}!")
    except Exception as e:
        print(f" Could not store glucose prediction to .csv!\n {e} \n")


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
        user_input = input(f"{prompt} {default} g/L or: ").strip()
        # Check if the input is empty (default value)
        if not user_input:
            print(f" --> Default value used: {default} g/L")
            return default
        # Validate the input as a number
        try:
            numeric_value = float(user_input)  # Change to int(user_input) if only integers are allowed
            print(f" --> New value set: {numeric_value} g/L")
            return numeric_value
        except ValueError:
            print("Invalid input! Please enter a number.")


#######################################################################################################################
# Set up the watchdog:
if __name__ == "__main__":
    # Set the path to the folder to be monitored and where the CSV with the measurements will appear
    #path = "C:/Users/basil/Documents/Studium/Master_Applied_Biosciences/X03 - Master's Thesis/Scripts and Code/Python_Control_Algorithm"
    path = "/home/natto/MSc/Python_Scripts/MPC_Algorithm"
    output_filename_string = 'MPC_MHE_Outputs_' + datetime.now().strftime("%d_%m_%Y") + '.csv'

    # Prompt the user for the density of the feed
    default_density = 1185  # [g_feed/L_feed]
    feed_density = get_input_with_default("\n Set density of feed solution to:", default_density)

    # Register the signal handler to protect the watchdog from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Create an event handler for the watchdog
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






















