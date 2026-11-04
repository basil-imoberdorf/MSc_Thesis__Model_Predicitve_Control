import os
import time
import signal
import socket
import numpy as np
import pandas as pd
import PIconnect as PI
from PIconnect.PIConsts import AuthenticationMode
from dotenv import load_dotenv

# PI-server credentials are read from a .env file (not part of the repository, see .env.example).
# load_dotenv() searches the folder of this script and its parent folders for the .env file.
load_dotenv()
from datetime import datetime

import matplotlib.pyplot as plt
from matplotlib.cm import plasma
from matplotlib.colors import Normalize


######################################################################################################################
# Safety measure to prevent accidental termination of the script (does not prevent termination if cmd window is closed!)
def signal_handler(sig, frame):
    """Handle the signal and ask for user confirmation."""
    user_input = input("\n Are you sure you want to exit? (y/n): ").strip().lower()
    if user_input == 'y':
        print("Exiting...")
        # Raise a KeyboardInterrupt to be caught by the except block
        send_command_to_pump(np.array([0]), 0)  # before stopping set the voltage to zero by commanding 0 RPM
        send_command_to_pump(np.array([420.69]), 420)  # send the special value to stop the pump
        raise KeyboardInterrupt
    else:
        print(" Continuing...\n")
        # Re-register the signal handler to continue
        signal.signal(signal.SIGINT, signal_handler)


#######################################################################################################################
def get_input_with_default(prompt, default):
    """Prompt the user for numeric input with a default value."""
    while True:
        user_input = input(f"{prompt} \033[93m\033[1m{default} RPM\033[0m or: ").strip()
        # Check if the input is empty (default value)
        if not user_input:
            print(f" --> Default value used: \033[92m\033[1m{default} RPM\033[0m\n")
            return default
        # Validate the input as a number
        try:
            numeric_value = float(user_input)  # Change to int(user_input) if only integers are allowed
            print(f" --> New value set: \033[92m\033[1m{numeric_value} RPM\033[0m\n")
            return numeric_value
        except ValueError:
            print("\033[91mInvalid input! Please enter a number.\033[0m\n")


#######################################################################################################################
def get_PI_server_value(tag_string):
    max_retries = 100  # I do not want to try indefinitely, because the rest of the algorithm must proceed
    retries = 0
    while retries < max_retries:  # this loop should make it more robust. If fetching the value fails, retry 100x
        try:
            with PI.PIServer(server=os.getenv('PI_SERVER', 'VMENPPIAFSERVER'),
                             username=os.getenv('PI_USERNAME'),
                             password=os.getenv('PI_PASSWORD'),
                             authentication_mode=AuthenticationMode.WINDOWS_AUTHENTICATION) as server:
                tag = server.search(f"*{tag_string}*")[0]
                # data = tag.recorded_value('-0m', retrieval_mode=RetrievalMode.AT_OR_BEFORE)  # get latest measurement
                # Extract the (interpolated) measurements from th current time and from two minutes ago, 2min = mpc delta_t
                data = tag.interpolated_values('-2m', '-0m',
                                               '2m')  # get two latest measurements, cant get a single one...
                measure_time = data.index[-1]  # extract the timestamp of the latest measurement
                measure_value = data.iloc[-1]  # extract the latest measurement
            return measure_time, measure_value

        except Exception as e:
            print(f" Could not retrieve data from PI server!\n Error: {e}\n Try again...\n")
            retries += 1
            time.sleep(0.01)  # Short delay between retries


#######################################################################################################################
def send_command_to_pump(data_to_send, RPM):
    """
    Sends a numpy array via TCP-IP to a client server
    :param data_to_send: ndarray
    :param RPM:
    :return: None
    """
    # Server address configuration, where to go with the data
    raspi_host = '192.168.0.200'  # <-- target IP address. For Raspberry Pi works also: "raspberrypi"
    raspi_port = 42069

    # Try to create a socket object and send the value to the Raspberry Pi
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as client_socket:
            # Convert numpy array to bytes
            array_data = data_to_send.tobytes()
            # Connect to the server and send data
            client_socket.connect((raspi_host, raspi_port))
            client_socket.sendall(array_data)
            time_date = datetime.now()
            print(f" {RPM} RPM command sent as voltage to Raspberry Pi at ", time_date.strftime("%H:%M:%S on %d.%m.%Y"))
            response = client_socket.recv(1024).decode()
            print(f" Raspberry Pi response: {response}\n")

    except Exception as e:
        print(f" Communication Error: Could not send the data!:\n {e} \n")


#######################################################################################################################
def linear_regression(x, y, forced_intercept):
    # Perform "normal" linear regression
    if forced_intercept is None:
        # Perform standard linear regression
        n = np.size(x)
        x_mean = np.mean(x)
        y_mean = np.mean(y)
        # Calculate slope and intercept
        slope = np.sum((x - x_mean) * (y - y_mean)) / np.sum((x - x_mean) ** 2)  # this is the same as Cov(x,y)/Var(x)
        intercept = y_mean - slope * x_mean

        # Calculate coefficient of determination to assess the quality of the fit
        residuals = y - (slope*x + intercept)
        R2 = 1 - np.var(residuals)/np.var(y)
        adj_R2 = 1 - (1 - R2) * (n - 1) / (n - 2 - 1)  # (n - 1) / (n - num_model_param - 1)

        # Group and return results as a list
        fit_results = [slope, intercept, R2, adj_R2, n]
        return fit_results

    # Linear regression where the intercept is forced through a constant
    else:
        # Compute the slope
        n = np.size(x)
        slope = (sum(y*x) - forced_intercept * n * np.mean(x))/sum(x**2)

        # Calculate coefficient of determination to assess the quality of the fit
        residuals = y - (slope*x + forced_intercept)
        R2 = 1 - np.var(residuals)/np.var(y)
        adj_R2 = 1 - (1 - R2) * (n - 1) / (n - 1 - 1)  # (n - 1) / (n - num_model_param - 1)

        # Group and return results as a list
        fit_results = [slope, R2, adj_R2, n]
        return fit_results


#######################################################################################################################
if __name__ == "__main__":
    # Register the signal handler to protect the script from accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Prompt the user for the feed calibration and the max rpm for the pump
    os.system('color')  # enables on Windows to print colored values to the command window
    default_max_rpm = 25  # [RPM]
    max_rpm = get_input_with_default("\n Set maximum pump RPM to:", default_max_rpm)

    # for now the points for which the pump is calibrated are hardcoded here
    rpm_calib_points_lib = np.array([1, 2, 4, 8, 12, 16, 20, 25, 30, 35, 40, 45, 50])
    rpm_calib_points = rpm_calib_points_lib[rpm_calib_points_lib <= max_rpm]

    # Presets
    pumped_mass_threshold = 6  # [g] with this threshold, low RPM settings read more points to get accurate estimates
    num_points_threshold = 40  # define here the number of values to be collected to get a good fit for the flow rates
    sample_interval = 20  # [s] define the interval at which the scale data is fetched from the PI server, not to fast!
    flow_fit_rows = []  # preallocate an empty list for the results of the flow rate fit
    flow_fit_forcedIntercept_rows = []  # preallocate list for the regression results where the intercept is forced
    index = 0
    plt.figure(figsize=(10, 6))
    start_time = pd.Timestamp.now()

    print(f"----------------------------------------------------")
    print(f" Start of calibration:\t {start_time}")
    print(f"----------------------------------------------------\n")
    # Loop over all calibration points and calculate the added mass over the corresponding time
    for rpm in rpm_calib_points:
        # (Re)set the counter and storage variables
        pumped_mass = 0
        calibration_time = 0
        time_array = np.array([0])
        weight_array = np.array([0])

        # Indicate when the data collection for a cliabration point has been finished
        if rpm > 1:
            print(f" Data collection completed for current calibration point!")
            print(f"----------------------------------------------------\n")

        # Wait a bit between switching form one RPM setting to the next to let the scale settle
        print(' Wait 10s for the scale to reach steady state...\n')
        time.sleep(10)
        print(f"----------------------------------------------------")
        print(f" Start of mass flow estimation for {rpm} RPM.")

        # Get the current absolute scale reading to subtract it form the coming ones to get only the mass added
        initial_time, initial_scale_reading = get_PI_server_value("BR07.read:A1_Scale_PV")

        # convert rpm to volt
        volts = np.array(rpm / max_rpm * 10)  # convert rpm to voltage
        #if 0 < rpm < max_rpm:
        #    # volts = 0.39721*rpm + 0.06329  # if 500g/L feed
        #    volts = 0.620998 * rpm + 0.055793  # if 700g/L feed
        #elif rpm > max_rpm:
        #    volts = 10
        #else:
        #    volts = 0

        send_command_to_pump(volts, rpm)  # Can't fetch initial time and start pump at same time but delay should be short

        print(f" DAC output voltage set to: {volts}V")
        print(" Start collecting data...\n")
        # Collect data until a minimum amount of feed has been pumped or a minimum number of points has been collected
        while pumped_mass < pumped_mass_threshold or np.size(weight_array) < num_points_threshold:
            # Wait before fetching the data
            time.sleep(sample_interval)

            # Fetch the current time and scale reading
            current_time, current_scale_reading = get_PI_server_value("BR07.read:A1_Scale_PV")

            # Calculate added mass up to now for the current calibration point and the time spend
            pumped_mass = current_scale_reading - initial_scale_reading  # [g]
            calibration_time = (current_time - initial_time).total_seconds() / 3600  # [h]

            # Append the data for the current flow rate estimation
            time_array = np.append(time_array, calibration_time)
            weight_array = np.append(weight_array, pumped_mass)

            if calibration_time > 120:  # if after 2h the pumped_mass_threshold has not been reached go to next point
                print(f" {120/60} h timeout reached at the {rpm} RPM calibration point!\n Moving to the next one...\n")
                break

        # Fit a linear model to the data to obtain the slope which corresponds to the mass flow rate at that RPM
        flow_fit_results = linear_regression(time_array, weight_array, forced_intercept=None)  # free intercept
        flow_fit_rows.append({"Slope":       flow_fit_results[0],
                              "Intercept":   flow_fit_results[1],
                              "R^2":         flow_fit_results[2],
                              "adj_R^2":     flow_fit_results[3],
                              "Nr. Samples": flow_fit_results[4]})

        # Fit a linear model with intercept forced through zero.
        flow_fit_forcedIntercept_results = linear_regression(time_array, weight_array, forced_intercept=0)
        flow_fit_forcedIntercept_rows.append({"Slope":       flow_fit_forcedIntercept_results[0],
                                              "R^2":         flow_fit_forcedIntercept_results[1],
                                              "adj_R^2":     flow_fit_forcedIntercept_results[2],
                                              "Nr. Samples": flow_fit_forcedIntercept_results[3]})

        try:  # safety try statement in case the plasma colormap fails not the whole calibration is lost
            # Plot the collected data and the fitted line (forced through zero) for each calibration point
            plt.scatter(time_array, weight_array, color='black', marker='o', edgecolor='black', facecolor='none', s=80)
            X = np.array([0, time_array[-1]])
            # Y = flow_fit_results[0] * X + flow_fit_results[1]
            Y = flow_fit_forcedIntercept_results[0] * X
            # Create a colormap and normalize the range
            norm = Normalize(vmin=0, vmax=len(rpm_calib_points)-1)
            plt.plot(X, Y, linewidth=3.5, color='white')  # this is for better contrast on the scatter markers
            plt.plot(X, Y, label=f"{rpm} RPM", linewidth=3, color=plasma(norm(index)))
            index += 1
        except Exception as e:
            print(f"Error while plotting the {rpm} RPM flow rate data and model!\n {e}")

    # Send the stop signal to the pump
    print(f" Data collection completed for last calibration point!")
    print(f"----------------------------------------------------\n\n")
    send_command_to_pump(np.array([420.69]), 420)

    # Show the collected data and linear models for the flow rate estimation
    plt.xlabel('Pumping Time [h]', fontsize=15)
    plt.ylabel('Pumped Mass [g]',  fontsize=15)
    plt.title('Flow Rate Estimation', fontsize=17)
    plt.legend(fontsize=15)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(f"Pump_Flow_Rates_{datetime.now().date()}.png", dpi=300, transparent=False, bbox_inches='tight')
    plt.show()

    # Now finally fit a linear model to the RPM and estimated mass flow rates
    flow_fit_df = pd.DataFrame(flow_fit_rows)
    forced_flow_fit_df = pd.DataFrame(flow_fit_forcedIntercept_rows)  # only the-forced-through-zero models are used
    rpm_fit_results = linear_regression(np.array(forced_flow_fit_df["Slope"]), rpm_calib_points, forced_intercept=None)
    rpm2flow_slope = rpm_fit_results[0]
    rpm2flow_intercept = rpm_fit_results[1]
    rpm2flow_R2 = rpm_fit_results[2]
    rpm2flow_adj_R2 = rpm_fit_results[3]

    # Fit again this time with intercept forced through zero
    rpm_fit_results_ForcedThroughZero = linear_regression(np.array(forced_flow_fit_df["Slope"]), rpm_calib_points, forced_intercept=0)
    rpm2flow_slope_ForcedThroughZero = rpm_fit_results_ForcedThroughZero[0]
    rpm2flow_R2_ForcedThroughZero = rpm_fit_results_ForcedThroughZero[1]
    rpm2flow_adj_R2_ForcedThroughZero = rpm_fit_results_ForcedThroughZero[2]

    print(f" Conversion factor WITHOUT intercept! (forced through zero):")
    print(f"--------------------------------------------------------------")
    print(f" RPM conversion factor is:\t  {rpm2flow_slope_ForcedThroughZero}")
    print(f" R^2 for the fit is:\t\t  {np.round(rpm2flow_R2_ForcedThroughZero, 6)}")
    print(f" Adjusted R^2 is:\t\t  {np.round(rpm2flow_adj_R2_ForcedThroughZero, 6)}\n")

    print(f" Conversion factor WITH intercept:")
    print(f"--------------------------------------------------------------")
    print(f" RPM conversion factor is:\t  {rpm2flow_slope}")
    print(f" Undesired intercept is:\t  {rpm2flow_intercept}")
    print(f" R^2 for the fit is:\t\t  {np.round(rpm2flow_R2, 6)}")
    print(f" Adjusted R^2 is:\t\t  {np.round(rpm2flow_adj_R2, 6)}\n")
    print(f" Calibrated on:\t\t\t\t  {len(rpm_calib_points)} points\n")
    print(f"--------------------------------------------------------------")
    print(f" End of calibration:\t  {pd.Timestamp.now()}")
    print(f"--------------------------------------------------------------\n")

    # Write the estimated values to a textfile
    with open(f"Pump_AutoCalibration_Report_{datetime.now().date()}.txt", "w") as file:
        file.write(f"===============================================================\n")
        file.write(f" Conversion factor WITHOUT intercept! (forced through zero):\n")
        file.write(f"---------------------------------------------------------------\n")
        file.write(f" RPM conversion factor:\t\t  {rpm2flow_slope_ForcedThroughZero}\n")
        file.write(f" R^2 for the fit:\t\t  {np.round(rpm2flow_R2_ForcedThroughZero, 6)}\n")
        file.write(f" Adjusted R^2 is:\t\t  {np.round(rpm2flow_adj_R2_ForcedThroughZero, 6)}\n")
        file.write(f"===============================================================\n\n\n")

        file.write(f"===============================================================\n")
        file.write(f" Conversion factor WITH intercept:\n")
        file.write(f"---------------------------------------------------------------\n")
        file.write(f" RPM conversion factor:\t\t {rpm2flow_slope}\n")
        file.write(f" Undesired intercept:\t\t {rpm2flow_intercept}\n")
        file.write(f" R^2 for the fit:\t\t {np.round(rpm2flow_R2,6)}\n")
        file.write(f" Adjusted R^2 is:\t\t {np.round(rpm2flow_adj_R2, 6)}\n")
        file.write(f"===============================================================\n\n\n")
        file.write(f"===============================================================\n")
        file.write(f" Calibrated on:\t\t\t {len(rpm_calib_points)} points\n")
        file.write(f" Start calibration:\t\t {start_time}\n")
        file.write(f" End of calibration:\t\t {pd.Timestamp.now()}\n")
        total_time_seconds = pd.Timestamp.now() - start_time
        hours = total_time_seconds.seconds // 3600  # Integer division to get the number of full hours
        minutes = (total_time_seconds.seconds % 3600) // 60  # Modulo operation to get the remaining minutes
        file.write(f" Calibration took:\t\t {hours} h and {minutes} min!\n")
        file.write(f"===============================================================\n\n\n")

        file.write(f"=================================================================================================\n")
        file.write(f" Results of the mass flow rate estimations with intercept forced through zero:\n")
        file.write(f"=================================================================================================")
        file.write(f"\n|\tRPM\t|\tSlope\t|\tR^2\t\t|\tadj_R2\t|\tNr. Samples\t|\n")
        file.write(f"-------------------------------------------------------------------------------------------------\n")
        for i in range(len(rpm_calib_points)):
            rpm = rpm_calib_points[i]
            slope = np.round(forced_flow_fit_df.iloc[i, 0], 3)
            R2 = np.round(forced_flow_fit_df.iloc[i, 1], 5)
            adj_R2 = np.round(forced_flow_fit_df.iloc[i, 2], 5)
            num_samples = forced_flow_fit_df.iloc[i, 3]
            file.write(f"|\t{rpm}\t|\t{slope}\t|\t{R2}\t\t|\t{adj_R2}\t|\t{num_samples}\t\t|\n")

        file.write(f"=================================================================================================\n")
        file.write(f" Results of the mass flow rate estimations with intercept:\n")
        file.write(f"========================================================================================================================")
        file.write(f"\n|\tRPM\t|\tSlope\t|\tIntercept\t|\tR^2\t\t|\tadj_R2\t|\tNr. Samples\t|\n")
        file.write(f"------------------------------------------------------------------------------------------------------------------------\n")
        for i in range(len(rpm_calib_points)):
            rpm = rpm_calib_points[i]
            slope = np.round(flow_fit_df.iloc[i, 0], 3)
            intercept = np.round(flow_fit_df.iloc[i, 1], 3)
            R2 = np.round(flow_fit_df.iloc[i, 2], 5)
            adj_R2 = np.round(flow_fit_df.iloc[i, 3], 5)
            num_samples = flow_fit_df.iloc[i, 4]
            file.write(f"|\t{rpm}\t|\t{slope}\t|\t{intercept}\t\t|\t{R2}\t\t|\t{adj_R2}\t|\t{num_samples}\t\t|\n")
        file.write(f"========================================================================================================================\n\n")

    # Display results of the calibration graphically
    plt.figure(figsize=(10, 6))
    plt.scatter(flow_fit_df["Slope"], rpm_calib_points, s=400, color='k', marker='+', linewidth=3)
    X = np.array([0, max(flow_fit_df["Slope"])])
    Y_WithIntercept = rpm2flow_slope*X + rpm2flow_intercept
    Y_ForcedThroughZero = rpm2flow_slope_ForcedThroughZero*X
    plt.plot(X, Y_WithIntercept, color='cyan', label='with intercept', linewidth=4)
    plt.plot(X, Y_ForcedThroughZero, color='red', label='zero intercept', linewidth=2)
    # plt.plot([-10,200], [0,0], color='k', linewidth=1)
    # plt.plot([0,0], [-1,25], color='k', linewidth=1)
    plt.plot([-4,4], [0,0], color='k', linewidth=1)
    plt.plot([0,0], [-1,1], color='k', linewidth=1)

    # Format plot
    plt.xlabel('Mass Flow Rate [g/h]', fontsize=15)
    plt.ylabel('Pump Speed [RPM]', fontsize=15)
    plt.title('Conversion Factor Estimation', fontsize=17)
    without_intercept_string = f'Slope = {rpm2flow_slope_ForcedThroughZero:.7f}\n$R^2 = {rpm2flow_R2_ForcedThroughZero:.7f}$'
    with_intercept_string = f"Slope = {rpm2flow_slope:.7f}\nIntercept = {rpm2flow_intercept:.7f}\n$R^2 = {rpm2flow_R2:.7f}$"
    plt.text(0.01 * max(flow_fit_df["Slope"]), 0.85 * max(rpm_calib_points), without_intercept_string, bbox=dict(facecolor='white', edgecolor='red', boxstyle='round,pad=0.5'), fontsize=14, color='black')  # Adjust position as needed
    plt.text(0.01 * max(flow_fit_df["Slope"]), 0.6 * max(rpm_calib_points), with_intercept_string, bbox=dict(facecolor='white', edgecolor='cyan', boxstyle='round,pad=0.5'), fontsize=14, color='black')  # Adjust position as needed
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=15, loc='lower right', facecolor='white')
    plt.savefig(f"Pump_AutoCalibration_Plot_{datetime.now().date()}.png", dpi=300, transparent=False, bbox_inches='tight')
    plt.show()






























