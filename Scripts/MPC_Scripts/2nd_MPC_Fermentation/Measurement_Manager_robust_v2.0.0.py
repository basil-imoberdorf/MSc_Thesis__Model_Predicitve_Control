# Standard Libraries
import os
import re
import time
import signal
import socket
from datetime import datetime

# Libraries that have to be installed
import serial  # to connect to the base and acid scale via the USB/RS-232 adapter
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


#######################################################################################################################
# Class to store conveniently the previous mass of the base and acid scale without the need for global variables
class MeasurementsProcessor(FileSystemEventHandler):
    def __init__(self, filename, base_density, acid_density, slope, intercept):
        super().__init__()  # Initialize the parent constructor, the watchdog
        # Initialize helper variables for the calculation of flow rates
        self.base_density = base_density  # [g/L] volumetric density used to convert mass flow to volume flow
        self.acid_density = acid_density  # [g/L] volumetric density used to convert mass flow to volume flow
        self.previous_time = pd.Timestamp.now()  # [h] variable to divide mass difference by the correct time difference
        # Define the reference values for the scale data. New values will be subtracted from these to get the added mass
        self.initial_base_bottle_mass = 0  # [g_Base_Bottle]
        self.initial_acid_bottle_mass = 0  # [g_Acid_Bottle]
        self.is_reference_set = False  # The reference values for the bottles are set with the first incoming values

        self.current_base_bottle_mass = 0  # [g_Base_Bottle]
        self.current_acid_bottle_mass = 0  # [g_Acid_Bottle]

        self.previous_base_mass = 0  # [g_Base_added_up_to_previous_time]
        self.previous_acid_mass = 0  # [g_Acid_added_up_to_previous_time]

        self.filename = filename
        self.numpy_data_array = np.zeros(6)  # array that will be stored to the CSV

        # Initialize the model parameters for CDW estimation based on a linear function
        self.cdw_slope = slope
        self.cdw_intercept = intercept
        self.current_glucose_prediction = 0

    def convert_measurement_data(self, current_time):
        # Check if this is the first time data arrives. If yes, set the data as reference for future data
        if not self.is_reference_set:
            self.initial_base_bottle_mass = self.current_base_bottle_mass  # [g_Base_Bottle]
            self.initial_acid_bottle_mass = self.current_acid_bottle_mass  # [g_Acid_Bottle]$
            self.is_reference_set = True

        # First, calculate how much feed has been added to the bioreactor up to the current time
        current_base_mass = self.initial_base_bottle_mass - self.current_base_bottle_mass  # [g_Base_added_up_to_now]
        current_acid_mass = self.initial_acid_bottle_mass - self.current_acid_bottle_mass  # [g_Acid_added_up_to_now]

        # Calculate mass flow rates based on the mass and the time difference of two consecutive measurements
        delta_time = (current_time - self.previous_time).total_seconds() / 3600  # time difference in hours
        intermediate_base_flow = (current_base_mass - self.previous_base_mass) / delta_time / self.base_density  # [L_Base/h]
        intermediate_acid_flow = (current_acid_mass - self.previous_acid_mass) / delta_time / self.acid_density  # [L_Acid/h]

        # Update helper variables for upcoming new data
        self.previous_base_mass = current_base_mass
        self.previous_acid_mass = current_acid_mass
        self.previous_time = current_time

        # Calculate CDW based on added amount of base using a linear calibration
        cdw_prediction = self.cdw_slope * current_base_mass + self.cdw_intercept

        # Arrange data in a numpy array to store it properly. Needs to have this form: [CDW, Glucose, Flow_Base, Flow_Acid, Vol_Base_added, Vol_Acid_added]
        self.numpy_data_array = np.array([cdw_prediction,
                                          self.current_glucose_prediction,
                                          intermediate_base_flow,
                                          intermediate_acid_flow,
                                          current_base_mass / self.base_density,
                                          current_acid_mass / self.acid_density])

        # Print the current measurements to the console for monitoring.
        print(f"  Current CDW Prediction:\t\033[94m\033[1m{round(cdw_prediction, 2)}\033[0m g/L")
        print(f"  Current Glucose Prediction:\t\033[33m{round(self.current_glucose_prediction, 2)}\033[0m g/L")
        print(f"  Current Base Flow:\t\t\033[95m\033[1m{round(intermediate_base_flow*1000, 2)}\033[0m mL/h")
        print(f"  Current Acid Flow:\t\t\033[92m\033[1m{round(intermediate_acid_flow*1000, 2)}\033[0m mL/h \n")
        # Also print the values obtained form the scales to compare with the reading on the actual scale
        print(f"  Current Base Scale Reading:\t\033[95m\033[1m{round(self.current_base_bottle_mass, 2)}\033[0m g")
        print(f"  Current Acid Scale Reading:\t\033[92m\033[1m{round(self.current_acid_bottle_mass, 2)}\033[0m g\n")

    # Function that will be triggered everytime a new file (or folder) is created within the monitored folder
    def on_created(self, event):
        # Get the path of the created file
        created_file_path = event.src_path
        filename = os.path.basename(created_file_path)  # Extract the filename from the path
        print("\n------------------------------------------------------------------")
        print(' New file detected:', filename)

        # Try to open and read the file, compute the glucose prediction, get scale values, compute CDW and send all
        time.sleep(0.2)  # wait a bit before trying to open the file, often the source is not yet finished writing to it
        try:
            # with open(created_file_path, mode='r') as file:
            new_spectral_df = pd.read_csv(f"{path}/{filename}")
            new_spectral = new_spectral_df.to_numpy()
            raman_spectrum = new_spectral[:, 1]  # extract the spectrum, we don't need the raman shift

            # Pre-Process the spectrum accordingly for the glucose PLS model
            processed_spectrum = pre_process_spectrum(raman_spectrum,
                                                      savgol_window=93,  # Savitzky-Golay filter width
                                                      savgol_order=2,    # Order of fitted polynomial
                                                      SNV=True,          # Standard Normal Variate
                                                      derivative=1,      # numerical differentiation
                                                      from_wave_num=1000, to_wave_num=1500)  # variable selection

            # Compute Glucose Prediction using the PLS loadings
            processed_spectrum = np.concatenate(([1], processed_spectrum))  # add the 1 for the bias term
            self.current_glucose_prediction = np.dot(processed_spectrum, glucose_loadings[:, 0])  # inference

            try:
                # Get the current process values of the base and acid scale
                self.current_base_bottle_mass = read_from_scale(port='COM4')
                self.current_acid_bottle_mass = read_from_scale(port='COM3')
            except Exception as e:
                print(f" Warning!: {e}\n")

            # Convert raw measurements to CDW. Then combine with the glucose prediction and flow rates to a nice array
            self.convert_measurement_data(pd.Timestamp.now())

            # Send the measurements and predictions via the network to ...
            send_data(self.numpy_data_array)

            # Finally append the data to a .csv file that is stored only on this computer
            save_to_csv(self.numpy_data_array)
            print("------------------------------------------------------------------\n")

        except Exception as e:
            filename = os.path.basename(created_file_path)
            file_extension = os.path.splitext(filename)[1]
            print(f" {file_extension}, file detected. Wrong file type! Could not open file!\n {e}")
            print("------------------------------------------------------------------\n")


#######################################################################################################################
def read_from_scale(port):
    """ Reads data over a serial connection and guarantees a valid response. """
    max_retries = 100  # I do not want to try indefinitely, because the rest of the algorithm must proceed
    retries = 0

    while retries < max_retries:
        try:
            # Attempt to open the serial connection
            with serial.Serial(port, 9600, timeout=2.1) as ser:  # 2.1sec because scale has 2s settling time
                ser.reset_input_buffer()
                # Read data
                line = ser.readline().decode('utf-8').strip()
                # Validate the data (check if line is non-empty)
                if line:
                    ## Regular expression to check the expected format (e.g., "+00353.41 g S")
                    #pattern = r"^[+-]?\d+\.\d+\s+[a-zA-Z]+\s*[a-zA-Z]*$"
                    #if not re.match(pattern, line):
                    #    raise ValueError(f"Invalid data format: {line}")

                    # Validate further using custom validation logic (can be modified as needed)
                    try:
                        # Extract the numeric value (e.g., "+00353.41")
                        match = re.search(r'[-+]?\d*\.\d+|\d+', line)  # Match the numeric value
                        numeric_value = np.array([float(match.group(0))], dtype=np.float64)[0]

                        if numeric_value > 200:  # sometimes the scale read like 0.5 despiet it being way more
                            # If valid, return the numeric value
                            return numeric_value  # this is where we want to leave the rabbit hole
                        else:
                            print(f"Process value below 200 g detected. Retrying...")
                    except (AttributeError, ValueError) as e:
                        print(f"Error extracting numeric value from data: {line}. Error: {str(e)}")

                # If data is invalid or empty, retry
                retries += 1
                time.sleep(0.01)  # Short delay between retries

        except (serial.SerialException, UnicodeDecodeError) as e:
            print(f"Warning: {e}. Retrying...\n")
            retries += 1
            time.sleep(0.5)  # Longer delay for connection issues
        except ValueError as ve:
            # Catch any value errors related to invalid data format or numeric extraction
            print(f"Error: {ve}. Retrying...\n")
            retries += 1
            time.sleep(0.01)  # Retry after a short delay

    # If all retries fail, raise a RuntimeError
    raise RuntimeError(f"No valid data obtained after {max_retries} retries!\n Previous value will be used!")


#######################################################################################################################
def clear_spikes(raw_Raman_spectrum):
    """
    This function removes systematic artifacts (spikes) in the Raman spectra caused by the sensor.
    The lower and upper bounds define the location of possible spikes. The algorithm loops over the spectrum
    and checks if a spike is present at each given location. If yes, the data within the bounds is removed and imputed
    with a linear interpolation of the values adjacent to the spike + some noise to keep the variance constant
    :param raw_Raman_spectrum:     type: numpy array
    :return: clean_spectrum:       type: numpy array
    """
    # Pre-allocate an array for the clean spectrum
    clean_spectrum = raw_Raman_spectrum

    # Define regions where spikes are located in terms of their Raman shift
    lower_bounds = np.array([3184, 3128, 3081, 3074, 3068, 3050, 3006, 2986, 2976, 2955, 2945, 2939, 2362, 2289, 2152,
                             1940, 1365, 1340, 1297, 1194, 1011, 768, 525, 243])
    upper_bounds = np.array([3188, 3132, 3084, 3078, 3072, 3053, 3009, 2991, 2978, 2960, 2950, 2943, 2366, 2292, 2156,
                             1944, 1373, 1345, 1302, 1203, 1016, 774, 533, 249])
    # Convert Raman-shifts to indices of the vector containing the spectrum
    lb_indices = 3200 - upper_bounds
    ub_indices = 3201 - lower_bounds

    # Check all possible locations for spikes
    for i in range(len(lower_bounds)):
        # Spikes that are too close need special treatment to properly detect them.
        # The number of neighbouring values is reduced to 3 to make it less likely to overlap with other spikes
        # However there is maybe also data imputed if there is no spike because the detection is then less reliable
        if lower_bounds[i] == 2945 or lower_bounds[i] == 2939 or lower_bounds[i] == 3074 or lower_bounds[i] == 3068:
            neighbour_range = 3
        else:
            neighbour_range = 5
        # Calculate standard deviation of values within the bounds and for those adjacent to the bounds
        peak_std = np.std(raw_Raman_spectrum[lb_indices[i]:ub_indices[i]], ddof=1)
        upstream_data = raw_Raman_spectrum[ub_indices[i]:ub_indices[i] + neighbour_range]
        downstream_data = raw_Raman_spectrum[lb_indices[i] - neighbour_range:lb_indices[i]]
        neighborhood_std = np.std(np.hstack([upstream_data, downstream_data]), ddof=1)
        # Detect if there is a peak and if yes impute the data with a linear interpolation + white noise
        threshold = 0.6
        if peak_std > (threshold * neighborhood_std):
            # Get parameters of linear function connecting the edge points of the peak
            slope = (raw_Raman_spectrum[ub_indices[i]] - raw_Raman_spectrum[lb_indices[i]]) / (
                    ub_indices[i] - lb_indices[i])
            intercept = raw_Raman_spectrum[lb_indices[i]] - lb_indices[i] * slope
            imputation_data = slope * np.arange(lb_indices[i], ub_indices[i]) + intercept
            clean_spectrum[lb_indices[i]:ub_indices[i]] = imputation_data + 0.1 * neighborhood_std * np.random.randn(
                len(imputation_data))
    return clean_spectrum


#######################################################################################################################
def pre_process_spectrum(raw_spectrum, savgol_window, savgol_order, SNV, derivative, from_wave_num, to_wave_num):
    """
    # Function to apply the pre-processing to the spectra to be compatible with the PLS loadings
    # Use the same parameters as in the training of the PLS model
    :param raw_spectrum: numpy array that holds the Raman spectrum, type: ndarray of floats
    :param savgol_window: set width of the Savitzky-Golay filter type: odd int
    :param savgol_order: set order of the fitted polynomial, type: int
    :param SNV:  enable or disable Standard Normal Variate,  type: bool
    :param derivative: set order of derivative (numerical),   type: int
    :param from_wave_num: for variable selection the start Raman Shift. Must be: 100 < start < to_wave_num, type: int
    :param to_wave_num:   for variable selection the end Raman Shift. Must be: to_wave_num < end < 3200 : int
    :return: final_spectrum the pre-processed spectrum,  type: ndarray of floats
    """
    # Clear Spikes
    inter_spectrum = clear_spikes(raw_spectrum)

    # Savitzky-Golay-Filter
    inter_spectrum = savgol_filter(inter_spectrum, savgol_window, savgol_order)

    # Standard Normal Variate (SNV)
    if SNV:  # degrees of freedom n-1 for sample standard deviation
        inter_spectrum = (inter_spectrum - np.mean(inter_spectrum)) / np.std(inter_spectrum, ddof=1)

    # Numerical differentiation, variables are exactly 1 cm-1 apart from the next so no need to divide diff(Raman_Shift)
    if derivative > 0:
        inter_spectrum = np.concatenate(([0], np.diff(inter_spectrum, derivative)))

    # Variable Selection
    inter_spectrum = inter_spectrum[3200 - to_wave_num:3201 - from_wave_num]
    final_spectrum = np.flip(inter_spectrum)  # flip array to have the wave numbers in increasing order

    return final_spectrum


#######################################################################################################################
def get_PLS_loadings(path_to_csv):
    """
    Sends a numpy array via TCP-IP to a client server
    :param path_to_csv: path to the CSV file where the loadings are stored for the PLS model
    :return: numpy array of the regression coefficients / loadings to perform inference on the glucose concentration
    """
    loadings_df = pd.read_csv(path_to_csv)
    loadings = loadings_df.to_numpy()
    return loadings


#######################################################################################################################
def send_data(data_to_send):
    """
    Sends a numpy array via TCP-IP to a client server
    :param data_to_send: ndarray
    :return: None
    """
    # Server address configuration, where to go with the data
    host = '192.168.0.150'  # <-- IP address of the internet laptop.
    port = 777

    # Try to create a socket object and send the glucose value to the Raspberry Pi
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as client_socket:
            # Convert numpy array to bytes
            array_data = data_to_send.tobytes()
            # Connect to the server and send data
            client_socket.connect((host, port))
            client_socket.sendall(array_data)
            time_date = datetime.now()
            print(" Data sent to the server.", time_date.strftime("%H:%M:%S on %d.%m.%Y"))
            response = client_socket.recv(1024).decode()
            print(f" Server response: {response}\n")

    except Exception as e:
        print(f" Communication Error: Could not send the data!:\n {e} \n")


#######################################################################################################################
def save_to_csv(numpy_data_array):
    """
    This function saves the aggregated measurements received from the BioSCADA platform and the ReactRaman to a .csv
    filename_string must be defined within the script
    :param numpy_data_array: The function needs a numpy array. One value for each measured quantity
    :return: None
    """
    try:  # in case something happens... better use try
        # Force data to be a row array to be able to add rows to a pandas df.
        column_data_array = numpy_data_array.reshape(1, -1)
        # If the file does not yet exist (or gets deleted) it will be (re)created here
        if not os.path.exists(filename_string):  # This is usually called only done on the first call of save_to_csv()
            # Create the CSV file and initialize with the column names and values for the first row
            column_names = ["CDW [g/L]", "Glucose [g/L]", "Base Flow [L/h]", "Acid Flow [L/h]", "Base Scale [g]",
                            "Acid Scale [g]"]
            new_data_df = pd.DataFrame(column_data_array, columns=column_names)
            new_data_df.insert(0, 'TimeStamp', pd.Timestamp.now())  # add datetime in ISO 8601 format
            new_data_df.to_csv(filename_string, index=False)  # create the CSV
            # Signal successful creation of the CSV file
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f"\n CSV file, {filename_string}, created and initialized!\n --> Task completed at: {time_date}\n")
        else:
            # Append values to the existing file
            new_data_df = pd.DataFrame(column_data_array)
            new_data_df.insert(0, 'TimeStamp', pd.Timestamp.now())
            with open(filename_string, mode='a', newline='') as file:
                new_data_df.to_csv(file, index=False, header=False)
            # Signal successful data treatment
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(
                f" Current measurements successfully appended to {filename_string}\n --> Task completed at {time_date}!")
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
        user_input = input(f"{prompt} \033[93m\033[1m{default}\033[0m or: ").strip()
        # Check if the input is empty (default value)
        if not user_input:
            print(f" --> Default value used: \033[92m\033[1m{default}\033[0m")
            return default
        # Validate the input as a number
        try:
            numeric_value = float(user_input)  # Change to int(user_input) if only integers are allowed
            print(f" --> New value set: \033[92m\033[1m{numeric_value}\033[0m")
            return numeric_value
        except ValueError:
            print("\033[91mInvalid input! Please enter a number.\033[0m")


#######################################################################################################################
# Set up the watchdog
if __name__ == "__main__":
    # enables on Windows to print colored values to the command window
    os.system('color')  # used to accentuate the glucose prediction

    # Get loadings/regression coefficients of the PLS model for glucose
    glucose_loadings = get_PLS_loadings(path_to_csv='PLS_Loadings_Glucose_4CalibFerm.csv')  # is in the same folder
    # Set the path to the folder to be monitored and where the Raman spectra will appear
    # path = "C:/Users/basil/Documents/Studium/Master_Applied_Biosciences/X03 - Master's Thesis/Scripts and Code/Python_Control_Algorithm/Raw_Raman_Spectra"
    path = "C:/Users/ReactRaman/Documents/iC Raman Experiments/Basil/MSc 2024/Raw_Raman_Spectra/Fermentations/6th_Fermentation_MPC"
    filename_string = 'Measurements_and_Predictions_' + datetime.now().strftime("%d_%m_%Y") + '.csv'

    # Get the density of the applied acid and base solutions
    default_base_density = 899  # [g/L] volumetric density of the 28-30% ammonia solution
    default_acid_density = 1146.2  # [g/L] volumetric density of the 25% phosphoric acid solution
    volumetric_base_density = get_input_with_default("\n Base density set to [g/L]:", default_base_density)
    volumetric_acid_density = get_input_with_default(" Acid density set to [g/L]:", default_acid_density)

    # Get the parameters to convert the mass of base added to a CDW prediction
    default_slope = 0.5328
    default_intercept = 0.8429
    calibration_slope = get_input_with_default("\n CDW calibration SLOPE set to:", default_slope)
    calibration_intercept = get_input_with_default(" CDW calibration INTERCEPT set to:", default_intercept)

    # Create the event handler for the watchdog. It is a for predictions extended class of the FileSystemEventHandler
    event_handler = MeasurementsProcessor(filename=filename_string,
                                          base_density=volumetric_base_density,
                                          acid_density=volumetric_acid_density,
                                          slope=calibration_slope,
                                          intercept=calibration_intercept)

    observer = Observer()
    observer.schedule(event_handler, path, recursive=False)  # No suborders have to be monitored so: recursive=False
    observer.start()
    # Print the time when the watchdog was activated
    start_time_date = datetime.now()
    print('\n ----\n ------')
    print(' Raman PLS Watchdog activated! Start at:', start_time_date.strftime("%H:%M:%S on %d-%m-%Y"))
    print(' ------\n ----\n')

    # Register the signal handler to protect the watchdog from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    try:
        while True:
            time.sleep(1)  # Define how long the watchdog "sleeps" before checking again for new files
    except KeyboardInterrupt:  # Press Ctrl + C in the command window to stop monitoring the folder
        print('\n Watchdog goes to sleep..zZz..')
        observer.stop()  # Stop monitoring the folder
        observer.join()  # Wait for the observer's thread to finish its execution before proceeding
        print(' Watchdog has fallen asleep! Target folder is no longer monitored!')

