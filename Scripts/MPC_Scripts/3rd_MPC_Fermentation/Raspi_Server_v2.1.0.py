import os
import time
import signal
import socket
import numpy as np
import pandas as pd
from datetime import datetime

start_timestamp = pd.Timestamp.now()  # keep track of the previous timestamp to calculate the time difference


#######################################################################################################################
def start_raspberry_server(HOST, PORT):
    """Start the Raspberry Pi server to receive data."""
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.bind((HOST, PORT))
    server_socket.listen(1)
    print(f"\n ---\n Raspberry Pi server listening on {HOST}:{PORT}\n ---")

    while True:
        try:
            print("\n------------------------------------------------------------------")
            print(" Waiting for new connection...\n")
            client_socket, client_address = server_socket.accept()
            print(f" Connection accepted from {client_address}")

            handle_client(client_socket)

        except KeyboardInterrupt:
            print("\n Shutting down the server gracefully.")
            server_socket.close()
            break
        except Exception as e:
            print(f" Error in server loop: {e}")
            continue


#######################################################################################################################
def handle_client(client_socket):
    """Handle incoming data from a connected client."""
    try:
        while True:
            data = client_socket.recv(1024)  # 1024 bytes should be more than enough to receive all values in one go
            if not data:
                print(" Data reception complete!")
                break

            # Deserialize the received NumPy array
            received_array = np.frombuffer(data, dtype=np.float64)
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" Received measurements: {np.round(received_array, 2)}\n")

            # Save the received data to CSV
            save_to_csv(received_array)

            # Send acknowledgment back to the client
            client_socket.send(b" Data received successfully!")
    except Exception as e:
        print(f" Error handling client: {e}")
    finally:
        client_socket.close()
        print(" Client connection closed.")
        print("------------------------------------------------------------------\n")


#######################################################################################################################
def save_to_csv(numpy_data_array):
    """
    This function saves the aggregated measurements received from the BioSCADA platform and the ReactRaman to a .csv
    'filename_string' must be defined within the script as global variable or inside 'if __name__ == "__main__":' block
    :param numpy_data_array: The function needs a numpy array of shape (n,). One value for each measured quantity
    :return: None
    """

    global start_timestamp  # this variable is needed to calculate the time relative to the first measurement

    try:  # in case something happens... better use try
        # Force data to be a row array to be able to add rows to a pandas dataframe.
        column_data_array = numpy_data_array.reshape(1, -1)
        # If the file does not yet exist (or gets deleted) it will be (re)created here
        if not os.path.exists(filename_string):  # This is usually called only done on the first call of save_to_csv()
            # Create the CSV file and initialize with the column names and values for the first row
            column_names = ["CDW [g/L]", "Glucose [g/L]", "Feed Flow [L/h]", "Base Flow [L/h]", "Acid Flow [L/h]",
                            "Vol. Feed added [L]", "Vol. Base added [L]", "Vol. Acid added [L]", "DOT [%]",
                            "xCO2 [%v/v]", "Integral xCO2 [%*h]", "CDW via xCO2 [g/L]"]
            new_data_df = pd.DataFrame(column_data_array, columns=column_names)
            current_timestamp = pd.Timestamp.now()
            start_timestamp = current_timestamp  # set the initial timestamp as reference for the relative time
            new_data_df.insert(0, 'TimeStamp', current_timestamp)  # add datetime in ISO 8601 format
            new_data_df.insert(1, 'Time [h]', 0)  # add relative time in hours
            new_data_df.to_csv(filename_string, index=False)  # create the CSV
            # Signal successful creation of the CSV file
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" CSV file, {filename_string}, created and initialized!\n --> Task completed at: {time_date}\n")
        else:
            # Create new data frame to append values to the existing file
            new_data_df = pd.DataFrame(column_data_array)
            # Calculate current time relative to measurement begin
            current_timestamp = pd.Timestamp.now()
            relative_time_h = (current_timestamp - start_timestamp).total_seconds()/3600  # get time difference in hours
            new_data_df.insert(0, 'TimeStamp', current_timestamp)
            new_data_df.insert(1, 'Time [h]', relative_time_h)
            with open(filename_string, mode='a', newline='') as file:
                new_data_df.to_csv(file, index=False, header=False)
            # Signal successful data treatment
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" Current measurements successfully appended to {filename_string}\n --> Task completed at {time_date}!\n")
    except Exception as e:
        print(f"Could not store received data to .csv!\n {e} \n")


######################################################################################################################
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
if __name__ == "__main__":
    # Define the name of the file where the incoming data will be stored to
    filename_string = 'Measurement_Inputs_for_MPC_' + datetime.now().strftime("%d_%m_%Y") + '.csv'

    # Register the signal handler to protect the watchdog from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Start the Raspberry Pi server
    host = "192.168.0.200"  # Static IP address of the Raspberry Pi
    port = 42069
    start_raspberry_server(HOST=host, PORT=port)

