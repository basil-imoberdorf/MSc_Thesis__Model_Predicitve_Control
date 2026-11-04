# Standard Libraries
import os
import time
import signal
import socket
import subprocess
from datetime import datetime

# Libraries that have to be installed
import numpy as np
import pandas as pd
import PIconnect as PI
from PIconnect.PIConsts import AuthenticationMode
from dotenv import load_dotenv

# PI-server credentials are read from a .env file (not part of the repository, see .env.example).
# load_dotenv() searches the folder of this script and its parent folders for the .env file.
load_dotenv()


class PIServerHandler:
    def __init__(self, server, username, password, feed_density):
        self.server = server
        self.username = username
        self.password = password

        # Initialize helper variables for the calculation of flow rates
        self.previous_time = pd.Timestamp.now()  # [h] variable to divide mass difference by the correct time difference
        self.previous_feed_bottle_mass = 0
        self.feed_density = feed_density

        # Define the reference value for the scale data. New values will be subtracted from these to get the added mass
        self.initial_feed_bottle_mass = 0  # [g_Feed_Bottle]
        self.previous_feed_mass = 0  # [g_Feed_added_up_to_previous_time]
        self.is_reference_set = False  # The reference values for the bottles are set with the first incoming values
        self.numpy_data_array = np.zeros(5)  # output array that will be stored to the CSV
        self.mpc_started = False  # bool used to signal the start of the mpc phase only once
        self.glucose_low = False  # bool used to signal glucose level below 5 g/L only once

    def process_data(self, current_time, current_feed_bottle_mass, current_DOT, received_data):
        # Check if this is the first time data arrives. If yes, set the data as reference for future data
        if not self.is_reference_set:
            self.initial_feed_bottle_mass = current_feed_bottle_mass  # [g_Feed_Bottle]
            self.previous_time = current_time
            self.is_reference_set = True

        # First, calculate how much feed has been added to the bioreactor up to the current time
        current_feed_mass = self.initial_feed_bottle_mass - current_feed_bottle_mass  # [g_Feed_added_up_to_now]

        # Calculate mass flow rates based on the mass and the time difference of two consecutive measurements
        delta_time = (current_time - self.previous_time).total_seconds() / 3600  # time difference in hours
        if delta_time > 0:  # same effect as: if delta_time > 0 self.is_reference_set
            intermediate_feed_flow = (current_feed_mass - self.previous_feed_mass) / delta_time / self.feed_density  # [L_Feed/h]
        else:  # if the reference is not set we would divide by zero therefore if it is not set we return 0 for the flow
            intermediate_feed_flow = 0

        # Update helper variables for upcoming new data
        self.previous_feed_mass = current_feed_mass
        self.previous_time = current_time

        # Arrange data in a numpy array to store it properly. Needs this form: [CDW, Glucose, Flow_F, Flow_B, Flow_A, Vol_Feed_added, Vol_Base_added, Vol_Acid_added, DOT]
        self.numpy_data_array = np.hstack([received_data[:2],  # [CDW, Glucose]
                                           intermediate_feed_flow,
                                           received_data[2:-2],  # [Flow_B, Flow_A]
                                           current_feed_mass / self.feed_density,
                                           received_data[4:],  # [Vol_Base_added, Vol_Acid_added]
                                           current_DOT])

        # Print the current measurements to the console for monitoring.
        print(f" Current Feed Scale Value:\t \033[93m{round(current_feed_bottle_mass, 2)} g\033[0m")
        print(f" Current Feed Flow:\t\t \033[91m{round(intermediate_feed_flow*1000, 2)} mL/h\033[0m")
        print(f" Current DOT:\t\t\t \033[94m{round(current_DOT, 2)} %\033[0m \n")

    def get_PI_server_value(self, tag_string):  # for this the user must be connected to the secure net of the HES-SO
        max_retries = 100  # I do not want to try indefinitely, because the rest of the algorithm must proceed
        retries = 0
        while retries < max_retries:  # this loop should make it more robust. If fetching the value fails, retry 100x
            try:
                with PI.PIServer(server=self.server,
                                 username=self.username,
                                 password=self.password,
                                 authentication_mode=AuthenticationMode.WINDOWS_AUTHENTICATION) as server:
                    tag = server.search(f"*{tag_string}*")[0]
                    # data = tag.recorded_value('-0m', retrieval_mode=RetrievalMode.AT_OR_BEFORE)  # get latest measurement
                    # Extract the (interpolated) measurements from th current time and from two minutes ago, 2min = mpc delta_t
                    data = tag.interpolated_values('-2m', '-0m', '2m')  # get two latest measurements, cant get a single one...
                    measure_time = data.index[-1]  # extract the timestamp of the latest measurement
                    measure_value = data.iloc[-1]  # extract the latest measurement
                return measure_time, measure_value

            except Exception as e:
                print(f" Could not retrieve data from PI server!\n Error: {e}\n Try again...\n")
                retries += 1
                time.sleep(0.01)  # Short delay between retries

    def start_socket_server(self, HOST, PORT):
        """Start the server to receive data."""
        server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        server_socket.bind((HOST, PORT))
        server_socket.setblocking(False)  # EFFECT NOT CLEAR WHEN TESTED. set socket to non-blocking mode to guarantee that recv() will never block indefinitely
        server_socket.listen(1)
        server_socket.settimeout(100)  # to be able to close the file using ctrl+c
        print(f"\n ----\n ------")
        print(f" Socket Server listening on {HOST} with Port: {PORT}")
        print(f" ------\n ----\n")
        while True:  # Listen indefinitely until the script is terminated
            try:
                print("\n-------------------------------------------------------------------------")
                print(" Waiting for new connection...\n")
                client_socket, client_address = server_socket.accept()
                print(f" Connection accepted from {client_address}")
                self.handle_client(client_socket)

            except KeyboardInterrupt:
                print("\n Shutting down the server gracefully.")
                server_socket.close()
                break
            except Exception as e:
                print(f" Error in server loop: {e}")
                continue

    def handle_client(self, client_socket):
        """Handle incoming data from a connected client."""
        try:
            while True:
                data = client_socket.recv(1024)
                if not data:
                    print(" Data transmission complete!")
                    break

                # Deserialize the received NumPy array
                received_array = np.frombuffer(data, dtype=np.float64)
                time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
                print(f" Received measurements: {np.round(received_array, 2)}\n")

                # Send an acknowledgment message back to the client
                client_socket.send(b" Data received successfully!")

                # Retrieve current feed data from the PI server
                print(f" Fetching data form PI server: {self.server}")
                feed_time, feed_bottle_value = self.get_PI_server_value("BR07.read:A1_Scale_PV")  # get feed bottle weight
                DOT_time, DOT_value = self.get_PI_server_value("BR07.read:A1_QIC680_PV")  # get dissolved oxygen tension

                if feed_time is None or feed_bottle_value is None:
                    print(" Error retrieving feed bottle data from PI server.\n")
                    continue
                if DOT_time is None or DOT_value is None:
                    print(" Error retrieving DOT value from PI server.\n")
                    continue

                # Process the PI server data and combine it with the received data to be able to send it off
                self.process_data(feed_time, feed_bottle_value, DOT_value, received_array)

                # # Send the data to the Raspberry Pi also via sockets:
                send_data(self.numpy_data_array)

                # If glucose concentration is below 5 g/L play a warning sound
                if received_array[1] < 5 and not self.glucose_low:
                    self.glucose_low = True  # this bool variable prevents the sound form being played over and over.
                    try:
                        subprocess.Popen(['python', 'play_glucose_low.py'])
                    except Exception as e:
                        print(f" Could not play Glucose Warning Beep anthem :(\n {e}\n")

                # if glucose concentration below 1 g/L the MPC will be started. Signal this by playing the epic anthem!
                if received_array[1] < 1 and not self.mpc_started:
                    self.mpc_started = True  # this bool variable prevents the sound form being played over and over.
                    try:
                        subprocess.Popen(['python', 'play_mpc_anthem.py'])
                    except Exception as e:
                        print(f" Could not play MPC anthem :(\n {e}\n")

        except Exception as e:
            print(f" Error handling client: {e}")
        finally:
            client_socket.close()
            print(" Client connection closed.")
            print("-------------------------------------------------------------------------\n")


#######################################################################################################################
def send_data(data_to_send):
    """
    Sends a numpy array via TCP-IP to a client server
    :param data_to_send: ndarray
    :return: None
    """
    # Server address configuration, where to go with the data
    raspi_host = '192.168.0.200'  # <-- target IP address. For Raspberry Pi works also: "raspberrypi"
    raspi_port = 42069

    # Try to create a socket object and send the glucose value to the Raspberry Pi
    try:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as client_socket:
            # Convert numpy array to bytes
            array_data = data_to_send.tobytes()
            # Connect to the server and send data
            client_socket.connect((raspi_host, raspi_port))
            client_socket.sendall(array_data)
            time_date = datetime.now()
            print(" Data sent to Raspberry Pi at ", time_date.strftime("%H:%M:%S on %d.%m.%Y"))
            response = client_socket.recv(1024).decode()
            print(f" Raspberry Pi response: {response}\n")

    except Exception as e:
        print(f" Communication Error: Could not send the data!:\n {e} \n")


#######################################################################################################################
# Safety measure to prevent accidental termination of the script (does not prevent termination if cmd window is closed!)
def signal_handler(sig, frame):
    """Handle the signal and ask for user confirmation."""
    user_input = input("\nAre you sure you want to exit? (y/n): ").strip().lower()
    if user_input == 'y':
        print("Exiting...")
        # Raise a KeyboardInterrupt to be caught by the except block
        raise KeyboardInterrupt
    else:
        print("Continuing...\n")
        # Re-register the signal handler to continue
        signal.signal(signal.SIGINT, signal_handler)


#######################################################################################################################
def get_input_with_default(prompt, default):
    """Prompt the user for numeric input with a default value."""
    while True:
        user_input = input(f"{prompt} \033[93m\033[1m{default} g/L\033[0m or: ").strip()
        # Check if the input is empty (default value)
        if not user_input:
            print(f" --> Default value used: \033[92m\033[1m{default} g/L\033[0m")
            return default
        # Validate the input as a number
        try:
            numeric_value = float(user_input)  # Change to int(user_input) if only integers are allowed
            print(f" --> New value set: \033[92m\033[1m{numeric_value} g/L\033[0m")
            return numeric_value
        except ValueError:
            print("\033[91mInvalid input! Please enter a number.\033[0m")


#######################################################################################################################
if __name__ == "__main__":
    # enables on Windows to print colored values to the command window
    os.system('color')  # used to accentuate the glucose prediction

    # get the feed density from the user to convert the mass flow to a volume flow
    default_density = 1185  # g_Feed/L_Feed of a 500 g_glucose/L_feed solution
    volumetric_feed_density = get_input_with_default("\n Feed Density set to:", default_density)

    # Create an instance of the handler with the necessary credentials
    pi_handler = PIServerHandler(server=os.getenv('PI_SERVER', 'VMENPPIAFSERVER'),
                                 username=os.getenv('PI_USERNAME'),
                                 password=os.getenv('PI_PASSWORD'),
                                 feed_density=volumetric_feed_density)

    # Register the signal handler to protect the watchdog from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Start the socket server to accept connections from the Raman laptop
    host = "192.168.0.150"  # IP address of the internet laptop
    port = 777

    pi_handler.start_socket_server(HOST=host, PORT=port)







