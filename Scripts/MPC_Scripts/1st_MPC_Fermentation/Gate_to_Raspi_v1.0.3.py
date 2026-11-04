# Standard Libraries
import os
import time
import signal
import socket
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
    def __init__(self, server, username, password):
        self.server = server
        self.username = username
        self.password = password

        # Initialize helper variables for the calculation of flow rates
        self.previous_time = pd.Timestamp.now()  # [h] variable to divide mass difference by the correct time difference
        self.previous_feed_bottle_mass = 0
        # Define the reference value for the scale data. New values will be subtracted from these to get the added mass
        self.initial_feed_bottle_mass = 0  # [g_Feed_Bottle]
        self.previous_feed_mass = 0  # [g_Feed_added_up_to_previous_time]
        self.is_reference_set = False  # The reference values for the bottles are set with the first incoming values
        self.numpy_data_array = np.zeros(5)  # output array that will be stored to the CSV

    def process_data(self, current_time, current_feed_bottle_mass, pi_server_data):
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
            intermediate_feed_flow = (current_feed_mass - self.previous_feed_mass) / delta_time / 1185  # [L_Feed/h]
        else:  # if the reference is not set we would divide by zero therefore if it is not set we return 0 for the flow
            intermediate_feed_flow = 0

        # Update helper variables for upcoming new data
        self.previous_feed_mass = current_feed_mass
        self.previous_time = current_time

        # Arrange data in a numpy array to store it properly. Needs this form: [CDW, Glucose, Flow_F, Flow_B, Flow_A]
        self.numpy_data_array = np.hstack([pi_server_data[:2], intermediate_feed_flow, pi_server_data[2:-2]])
        print(f" Current Feed Scale Value:\t \033[93m{round(current_feed_bottle_mass, 2)} g\033[0m")
        print(f" Current Feed Flow:\t\t \033[91m{round(intermediate_feed_flow, 2)*1000} mL/h\033[0m \n")

    def get_value(self, tag_string):
        max_retries = 100  # I do not want to try indefinitely, because the rest of the algorithm must proceed
        retries = 0
        while retries <= max_retries:  # this loop should make it more robust. If fetching the value fails, retry 100x
            try:
                with PI.PIServer(server=self.server,
                                 username=self.username,
                                 password=self.password,
                                 authentication_mode=AuthenticationMode.WINDOWS_AUTHENTICATION) as server:
                    tag = server.search(f"*{tag_string}*")[0]
                    # data = tag.recorded_value('-0m', retrieval_mode=RetrievalMode.AT_OR_BEFORE)  # get latest measurement
                    # Extract the (interpolated) measurements from th current time and from two minutes ago, 2min = mpc delta_t
                    data = tag.interpolated_values('-2m', '-0m', '2m')  # get two latest measurements
                    measure_time = data.index[-1]  # extract the timestamp of the latest measurement
                    measure_value = data.iloc[-1]  # extract the latest measurement
                    print(f" Feed scale read via {server.server_name}")
                return measure_time, measure_value

            except Exception as e:
                print(f" Could not retrieve data from PI server!\n Error: {e}\n Try again...\n")
                retries += 1
                time.sleep(0.01)  # Short delay between retries

    def start_socket_server(self, HOST, PORT):
        """Start the server to receive data."""
        server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        server_socket.bind((HOST, PORT))
        server_socket.listen(1)  # listen to only 1 device
        server_socket.settimeout(100)  # to be able to close the file using ctrl+c the defined timeout
        print("\n ----\n ------")
        print(f" Socket Server listening on {HOST} with Port: {PORT}")
        print(" ------\n ----\n")
        while True:
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
            while True:  # loop to make sure all data is received
                data = client_socket.recv(1024)  # 1024 bytes is more than enough buffer to receive the data in one go
                if not data:  # No more data is available when client_socket.recv() returns an empty byte string
                    print(" Data reception complete!")
                    break

                # Deserialize the received NumPy array
                received_array = np.frombuffer(data, dtype=np.float64)
                print(f" Received measurements: {np.round(received_array, 2)}\n")

                # Retrieve current feed bottle weight from the PI server
                feed_time, feed_bottle_value = self.get_value("BR07.read:A1_Scale_PV")  # PI tag for the feed scale

                if feed_time is None or feed_bottle_value is None:
                    print(" Error retrieving data from PI server.\n")
                    continue

                # Process the PI server data and combine it with the received data to send it off
                self.process_data(feed_time, feed_bottle_value, received_array)

                # Send the data to the Raspberry Pi also via sockets:
                send_data(self.numpy_data_array)

                # Send an acknowledgment message back to the client
                client_socket.send(b" Data received successfully")

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
            print(" Data sent to the server.", time_date.strftime("%H:%M:%S on %d.%m.%Y"))
            response = client_socket.recv(1024).decode()
            print(f" Server response: {response}\n")

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
if __name__ == "__main__":
    # enables on Windows to print colored values to the command window
    os.system('color')  # used to accentuate the glucose prediction

    # Create an instance of a PI server handler with the necessary credentials
    pi_handler = PIServerHandler(server=os.getenv('PI_SERVER', 'VMENPPIAFSERVER'),
                                 username=os.getenv('PI_USERNAME'),
                                 password=os.getenv('PI_PASSWORD'))

    # Register the signal handler to protect the watchdog from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Start the socket server to accept connections from the Raman laptop
    host = "192.168.0.150"  # IP address of the internet laptop
    port = 777
    pi_handler.start_socket_server(HOST=host, PORT=port)







