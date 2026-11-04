import os
import time
import gpiod
import signal
import socket
import megaind
import numpy as np
import pandas as pd
from datetime import datetime


finish_flag = False
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


######################################################################################################################
def get_auto_mode():
    # $ sudo gpiodetect -> list all gpiochips present on the system, their names, labels and number of GPIO lines

    # Define the GPIO line
    chip = gpiod.Chip('gpiochip4')  # Use the appropriate gpiochip number
    line_num = 27  # GPIO pin number
    line = chip.get_line(line_num)

    # Request the line and turn it on
    line.request(consumer='blink', type=gpiod.LINE_REQ_DIR_OUT)

    # Set the line to 'high' (turn on)
    line.set_value(1)
    print(f"\n GPIO {line_num} to Pin 7 is now ON. \nPump now in AUTO mode!\n")

    # Release the line and close the chip
    line.release()
    chip.close()


######################################################################################################################
def get_manual_mode():
    # $ sudo gpiodetect -> list all gpiochips present on the system, their names, labels and number of GPIO lines

    # Define the GPIO line
    chip = gpiod.Chip('gpiochip4')  # Use the appropriate gpiochip number
    line_num = 27  # GPIO pin number
    line = chip.get_line(line_num)

    # Request the line and turn it on
    line.request(consumer='blink', type=gpiod.LINE_REQ_DIR_OUT)

    # Set the line to 'high' (turn on)
    line.set_value(0)
    print(f"\nGPIO {line_num} to Pin 7 is now OFF. \nPump now in MANUAL mode!\n")

    # Release the line and close the chip
    line.release()
    chip.close()
    print(f"GPIO {line_num} control released.")

    
#######################################################################################################################
def send_command_to_pump(volts):
    try:
        volts = volts.item()  # convert np array to a simple float
        megaind.set0_10Out(0, 1, volts)  # Send the voltage to the pump via the DAC      
        print(f" Current DAC Output:\t\t {round(volts, 2)} V")
        time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
        print(f" New command send! --> Task completed at: {time_date}!\n")
    except Exception as e:
        print(f"\n Error! Could not send desired flow rate command to the pump!\n {e}")


#######################################################################################################################
def start_raspberry_server(HOST, PORT):
    """Start the Raspberry Pi server to receive data."""
    server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    server_socket.bind((HOST, PORT))
    server_socket.listen(1)
    get_auto_mode()  # turn the pump into auto mode to avoid stopping it
    print(f"\n ---\n Raspberry Pi server listening on {HOST}:{PORT}\n ---")

    while True:
        try:
            print("\n------------------------------------------------------------------")
            print(" Waiting for new connection...\n")
            client_socket, client_address = server_socket.accept()
            print(f" Connection accepted from {client_address}")
            handle_client(client_socket)
            if finish_flag:
                print("Pump calibration completed!")
                raise KeyboardInterrupt
        except KeyboardInterrupt:
            send_command_to_pump(np.array([0]))  # set voltage to 0 before cutting the connection so that if connected agian the pump doesnt start outta nowhere
            get_manual_mode()  # release the pump
            print("\n Shutting down the server gracefully.")
            server_socket.close()
            break
        except Exception as e:
            print(f" Error in server loop: {e}")
            continue


#######################################################################################################################
def handle_client(client_socket):
    """Handle incoming data from a connected client."""
    global finish_flag
    try:
        while True:
            data = client_socket.recv(1024)
            if not data:
                print(" Data reception complete!")
                break

            # Deserialize the received NumPy array, the volts command for the pump 
            received_array = np.frombuffer(data, dtype=np.float64)
            print(f" Received voltage command: {np.round(received_array, 2)}\n")
            
            # Check if the abort criterion has been sent
            if received_array == 420.69:
                finish_flag = True
                break

            # The received data is the command signal in volts to be sent to the pump
            send_command_to_pump(received_array)

            # Send acknowledgment back to the client
            client_socket.send(b" Data received successfully!")
    except Exception as e:
        print(f" Error handling client: {e}")
    finally:
        client_socket.close()
        print(" Client connection closed.")
        print("------------------------------------------------------------------\n")


#######################################################################################################################
if __name__ == "__main__":
    # Register the signal handler to protect the watchdog from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)
    
    # Set the finishing criterion
    finish_flag = False

    # Start the Raspberry Pi server
    host = "192.168.0.200"
    port = 42069
    start_raspberry_server(HOST=host, PORT=port)
