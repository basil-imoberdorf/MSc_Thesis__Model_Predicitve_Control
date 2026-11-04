import time
import signal
import megaind

from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from datetime import datetime


class TextFileHandler(FileSystemEventHandler):
    def __init__(self, filename, max_rpm, calibration_slope):
        self.flow_rate = 0
        self.filename = filename
        self.max_rpm = max_rpm  # to map the RPM to the analog signal
        self.calibration_slope = calibration_slope  # to map the desired flow rate to the necessary RPM
        self.last_modified_time = 0
        self.debounce_delay = 2  # Time in seconds to wait before processing the file change


    def on_modified(self, event):
        current_time = time.time()
        # Both writing and saving the file causes this function to trigger. We only need one activation
        if current_time - self.last_modified_time > self.debounce_delay:  # wait to prevent triggering twice the same
            self.last_modified_time = current_time
            print(f"\n-------------------------------------------------------------------------")  # indicate step start
            print(f" {self.filename} has been modified.")
            self.read_value()
            
            
    def read_value(self):
        try:
            with open(self.filename, 'r') as file:
                self.flow_rate = float(file.read().strip())
        except Exception as e:
            print(f" !!!Could not read flow rate string!!!\n Error: {e}\n")

        rpm = calibration_slope * self.flow_rate  # convert flow rate to rpm 
        
        if rpm < 0:  # handle the case where a negative flow rate is written to prevent a crash from the DAC
            print(f" Warning! Cannot set negative flow rate: {round(self.flow_rate,2)} g/h")
            print(f" Current RPM set to:\t\t 0")
            volts = 0
            print(f" Current Voltage set to:\t {round(volts,2)} V")
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" New command send! --> Task completed at: {time_date}!")
            print(f"-------------------------------------------------------------------------\n")  # indicate ste
        
        elif rpm <= self.max_rpm:
            print(f" Current Flow Rate:\t {round(self.flow_rate,2)} g/h")
            print(f" Current RPM:\t\t {round(rpm,2)}")
            volts = rpm/self.max_rpm*10
            print(f" Current Voltage:\t {round(volts,2)} V")
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" New command send! --> Task completed at: {time_date}!")
            print(f"-------------------------------------------------------------------------\n")  # indicate step
        else:
            print(f" Current Flow Rate:\t {round(self.flow_rate,2)} mL/h")
            print(f" Warning! Desired flow rate cannot be met with current max RPM!")
            volts = 10
            print(f" Current Voltage:\t {round(volts, 2)} V")
            time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
            print(f" New command send! --> Task completed at: {time_date}!")
            print(f"-------------------------------------------------------------------------\n")  # indicate step start

    
        # Send the voltage to the pump via the DAC
        megaind.set0_10Out(0, 1, volts)

        # amps = rpm/maxRPM*16 + 4
        # print(f"Current Voltage:\t {round(amps, 2)}\n")
        # megaind.set4_20Out(0, 1, rpm2amps(RPM))


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


if __name__ == "__main__":
    # First prompt the user for the maximum RPM set on the Pump and the factor to convert from Flow to RPM
    while True:
        try:
            max_rpm = float(input("\n\nEnter the max RPM: "))
            calibration_slope = float(input("Enter the calibration slope: "))
            break
        except ValueError:
            print("Error: Invalid Input! Please enter numeric values.\n")

    path_to_watch = "."  # The directory to watch, '.' is the current folder
    file_to_watch = "FlowRateCommand.txt"  # Replace with your actual file name
    event_handler = TextFileHandler(file_to_watch, max_rpm, calibration_slope)
    observer = Observer()
    observer.schedule(event_handler, path=path_to_watch, recursive=False)

    print(f"\nStart monitoring changes to {file_to_watch}...\n")
    observer.start()

    # Register the signal handler
    signal.signal(signal.SIGINT, signal_handler)

    try:
        while True:
            time.sleep(1)  # Define how long the watchdog "sleeps" before checking again for changes in the .txt file
    except KeyboardInterrupt:
        observer.stop()
        observer.join()
