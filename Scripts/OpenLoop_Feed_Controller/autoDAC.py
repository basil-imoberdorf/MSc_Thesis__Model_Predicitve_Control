import time
import signal
import megaind
from datetime import datetime
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


def main():
    # First prompt the user for the maximum RPM set-point of the pump and the factor to convert from mass flow to RPM
    while True:
        try:
            max_rpm = float(input("\n Enter the max RPM: "))
            calibration_slope = float(input(" Enter the calibration slope: "))
            break
        except ValueError:
            print(" Error: Invalid Input! Please enter a numeric value.\n")

    # Initialize the watchdog
    path_to_watch = "."  # '.' is the current folder, the directory to monitor,
    file_to_watch = "FlowRateCommand.txt"  # name of the textfile that must be modified by another script
    event_handler = TextFileHandler(file_to_watch, max_rpm, calibration_slope)
    observer = Observer()
    observer.schedule(event_handler, path=path_to_watch, recursive=False)

    print(f"\n Start monitoring changes to {file_to_watch}...\n")
    observer.start()

    # Register the signal handler which will protect the script from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    try:
        while True:
            # Define how long the watchdog "sleeps" before checking again for changes in the .txt file
            time.sleep(1)
    except KeyboardInterrupt:
        observer.stop()
        observer.join()


class TextFileHandler(FileSystemEventHandler):
    def __init__(self, filename, max_rpm, calibration_slope):
        self.flow_rate = 0
        self.filename = filename
        self.max_rpm = max_rpm  # to map the RPM to the analog signal
        self.calibration_slope = calibration_slope  # to map the desired flow rate to the necessary RPM
        self.last_modified_time = 0
        self.debounce_delay = 2  # Time in seconds before the textfile can be accessed again. To prevent double trigger

    def on_modified(self, event):
        current_time = time.time()
        # Both writing and saving the file causes this function to trigger. We only need one activation
        if current_time - self.last_modified_time > self.debounce_delay:  # wait to prevent triggering twice the same
            self.last_modified_time = current_time
            print(f"\n-------------------------------------------------------------------------")  # indicate step start
            print(f" {self.filename} has been modified.")
            self.read_and_send_command()

    def read_and_send_command(self):
        try:
            with open(self.filename, 'r') as file:
                self.flow_rate = float(file.read().strip())
        except Exception as e:
            print(f" !!!Could not read flow rate string!!!\n Error: {e}\n")

        # convert flow rate to rpm
        rpm = self.calibration_slope * self.flow_rate
        time_date = datetime.now().strftime("%H:%M:%S on %d.%m.%y")
        if rpm < 0:  # handle the case where a negative flow rate is written to prevent a crash from the DAC
            print(f" Warning! Cannot set negative flow rate: {round(self.flow_rate, 2)} g/h")
            print(f" Current RPM set to:\t 0")
            volts = 0
        
        elif rpm <= self.max_rpm:
            print(f" Current Flow Rate:\t {round(self.flow_rate, 2)} g/h")
            print(f" Current RPM set to:\t {round(rpm, 2)}")
            volts = rpm/self.max_rpm*10

        else:
            print(f" Current Flow Rate:\t {round(self.flow_rate, 2)} g/h")
            print(f" Warning! Desired flow rate cannot be met with current max RPM!")
            print(f" Current RPM set to:\t {round(self.max_rpm, 2)}")
            volts = 10

        # Send the voltage to the pump via the DAC
        megaind.set0_10Out(0, 1, volts)
        # megaind.set4_20Out(0, 1, rpm/maxRPM*16 + 4)
        print(f" Current Voltage:\t {round(volts, 2)} V")
        print(f" New command send! --> Task completed at: {time_date}!")
        print(f"-------------------------------------------------------------------------\n")  # indicate step end


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
    main()
