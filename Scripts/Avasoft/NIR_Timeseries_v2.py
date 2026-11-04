# Standard Python libraries
import sys
import csv
import time
import signal
from datetime import datetime

# Third-party libraries
from avaspec import *
from apscheduler.schedulers.background import BackgroundScheduler
from apscheduler.triggers.interval import IntervalTrigger


def main():
    # Initialize the connection to the NIR spectrometer
    AVS_Init(0)  # Opens the communication with the spectrometer via USB. No try statement, because if this fails...
    mylist = AVS_GetList(1)  # Returns device information
    print('Serial Number:', mylist[0].SerialNumber)
    print('Spectrometer Status:', mylist[0].Status, '\n')

    # handle = AVS_Activate(mylist[0])  # get handle to manipulate the device
    handle = AVS_GetHandleFromSerial(mylist[0].SerialNumber)  # alternative

    # Switch from default 14-bit ADC to the 16-bit ADC which has higher resolution
    AVS_UseHighResAdc(handle, True)

    # Prompt the user to name the file where the spectra will be stored to
    Filename = input(" Name the output CSV file for the spectra [text without .csv]: \n  --> ")
    Filename = '../../Sensor_Data/NIR/' + Filename + '.csv'

    # Prompt the user to for the spectrometer configuration parameters
    IntegrationTime = ask_integer(instruction_message="\n Enter the integration time [ms] :")
    Xenon_FlashesPerScan = ask_integer(instruction_message="Enter the number of Xenon flashes per scan [int] :")
    NrAverages = ask_integer(instruction_message="Enter the number of scans to be averaged [int] :")

    # The first entry into the CSV is always the wavelengths so write that row already now
    wavelengths = (AVS_GetLambda(handle)[0:2048])  # get wavelengths, this one has only 2048 pixel but 4096 are possible
    save_spectrum(Filename, wavelengths)  # store vector to CSV

    # Initialize the scheduler which will call the measure function at the defined sampling interval
    scheduler = BackgroundScheduler()

    # Get desired start time and sampling interval from the user
    start_time = ask_start_datetime()
    interval_minutes = int(input("Enter the measurement time interval in minutes [integer]: "))

    # Add the job to the scheduler
    scheduler.add_job(measure,  # name of the function that will be executed in the specified interval
                      IntervalTrigger(minutes=interval_minutes, start_date=start_time),  # Create an interval trigger
                      id='interval_job',
                      args=[handle, IntegrationTime, NrAverages, Xenon_FlashesPerScan, Filename])  # args for measure()
    scheduler.start()   # Start the scheduler
    print('\n\n*** Time series successfully started! ***\n')

    # Register the signal handler to prevent the script from accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Keep the script running
    try:
        while True:
            time.sleep(1)
    except (KeyboardInterrupt, SystemExit):
        kill = AVS_Deactivate(handle)  # AVS_Deactivate returns True when device successfully closed
        print('\nConnection to spectrometer terminated: ', kill)
        AVS_Done()  # Closes the communication and releases internal storage.
        scheduler.shutdown()  # Shut down the scheduler gracefully
        print('Job scheduler deactivated!')


def ask_start_datetime():
    """
    Prompt user for the time (of today) at which the time series should start
    Input string must be of the form "HH:MM" where HH stands for hours and MM for minutes
    """
    while True:
        try:
            user_time_str = input("\n Enter the start time [HH:MM]: ")
            user_time = datetime.strptime(user_time_str, "%H:%M").time()
            user_datetime = datetime.combine(datetime.now().date(), user_time)
            return user_datetime
        except ValueError:
            print(" Invalid time format! Please enter in HH:MM format\n")


def ask_integer(instruction_message):
    """ Prompt user for an integer. If no integer is returned prompt again until an integer is received """
    while True:
        try:
            user_integer = int(input(f" {instruction_message} "))
            return user_integer
        except ValueError:
            print(f" Invalid input! Please enter an integer > 0\n")


def measure(handle, IntegrationTime, NrAverages, Xenon_FlashesPerScan, Filename):
    """ Initialize the spectrometer with the desired settings and perform the measurement. """
    print('\nNew measurement started at:\t', datetime.now())
    measure_config = MeasConfigType()  # get the object with all the measurement parameters
    measure_config.m_StartPixel = 0
    measure_config.m_StopPixel = 2047  # not 2048 because python starts counting at 0
    measure_config.m_IntegrationTime = IntegrationTime  # in milliseconds
    measure_config.m_NrAverages = NrAverages
    measure_config.m_Control_m_StrobeControl = Xenon_FlashesPerScan  # Determines the number of pulses the spectrometer sends out during one integration time cycle

    AVS_PrepareMeasure(handle, measure_config)  # Configure the spectrometer for the measurement
    AVS_Measure(handle, 0, 1)  # actual command to trigger the measurement

    spectrum_ready = False
    while not spectrum_ready:
        spectrum_ready = AVS_PollScan(handle)  # check if measurement is ready
        time.sleep(0.2)  # checking each 200ms for a new spectrum

    spectrum = AVS_GetScopeData(handle)[1][0:2048]  # [0] is timestamp in clock cycles
    # Append the new spectrum as row to the CSV file
    save_spectrum(Filename, spectrum)
    print('New spectrum appended at:\t', datetime.now())


def save_spectrum(filename, data_vector):
    """ Append data as row to the CSV file. """
    with open(filename, 'a', newline='') as csvfile:  # 'a' stands for append mode
        writer = csv.writer(csvfile)
        writer.writerow(data_vector)


def signal_handler(sig, frame):
    """ Handle the signal and ask for user confirmation. """
    user_input = input("\n Are you sure you want to exit? (y/n): ").strip().lower()
    if user_input == 'y':
        print(" Exiting...")
        # Raise a KeyboardInterrupt to be caught by the except block
        raise KeyboardInterrupt
    else:
        print(" Continuing...")
        # Re-register the signal handler to continue
        signal.signal(signal.SIGINT, signal_handler)


if __name__ == "__main__":
    main()





