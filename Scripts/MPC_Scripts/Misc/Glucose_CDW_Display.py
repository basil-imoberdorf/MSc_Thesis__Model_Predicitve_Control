import time
import glob
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler


#######################################################################################################################
# Update the plot with new data
def update_plot(x_data, glucose_data, CDW_data):
    try:
        plt.plot(x_data, 0*CDW_data, color='black', linewidth=2)  # indicate where zero is
        plt.plot(x_data, glucose_data, label="Glucose via PLS", color='orange', linewidth=3)
        plt.plot(x_data, CDW_data, label="Biomass via Base", color='blue', linewidth=3)
        ax.set_title(f"Current Glucose Prediction: {round(glucose_data.iloc[-1], 2)} g/L", fontsize=18)
        # plots will be updated in the while loop of the watchdog
    except Exception as exe:
        print(f"Error updating axes: {exe}")


#######################################################################################################################
# This function is triggered each time a file is modified inside the monitored directory
def on_modified(event):
    if event.is_directory:  # Ignore directories
        return  # the function will not trigger when a directory is modified

    try:  # looking for the filename inside here allows to start this script before the necessary file is even created
        Filename = glob.glob(f"{'Measurements_and_Predictions'}*")[0]  # find the file that starts with 'Glucose_Predi'
        # Read the updated CSV file
        glucose_df = pd.read_csv(Filename)  # pandas handles the file opening and closing internally, no need for with
        x_data = pd.to_datetime(glucose_df['TimeStamp'], format='ISO8601')
        gluc_data = glucose_df['Glucose [g/L]']
        cdw_data = glucose_df['CDW [g/L]']

        # Update the plot with new data
        update_plot(x_data, gluc_data, cdw_data)

        current_datetime = datetime.now().strftime("%H:%M on %d.%m.%Y")
        print(f"Plot updated at {current_datetime}")

    except Exception as exe:
        print(f"Error reading file: {exe}")


#######################################################################################################################
# Set up the watchdog
if __name__ == "__main__":

    # Set the path to the folder to be monitored and where the Raman spectra will appear
    #path = "C:/Users/basil/Documents/Studium/Master_Applied_Biosciences/X03 - Master's Thesis/Scripts and Code/Python_Control_Algorithm"
    path = "C:/Users/ReactRaman/Documents/iC Raman Experiments/Basil/MSc 2024/Python"

    # Initialize and configure the figure
    fig, ax = plt.subplots()
    # ax.set_xlabel('Datetime', fontsize=14)
    ax.set_ylabel('Concentration [g/L]', fontsize=16)
    ax.grid(True)
    # ax.set_ylim([-5, 26])  # Set y-axis limits
    plt.subplots_adjust(left=0.1, right=0.9, top=0.9, bottom=0.25)
    plt.gca().xaxis.set_major_formatter(plt.matplotlib.dates.DateFormatter('%d.%m %H:%M'))
    plt.xticks(rotation=50)

    # Load the currently available data and plot it, Glucose_PLS_Manager.py creates and updates the target CSV
    try:  # Try to load the file called Glucose_Prediction_dd_mm_yyyy...
        filename = glob.glob(f"{'Measurements_and_Predictions'}*")[0]  # find the file that starts with 'Glucose_Predi'
        with open(filename, mode='r') as file:
            Glucose_df = pd.read_csv(file)
            time_stamps = pd.to_datetime(Glucose_df['TimeStamp'], format='ISO8601')
            glucose_values = Glucose_df['Glucose [g/L]']
            CDW_values = Glucose_df['CDW [g/L]']

        # Update the plot with initial data
        update_plot(time_stamps, glucose_values, CDW_values)
        plt.legend(loc="upper center",  # Place the legend at the top center of the bbox_to_anchor area
                   bbox_to_anchor=(0.5, -0.2),  # Move legend to the bottom center outside the plot
                   ncol=2,  # Arrange legend entries in 2 columns
                   fontsize=12)  # Adjust font size for better readability
        plt.tight_layout()
    except Exception as e:
        print(f"Error: Could not update axes {e}")

    # Create an event handler for the watchdog
    event_handler = FileSystemEventHandler()
    event_handler.on_modified = on_modified  # callback to the function that is triggered each time a file is modified
    observer = Observer()
    observer.schedule(event_handler, path, recursive=False)  # No suborders have to be monitored so: recursive=False
    observer.start()
    # Print the time when the watchdog was activated
    start_time_date = datetime.now()
    print('\n Glucose Display Watchdog activated:', start_time_date.strftime("%H:%M:%S  %d-%m-%Y"), '\n')

    try:
        while True:
            plt.draw()
            plt.pause(0.1)  # Allow for interactive updates
            time.sleep(1)  # Define how long the watchdog "sleeps" before checking again for new files
    except KeyboardInterrupt:  # Press Ctrl + C in the command window to stop monitoring the folder
        print('\n Watchdog goes to sleep..zZz..')
        observer.stop()  # Stop monitoring the folder
        observer.join()  # Wait for the observer's thread to finish its execution before proceeding
        print(' Watchdog has fallen asleep! Target folder is no longer monitored!')

