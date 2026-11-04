import time
import signal
import numpy as np
from datetime import datetime
from datetime import timedelta


def main():
    # Register the signal handler which will protect the script from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Prompt the user for the model parameters
    mu_set = ask_float("\n Enter the desired growth rate, mu_set [h-1]: ")
    Yxs = ask_float("Enter the substrate yield coefficient, Yxs [gg-1]: ")
    ms = ask_float("Enter the maintenance coefficient, ms [gg-1h-1]: ")
    x0 = ask_float("Enter the initial biomass concentration, x0 [gL-1]: ")
    V0 = ask_float("Enter the initial working volume, V0 [L]: ")

    feed_density = ask_float("Enter the density of the feed solution, roh [gL-1]: ")
    glucose_conc = ask_float("Enter the glucose concentration in the feed, c_gluc [gL-1]: ")

    # Construct the flow rate reference trajectory
    relative_time = np.linspace(0, 24, 8641)  # design profile for max 24h with a time step of 10 s
    feed_mass_flow_rate = (feed_density/glucose_conc) * (mu_set / Yxs + ms) * x0 * V0 * np.exp(mu_set*relative_time)

    # Prompt the user to specify the start time of the feed profile. Must be of the form
    start_time = get_start_time()
    print('\n Feed profile successfully initialized!')
    print(' Waiting for start time...\n')
    wait_until_target_time(start_time)
    print('\n Start time reached! Feed values are being written...\n')
    start_time = datetime.now()

    # Enter the loop and write
    for t, f in zip(relative_time, feed_mass_flow_rate):
        target_time = start_time + timedelta(hours=t)
        while True:
            now = datetime.now()
            if now >= target_time:
                write_function_value(f)
                break
            time.sleep(1)  # Check every second

    print('Last value written! Feed is know held constant!')


def ask_float(instruction_message):
    """ Prompt user for a number and convert it to a float. If no integer is returned prompt again until a number is received """
    while True:
        try:
            user_float = float(input(f" {instruction_message} "))
            return user_float
        except ValueError:
            print(f" Invalid input! Please enter a number > 0\n")


def get_start_time():
    """Prompt user for a target time in HH:MM format."""
    while True:
        try:
            # The start time must be a time of the current day that has not yet passed
            user_time_str = input("\n Enter the feed start time [HH:MM]: ")
            user_time = datetime.strptime(user_time_str, "%H:%M").time()

            # Make sure the desired time has not already elapsed today
            if user_time < datetime.now().time():
                raise TypeError
            return user_time

        except ValueError:
            print(" Invalid time format. Please enter in HH:MM format.")
        except TypeError:
            now_string = datetime.now().strftime("%H:%M")
            print(f" The desired start time has already elapsed!\n Please enter a time between {now_string} and 24:00 ")


def wait_until_target_time(target_time):
    """Wait until the current time matches the target time."""
    while True:
        now = datetime.now().time()
        if now >= target_time:
            return  # leave loop as soon as the current time is "bigger" than the target time
        time.sleep(5)  # Check every 5 seconds


def write_function_value(function_value):
    """Write a single function value to the file, opening and closing the file each time."""
    with open("FlowRateCommand.txt", "w") as file:
        file.write(f"{function_value}")
    print(' Current mass flow rate command:', round(function_value, 3), ' g/h')


def signal_handler(sig, frame):
    """Handle the signal and ask for user confirmation."""
    user_input = input("\nAre you sure you want to exit? (y/n): ").strip().lower()
    if user_input == 'y':
        print(" Exiting...")
        # Raise a KeyboardInterrupt to be caught by the except block
        raise KeyboardInterrupt
    else:
        print(" Continuing...\n")
        # Re-register the signal handler to continue
        signal.signal(signal.SIGINT, signal_handler)


if __name__ == "__main__":
    main()
