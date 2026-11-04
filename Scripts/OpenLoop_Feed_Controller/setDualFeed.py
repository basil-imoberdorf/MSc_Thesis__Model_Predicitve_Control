import time
import signal
import numpy as np
from datetime import datetime
from datetime import timedelta


def main():
    # Register the signal handler which will protect the script from (one) accidental keyboard interrupt
    signal.signal(signal.SIGINT, signal_handler)

    # Prompt the user for the model parameters
    mu_set1 = ask_float("\n Enter the first desired growth rate, mu_set1 [h-1]: ")
    t_feed1 = ask_float("Enter the duration of the first growth phase, t_feed1 [h]: ")
    mu_set2 = ask_float("Enter the second desired growth rate, mu_set2 [h-1]: ")
    t_feed2 = ask_float("Enter the duration of the second growth phase, t_feed2 [h]: ")
    Yxs = ask_float("Enter the substrate yield coefficient, Yxs [gg-1]: ")
    ms = ask_float("Enter the maintenance coefficient, ms [gg-1h-1]: ")
    x0 = ask_float("Enter the initial biomass concentration, x0 [gL-1]: ")
    V0 = ask_float("Enter the initial working volume, V0 [L]: ")
    feed_density = ask_float("Enter the density of the feed solution, roh [gL-1]: ")
    glucose_conc = ask_float("Enter the glucose concentration in the feed, c_gluc [gL-1]: ")

    ################################
    # make first exp. feed
    relative_times1 = np.linspace(0, t_feed1, int(t_feed1*360+1))  # time in hours with 10s interval for first feed
    feed_mass_flow_rate_1 = np.round((feed_density/glucose_conc) * (mu_set1/Yxs + ms) * x0 * V0 * np.exp(mu_set1 * relative_times1), decimals=4)

    # calculate volume and biomass concentration after the first exp. feed based on model
    V_feed = np.trapz(feed_mass_flow_rate_1, relative_times1)/feed_density
    V1 = V_feed + V0
    x1 = (x0 * V0 * np.exp(mu_set1 * relative_times1[-1])) / V1

    # make second exp. feed
    relative_times2 = np.linspace(0, t_feed2, int(t_feed2*360+1))  # *360 to have an entry every 10s
    feed_mass_flow_rate_2 = np.round((feed_density/glucose_conc) * (mu_set2/Yxs + ms) * x1 * V1 * np.exp(mu_set2 * relative_times2), decimals=4)

    # Stack the two feeds
    relative_times = np.hstack((relative_times1, np.linspace(t_feed1+1/360, t_feed1+t_feed2, int(t_feed2*360)))) # +1/360 to shift 10s
    feed_mass_flow_rate = np.hstack((feed_mass_flow_rate_1, feed_mass_flow_rate_2[1:]))  # omit the first value of the second phase as it is the same as the last of the first phase
    ################################

    # Prompt the user to specify the start time of the feed profile. Must be of the form
    start_time = get_start_time()
    print('\n Bi-feed profile successfully initialized!')
    print('\n Waiting for starting time...')
    wait_until_target_time(start_time)
    print(' Starting time reached! Values are being written\n')
    start_time = datetime.now()

    # Enter the loop and write each 10s the feed mass flow rate command
    for t, f in zip(relative_times, feed_mass_flow_rate):
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
            return
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
        print("Exiting...")
        # Raise a KeyboardInterrupt to be caught by the except block
        raise KeyboardInterrupt
    else:
        print("Continuing...\n")
        # Re-register the signal handler to continue
        signal.signal(signal.SIGINT, signal_handler)


if __name__ == "__main__":
    main()