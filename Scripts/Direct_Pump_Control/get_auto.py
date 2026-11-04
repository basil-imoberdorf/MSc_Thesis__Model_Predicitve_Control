import gpiod

# $ sudo gpiodetect -> list all gpiochips present on the system, their names, labels and number of GPIO lines

# Define the GPIO line
chip = gpiod.Chip('gpiochip4')  # Use the appropriate gpiochip number
line_num = 27  # GPIO pin number
line = chip.get_line(line_num)

# Request the line and turn it on
line.request(consumer='blink', type=gpiod.LINE_REQ_DIR_OUT)

# Set the line to 'high' (turn on)
line.set_value(1)
print(f"\nGPIO {line_num} to Pin 7 is now ON. \nPump now in AUTO mode!\n")

# Release the line and close the chip
line.release()
chip.close()
print(f"GPIO {line_num} control released.")

