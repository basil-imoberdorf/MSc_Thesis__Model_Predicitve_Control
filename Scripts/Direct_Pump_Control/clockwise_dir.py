import gpiod

# $ sudo gpiodetect -> list all gpiochips present on the system, their names, labels and number of GPIO lines


# Define the GPIO line
chip = gpiod.Chip('gpiochip4')  # Use the appropriate gpiochip number
line_num1 = 22  # GPIO pin number
line_num2 = 10  # GPIO pin number
line1 = chip.get_line(line_num1)
line2 = chip.get_line(line_num2)

# Request the line and turn it on
line1.request(consumer='blink', type=gpiod.LINE_REQ_DIR_OUT)
line2.request(consumer='blink', type=gpiod.LINE_REQ_DIR_OUT)

# Set the line to 'high' (turn on)
line2.set_value(1) # enable direction control
line1.set_value(0) # change direction
line2.set_value(0) # disable direction control
print(f"\nGPIO {line_num1} to Pin 5 is now OFF! \n Pump rotates CLOCKWISE!\n")


# Release the line and close the chip
line1.release()
line2.release()
chip.close()
print(f"GPIO {line_num1} control released.")
