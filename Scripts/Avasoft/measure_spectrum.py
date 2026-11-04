# Note: This file requires avaspecx64.dll and avaspec.py to be in the same folder!

from avaspec import *
import matplotlib.pyplot as plt
import time

########################################################################################################################
# Connecting to the Spectrometer
########################################################################################################################
# First, AVS_Init() initializes the library with the a_Port parameter by looking for USB devices (0),
# Ethernet devices (256), or both (-1).
AVS_Init(0)


mylist = AVS_GetList(1)  # input is an integer of the number of connected devices

print(mylist[0].Status)
print(type(mylist[0].Status))
print(mylist[0].SerialNumber)
print(type(mylist[0].SerialNumber))
handle = AVS_Activate(mylist[0])
print(type(handle))
########################################################################################################################
# Configuring the Spectrometer
########################################################################################################################
AVS_UseHighResAdc(handle, True)  # Switch from default 14-bit ADC to the 16-bit ADC

# With the handle, we can now prepare the spectrometer to collect data. We call AVS_GetParameter() to retrieve
# configuration information about the spectrometer in the form of a DeviceConfigType. From here we retrieve the number
# of pixels, which we’ll use to get the wavelength data and later plot the spectrum. Next, we use AVS_GetLambda() to
# retrieve the wavelength calibration for each pixel – we’ll iterate through and store the values for later.

devcon = AVS_GetParameter(handle, 63484)  # "size" is the size in bytes allocated to store DeviceConfigType
pixels = devcon.m_Detector_m_NrPixels
wavelengths = AVS_GetLambda(handle)

# Measurement Parameters
measconfig = MeasConfigType()
measconfig.m_StartPixel = 0
measconfig.m_StopPixel = pixels - 1
measconfig.m_IntegrationTime = 100
measconfig.m_IntegrationDelay = 0
measconfig.m_NrAverages = 1
measconfig.m_CorDynDark_m_Enable = 0
measconfig.m_CorDynDark_m_ForgetPercentage = 100
measconfig.m_Smoothing_m_SmoothPix = 0
measconfig.m_Smoothing_m_SmoothModel = 0
measconfig.m_SaturationDetection = 0
measconfig.m_Trigger_m_Mode = 0
measconfig.m_Trigger_m_Source = 0
measconfig.m_Trigger_m_SourceType = 0
measconfig.m_Control_m_StrobeControl = 0
measconfig.m_Control_m_LaserDelay = 0
measconfig.m_Control_m_LaserWidth = 0
measconfig.m_Control_m_LaserWaveLength = 0.0
measconfig.m_Control_m_StoreToRam = 0

########################################################################################################################
# Collecting and displaying spectral data
########################################################################################################################
# Now we can have the spectrometer begin taking data.
AVS_PrepareMeasure(handle, measconfig)  # send the measurement parameters to the spectrometer
AVS_Measure(handle, 0, 1)

# A 50 ms exposure time and single scan will complete quickly.
# Longer experiments may require us to wait for the collection to finish
# so we’ll add a small check to see if the spectral data is ready.

data_ready = False
while not data_ready:
    data_ready = AVS_PollScan(handle) # AVS_PollScan() checks if the measurement has finished
    time.sleep(0.1)

# And now it's time to grab the data!
timestamp, spectral_data = AVS_GetScopeData(handle)

# Finally, we can plot the spectrum using matplotlib
plt.clf()
plt.plot(wavelengths[0:2048], spectral_data[0:2048])
plt.show()

print('Finito')


