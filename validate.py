import time
import allantools
import numpy as np
import matplotlib.pyplot as plt
from struct import *

f = open('FreqBTB-Tstamps_Tau1e-4s_Samples5e6_VCO_Binary', 'r+b')
BigBytesString = f.read()
f.close()

# Convert (unpack) binary string into a ASCII tuple

Bytes = len(BigBytesString)                                    # Number of bytes of the binary string
freqtimeList = unpack('>'+'dQ'*int(Bytes/16), BigBytesString)  # Bytes of the single (freq/tstamp) sample = 16 -> int(Bytes/16) is the total number of samples


# In[ ]:


# Overlapping ADEV calculation

start = time.time()

freqtimearray = np.array(freqtimeList)                     # Convert list to numpy array
freqarray=(freqtimearray[::2])                             # Take the elements of the array at even steps
timearray=(freqtimearray[1::2])                            # Take the elements of the array at odd steps
 
StepSize = (timearray[2:]-timearray[1:-1]).mean()          # Calculate the mean time delay between two consecutive frequency samples
StepSize *= 1e-12                                          # From psec -> sec (only if data are taken in PAKED form)
#print(StepSize)

freqmean = freqarray.mean()                                # Calculate the mean frequency over all the frequency samples
#print(freqmean)
                                              
Adevlist = []                                              # Define list: list=[] ; Define numpy array:  arr=np.array([])
Sampleslist = []

Taulist = [1e-4 * 10**(n/5) for n in range(20)]            # Define list of gate times (sec) - NotEvenlySpaced: 1e-5 * 2n - EvenlySpaced: 1e-5 * 10**(n/5) - ...

for tau in Taulist:
    
    # Make sure we follow through with Tratio instead of tau later, in case rounding matters
    Tratio = int(round(tau/StepSize))                      # Take the nearest integer to the Tratio (i.e. Number of freq samples inside a specific Gate time)
    #print(Tratio)
        
    # This first line is meaningless because of the following line, right?
    # The first is a simple decimation (not used) but the second averages in gate time
    # Not clear that the first is not too naive or the second too clever
    freqarr = freqarray[::Tratio]                                                                   # Define an array with number of elements = Tratio
    freqarr = np.array( [freqarray[Tratio*i:Tratio*(i+1)].mean() for i in range(len(freqarr))] )    # Average all the freq samples (# = Tratio) inside a specific Gate time
    #print(len(freqarr))
        
    print(f'Gate time = {tau} s, Samples = {len(freqarr)}') 
                
    deltafreqarray = (freqarr[1:]-freqarr[:-1])**2         # Calculate the square of all the differences of consecutive pairs of frequencies
    #print(len(deltafreqarray))
    adev = (np.mean(deltafreqarray)/2)**(1/2)              # Calculate Adev
    print(adev)
    #print(adev/freqmean)
        
    Adevlist.append(adev)                                                 # Append a value to the list at each iteration
    Sampleslist.append(len(freqarr))                                      # Creating the Samples array for each Gate time (append the number of samples value at each iteration)
         
    #end = time.time()     
    #print(end-start)                                                      # Acquisition time for each Gate time (sec)

end = time.time()
print("Time to run our own calculation:")
print(end-start)
print()

# fractional frequency y
y1 = np.array(freqarray)/freqmean

start = time.time()
# allantools.adev gives tuple of (taus, adevs, adev_errors, ns)
libans = allantools.adev(y1, rate=(1/StepSize), data_type='freq', taus=Taulist)
end = time.time()
print("Time to run allantime library calculation:")
print(end-start)

# convert our calc to fractional frequency

xl = libans[0]
yl = libans[1]
xr = Taulist
yr = Adevlist/freqmean

#plt.plot(xl, yl, 'ko')
#plt.plot(xr, yr, 'ro')
plt.errorbar(xl, yl, yerr=libans[2], fmt='ko')
plt.plot(xr, yr, 'ro')
plt.xscale('log')
plt.yscale('log')
plt.title('ADEV vs. gate time, two calculations')
plt.xlabel('sec')
plt.show()
#plt.savefig('adev.png')
