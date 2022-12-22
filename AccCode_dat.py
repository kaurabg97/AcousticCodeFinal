# Title : Code for Acoustics
#Coded by: Kaurab Gautam 
#%% Pacakages installation
#get_ipython().system('pip install npTDMS') # Installating nptdms package

import pandas as pd # pandas import
import numpy as np # numpy import
import matplotlib.pyplot as plt # matplotlib import
from nptdms import TdmsFile # tdmsfile import
import math # maths import
from sklearn.metrics import mean_squared_error #scikit learn import
from numpy.fft import fft, ifft # FFT import

#%% Input required
T_aq = 2 #sec
F_aq = 204800 ## Acquizition frequency
N = T_aq*F_aq 
n_mics = 14 #no of mics

#%%  Import and reading of the tdms file
#function to calculate RMS values
def RMS(X):
  N = len(X)
  ms = 0
  for i in range(N):
    ms = ms+X[i]**2
  ms = ms/N
  return ms**0.5

data = np.loadtxt('C:/Users/Kaurab Gautam/OneDrive - University of Cincinnati/Desktop/ResearchFi/TestsAndResults/Baseline/Version2/V2/v3.dat', unpack = True)
pressure = np.delete(data, 0, 0)
#%%%
pre_fluc=np.zeros((N,n_mics))   # 
OASPL_rms = np.zeros((1,n_mics)) # creating array to store OASPL_rms for each microphons
OASPL=np.zeros((1,n_mics))  # creating array to store OASPL_rms for each microphons
#P.shape # P is an array of shape 1024000*12 and has pressure values

## Normalizing the value of P and assigning new vwctor Pp
sum = 0
for j in range(n_mics):
  p_mean=np.mean(pressure[:,j]) # makes average of each colums ( which has 1024000 data vales)
  for i in range (N):
    pre_fluc[i,j]=pressure[i,j]-p_mean # Pp is pressure fluctuation
  
  OASPL_rms[0,j] = RMS(pre_fluc[:,j])
OASPL_rms

# this gives OASPL in pa
#%% Calculation OASPL in db
P_ref = 20e-6 # referecne pressure
OASPL_rms = 20*np.log10(OASPL_rms/P_ref)
OASPL_rms
OASPL_rms = OASPL_rms.flatten()
angles = [152, 148, 144, 140, 136, 132, 126, 120, 116, 110, 105, 90, 70, 45]
#Angles = [45, 70, 90, 105, 110, 116, 120, 126, 132, 136, 140, 144, 148, 152]
#%% Renaming the files
OASPL_DD3p0= OASPL_rms
DF = pd.DataFrame(OASPL_rms)

DF.to_csv("OASPL.csv")
#%% Plotting results
plt.plot(angles,OASPL_DD3p0,"-o",label='DA type D')
plt.ylabel("OASPL, dB")
plt.xlabel("Angle, degree")
plt.title("Performance of double array MVGs")
plt.legend( loc = (1.01,0))
plt.savefig('OASPL for nozzle A')
plt.grid()
#%% OASPL calculation from another technique and and SPL plots

from matplotlib.ticker import FuncFormatter
## Generating SPL graph
blk_size = 5120  ## this is block size whicah is defined as bs = N/nb where nb is number of block = 200
F_res = F_aq/blk_size # Frequency data resolution 
ffi = np.linspace(F_res, F_aq, blk_size) # integration domain
Fb = np.zeros((np.int_(blk_size), np.int_(N*F_res/F_aq)))   # Fb = np.zeros(bs,nb) Coz we know nb = N*(df/fs))
SPL = np.zeros((np.int_(blk_size/2), np.int_(n_mics))) # Calculating SPL for each frequency
Fdom = ffi[0:np.int_(blk_size/2)] # Frequency domain 
for u in range(n_mics): 
  X = pre_fluc[:,u]
  b = 0
  for k in range(np.int_((N*F_res)/F_aq)):
    a = np.int_(b)
    b = np.int((k+1)*blk_size)
    Fb[:,k]=np.abs(fft(X[a:b])/blk_size)
  Fb = Fb**2
  F = np.mean(Fb,axis=1)
  PSDB = (1/F_res)*F
  q = np.trapz(PSDB,ffi)
  OASPL[0,u] = 10*np.log10(q/P_ref**2)
  OASPL
  SPL[:,u] = 10*np.log10((2*F_res*PSDB[0:np.int_(blk_size/2)])/(P_ref**2))

#%% Renaming the SPL calculated here
SPLDD= SPL

#%% To plot all three graphs

fig, axs = plt.subplots(3, figsize= (4,5))
plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=1.3, 
                    wspace=0.6, 
                    hspace=0.6)

axs[0].semilogx(Fdom,SPLDD[:,0], color = 'k',label = '152 degree')
axs[0].grid(b=True, which='major', color='#666666', linestyle='-')
axs[0].minorticks_on()
axs[0].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
axs[0].set_title("152 degree microphone")
axs[0].set_xlabel('Frequency, Hz')
axs[0].set_ylabel('SPL, db')
axs[0].legend(loc='upper left')
plt.savefig('SPL for for nozzle for 152')

axs[1].semilogx(Fdom,SPLDD[:,12], color = 'b',label = '90 degree')
axs[1].grid(b=True, which='major', color='#666666', linestyle='-')
axs[1].minorticks_on()
axs[1].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
axs[1].set_xlabel('Frequency, Hz')
axs[1].set_ylabel('SPL, db')
axs[1].legend(loc='upper left')
axs[1].set_title("90 degree microphone")
plt.savefig('SPL for for nozzle for 90')

axs[2].semilogx(Fdom,SPLDD[:,13], color = 'r',label = '70 degree')
axs[2].grid(b=True, which='major', color='#666666', linestyle='-')
axs[2].minorticks_on()
axs[2].grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)
axs[2].set_xlabel('Frequency, Hz')
axs[2].set_ylabel('SPL, db')
axs[2].legend(loc='upper left')
axs[2].set_title("70 degree microphone")

plt.savefig('SPL for for nozzle for 70 degree')


