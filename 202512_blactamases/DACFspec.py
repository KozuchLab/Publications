# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import math
from scipy import signal

####################################################
# INPUT:
# Loading data from GROMACS or TINKER file
input_file          = 'ProSH_md_ALL.dipole'
dt_in_fs            = 1
square_window_in_fs = 10000
# Windows for FFT: 'hann' 'blackmannharris' 'hamming' 'gauss' 'lorentz' 'square'
window_type         = 'hann' 

Debye_to_Cm = 3.33564e-30
c = 2.99792458e10

####################################################

def spec_axes():
    
    '''
    Template spectrum with 1 cm-1 datapoint spacing and zero intensites.
    '''
    
    v_min = 0
    v_max = 4000
    v_axis = np.arange(v_min,v_max+1)
    I_axis = np.zeros(len(v_axis))
    
    return v_axis, I_axis

def fit_to_spec_axes(in_v_axis,in_I_axis):
    
    '''
    Spectra obtained from FFA are formatted to uniform 1 cm-1 data point spacing, and
    to spectral resolution specified in the vibspek.json file.
    '''
    spec_res_in_cm1 = 2
    
    out_v_axis, out_I_axis = spec_axes()
    
    for i in range(len(out_v_axis)):
        temp_I = []
        lower_v = out_v_axis[i] - spec_res_in_cm1
        upper_v = out_v_axis[i] + spec_res_in_cm1
        for v,I in zip(in_v_axis,in_I_axis):
            if abs(v) >= lower_v and abs(v) <= upper_v:
                temp_I.append(I)
        if temp_I != []:        
            out_I_axis[i] = np.mean(temp_I)
        
    return out_v_axis, out_I_axis

dipole_data = np.loadtxt(input_file)#  , skiprows= 27)
if len(dipole_data[1,:]) > 3:
    total_dipole_moment = dipole_data[:,1:4]
elif len(dipole_data[1,:]) == 3:
    total_dipole_moment = dipole_data[:,0:3]

deriv_dipole_moment_x = np.gradient(total_dipole_moment[:,0])
deriv_dipole_moment_y = np.gradient(total_dipole_moment[:,1])
deriv_dipole_moment_z = np.gradient(total_dipole_moment[:,2])
deriv_dipole_moment_x = deriv_dipole_moment_x/dt_in_fs
deriv_dipole_moment_y = deriv_dipole_moment_y/dt_in_fs
deriv_dipole_moment_z = deriv_dipole_moment_z/dt_in_fs


#Calc. Dipole Autocorrelation Function
def calc_autocorrelation_function (u1, u2):
    autocorrelation_function = np.correlate(u1, u2, 'full')
    autocorrelation_function = autocorrelation_function[autocorrelation_function.size // 2:]
    n= np.arange (len(autocorrelation_function), 0, -1)
    return autocorrelation_function/n


dacf_x = calc_autocorrelation_function(deriv_dipole_moment_x, deriv_dipole_moment_x)
dacf_y = calc_autocorrelation_function(deriv_dipole_moment_y, deriv_dipole_moment_y)
dacf_z = calc_autocorrelation_function(deriv_dipole_moment_z, deriv_dipole_moment_z)
dacf_full = dacf_x + dacf_y + dacf_z

# Square window:
correct_time = int(1/dt_in_fs)
square_window = square_window_in_fs*correct_time
if len(dacf_full) > square_window:
    crop = square_window
else:
    crop = len(dacf_full)
dacf = dacf_full[0:crop]

# # Gaussian window 
if window_type == 'gauss':
    window = signal.windows.gaussian(crop*2,1,0.4 * crop/2)[crop:]

# Lorentzian/exponential window 
elif window_type == 'lorentz':
    window = signal.windows.exponential(crop,center=0,tau=1000*correct_time)
    
# Hamming window 
elif window_type == 'hamming':
    window = signal.windows.hamming(crop*2)[crop:]
    
# Hann window 
elif window_type == 'hann':
    window = signal.windows.hann(crop*2)[crop:]
    
# Blackmann-Harris window 
elif window_type == 'blackmannharris':
    window = signal.windows.blackmanharris(crop*2)[crop:]
    
elif window_type == 'square':
    window = 1

dacf = dacf * window

plt.plot(dacf)
plt.plot(window*max(dacf))
plt.ylabel('DACF')
plt.xlabel('Time point')
plt.show()

#Zero padding and shift
length_dacf = len(dacf)
desired_length = math.log(len(dacf),2)
desired_length = np.ceil(desired_length) + 4
desired_length = int(2**desired_length)

dacf_padded = np.pad(dacf, (0, desired_length - length_dacf))
ir_spectra = np.fft.fft(dacf_padded)
shifted_fft = np.fft.fftshift(ir_spectra)

#x-axis
N = len(dacf_padded)
dt = dt_in_fs*1e-15
frequency = np.fft.fftfreq(N, dt)
shifted_fft_freq = np.fft.fftshift(frequency)
wavenumber = shifted_fft_freq/c

#plotting
x = wavenumber
y = np.real(shifted_fft)*(frequency**2)
y = y*Debye_to_Cm**2/(1e-15)**2
x, y = fit_to_spec_axes(x,y)
plt.plot(x,y)
plt.ylabel('Intensity / a.u.')
plt.xlabel('Wavenumber / cm-1')
plt.show()

output = []
for v,I in zip(x,y):
    output.append([v,I])
np.savetxt(input_file.split('.')[0]+'.spec',output)

