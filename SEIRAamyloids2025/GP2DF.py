#############################################################
#                                                           #
#       Global Peak and 2nd Derivative Fitting (GP2DF)      #
#                      (2025/02/15)                         #
#                                                           #
#       please cite the following reference:                #
#       Heermant et al., 2025, ChemRxiv, DOI:               #
#       xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx                     #
#                                                           #
#############################################################

'''This script performs a global fitting to spectra AND their 2nd derivatives. 
This enables to extract narrow bands from a background of low curvature or broad
bands. The concept is described in the SI to XXXXXXXXXXXXXXXXXXXXXXXXX.
Feel free to reach out to the authors for advice in proper use of this tool with
all its options.

This script goes as follows:
    1. Copy to a folder with spectra with same naming and a index specifying the sequence.
    2. Set up your fitting model, i.e. species with narrow (here sharp) bands or broad ones.
    3. Provide boundaries if species can appear only as positives or also negatives (e.g. as in difference spectra).
    4. Automatized baselining can be tricky. There are few options.

The routine goes as follows:
    1. Coarse baselining and normalization of spectra - this necessary so that fitting is not biased towards large intensities.
       The baseline is discarded and normalization factors are stored.
    2. First proper baselining using a polynom through specified regions.
    3. Intensities of species are refined (guessed) for fitting parameters as provided.
    4. First true fitting step, where intensities, single bands within species, and baseline is refined - some constraints are present.
    5. Unconstraint fitting step, where intensities, single bands within species, and baseline is refined.
    6. Normalization factors are reapplied and results are printed.'''

import numpy as np
import matplotlib as mplt
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter
from datetime import datetime
from sklearn.metrics import auc

# Global variables will be replaced in a newer version
global input_sharp
global input_broad
global input_combined
global num_of_specs
global factor_2ndderiv
global polynom_order
global input_fraction_bounds
global sg_window 
global polynom_midpoint
global name_species
global factors_to_1
global new_factors_to_1
global output_folder
global sharp_max
global broad_min

#################
####  Input  ####

# Input info on the spectra that are supposed to be fit:
input_pref = 'out_Spec_avg'
first_index = 1
last_index = 50
use_every_nth = 1             #this must fit with first_index and last_index, i.e. (last - first)/nth must be an integer!
wavenum_top = 1850
wavenum_bot = 1300
output_folder = 'GP2DF_out'
species_for_ratio = [[1,2,3],[4,5]]

#############################
#### Spectral Components ####
name_species = ['Lipid',\
                'Amide I',\
                'Amide I',\
                'Amide I',\
                'Amide II',\
                'Amide II',\
                'Lip',\
                'Lip',\
                'Lip',\
                'Lip',\
                'Lip',\
                'Water']

# Species with sharp peaks < 20 cm-1
input_sharp = [
                [[1741.87, 0.136462, 6.77, 1]], 
                [[1700.69, 0.003636, 10.81, 1]], 
                [[1667.29, 0.096147, 19.96, 1]], 
                [[1650.3, 0.047986, 8.92, 1], [1634.09, 0.036569, 7.05, 1]], 
                [[1556.13, 1.09, 17.38, 1]], 
                [[1526.87, 1.04, 20.0, 1]], 
                [[1473.99, 0.510169, 11.86, 1], [1463.46, 0.334115, 5.93, 1]], 
                [[1406.44, 0.702937, 19.8, 1]], 
                [[1401.2, 0.078052, 7.27, 1]], 
                [[1349.91, 0.019195, 5.97, 1]], 
                [[1342.93, 0.108322, 18.69, 1]]
                ]

# Species with wide peaks > 35 cm-1
input_broad = [ [[1644.86, 10.06, 41.08, 1]] ]

# Specify if the 'concentration' can be only positiv, negative or both. First and second row are lower and upper boundary, respectively. 
# Number of boundaries has to match the number of overall species (sharp and broad):
input_fraction_bounds = [[      0,      0,      0,      0,      0,-np.inf,-np.inf,      0,      0,      0,      0,-np.inf],
                         [ np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf, np.inf]]

#########################
#### Polynom Guess ######
polynom_guess_range = [x for x in range(wavenum_bot, wavenum_bot+1)] +\
                      [x for x in range(1504, 1505)] +\
                      [x for x in range(1780, wavenum_top+1)]

# polynom_guess_range = [x for x in range(wavenum_bot, wavenum_bot+3)] + [x for x in range(1780, wavenum_top+1)]

#########################
#### Fitting Options ####

# fit options:
band_shape = 'gauss'
peak_fix = 'medium' # 'tight' or 'loose'
factor_2ndderiv = 100
sg_window = 21
step_size_fit = 0.1
sharp_max = 20
broad_min = 35
init_fraction_guess = 0.1
# init_fraction_guess = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.0,0.0]

# Polynomial Baseline:
polynom_order = 3
polynom_midpoint = (wavenum_top+wavenum_bot)/2
#polynom_midpoint = wavenum_top

# Other stuff:
np.set_printoptions(precision=4)
num_of_specs = int((last_index - first_index)/use_every_nth + 1)
try: os.mkdir(output_folder) 
except: print('Folder %s already exists' %(output_folder))
#################
fontsize = 7
mplt.rcParams['svg.fonttype'] = 'none'
mplt.rcParams['font.family'] = ['Arial', 'sans-serif']
mplt.rcParams['mathtext.fontset'] = 'custom'
mplt.rcParams['mathtext.it'] = 'Arial:italic'
mplt.rcParams['mathtext.rm'] = 'Arial'
mplt.rcParams['figure.dpi'] = 600
mplt.rcParams['savefig.dpi'] = 600
mplt.rcParams['lines.markersize'] = 2
mplt.rcParams['lines.linewidth'] = 0.7
mplt.rc('xtick', labelsize=fontsize) 
mplt.rc('ytick', labelsize=fontsize) 

def load_specs(input_pref,index_range,use_every_nth, man_factors_to_1):

    in_set_wavenum = []
    in_set_spectra = []
    
    first = index_range[0]
    last  = index_range[1]
    
    input_filenames = []
    factors_to_1 = []
    if man_factors_to_1 == []:
        for a in range(first, int(last+1), use_every_nth):
            suffix = str(a).zfill(3)+'.txt'
            filename = input_pref + suffix
            x_temp,y_temp = np.loadtxt(filename).transpose()
            input_filenames.append(filename)
            sum_of_y = sum(abs(y_temp))
            norm_y = y_temp/sum_of_y
            factors_to_1.append(sum_of_y)
            in_set_wavenum.append(x_temp)
            in_set_spectra.append(norm_y)
    else:
        i = 0
        for a in range(first, int(last+1), use_every_nth):
            suffix = str(a).zfill(3)+'.txt'
            filename = input_pref + suffix
            x_temp,y_temp = np.loadtxt(filename).transpose()
            input_filenames.append(filename)
            norm_y = y_temp/man_factors_to_1[i]
            in_set_wavenum.append(x_temp)
            in_set_spectra.append(norm_y)
            i = i + 1
        factors_to_1 = man_factors_to_1
            
        
        # print(sum_of_y)
    return in_set_wavenum, in_set_spectra, factors_to_1, input_filenames

def cut_spectra(in_wavenum,in_inten,maxim,minim):
    new_wavenum = []
    new_inten = []
    for wavenum, inten in zip(in_wavenum, in_inten):
        temp_wavenum = []
        temp_inten = []
        for x, y in zip(wavenum,inten):
            if x > maxim:
                continue
            elif x < minim:
                continue
            else:
                temp_wavenum.append(x)
                temp_inten.append(y)
        new_wavenum.append(temp_wavenum)
        new_inten.append(temp_inten)
            
    return new_wavenum, new_inten

def cut_spectra_for_ploynom_guess(in_wavenum,in_inten,maxim,minim,polynom_guess_range):
    new_wavenum = []
    new_inten = []
    for wavenum, inten in zip(in_wavenum, in_inten):
        temp_wavenum = []
        temp_inten = []
        for x, y in zip(wavenum,inten):
            if x > maxim:
                continue
            elif x < minim:
                continue
            elif round(x) not in polynom_guess_range:
                continue
            else:
                temp_wavenum.append(x)
                temp_inten.append(y)
        new_wavenum.append(temp_wavenum)
        new_inten.append(temp_inten)
            
    return new_wavenum, new_inten

def combine_spectra(set_wavenum, set_spectra):

    x = np.array(set_wavenum)
    y = np.array(set_spectra)

    return x.flatten().tolist(), y.flatten().tolist()

def combine_inten_only(set_spectra):

    y = np.array(set_spectra)

    return y.flatten().tolist()

def sepate_spectra(lin_wavenum, lin_spectra, num):
    
    length = int(len(lin_spectra)/num)
    # print(length)
    
    x = []
    y = []
    
    for i in range(0,num):
        x.append(lin_wavenum[i*length : (i+1)*length])
        y.append(lin_spectra[i*length : (i+1)*length])
        
    return x, y

def second_derivative(inten):
    global factor_2ndderiv
    global sg_window
    
    new_inten = savgol_filter(inten, sg_window, 3)
    
    sec_inten = np.gradient(np.gradient(new_inten))
    
    return sec_inten*factor_2ndderiv

def deriv_spectra(set_inten):
    
    diff_spec = []
    for spec in set_inten:
        diff_spec.append(second_derivative(spec))
    
    return diff_spec

def polynom(x, p):
    global polynom_midpoint
    
    order = len(p)
    
    x = np.array(x)  
    y = np.zeros(len(x))
    
    mid_x = polynom_midpoint
    
    for i in range(0, order):
        y = y + (p[i]/(i+1)**(2*i))*(x-mid_x)**i

    return y

def gaussian(x, sigma, x_0):
    y = (1/sigma/np.sqrt(2*np.pi))*np.exp(-(x-x_0)**2/2/sigma**2)     
    return y

def lorentzian(x, sigma, x_0):
    y = (sigma/np.pi)*(1/((x - x_0)**2 + sigma**2))
    return y     

def pseudovoigt(x, frac, sigma, x_0):
    sigma_G = sigma
    sigma_L = sigma * np.sqrt(2*np.log(2))
    y = frac * gaussian(x, sigma_G, x_0) + (1 - frac) * lorentzian(x, sigma_L, x_0)
    return y  

def gauss_band(x, x_0, A, sigma):
    y = A * gaussian(x, sigma, x_0) 
    return y

def pseudovoigt_band(x, x_0, A, sigma, frac):
    x = np.array(x)
    y = A * pseudovoigt(x, frac, sigma, x_0)
    return y

def linearized_params(input_sharp, input_broad):
    input_combi = [] 
    input_linear = []
    for el1 in input_sharp:
        input_combi.append(el1)
        for el2 in el1:
            for p in el2:
                input_linear.append(p)
    for el1 in input_broad:
        input_combi.append(el1)
        for el2 in el1:
            for p in el2:
                input_linear.append(p)
                
    return input_combi, input_linear

def unlinearized_params(input_linear):
    global input_combined
    
    counter = 0
    level_0 = []
    for el1 in input_combined:
        level_1 = []
        for el2 in el1:
            level_2 = []
            for p in el2:
                level_2.append(input_linear[counter])
                counter += 1
            level_1.append(level_2)
        level_0.append(level_1)
        
    return level_0

def unlinearized_params_and_round(input_linear):
    global input_combined
    
    counter = 0
    level_0 = []
    for el1 in input_combined:
        level_1 = []
        for el2 in el1:
            level_2 = []
            for p in el2:
                if input_linear[counter] > 1:
                    val = round(input_linear[counter],2)
                elif input_linear[counter] > 0.99:
                    val = 1
                elif input_linear[counter] > 0.001:
                    val = round(input_linear[counter],6)
                else:
                    val = round(input_linear[counter],8)
                level_2.append(val)
                counter += 1
            level_1.append(level_2)
        level_0.append(level_1)
        
    return level_0

def number_of():
    global input_combined
    
    number_of_species = len(input_combined)
    counter_param = 0
    counter_bands = 0
    for el1 in input_combined:
        for el2 in el1:
            counter_bands += 1
            for p in el2:
                counter_param += 1
    number_of_bands = counter_bands
    number_of_param = counter_param
    
    return number_of_species, number_of_bands, number_of_param

def set_polynom_input():
    global num_of_specs
    global polynom_order

    polynom_input = []
    for i in range(0, num_of_specs):
        for j in range(0, int(polynom_order+1)):
            polynom_input.append(0)
            
    return polynom_input

def set_fraction_input(frac_init):
    global num_of_specs
    global input_combined
    
    def isfloat(x):
        try:
            float(x)
            return True
        except:
            pass
        
    def islist(x):
        try:
            len(x)
            return True
        except:
            pass
        
    if isfloat(frac_init):
        number_of_species, number_of_bands, number_of_param = number_of()
        fraction_input = []
        for j in range(0, int(number_of_species)):
            for i in range(0, num_of_specs):
                fraction_input.append(frac_init)    
    elif islist(frac_init):
        number_of_species, number_of_bands, number_of_param = number_of()
        fraction_input = []
        for j in range(0, int(number_of_species)):
            for i in range(0, num_of_specs):
                fraction_input.append(frac_init[j])    
    
    return fraction_input

def bounds_param(info):
    global input_sharp
    global input_broad              
                
    if 'fix_all' in info:
        bound_low = []
        bound_up  = []
        for el1 in input_sharp:
            for el2 in el1:
                
                bound_low.append(el2[0]-0.000001)
                bound_up.append(el2[0]+0.000001)
                
                bound_low.append(el2[1]-0.000001)
                bound_up.append(el2[1]+0.000001)
                if el2[2]-0.000001 < 2:
                    bound_low.append(2)
                else:
                    bound_low.append(el2[2]-0.000001)
                if el2[2]+0.000001 > sharp_max:
                    bound_up.append(sharp_max)
                else:
                    bound_up.append(el2[2]+0.000001)
                    
                if 'gauss' in info:
                    bound_low.append(el2[3]-0.000001)
                    bound_up.append(1)
                else:
                    bound_low.append(el2[3]-0.000001)
                    bound_up.append(el2[3]+0.000001)
                    

                    
    
        for el1 in input_broad:
            for el2 in el1:
                bound_low.append(el2[0]-0.000001)
                bound_up.append(el2[0]+0.000001)
                
                bound_low.append(el2[1]-0.000001)
                bound_up.append(el2[1]+0.000001)
                if el2[2]-0.000001 < broad_min:
                    bound_low.append(broad_min)
                else:
                    bound_low.append(el2[2]-0.000001)
                if el2[2]+0.000001 > 100:
                    bound_up.append(100)
                else:
                    bound_up.append(el2[2]+0.000001)
                if 'gauss' in info:
                    bound_low.append(el2[3]-0.000001)
                    bound_up.append(1)
                else:
                    bound_low.append(el2[3]-0.000001)
                    bound_up.append(el2[3]+0.000001)

                    
    else:
        
        bound_low = []
        bound_up  = []
        
        if 'gauss' in info:
            bl_frac = 0.99999999
            bu_frac = 1.00000000
        if 'lorentz' in info:
            bl_frac = 0.00000000
            bu_frac = 0.00000001
        if 'voigt' in info:
            bl_frac = 0.00000000
            bu_frac = 1.00000000
            
        if 'tight' in info:
            bl_peak = -1
            bu_peak = +1
        if 'medium' in info:
            bl_peak = -5
            bu_peak = +5
        if 'loose' in info:
            bl_peak = -10
            bu_peak = +10
    
        for el1 in input_sharp:
            for el2 in el1:
                # for p in el2:
                    bound_low.append(el2[0] + bl_peak)
                    bound_up.append(el2[0] + bu_peak)
                    bound_low.append(0)
                    bound_up.append(np.inf)
                    bound_low.append(3)
                    bound_up.append(sharp_max)
                    bound_low.append(bl_frac)
                    bound_up.append(bu_frac)
                    
    
        for el1 in input_broad:
            for el2 in el1:
                # for p in el2:
                    bound_low.append(el2[0] + bl_peak)
                    bound_up.append(el2[0] + bu_peak)
                    bound_low.append(0)
                    bound_up.append(np.inf)
                    bound_low.append(broad_min)
                    bound_up.append(100)
                    bound_low.append(bl_frac)
                    bound_up.append(bu_frac)
                
    return [bound_low, bound_up]
    
def bounds_polynom(info):
    global num_of_specs
    global polynom_order 

    bound_low = []
    bound_up  = []
    for i in range(0, num_of_specs):
        # for j in range(0, int(polynom_order+1)):

            if 'fix' in info:               
                # bound_low.append(-np.inf)
                # bound_up.append (np.inf)   
                # bound_low.append(-np.inf)
                # bound_up.append (np.inf)   
        
                # for i in range(2, int(polynom_order+1)):
                #     bound_low.append(-(10**(-5*i)))
                #     bound_up.append ( 10**(-5*i))   
                for i in range(0, int(polynom_order+1)):
                    bound_low.append(-(10**(-4*i)))
                    bound_up.append ( (10**(-4*i)))   
                    
            elif 'offset' in info:
                bound_low.append(-np.inf)
                bound_up.append (np.inf)   
        
                for i in range(1, int(polynom_order+1)):
                    bound_low.append(-(10**(-4*i)))
                    bound_up.append ( 10**(-4*i))  
         
            else:
                for i in range(0, int(polynom_order+1)):
                    bound_low.append(-np.inf)
                    bound_up.append (np.inf)    

    return [bound_low, bound_up]

def bounds_fractions(info, fracs):
    global input_fraction_bounds
    global num_of_specs
    
    number_of_species, number_of_bands, number_of_param = number_of()
            
    if 'fix_all' in info:
        bound_low = []
        bound_up  = []
        for b in fracs:
                bound_low.append(b-0.000001)    
                bound_up.append(b+0.000001)  

    elif 'zero' in info:
        bound_low = []
        bound_up  = []
        for j, bl, bu in zip(range(0, int(number_of_species)), input_fraction_bounds[0], input_fraction_bounds[1]):
            for i in range(0, num_of_specs):
                bound_low.append(-0.00000000001)    
                bound_up.append(0.00000000001)   
                
    else:
        bound_low = []
        bound_up  = []
        for j, bl, bu in zip(range(0, int(number_of_species)), input_fraction_bounds[0], input_fraction_bounds[1]):
            for i in range(0, num_of_specs):
                bound_low.append(bl)    
                bound_up.append(bu)   

    return [bound_low, bound_up]    

def p0_and_bounds(p0_param, p0_polynom, p0_fractions, b_param, b_polynom, b_fractions):
    p0 = []
    bounds_low = []
    bounds_up  = []
    
    # print(b_polynom)
    # print(b_fractions)
    
    for a, b, c in zip(p0_param, b_param[0], b_param[1]):
        p0.append(a)
        bounds_low.append(b)
        bounds_up.append(c)
    for a, b, c in zip(p0_polynom,  b_polynom[0],  b_polynom[1]):
        p0.append(a)
        bounds_low.append(b)
        bounds_up.append(c)
    for a, b, c in zip(p0_fractions,  b_fractions[0],  b_fractions[1]):
        p0.append(a)
        bounds_low.append(b)
        bounds_up.append(c)
        
    return p0, [bounds_low, bounds_up]

def polynom_over_time(lin_waven, *p):
    global num_of_specs
    global polynom_order

    # print('######')   


    if len(p) == 1:
        p = p[0]
        
       
    number_of_species, number_of_bands, number_of_param = number_of()   
    waven = lin_waven[0:int(len(lin_waven)/num_of_specs/2)]
             
    p_polynom = p
    p_polynom_use = []
    for i in range(0, num_of_specs):
        p_polynom_use.append([x for x in p_polynom[i*(polynom_order+1) : (i+1)*(polynom_order+1)]])
        # print([i*(polynom_order+1), (i+1)*(polynom_order+1)])
    # print(p_polynom_use)
            
    spectra = []
    poly = []
    spectra_with_poly = []
    for i in range(0, num_of_specs):
        sub_poly = polynom(waven, p_polynom_use[i])
        sub_spectra = np.zeros(len(waven))
        spectra_with_poly.append(sub_spectra + sub_poly)    
        spectra.append(sub_spectra)
        poly.append(sub_poly)

    deriv = deriv_spectra(spectra)
    lin_inten = combine_inten_only(spectra_with_poly)
    lin_2ndin = combine_inten_only(deriv)
    lin_i = combine_inten_only([lin_inten,lin_2ndin])
    
    # fig, axs = plt.subplots(2, sharex=True)
    # for s, p, d in zip(spectra_with_poly, poly, deriv):
    #     axs[0].plot(waven, s)
    #     axs[0].plot(waven, p)
    #     axs[1].plot(waven, d)
    # plt.show()
    
    # plt.plot(lin_i)
    # plt.show()
       
    return lin_i

def spectra_over_time(lin_waven, *p):
    global num_of_specs
    global polynom_order

    # print('######')   


    if len(p) == 1:
        p = p[0]
        
       
    number_of_species, number_of_bands, number_of_param = number_of()   
    waven = lin_waven[0:int(len(lin_waven)/num_of_specs/2)]
        
    p_param = p[0 : number_of_param]
    p_param_use = unlinearized_params(p_param)
    # print(p_param_use)
    # print('Print for debug:')
    # print([0,int(len(lin_waven)/num_of_specs/2)])
    # print([0,number_of_param])  
      
    p_polynom = p[number_of_param : number_of_param+(polynom_order+1)*num_of_specs]
    # print([number_of_param,number_of_param+(polynom_order+1)*num_of_specs])
    p_polynom_use = []
    for i in range(0, num_of_specs):
        p_polynom_use.append([x for x in p_polynom[i*(polynom_order+1) : (i+1)*(polynom_order+1)]])
        # print([i*(polynom_order+1), (i+1)*(polynom_order+1)])
    # print(p_polynom_use)
        
    p_fractions = p[-number_of_species*num_of_specs : ]
    # print([-number_of_species*num_of_specs, len(p)])
    p_fractions_use = []
    for i in range(0, number_of_species):
        p_fractions_use.append([x for x in p_fractions[i*(num_of_specs) : (i+1)*(num_of_specs)]])
        # print([i*(num_of_specs), (i+1)*(num_of_specs)])
    # print(p_fractions_use)
    
    

    
    species_spectra = []
    for species in p_param_use:
        comp_spec = np.zeros(len(waven))
        for comp in species:
            # print(comp)
            comp_spec = comp_spec + pseudovoigt_band(waven, comp[0], comp[1], comp[2] , comp[3])
        species_spectra.append(comp_spec)
        
    # for s in species_spectra:
    #     plt.plot(waven, s)
    # plt.show()
    
    # print(p_param_use)
    # print(p_polynom_use)
    # print(p_fractions_use)


    spectra = []
    poly = []
    spectra_with_poly = []
    for i in range(0, num_of_specs):
        sub_poly = polynom(waven, p_polynom_use[i])
        sub_spectra = np.zeros(len(waven))
        for j in range(0, number_of_species):
            sub_spectra = sub_spectra + p_fractions_use[j][i] * species_spectra[j]
        spectra_with_poly.append(sub_spectra + sub_poly)    
        spectra.append(sub_spectra)
        poly.append(sub_poly)

    deriv = deriv_spectra(spectra)
    lin_inten = combine_inten_only(spectra_with_poly)
    lin_2ndin = combine_inten_only(deriv)
    lin_i = combine_inten_only([lin_inten,lin_2ndin])
    
    # fig, axs = plt.subplots(2, sharex=True)
    # for s, p, d in zip(spectra_with_poly, poly, deriv):
    #     axs[0].plot(waven, s)
    #     axs[0].plot(waven, p)
    #     axs[1].plot(waven, d)
    # plt.show()
    
    # plt.plot(lin_i)
    # plt.show()
       
    return lin_i

def output_fit_spectra(lin_waven, *p):
    global num_of_specs
    global polynom_order
    global input_sharp
    global input_broad

    # print('######')    

    if len(p) == 1:
        p = p[0]
        
       
    number_of_species, number_of_bands, number_of_param = number_of()   
    waven = lin_waven[0:int(len(lin_waven)/num_of_specs/2)]
        
    p_param = p[0 : number_of_param]
    p_param_use = unlinearized_params(p_param)
      
    p_polynom = p[number_of_param : number_of_param+(polynom_order+1)*num_of_specs]
    p_polynom_use = []
    for i in range(0, num_of_specs):
        p_polynom_use.append([x for x in p_polynom[i*(polynom_order+1) : (i+1)*(polynom_order+1)]])
        
    p_fractions = p[-number_of_species*num_of_specs : ]
    p_fractions_use = []
    for i in range(0, number_of_species):
        p_fractions_use.append([x for x in p_fractions[i*(num_of_specs) : (i+1)*(num_of_specs)]])    

    species_waven = []
    for species in p_param_use:
        species_waven.append(waven)
    
    species_spectra = []
    comp_spec = []
    for species in p_param_use:
        comp_spec = np.zeros(len(waven))
        for comp in species:
            comp_spec = comp_spec + pseudovoigt_band(waven, comp[0], comp[1], comp[2] , comp[3])
        species_spectra.append(comp_spec)
                   
    sharp_species_spectra = []
    comp_spec = []
    for species in p_param_use[0:len(input_sharp)]:
        comp_spec = np.zeros(len(waven))
        for comp in species:
            comp_spec = comp_spec + pseudovoigt_band(waven, comp[0], comp[1], comp[2] , comp[3])
        sharp_species_spectra.append(comp_spec)
        
    broad_species_spectra = []
    comp_spec = []    
    for species in p_param_use[-len(input_broad):]:
        comp_spec = np.zeros(len(waven))
        for comp in species:
            comp_spec = comp_spec + pseudovoigt_band(waven, comp[0], comp[1], comp[2] , comp[3])
        broad_species_spectra.append(comp_spec)

    spectra = []
    poly = []
    spectra_with_poly = []
    spectra_broad_only = []
    spectra_sharp_only = []
    for i in range(0, num_of_specs):
        sub_poly = polynom(waven, p_polynom_use[i])
        sub_spectra = np.zeros(len(waven))
        for j in range(0, number_of_species):
            sub_spectra = sub_spectra + p_fractions_use[j][i] * species_spectra[j]
        spectra_with_poly.append(sub_spectra + sub_poly)    
        spectra.append(sub_spectra)
        poly.append(sub_poly)
        
    sub_spectra = []
    for i in range(0, num_of_specs):
        sub_spectra = np.zeros(len(waven))
        for j in range(0, len(sharp_species_spectra)):
            sub_spectra = sub_spectra + p_fractions_use[j][i] * species_spectra[j]
        spectra_sharp_only.append(sub_spectra)   
        


    sub_spectra = []        
    for i in range(0, num_of_specs):
        sub_spectra = np.zeros(len(waven))
        for j in range(len(sharp_species_spectra), number_of_species):
            sub_spectra = sub_spectra + p_fractions_use[j][i] * species_spectra[j]
        spectra_broad_only.append(sub_spectra)    

    # for s in spectra_sharp_only:
    #     plt.plot(waven, s)
    # plt.show()
    # for s in spectra_broad_only:
    #     plt.plot(waven, s)
    # plt.show()

    deriv = deriv_spectra(spectra)
       
    
    return species_waven, spectra_with_poly, spectra, poly, deriv, species_spectra, spectra_sharp_only, spectra_broad_only



def R2_RMSD(v,v_fit):
    SumSqr_res = np.sum((np.subtract(v,v_fit))**2)
    SumSqr_tot = np.sum((np.subtract(v,np.mean(v)))**2)   
    R2         = 1 - SumSqr_res/SumSqr_tot
    RMSD       = np.sqrt(SumSqr_res)/len(v_fit)
    return R2, RMSD


def show_polynoms(lin_w, set_wavenum, set_spectra, set_deriv, popt, info):
    global num_of_specs    
    global output_folder
    global factors_to_1
  
    
    result = polynom_over_time(lin_w, popt)
    fit_w, fit_result = sepate_spectra(lin_w, result, num_of_specs*2)
    fit_lin_i = fit_result[0:int(len(fit_result)/2)]
    fit_lin_deriv = fit_result[int(len(fit_result)/2):]
    R2, RMSD = R2_RMSD(lin_i, result)
    
    plt.clf()
    new_factors_to_1 = []
    fig, axs = plt.subplots(2,2, sharex=True)
    fig.suptitle("%s - R2 = %s and RMSD = %s" %(info,round(R2,3), round(RMSD, 10)))
    for w, i, i_fit, di, di_fit in zip(set_wavenum, set_spectra, fit_lin_i, set_deriv, fit_lin_deriv):
        axs[0,0].plot(w, 1000*np.array(i))
        axs[0,0].plot(w, 1000*np.array(i_fit),'--')
        axs[0,1].plot(w, 1000*(np.array(i)-np.array(i_fit)))
        axs[1,0].plot(w, 1000*np.array(di))
        axs[1,0].plot(w, 1000*np.array(di_fit),'--')
        axs[1,1].plot(w, 1000*(np.array(di)-np.array(di_fit)))
        new_factors_to_1.append(sum(abs(np.array(i)-np.array(i_fit))))
    plt.show()
    plt.savefig('zzz_00_PolynomGuess.png')
    
    factors_to_1 = np.multiply(factors_to_1,new_factors_to_1)
    
    
    print("INFO: %s - R2 = %s and RMSD = %s" %(info,round(R2,3), round(RMSD, 10)))
    
def show_polynoms_keepfactorsto1(lin_w, set_wavenum, set_spectra, set_deriv, popt, info):
    global num_of_specs    
    global output_folder
    global factors_to_1
  
    
    result = polynom_over_time(lin_w, popt)
    fit_w, fit_result = sepate_spectra(lin_w, result, num_of_specs*2)
    fit_lin_i = fit_result[0:int(len(fit_result)/2)]
    fit_lin_deriv = fit_result[int(len(fit_result)/2):]
    R2, RMSD = R2_RMSD(lin_i, result)
    
    plt.clf()
    fig, axs = plt.subplots(2,2, sharex=True)
    fig.suptitle("%s - R2 = %s and RMSD = %s" %(info,round(R2,3), round(RMSD, 10)))
    for w, i, i_fit, di, di_fit in zip(set_wavenum, set_spectra, fit_lin_i, set_deriv, fit_lin_deriv):
        axs[0,0].plot(w, 1000*np.array(i))
        axs[0,0].plot(w, 1000*np.array(i_fit),'--')
        axs[0,1].plot(w, 1000*(np.array(i)-np.array(i_fit)))
        axs[1,0].plot(w, 1000*np.array(di))
        axs[1,0].plot(w, 1000*np.array(di_fit),'--')
        axs[1,1].plot(w, 1000*(np.array(di)-np.array(di_fit)))
    plt.show()
    plt.savefig('zzz_01_PolynomGuess.png')
    
    print("INFO: %s - R2 = %s and RMSD = %s" %(info,round(R2,3), round(RMSD, 10)))
    

def write_specs(input_filenames, lin_w, set_wavenum, set_spectra, set_deriv, popt, info):
    global num_of_specs    
    
    global output_folder
    
    homepath = os.getcwd()
    os.chdir(homepath + '/' + output_folder)
    
    result = spectra_over_time(lin_w, popt)
    fit_w, fit_result = sepate_spectra(lin_w, result, num_of_specs*2)
    fit_lin_i = fit_result[0:int(len(fit_result)/2)]
    fit_lin_deriv = fit_result[int(len(fit_result)/2):]
    R2, RMSD = R2_RMSD(lin_i, result)
    
    plt.clf()
    fig, axs = plt.subplots(2,2, sharex=True)
    fig.suptitle("%s - R2 = %s and RMSD = %s" %(info,round(R2,3), round(RMSD, 10)))
    for w, i, i_fit, di, di_fit in zip(set_wavenum, set_spectra, fit_lin_i, set_deriv, fit_lin_deriv):
        axs[0,0].plot(w, 1000*np.array(i))
        axs[0,0].plot(w, 1000*np.array(i_fit),'--')
        axs[0,1].plot(w, 1000*(np.array(i)-np.array(i_fit)))
        axs[1,0].plot(w, 1000*np.array(di))
        axs[1,0].plot(w, 1000*np.array(di_fit),'--')
        axs[1,1].plot(w, 1000*(np.array(di)-np.array(di_fit)))
    plt.show()
    plt.savefig('zzz_02_Fit.png')
    
    print("INFO: %s - R2 = %s and RMSD = %s" %(info,round(R2,3), round(RMSD, 10)))

    
    waven_i, spectra_with_poly, spectra, poly, fit_deriv, species_i, spec_sharp, spec_broad = output_fit_spectra(lin_w, popt)
    
   
    for_output = []
    for_output.append(set_wavenum[0])
    for i in species_i:
        plt.plot(set_wavenum[0], i)
        for_output.append(i)
    plt.xlabel('Wavenumber / cm-1')
    plt.ylabel('Absorbance / a.u.')
    plt.legend(name_species)
    plt.savefig('zzz_species_spectra.png')
    # plt.show()
    plt.clf()
      
    np.savetxt('zzz_species_spectra.txt', np.array(for_output).transpose())
      
    for i in range(0, len(species_i)):
        np.savetxt('zzz_species_spectra_%s.txt' %(i), np.array([waven_i[i], species_i[i]]).transpose())
        
        
    for f, filename, w, i, i_fit, i_fit_bsled, bsl, d, d_fit, i_fit_sharp, i_fit_broad in zip(factors_to_1, input_filenames, set_wavenum, set_spectra, spectra_with_poly, spectra, poly, set_deriv, fit_deriv, spec_sharp, spec_broad):
        fig, axs = plt.subplots(2,sharex=True)
        axs[0].plot(w, np.multiply(f,i), 'o')
        axs[0].plot(w, np.multiply(f,i_fit), '-')
        axs[0].plot(w, np.multiply(f,bsl), '-')
        axs[0].plot(w, np.multiply(f,i_fit_bsled), '-')
        axs[1].plot(w, d, 'o')
        axs[1].plot(w, d_fit, '-')
        plt.savefig('zzz_%s.png' %(filename.split('.')[0]))
        plt.clf()
        
        np.savetxt('bsled_%s.txt' %(filename), np.array([w, np.multiply(f,i - bsl)]).transpose())
        np.savetxt('sharp_%s.txt' %(filename), np.array([w, np.multiply(f,i - bsl - i_fit_broad)]).transpose())
        np.savetxt('broad_%s.txt' %(filename), np.array([w, np.multiply(f,i - bsl - i_fit_sharp)]).transpose())
        
        
        plt.close()
    os.chdir(homepath)        
    
    return np.array(for_output).transpose()
     


def write_conc(input_pref,index_range,use_every_nth, popt, species_for_ratio):
    global name_species
    global factors_to_1
    global output_folder
    
    homepath = os.getcwd()
    os.chdir(homepath + '/' + output_folder)
    
    waven_i, spectra_with_poly, spectra, poly, fit_deriv, species_i, spec_sharp, spec_broad = output_fit_spectra(lin_w, popt)
    
    first = index_range[0]
    last  = index_range[1]
    spec_indices = []
    for a in range(first, int(last+1), use_every_nth):
        spec_indices.append(a)
    
    fractions = popt[-num_of_specs*number_of_species:]
    fracs = []
    for i in range(0, number_of_species):
        fracs.append(fractions[i*num_of_specs:(i+1)*num_of_specs])
      
    for_output = []
    for_output.append(spec_indices)
    for f in fracs:
        plt.plot(spec_indices, f*factors_to_1)
        for_output.append(f*factors_to_1)

    plt.legend(name_species)
    plt.xlabel('Spectrum number')
    plt.ylabel('Fraction of species / a.u.')
    plt.savefig('zzz_fraction_traces.png')
    plt.show()
    
    np.savetxt('zzz_fraction_traces.txt', np.array(for_output).transpose())
    
           
    # calculate the area under the the curves of each specie and multiply with the fraction of each specie
    fraction_traces = np.array(for_output).transpose()    
    species_spectra = []
    species_spectra.append(set_wavenum[0])
    for i in species_i:
        species_spectra.append(i)    
    # print(np.array(species_spectra).transpose())
    species_spectra = np.array(species_spectra).transpose()
    real_areas=[]
    real_areas.append(fraction_traces[:, 0])
    real_areas_per_species = []
    for i in range(1, number_of_species+1):
        area_species = auc(species_spectra[:, 0], species_spectra[:, i])
        real_areas_per_species.append(area_species * fraction_traces[:, i])
        plt.plot(fraction_traces[:, 0], area_species * fraction_traces[:, i], '-')
        real_areas.append(area_species * fraction_traces[:, i])
    plt.legend(name_species)
    plt.xlabel('Spectrum number')
    plt.ylabel('Relative areas / a.u.')
    plt.savefig('zzz_relative_areas.png')
    plt.show()
    #print(area_species)
    np.savetxt('zzz_relative_areas.txt', np.array(real_areas).transpose())
    plt.clf()
    # calculate the ratios of the areas of amide I to amide II

    print(real_areas_per_species)
    real_areas_per_species = np.array(real_areas_per_species).transpose()
    ratio_per_timepoint = []
    for t in range(len(spectra_with_poly)):
        temp_num = 0
        temp_den = 0
        for i in species_for_ratio[0]:
          temp_num = temp_num + real_areas_per_species[t][i]
        for i in species_for_ratio[1]:
                temp_den = temp_den + real_areas_per_species[t][i]
        ratio_per_timepoint.append(temp_num/temp_den)
    
    ratio_per_timepoint = ratio_per_timepoint
    plt.plot(fraction_traces[:, 0], ratio_per_timepoint, '-')
    # print(for_ratios[:, 2]/for_ratios[:, 3])
    plt.xlabel('Spectrum number')
    plt.ylabel('Amide I / Amide II ratio')
    plt.savefig('zzz_amide_ratios.png')
    plt.show()
    
    np.savetxt('zzz_ratio.txt', np.array(ratio_per_timepoint).transpose())
    
    os.chdir(homepath)
   
    return np.array(for_output).transpose()

def write_params(params,name_species):
    global output_folder
    
    homepath = os.getcwd()
    os.chdir(homepath + '/' + output_folder)
    
    with open('zzz_params.txt', 'w') as out_file:
        for ps,ns in zip(params,name_species):
            p_string = ns.rjust(16)
            for p in ps:
                p_string = p_string + ' | ' + str(p[0]).rjust(10) + str(p[1]).rjust(10) + str(p[2]).rjust(8) + str(p[3]).rjust(8)
            out_file.write(p_string+'\n')
        out_file.write('\n')
        out_file.write('\n')
        out_file.write(str(name_species)+'\n')
        out_file.write(str(params)+'\n')
        
            
    os.chdir(homepath)        
        
     

##############
##############
#### MAIN ####
##############
##############

now    = datetime.now()
dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_str + ' INFO: Starting.')


#####################
# Load set of spectra
init_set_wavenum, init_set_spectra, factors_to_1, input_filenames = load_specs(input_pref,[first_index, last_index],use_every_nth, [])

###########################
# Make generally compatible
init_set_wavenum, init_set_spectra = cut_spectra(init_set_wavenum, init_set_spectra, wavenum_top, wavenum_bot)
set_wavenum, set_spectra = init_set_wavenum, init_set_spectra
lin_wavenum, lin_inten = combine_spectra(set_wavenum, set_spectra)
set_deriv = deriv_spectra(set_spectra)
lin_wavenum, lin_2ndin = combine_spectra(set_wavenum, deriv_spectra(set_spectra))
lin_w, lin_i = combine_spectra([lin_wavenum,lin_wavenum], [lin_inten,lin_2ndin])

init_set_wavenum_for_polynom, init_set_spectra_for_polynom = cut_spectra_for_ploynom_guess(init_set_wavenum, init_set_spectra, wavenum_top, wavenum_bot,polynom_guess_range)
set_wavenum_for_polynom, set_spectra_for_polynom = init_set_wavenum_for_polynom, init_set_spectra_for_polynom
lin_wavenum_for_polynom, lin_inten_for_polynom = combine_spectra(set_wavenum_for_polynom, set_spectra_for_polynom)
set_deriv_for_polynom = deriv_spectra(set_spectra_for_polynom)
lin_wavenum_for_polynom, lin_2ndin_for_polynom = combine_spectra(set_wavenum_for_polynom, deriv_spectra(set_spectra_for_polynom))
lin_w_for_polynom, lin_i_for_polynom = combine_spectra([lin_wavenum_for_polynom,lin_wavenum_for_polynom], [lin_inten_for_polynom,lin_2ndin_for_polynom])



###############
#Polynom guess:
now    = datetime.now()
dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_str + ' INFO: Polynom guess...') 
    
p0_param_combined, p0_param_linear = linearized_params(input_sharp, input_broad)
input_combined = p0_param_combined
number_of_species, number_of_bands, number_of_param = number_of() 

p0_polynom = set_polynom_input()
b_polynom = bounds_polynom([''])
p0, b0 = p0_polynom, b_polynom
    
popt, pcov = curve_fit(polynom_over_time, lin_w_for_polynom, lin_i_for_polynom, p0, bounds = b0, method = 'trf', loss = 'soft_l1', diff_step=step_size_fit, x_scale='jac')

show_polynoms(lin_w, set_wavenum, set_spectra, set_deriv, popt, 'Polynom Guess')



###########################
# Re-make generally compatible
init_set_wavenum, init_set_spectra, factors_to_1, input_filenames = load_specs(input_pref,[first_index, last_index],use_every_nth, factors_to_1)
init_set_wavenum, init_set_spectra = cut_spectra(init_set_wavenum, init_set_spectra, wavenum_top, wavenum_bot)
set_wavenum, set_spectra = init_set_wavenum, init_set_spectra
lin_wavenum, lin_inten = combine_spectra(set_wavenum, set_spectra)
set_deriv = deriv_spectra(set_spectra)
lin_wavenum, lin_2ndin = combine_spectra(set_wavenum, deriv_spectra(set_spectra))
lin_w, lin_i = combine_spectra([lin_wavenum,lin_wavenum], [lin_inten,lin_2ndin])

init_set_wavenum_for_polynom, init_set_spectra_for_polynom = cut_spectra_for_ploynom_guess(init_set_wavenum, init_set_spectra, wavenum_top, wavenum_bot,polynom_guess_range)
set_wavenum_for_polynom, set_spectra_for_polynom = init_set_wavenum_for_polynom, init_set_spectra_for_polynom
lin_wavenum_for_polynom, lin_inten_for_polynom = combine_spectra(set_wavenum_for_polynom, set_spectra_for_polynom)
set_deriv_for_polynom = deriv_spectra(set_spectra_for_polynom)
lin_wavenum_for_polynom, lin_2ndin_for_polynom = combine_spectra(set_wavenum_for_polynom, deriv_spectra(set_spectra_for_polynom))
lin_w_for_polynom, lin_i_for_polynom = combine_spectra([lin_wavenum_for_polynom,lin_wavenum_for_polynom], [lin_inten_for_polynom,lin_2ndin_for_polynom])

#Polynom guess II:
now    = datetime.now()
dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_str + ' INFO: Polynom guess II...') 
    
p0_param_combined, p0_param_linear = linearized_params(input_sharp, input_broad)
input_combined = p0_param_combined
number_of_species, number_of_bands, number_of_param = number_of() 

p0_polynom = set_polynom_input()
b_polynom = bounds_polynom([''])
p0, b0 = p0_polynom, b_polynom
    
popt, pcov = curve_fit(polynom_over_time, lin_w_for_polynom, lin_i_for_polynom, p0, bounds = b0, method = 'trf', loss = 'soft_l1', diff_step=step_size_fit, x_scale='jac')

show_polynoms_keepfactorsto1(lin_w, set_wavenum, set_spectra, set_deriv, popt, 'Polynom Guess')


###############
#Fit fractions I:   
now    = datetime.now()
dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_str + ' INFO: Guess fractions...') 
    
p0_param_combined, p0_param_linear = linearized_params(input_sharp, input_broad)
input_combined = p0_param_combined
number_of_species, number_of_bands, number_of_param = number_of() 

p0_param = p0_param_linear
p0_polynom = popt
p0_fractions = set_fraction_input(init_fraction_guess)

b_param = bounds_param([band_shape,'tight']) # with 'fix_all' instead of 'tight' species won't change
b_polynom = bounds_polynom(['fix'])
b_polynom[0] = np.add(b_polynom[0],p0_polynom)
b_polynom[1] = np.add(b_polynom[1],p0_polynom)
b_fractions = bounds_fractions([''],[])
p0, b0 = p0_and_bounds(p0_param, p0_polynom, p0_fractions, b_param, b_polynom, b_fractions)

popt, pcov = curve_fit(spectra_over_time, lin_w, lin_i, p0, bounds = b0, method = 'trf', loss = 'soft_l1', diff_step=step_size_fit, x_scale='jac')

write_specs(input_filenames, lin_w, set_wavenum, set_spectra, set_deriv, popt, 'Guess fractions')



##################
#Fit fractions II:
now    = datetime.now()
dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_str + ' INFO: Refine fractions and species...') 
    
p0_param = popt[0:number_of_param]
p0_polynom = popt[number_of_param : number_of_param+(polynom_order+1)*num_of_specs]
p0_fractions = set_fraction_input(init_fraction_guess)

b_param = bounds_param([band_shape,'tight'])
b_polynom = bounds_polynom([''])
b_fractions = bounds_fractions([''],popt[-num_of_specs*number_of_species:])
p0, b0 = p0_and_bounds(p0_param, p0_polynom, p0_fractions, b_param, b_polynom, b_fractions)
p0 = popt

popt, pcov = curve_fit(spectra_over_time, lin_w, lin_i, p0, bounds = b0, method = 'trf', loss = 'soft_l1', diff_step=step_size_fit, x_scale='jac') 

write_specs(input_filenames, lin_w, set_wavenum, set_spectra, set_deriv, popt, 'Refine fractions and species')



# ######################
# #Fit spectral species: 
# b_param = bounds_param([peak_fix,band_shape])
# b_polynom = bounds_polynom([''])
# b_fractions = bounds_fractions(['fix_all'],popt[-num_of_specs*number_of_species:])
# p0, b0 = p0_and_bounds(p0_param, p0_polynom, p0_fractions, b_param, b_polynom, b_fractions)
# p0 = popt

# # for p, bl, bu in zip(p0, b0[0], b0[1]):
# #     print(' %s %s %s' %(str(bl), str(p), str(bu)))

# popt, pcov = curve_fit(spectra_over_time, lin_w, lin_i, p0, bounds = b0, method = 'trf', loss = 'soft_l1', diff_step=step_size_fit, x_scale='jac')

# write_specs(input_filenames, lin_w, set_wavenum, set_spectra, popt, 'Fit species')

# print("INFO: Fit species - R2 = %s and RMSD = %s" %(round(R2,3), round(RMSD, 10)))

######################
#Fit all: 
now    = datetime.now()
dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_str + ' INFO: Fit all...')

b_param = bounds_param([peak_fix,band_shape])
b_polynom = bounds_polynom([''])
b_fractions = bounds_fractions([''],popt[-num_of_specs*number_of_species:])
p0, b0 = p0_and_bounds(p0_param, p0_polynom, p0_fractions, b_param, b_polynom, b_fractions)
p0 = popt

popt, pcov = curve_fit(spectra_over_time, lin_w, lin_i, p0, bounds = b0, method = 'trf', loss = 'soft_l1', diff_step=step_size_fit, x_scale='jac')

write_specs(input_filenames, lin_w, set_wavenum, set_spectra, set_deriv, popt, 'Fit all')

out_params = unlinearized_params_and_round(popt[0:number_of_param])
print(out_params)
write_params(out_params,name_species)

write_conc(input_pref,[first_index, last_index],use_every_nth, popt, species_for_ratio)
now    = datetime.now()
dt_str = now.strftime("%d/%m/%Y %H:%M:%S")
print(dt_str + ' INFO: Done... Feierabend')

