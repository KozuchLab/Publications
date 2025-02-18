import numpy as np
import copy

'''This script averages spectra based on ranges or logarithmically'''


#################
####  Input  ####

# Input info on the spectra that are supposed to be fit:
input_pref = 'out_Spec'
output_pref = input_pref + '_avg'
total_num_of_specs = 502
# Logarythmic averaging? If yes, specify how many spectra you want as result:
log_avg = True
log_dp = 100
# Averaging according to ranges? Specify First, last, steps:
if log_avg != True:
    ranges = [[ 1, 101,  1],
              [121, 501, 20]]

if log_avg == True:
    ranges = []
    dp = 0
    org_log_dp = copy.deepcopy(log_dp)
    while dp <= org_log_dp:
        log_idx = []
        for i in range(1, total_num_of_specs+1):
            log_idx.append(round(log_dp*np.log(i)/np.log(total_num_of_specs) + 1))
        dp = len(set(log_idx)) 
        log_dp = log_dp + 1
        
    temp_log_idx = log_idx[0]
    first_idx = 1
    last_idx = 1
    for i in range(0, total_num_of_specs):
        temp_range = []
        if log_idx[i] > temp_log_idx:
            # print('#########')
            last_idx = i + 1            
            mid_idx = int(np.floor((last_idx + first_idx)/2))
            delta = int(last_idx - first_idx)
            # print(i)
            # print(first_idx)
            # print(last_idx)
            temp_range = [mid_idx, mid_idx, delta]
            first_idx = i + 1
            temp_log_idx = log_idx[i]
        if temp_range != []:
            ranges.append(temp_range)
            temp_range = []
            
# Example:
# ranges = [[start_1, end_1, delta_1],
#           [start_2, end_2, delta_2]]
#
# ranges = [[ 1,  21,  1],
#           [31, 501, 20]]
#
# This make a first range with spectra 
# - 1,2,3,4, ..., 21      - without any averaging, because floor(delta/2 = 0.5) = 0
# - 31,51,71, ..., 511    - with averaging of 10 before and 10 after, because floor(delta/2 = 10) = 10

########################
# Main:
    
    
def load_specs(input_pref, num_specs):

    in_set_wavenum = []
    in_set_spectra = []
       
    for a in range(1, num_specs+1):
        
        suffix = str(a).zfill(3)+'.txt'
        filename = input_pref + suffix
        x,y = np.loadtxt(filename).transpose()
        in_set_wavenum.append(x)
        in_set_spectra.append(y)
        
    return in_set_wavenum, in_set_spectra

def write_specs(output_pref, idx, waven, spectrum):

    suffix = str(idx).zfill(3)+'.txt'
    filename = output_pref + suffix
    np.savetxt(filename, np.array([waven,spectrum]).transpose())
        

set_waven, set_spectra = load_specs(input_pref, total_num_of_specs)

num_spec = total_num_of_specs
    
# creating array with info for averaging:
idxs_for_avg = []
for seq in ranges:
    

    print(seq)
    for idx in range(seq[0], seq[1]+seq[2], seq[2]):
        idxs_temp = []
        bottom = int(idx - np.floor(seq[2]/2))
        top = int(idx + np.floor(seq[2]/2))
        if bottom < 1:
            bottom = 1
        if top > num_spec:
            top = num_spec
        for i in range(bottom, top+1):
            idxs_temp.append(i-1)
        idxs_for_avg.append(idxs_temp)
    

with open(output_pref+'_LOG.txt', 'w') as outfile:
    counter = 1
    for idxs in idxs_for_avg:
        collected_specs = []
        collected_waven = []
        string_idx = ''
        for idx in idxs:
            string_idx = string_idx + ' ' + str(idx+1).rjust(4)
            collected_specs.append(set_spectra[idx])
            collected_waven.append(set_waven[idx])
        if len(collected_specs) == 1:
            avg_waven = collected_waven[0]
            std_waven = np.zeros(len(avg_waven))
            avg_spec  = collected_specs[0]
            std_spec  = np.zeros(len(avg_spec))
        elif len(collected_specs) > 1:
            avg_waven = np.mean(collected_waven, axis=0)
            std_waven = np.std(collected_waven, axis=0)
            avg_spec  = np.mean(collected_specs, axis=0)
            std_spec  = np.std(collected_specs, axis=0)
                       
        write_specs(output_pref, counter, avg_waven, avg_spec)
        outfile.write(str(counter).rjust(5) + str(int(np.mean(idxs)+1)).rjust(5) + '  ' + string_idx + '\n')
        
        counter = counter + 1
        
    
        