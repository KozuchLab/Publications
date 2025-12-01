from auxfuncs import openjson,xyz_linear,NMA_watch,int_coords,invert_vectors,assign_by_displ
from qmandmm import readkey,writekey,submit_tink_min,align_tink_xyz,mm_data
from manual_nma import manual_nma
import os

def fit_MM_xyz(fake_x_axis, *p):
    ini = openjson('ini.json')
    number_of_atoms = ini['number_of_atoms'][0]
    number_of_NMs = ini['number_of_NMs'][0]
    scaling_factor = ini['scaling_factor']
    all_MM_data = []
    for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])        

        MM_params_init = readkey(filenames['tink_temp_key'])
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params_init, p, ['fit_xyz', 'minimize'])    
        submit_tink_min(filenames)
        align_tink_xyz(filenames,number_of_atoms)
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
        all_tink_data = openjson(filenames['tink_json'])
        all_gau_data = openjson(filenames['gau_json'])
        MM_data = xyz_linear(all_tink_data['coords'])
        all_MM_data.append(MM_data)
    MM_data = xyz_linear(all_MM_data)    
                
    os.remove('zzz_NMA_watch.txt')
    for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])
        all_tink_data = openjson(filenames['tink_json'])
        NMA_watch(all_gau_data, all_tink_data)
    return MM_data

def fit_MM_xyz_int(fake_x_axis, *p):
    ini = openjson('ini.json')
    number_of_atoms = ini['number_of_atoms'][0]
    number_of_NMs = ini['number_of_NMs'][0]
    scaling_factor = ini['scaling_factor']
    all_MM_data = []
    for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])

        MM_params_init = readkey(filenames['tink_temp_key'])
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params_init, p, ['fit_xyz', 'minimize'])    
        submit_tink_min(filenames)
        align_tink_xyz(filenames,number_of_atoms)
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
        all_tink_data = openjson(filenames['tink_json'])
        MM_data_1 = xyz_linear(all_tink_data['coords'])
        MM_data_2 = int_coords(filenames,all_tink_data['coords'])
        MM_data = MM_data_1 + MM_data_2
        all_MM_data.append(MM_data)
        
    MM_data = xyz_linear(all_MM_data)  
    
    os.remove('zzz_NMA_watch.txt')
    for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])
        all_tink_data = openjson(filenames['tink_json'])
        NMA_watch(all_gau_data, all_tink_data)

    return MM_data

def fit_MM_xyz_freq_dipder(fake_x_axis, *p):

    ini = openjson('ini.json')
    number_of_atoms = ini['number_of_atoms'][0]
    number_of_NMs = ini['number_of_NMs'][0]
    scaling_factor = ini['scaling_factor']
    NM_low = ini['lowest_normalmode']
    all_MM_data = []
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])

        MM_params_init = readkey(filenames['tink_temp_key'])
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params_init, p, ['fit_freq', 'minimize', filenames['tink_min_xyz']])    
        submit_tink_min(filenames)
        align_tink_xyz(filenames,number_of_atoms)
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params_init, p, ['fit_freq', 'vibrate', filenames['tink_min_xyz']])   
        manual_nma(filenames)
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
        assign_by_displ(filenames)
        all_tink_data = openjson(filenames['tink_json'])   
        all_gau_data = openjson(filenames['gau_json'])   
        QM_data_3 = all_gau_data['dipders'][NM_low:]
        
        MM_data_1 = xyz_linear(all_tink_data['coords'])
        MM_data_2 = all_tink_data['freqs'][NM_low:]
        MM_data_3 = all_tink_data['dipders'][NM_low:]
        MM_data_3 = invert_vectors(QM_data_3,MM_data_3)
        MM_data_3 = xyz_linear(MM_data_3)
        MM_data = MM_data_1 + MM_data_2 + MM_data_3
        all_MM_data.append(MM_data)
        
    MM_data = xyz_linear(all_MM_data)  
    
    os.remove('zzz_NMA_watch.txt')
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])
        all_tink_data = openjson(filenames['tink_json'])
        NMA_watch(all_gau_data, all_tink_data)

    return MM_data

def fit_MM_xyz_freq(fake_x_axis, *p):

    ini = openjson('ini.json')
    number_of_atoms = ini['number_of_atoms'][0]
    number_of_NMs = ini['number_of_NMs'][0]
    scaling_factor = ini['scaling_factor']
    NM_low = ini['lowest_normalmode']
    all_MM_data = []
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])

        MM_params_init = readkey(filenames['tink_temp_key'])
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params_init, p, ['fit_freq', 'minimize', filenames['tink_min_xyz']])    
        submit_tink_min(filenames)
        align_tink_xyz(filenames,number_of_atoms)
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params_init, p, ['fit_freq', 'vibrate', filenames['tink_min_xyz']])   
        manual_nma(filenames)
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
        assign_by_displ(filenames)
        all_tink_data = openjson(filenames['tink_json'])   
        all_gau_data = openjson(filenames['gau_json'])   

        
        MM_data_1 = xyz_linear(all_tink_data['coords'])
        MM_data_2 = all_tink_data['freqs'][NM_low:]
        MM_data = MM_data_1 + MM_data_2
        all_MM_data.append(MM_data)
        
    MM_data = xyz_linear(all_MM_data)  
    
    os.remove('zzz_NMA_watch.txt')
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])
        all_tink_data = openjson(filenames['tink_json'])
        NMA_watch(all_gau_data, all_tink_data)

    return MM_data



