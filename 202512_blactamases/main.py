# ReValAMOEBA - Python script to refine valence parameters in AMOEBA against
#               quantum mechanical normal mode analysis. For help contact
#               Jacek Kozuch (jacek.kozuch@tu-braunschweig.de). If you use
#               this script or parts of it, please cite:
#               10.26434/chemrxiv-2025-9tlxj



import json
import numpy as np
from subprocess import run
import os
from scipy.optimize import curve_fit
import scipy.constants as sciconst
from scipy.spatial.transform import Rotation as R
import sys
import shutil
import time
import platform
from subprocess import run
import matplotlib.pyplot as plt
import itertools

from init import gen_filenames
from auxfuncs import number_of,\
                     openjson,\
                     NMA_watch,\
                     write_pdb,\
                     write_pdb_gau,\
                     xyz_linear,\
                     xyz_Nby3,\
                     int_coords,\
                     delete_annoying_files,\
                     NMA_plot
from qmandmm import generate_gau_init_files,\
                    generate_gau_opt_files,\
                    generate_gau_harm_files,\
                    generate_gau_cart_ef_files,\
                    submit_gau_init,\
                    submit_gau_opt,\
                    submit_gau_harm,\
                    qm_data,\
                    generate_tink_min_files,\
                    submit_tink_min,\
                    align_tink_xyz,\
                    mm_data,\
                    readkey,\
                    init_b0,\
                    k_params,\
                    writekey,\
                    add_restrainplane,\
                    b0_params,\
                    rewritekey
from manual_nma import manual_nma
from fittingfunctions import fit_MM_xyz, fit_MM_xyz_int, fit_MM_xyz_freq_dipder, fit_MM_xyz_freq


if __name__ == '__main__':
    
    '''Starting with relevant QM calculations'''
    with open('ReValAMOEBA.json','r') as file:
        ini = json.load(file)     
    input_xyz = ini['main_input_xyz']
    input_key = ini['input_key']
    input_suffs = ini['input_suffices']
    ef_suffs = ini['efield_suffices']
    scaling_factor = ini['scaling_factor']
        
    number_of_atoms, number_of_NMs = number_of(input_xyz)
    ini['number_of_atoms'] = number_of_atoms,
    ini['number_of_NMs'] = number_of_NMs,
    ini['displacements'] = 'cartesians'
    NM_low = 0

    try:
        run('rm zzz_*.*', shell=True)
    except:
        pass    
    if ini['redosteps']:
        try:
            run('rm tink*.*', shell=True)
        except:
            pass

    ini['filenamejsons'] = []
    for suff in input_suffs:
        gen_filenames(input_xyz, input_key, suff)
        ini['filenamejsons'].append('filenames'+suff+'.json')
    filenames = openjson(ini['filenamejsons'][0])  

    with open('ini.json', 'w') as file:
        json.dump(ini, file, indent = 4, ensure_ascii = False)
    if os.path.isfile(filenames['tink_temp_key']) == True:
        shutil.copy(filenames['tink_temp_key'],filenames['input_key'])
    elif os.path.isfile(filenames['tink_DONE_key']) == True:
        shutil.copy(filenames['tink_DONE_key'],filenames['input_key'])
    else:
        shutil.copy(filenames['org_input_key'],filenames['input_key'])
    shutil.copy(filenames['input_key'],filenames['tink_temp_key'])
                
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        generate_gau_init_files(filenames,ini)
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)        
        submit_gau_init(filenames)
        generate_gau_opt_files(filenames, ini)
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        submit_gau_opt(filenames)
        generate_gau_harm_files(filenames, ini)
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        submit_gau_harm(filenames)  

    ini['filenamejsons'] = []
    for suff in input_suffs+ef_suffs:
        gen_filenames(input_xyz, input_key, suff)
        ini['filenamejsons'].append('filenames'+suff+'.json')
    filenames = openjson(ini['filenamejsons'][0])  
    generate_gau_cart_ef_files(filenames, ef_suffs)
            
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        generate_gau_opt_files(filenames, ini)
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        submit_gau_opt(filenames)
        generate_gau_harm_files(filenames, ini)
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        submit_gau_harm(filenames)  
        
    '''Setting up everything for MM'''
    fit_min_frequency = ini["min_freq_nma"]
    NM_nums = []
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        qm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor)
        all_gau_data = openjson(filenames['gau_json'])
        i = 0
        for v in all_gau_data['freqs']:
            if v > fit_min_frequency:
                NM_low = i
                break
            i = i + 1
        NM_nums.append(NM_low)
            
    ini['lowest_normalmode'] = min(NM_nums)
    with open('ini.json', 'w') as file:
        json.dump(ini, file, indent = 4, ensure_ascii = False)
        
    '''Initialize everything from MM side'''
    suffix = '_00_pre'
    try:
        os.remove('zzz_NMA_watch.txt')
    except:
        pass
    
    for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
        filenames = openjson(filenamejson)
        generate_tink_min_files(filenames, number_of_atoms)
        submit_tink_min(filenames)
        align_tink_xyz(filenames,number_of_atoms)
        manual_nma(filenames)
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
        all_tink_data = openjson(filenames['tink_json'])
        all_gau_data = openjson(filenames['gau_json'])
        NMA_watch(all_gau_data, all_tink_data)
        shutil.copy(filenames['tink_json'], filenames['tink_json']+suffix)
        shutil.copy('zzz_NMA_watch.txt', 'zzz_NMA_watch.txt'+suffix)
        write_pdb_gau(filenames)
        shutil.copy(filenames['tink_min_xyz'].replace('.xyz','.pdb'), 'PDB_'+filenames['tink_min_xyz'].replace('_min',suffix).replace('.xyz','.pdb'))
        os.remove(filenames['tink_min_xyz'].replace('.xyz','.pdb'))
    
    # try:
    #     os.remove('zzz_NMA_watch.png')
    # except:
    #     pass    
    # NMA_plot(ini)
    # shutil.copy('zzz_NMA_watch.png', 'zzz_NMA_watch.png'+suffix)
    shutil.copy(filenames['tink_temp_key'], filenames['tink_temp_key']+suffix)
    add_restrainplane(filenames['tink_temp_key'],ini['add_restrains'])

    print("###############")
    ''''Run 1 - Optimize coords'''
    print('Run 1 - Optimize coords')
    if ini['adapt_structure'] == True:
        suffix = '_01_initxyz'
        if os.path.isfile(filenames['tink_json']+suffix) == False or ini['redosteps'] == True:
            all_popt = []
            for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
                filenames = openjson(filenamejson)
                all_gau_data = openjson(filenames['gau_json'])
                mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'tinker_vibrate')
                MM_params = readkey(filenames['tink_temp_key'])
                popt = init_b0(filenames,all_gau_data['coords'],number_of_atoms)
                popt = np.around(popt,4)
                all_popt.append(popt)
            popt = np.mean(all_popt,axis=0)
            os.remove('zzz_NMA_watch.txt')
            for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
                filenames = openjson(filenamejson)
                writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_xyz', 'minimize'])    
                submit_tink_min(filenames)
                align_tink_xyz(filenames,number_of_atoms)
                writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_xyz', 'vibrate']) 
                manual_nma(filenames)
                mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
                all_tink_data = openjson(filenames['tink_json'])
                NMA_watch(all_gau_data, all_tink_data)
                shutil.copy(filenames['tink_json'], filenames['tink_json']+suffix)
                shutil.copy('zzz_NMA_watch.txt', 'zzz_NMA_watch.txt'+suffix)
                write_pdb(filenames)
                shutil.copy(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'), 'PDB_'+filenames['tink_min_xyz'].replace('_min',suffix).replace('.xyz','.pdb'))
                os.remove(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'))
                shutil.copy(filenames['tink_temp_key'], filenames['tink_temp_key']+suffix)
            # try:
            #     os.remove('zzz_NMA_watch.png')
            # except:
            #     pass    
            # NMA_plot(ini)
            # shutil.copy('zzz_NMA_watch.png', 'zzz_NMA_watch.png'+suffix)
                
    if ini['fit_structure'] == True:
        suffix = '_02_optxyz'
        if os.path.isfile(filenames['tink_json']+suffix) == False or ini['redosteps'] == True:
            all_MM_p = []
            all_QM_data = []
            all_sigma = []
            for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
                filenames = openjson(filenamejson)
                all_gau_data = openjson(filenames['gau_json'])
                mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'tinker_vibrate')
                MM_params = readkey(filenames['tink_temp_key'])
                MM_p_1, bound0 = b0_params(MM_params)
                MM_p_2 = init_b0(filenames,all_gau_data['coords'],number_of_atoms)
                if len(MM_p_1) == len(MM_p_2):
                    MM_p = MM_p_2
                else:
                    MM_p = MM_p_2 + MM_p_1[len(MM_p_2):]
                all_MM_p.append(MM_p)
                if ini['fit_carts_and_ints'] == False:
                    QM_data = xyz_linear(all_gau_data['coords'])
                    fake_x_axis = np.arange(0, len(QM_data))  
                    sigma = []
                    
                elif ini['fit_carts_and_ints'] == True:
                    QM_data_1 = xyz_linear(all_gau_data['coords'])
                    QM_data_2 = int_coords(filenames,all_gau_data['coords'],number_of_atoms)
                    QM_data = QM_data_1 + QM_data_2
                    fake_x_axis = np.arange(0, len(QM_data))  
                    sigma_1 = [0.001]*len(QM_data_1)
                    sigma_2 = [0.1]*len(QM_data_2)
                    sigma = sigma_1 + sigma_2
                all_QM_data.append(QM_data)
                if sigma != []:
                    all_sigma.append(sigma)
                
            MM_p = np.mean(all_MM_p,axis=0)
            QM_data = xyz_linear(all_QM_data)
            if all_sigma != []:
                sigma = xyz_linear(all_sigma)
            if ini['fit_carts_and_ints'] == False:
                [popt, pcov] = curve_fit(fit_MM_xyz, fake_x_axis, QM_data, MM_p)
            elif ini['fit_carts_and_ints'] == True:
                [popt, pcov] = curve_fit(fit_MM_xyz_int, fake_x_axis, QM_data, MM_p, bounds= bound0, method='trf', diff_step= 0.0001)
                popt = np.around(popt,4)
    
            os.remove('zzz_NMA_watch.txt')
            for filenamejson in ini['filenamejsons'][0:len(ini["input_suffices"])]:
                filenames = openjson(filenamejson)
                all_gau_data = openjson(filenames['gau_json'])
                writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_xyz', 'minimize'])    
                submit_tink_min(filenames)
                align_tink_xyz(filenames,number_of_atoms)
                writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_xyz', 'vibrate']) 
                manual_nma(filenames)
                mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
                all_tink_data = openjson(filenames['tink_json'])
                NMA_watch(all_gau_data, all_tink_data)
                shutil.copy(filenames['tink_json'], filenames['tink_json']+suffix)
                shutil.copy('zzz_NMA_watch.txt', 'zzz_NMA_watch.txt'+suffix)
                write_pdb(filenames)
                shutil.copy(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'), 'PDB_'+filenames['tink_min_xyz'].replace('_min',suffix).replace('.xyz','.pdb'))
                os.remove(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'))
                shutil.copy(filenames['tink_temp_key'], filenames['tink_temp_key']+suffix)

            # try:
            #     os.remove('zzz_NMA_watch.png')
            # except:
            #     pass    
            # NMA_plot(ini)
            # shutil.copy('zzz_NMA_watch.png', 'zzz_NMA_watch.png'+suffix)

    '''Initialize everything from MM side'''
    suffix = '_03_withef'
    try:
        os.remove('zzz_NMA_watch.txt')
    except:
        pass
    
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        generate_tink_min_files(filenames, number_of_atoms)
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_xyz', 'minimize',filenames['tink_min_xyz']])    
        submit_tink_min(filenames)
        align_tink_xyz(filenames,number_of_atoms)
        writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_xyz', 'vibrate',filenames['tink_min_xyz']]) 
        manual_nma(filenames)
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
        all_tink_data = openjson(filenames['tink_json'])
        all_gau_data = openjson(filenames['gau_json'])
        NMA_watch(all_gau_data, all_tink_data)
        shutil.copy(filenames['tink_json'], filenames['tink_json']+suffix)
        shutil.copy('zzz_NMA_watch.txt', 'zzz_NMA_watch.txt'+suffix)
        write_pdb(filenames)
        shutil.copy(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'), 'PDB_'+filenames['tink_min_xyz'].replace('_min',suffix).replace('.xyz','.pdb'))
        os.remove(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'))
        shutil.copy(filenames['tink_temp_key'], filenames['tink_temp_key']+suffix)
    try:
        os.remove('zzz_NMA_watch.png')
    except:
        pass    
    NMA_plot(ini)
    shutil.copy('zzz_NMA_watch.png', 'zzz_NMA_watch.png'.replace('.png',suffix+'.png'))


    '''Based on how similar QM and MM structures are...'''
    RMSDs = []
    for filenamejson in ini['filenamejsons']:
        filenames = openjson(filenamejson)
        all_gau_data = openjson(filenames['gau_json'])
        all_tink_data = openjson(filenames['tink_json'])
        QM_data = xyz_linear(all_gau_data['coords'])
        MM_data = xyz_linear(all_tink_data['coords'])
        RMSD = np.sqrt(np.mean(np.subtract(QM_data,MM_data)**2))
        RMSDs.append(RMSD)
    RMSD = max(RMSDs) 
    if RMSD > ini['int_displ_threshold']:
        ini['displacements'] = 'internals'
        print('RSMD = %s between QM and MM structures. Pre-optimizing force constants in internal coordinates' %(str(round(RMSD,4))))
    else:
        ini['displacements'] = 'cartesians'
        print('RSMD = %s between QM and MM structures. Pre-optimizing force constants in cartesians'%(str(round(RMSD,4))))
    with open('ini.json', 'w') as file:
        json.dump(ini, file, indent = 4, ensure_ascii = False)

    '''Continue with fit'''
    if ini['fit_torsion_nma'] == False:
        rewritekey(filenames['tink_temp_key'], ['bond','angle'], ['strbnd','torsion']) 
        #rewritekey(filenames['tink_temp_key'], ['bond','angle','strbnd'], ['torsion']) 
        
    print('Run 3 - Optimize force constants using NMA')
    '''Run 3 - Optimize force constants using NMA'''
    rewritekey(filenames['tink_temp_key'], ['bond','angle'], ['strbnd','torsion']) 
    suffix = '_21_optbndang'
    if os.path.isfile(filenames['tink_json']+suffix) == False or ini['redosteps'] == True:
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'tinker_vibrate')
        MM_params = readkey(filenames['tink_temp_key'])
        MM_p, b0 = k_params(MM_params)
        all_MM_p = []
        all_QM_data = []
        all_sigma = []
        for filenamejson in ini['filenamejsons']:
            filenames = openjson(filenamejson)
            all_gau_data = openjson(filenames['gau_json'])            
            QM_data_1 = xyz_linear(all_gau_data['coords'])
            QM_data_2 = all_gau_data['freqs'][NM_low:]
            QM_data_3 = xyz_linear(all_gau_data['dipders'][NM_low:])
            QM_data = QM_data_1 + QM_data_2 + QM_data_3
            fake_x_axis = np.arange(0, len(QM_data))    
            sigma_1 = [0.01]*len(QM_data_1)
            sigma_2 = [1]*len(QM_data_2)
            sigma_3 = [1]*len(QM_data_3)
            for nm in ini['importantmodes']:
                if nm >= 0:
                    nm_idx = nm - number_of_NMs -1
                    sigma_2[nm_idx] = 0.01
                    
                    # nm_idx = 3*(nm_idx + 1) - 1
                    # sigma_3[nm_idx-2] = 0.1
                    # sigma_3[nm_idx-1] = 0.1
                    # sigma_3[nm_idx-0] = 0.1
                elif nm < 0:
                    nm_idx = nm
                    sigma_2[nm_idx] = 0.01
                    
                    # nm_idx = 3*(nm_idx + 1) - 1
                    # sigma_3[nm_idx-2] = 0.1
                    # sigma_3[nm_idx-1] = 0.1
                    # sigma_3[nm_idx-0] = 0.1
                    
            sigma = sigma_1 + sigma_2 + sigma_3
            all_QM_data.append(QM_data)
            all_sigma.append(sigma)
        QM_data = xyz_linear(all_QM_data)
        sigma = xyz_linear(all_sigma)        
        
        [popt, pcov] = curve_fit(fit_MM_xyz_freq_dipder, fake_x_axis, QM_data, p0=MM_p, bounds=b0, sigma=sigma, method='trf', diff_step= 0.01)
        popt = np.around(popt,4)
         
        os.remove('zzz_NMA_watch.txt')        
        for filenamejson in ini['filenamejsons']:
            filenames = openjson(filenamejson)
            generate_tink_min_files(filenames, number_of_atoms)
            writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_freq', 'minimize',filenames['tink_min_xyz']])    
            submit_tink_min(filenames)
            align_tink_xyz(filenames,number_of_atoms)
            writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_freq', 'vibrate',filenames['tink_min_xyz']]) 
            manual_nma(filenames)
            mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
            all_tink_data = openjson(filenames['tink_json'])
            all_gau_data = openjson(filenames['gau_json'])
            NMA_watch(all_gau_data, all_tink_data)
            shutil.copy(filenames['tink_json'], filenames['tink_json']+suffix)
            shutil.copy('zzz_NMA_watch.txt', 'zzz_NMA_watch.txt'+suffix)
            write_pdb(filenames)
            shutil.copy(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'), 'PDB_'+filenames['tink_min_xyz'].replace('_min',suffix).replace('.xyz','.pdb'))
            os.remove(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'))
            shutil.copy(filenames['tink_temp_key'], filenames['tink_temp_key']+suffix)
            
        try:
            os.remove('zzz_NMA_watch.png')
        except:
            pass    
        NMA_plot(ini)
        shutil.copy('zzz_NMA_watch.png', 'zzz_NMA_watch.png'.replace('.png',suffix+'.png'))
    
    if ini['fit_strbnd_nma']:
        rewritekey(filenames['tink_temp_key'], ['bond','angle','strbnd'], ['torsion']) 
        suffix = '_22_optbndangstrnd'
        if os.path.isfile(filenames['tink_json']+suffix) == False or ini['redosteps'] == True:
            mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'tinker_vibrate')
            MM_params = readkey(filenames['tink_temp_key'])
            MM_p, b0 = k_params(MM_params)
            all_MM_p = []
            all_QM_data = []
            all_sigma = []
            for filenamejson in ini['filenamejsons']:
                filenames = openjson(filenamejson)
                all_gau_data = openjson(filenames['gau_json'])            
                QM_data_1 = xyz_linear(all_gau_data['coords'])
                QM_data_2 = all_gau_data['freqs'][NM_low:]
                QM_data_3 = xyz_linear(all_gau_data['dipders'][NM_low:])
                QM_data = QM_data_1 + QM_data_2 + QM_data_3
                fake_x_axis = np.arange(0, len(QM_data))    
                sigma_1 = [0.01]*len(QM_data_1)
                sigma_2 = [1]*len(QM_data_2)
                sigma_3 = [1]*len(QM_data_3)
                for nm in ini['importantmodes']:
                    if nm >= 0:
                        nm_idx = nm - number_of_NMs -1
                        sigma_2[nm_idx] = 0.01
                        
                        # nm_idx = 3*(nm_idx + 1) - 1
                        # sigma_3[nm_idx-2] = 0.1
                        # sigma_3[nm_idx-1] = 0.1
                        # sigma_3[nm_idx-0] = 0.1
                    elif nm < 0:
                        nm_idx = nm
                        sigma_2[nm_idx] = 0.01
                        
                        # nm_idx = 3*(nm_idx + 1) - 1
                        # sigma_3[nm_idx-2] = 0.1
                        # sigma_3[nm_idx-1] = 0.1
                        # sigma_3[nm_idx-0] = 0.1
                sigma = sigma_1 + sigma_2 + sigma_3
                all_QM_data.append(QM_data)
                all_sigma.append(sigma)
            QM_data = xyz_linear(all_QM_data)
            sigma = xyz_linear(all_sigma)        
            
            [popt, pcov] = curve_fit(fit_MM_xyz_freq_dipder, fake_x_axis, QM_data, p0=MM_p, bounds=b0, sigma=sigma, method='trf', diff_step= 0.01)
            popt = np.around(popt,4)
             
            os.remove('zzz_NMA_watch.txt')        
            for filenamejson in ini['filenamejsons']:
                filenames = openjson(filenamejson)
                generate_tink_min_files(filenames, number_of_atoms)
                writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_freq', 'minimize',filenames['tink_min_xyz']])    
                submit_tink_min(filenames)
                align_tink_xyz(filenames,number_of_atoms)
                writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_freq', 'vibrate',filenames['tink_min_xyz']]) 
                manual_nma(filenames)
                mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
                all_tink_data = openjson(filenames['tink_json'])
                all_gau_data = openjson(filenames['gau_json'])
                NMA_watch(all_gau_data, all_tink_data)
                shutil.copy(filenames['tink_json'], filenames['tink_json']+suffix)
                shutil.copy('zzz_NMA_watch.txt', 'zzz_NMA_watch.txt'+suffix)
                write_pdb(filenames)
                shutil.copy(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'), 'PDB_'+filenames['tink_min_xyz'].replace('_min',suffix).replace('.xyz','.pdb'))
                os.remove(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'))
                shutil.copy(filenames['tink_temp_key'], filenames['tink_temp_key']+suffix)
                
            try:
                os.remove('zzz_NMA_watch.png')
            except:
                pass    
            NMA_plot(ini)
            shutil.copy('zzz_NMA_watch.png', 'zzz_NMA_watch.png'.replace('.png',suffix+'.png'))
            
    '''Run 4 - Optimize force constants using NMA'''
    suffix = '_31'
    if os.path.isfile(filenames['tink_json']+suffix) == False or ini['redosteps'] == True:
        mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'tinker_vibrate')
        MM_params = readkey(filenames['tink_temp_key'])
        MM_p, b0 = k_params(MM_params)
        all_MM_p = []
        all_QM_data = []
        all_sigma = []
        for filenamejson in ini['filenamejsons']:
            filenames = openjson(filenamejson)
            all_gau_data = openjson(filenames['gau_json'])            
            QM_data_1 = xyz_linear(all_gau_data['coords'])
            QM_data_2 = all_gau_data['freqs'][NM_low:]
            QM_data = QM_data_1 + QM_data_2
            fake_x_axis = np.arange(0, len(QM_data))    
            sigma_1 = [0.01]*len(QM_data_1)
            sigma_2 = [1]*len(QM_data_2)
            for nm in ini['importantmodes']:
                if nm >= 0:
                    nm_idx = nm - number_of_NMs -1
                    sigma_2[nm_idx] = 0.01
                    
                    # nm_idx = 3*(nm_idx + 1) - 1
                    # sigma_3[nm_idx-2] = 0.1
                    # sigma_3[nm_idx-1] = 0.1
                    # sigma_3[nm_idx-0] = 0.1
                elif nm < 0:
                    nm_idx = nm
                    sigma_2[nm_idx] = 0.01
                    
                    # nm_idx = 3*(nm_idx + 1) - 1
                    # sigma_3[nm_idx-2] = 0.1
                    # sigma_3[nm_idx-1] = 0.1
                    # sigma_3[nm_idx-0] = 0.1
                    
            sigma = sigma_1 + sigma_2
            all_QM_data.append(QM_data)
            all_sigma.append(sigma)
        QM_data = xyz_linear(all_QM_data)
        sigma = xyz_linear(all_sigma)        
        
        [popt, pcov] = curve_fit(fit_MM_xyz_freq, fake_x_axis, QM_data, p0=MM_p, bounds=b0, sigma=sigma, method='trf', diff_step= 0.01)
        popt = np.around(popt,4)
         
        os.remove('zzz_NMA_watch.txt')        
        for filenamejson in ini['filenamejsons']:
            filenames = openjson(filenamejson)
            generate_tink_min_files(filenames, number_of_atoms)
            writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_freq', 'minimize',filenames['tink_min_xyz']])    
            submit_tink_min(filenames)
            align_tink_xyz(filenames,number_of_atoms)
            writekey(filenames['input_key'], filenames['tink_temp_key'], MM_params, popt, ['fit_freq', 'vibrate',filenames['tink_min_xyz']]) 
            manual_nma(filenames)
            mm_data(filenames,number_of_atoms,number_of_NMs,scaling_factor,'manual_nma')
            all_tink_data = openjson(filenames['tink_json'])
            all_gau_data = openjson(filenames['gau_json'])
            NMA_watch(all_gau_data, all_tink_data)
            shutil.copy(filenames['tink_json'], filenames['tink_json']+suffix)
            shutil.copy('zzz_NMA_watch.txt', 'zzz_NMA_watch.txt'+suffix)
            write_pdb(filenames)
            shutil.copy(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'), 'PDB_'+filenames['tink_min_xyz'].replace('_min',suffix).replace('.xyz','.pdb'))
            os.remove(filenames['tink_min_aligned_xyz'].replace('.xyz','.pdb'))
            shutil.copy(filenames['tink_temp_key'], filenames['tink_temp_key']+suffix)
            
        try:
            os.remove('zzz_NMA_watch.png')
        except:
            pass    
        NMA_plot(ini)
        shutil.copy('zzz_NMA_watch.png', 'zzz_NMA_watch.png'.replace('.png',suffix+'.png'))

    shutil.copy(filenames['tink_temp_key'],filenames['tink_DONE_key'])
    delete_annoying_files(filenames['tink_dip_out'],number_of_NMs+7)
    delete_annoying_files(filenames['tink_nrg_out'],number_of_NMs+7)
    delete_annoying_files(filenames['tink_vib_xyz'],number_of_NMs+7)


