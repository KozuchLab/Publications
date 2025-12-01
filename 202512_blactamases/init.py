import json

def gen_filenames(input_xyz, input_key, suff):
    gen_prefix = input_xyz.replace('.xyz','')
    prefix = input_xyz.replace('.xyz',suff)
    
    if any(check in ['x','y','z'] for check in suff):
        new_input_xyz = input_xyz
    else:
        new_input_xyz = input_xyz.replace('.xyz',suff+'.xyz')
    
    filenames = {'prefix' : prefix,
                 'gen_prefix' : gen_prefix,
                 'input_xyz' : new_input_xyz,
                 'org_input_key' : input_key,
                 'input_key' : 'tink_'+input_key,
                 'gau_json' : 'gau'+suff+'.json',
                 'gau_init_com' : 'gau_'+prefix+'_init.com',
                 'gau_init_out' : 'gau_'+prefix+'_init.out',
                 'gau_init_chk' : 'gau_'+prefix+'_init.chk',
                 'gau_init_fchk': 'gau_'+prefix+'_init.fchk',
                 'gau_init_cart': 'gau_'+prefix+'_init.cart',
                 'gau_init_zmat': 'gau_'+prefix+'_init.zmat',
                 'gau_opt_com' : 'gau_'+prefix+'.com',
                 'gau_opt_out' : 'gau_'+prefix+'.out',
                 'gau_opt_chk' : 'gau_'+prefix+'.chk',
                 'gau_opt_fchk': 'gau_'+prefix+'.fchk',
                 'gau_opt_cart': 'gau_'+prefix+'.cart',
                 'gau_opt_zmat': 'gau_'+prefix+'.zmat',
                 'gau_opt_pdb': 'PDB_gau_'+prefix+'.pdb',
                 'gau_harm_com' : 'gau_'+prefix+'_harm.com',
                 'gau_harm_out' : 'gau_'+prefix+'_harm.out',
                 'gau_harm_chk' : 'gau_'+prefix+'_harm.chk',
                 'gau_harm_fchk': 'gau_'+prefix+'_harm.fchk',
                 'gau_harm_cart': 'gau_'+prefix+'_harm.cart',
                 'gau_harm_zmat': 'gau_'+prefix+'_harm.zmat',
                 'tink_json'    : 'tink'+suff+'.json',
                 'tink_min_xyz' : 'tink_'+prefix+'_min.xyz',
                 'tink_min_aligned_xyz' : 'tink_'+prefix+'_min_aligned.xyz',
                 'tink_min_out' : 'tink_'+prefix+'_min.out',
                 'tink_vib_xyz' : 'tink_'+prefix+'_vib.xyz',
                 'tink_vib_out' : 'tink_'+prefix+'_vib.out',
                 'tink_nma_xyz'   : 'tink_'+prefix+'_nma.xyz',
                 'tink_nrg_out'   : 'tink_'+prefix+'_nrg.out',
                 'tink_dip_out'   : 'tink_'+prefix+'_dip.out',
                 'tink_temp_key' : 'tink_'+gen_prefix+'_temp.key',
                 'tink_preopt_key' : 'tink_'+gen_prefix+'_preopt.key',
                 'tink_opt_key' : 'tink_'+gen_prefix+'_opt.key',
                 'tink_DONE_key' : 'DONE_'+gen_prefix+'.key'}
    
    with open('filenames'+suff+'.json', 'w') as file:
        json.dump(filenames, file, indent = 4, ensure_ascii = False)
    
    return filenames
