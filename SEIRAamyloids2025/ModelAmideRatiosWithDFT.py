import numpy as np
import os
from scipy.optimize import curve_fit

##################################
# INPUT:
filename = '20230503_WK3SH-6MH_POPC-POPG_Intensities_AmideBands.txt'
verbose = True
R2_threshold = 0.99
##################################

OD_to_mOD = 1000
data = np.loadtxt(filename, skiprows=1)


def Amide_Inten(theta, A, phi, frac):    
    '''The profiles of orientation dependent bands are modelled as described in
    Forbrig et al., Langmuir 2018, 34, 6, 2373–2385, https://doi.org/10.1021/acs.langmuir.7b04265
    as an angular dependent part and an offset that is not dependent of angle. 
    See def p_bSheet(theta) and def ap_bSheet(theta).'''

    Inten = A*((1-frac)*(np.cos((theta-phi)*np.pi/180))**2 + frac)
    
    return Inten


def p_bSheet(theta):
    '''Orientation dependent intensities of parallel b-sheets'''
    
    A_1690 = 2133.31429
    frac_1690 = 0
    phi_1690 = 97.81381
    A_1630 = 60257.17254
    frac_1630 = 0
    phi_1630 = 92.46796
    A_1530 = 40170.80263
    frac_1530 = 0.10926
    phi_1530 = -41.42413
    
    Inten_1690 = Amide_Inten(theta, A_1690, phi_1690, frac_1690)
    Inten_1630 = Amide_Inten(theta, A_1630, phi_1630, frac_1630)
    Inten_1530 = Amide_Inten(theta, A_1530, phi_1530, frac_1530)

    return np.array([         0, Inten_1630, Inten_1530])
    # return np.array([Inten_1690, Inten_1630, Inten_1530])

def ap_bSheet(theta):
    '''Orientation dependent intensities of antiparallel b-sheets'''
    
    A_1690 = 3723.4012
    frac_1690 = 0.01429
    phi_1690 = 13.7061
    A_1630 = 59823.22511
    frac_1630 = 0
    phi_1630 = 90.37562
    A_1530 = 20919.61169
    frac_1530 = 0.43037
    phi_1530 = -9.21653
    
    Inten_1690 = Amide_Inten(theta, A_1690, phi_1690, frac_1690)
    Inten_1630 = Amide_Inten(theta, A_1630, phi_1630, frac_1630)
    Inten_1530 = Amide_Inten(theta, A_1530, phi_1530, frac_1530)

    return np.array([Inten_1690, Inten_1630, Inten_1530])


def DFT_based_Amide_Intensities(theta_p, theta_ap, frac_ap):
    '''Calculating intensities for a given mode, at specified fractions of 
    parallel to antiparallel b-sheets and corresponding angles'''
    
    Intens_p  = (1-frac_ap)*  p_bSheet(theta_p)
    Intens_ap =    frac_ap * ap_bSheet(theta_ap)

    return np.array(Intens_p) + np.array(Intens_ap)

def fitted_expression(fake_x, theta_p, theta_ap, frac_ap, scaling):
    '''Function that is fit against experimental data. The scaling is applied as normalization 
    to experimental data, but is irrelevant (discarded)'''
    
    Intens = scaling * DFT_based_Amide_Intensities(theta_p, theta_ap, frac_ap)
    Intens = np.concatenate((Intens, Intens), axis=0)
    
    return Intens

def stardard_errors(popt,exp_amide):
    '''As bounds are used, errors can be calculated from curve_fit in every case. 
    That's a manual calculation of errors from a manually set up Jacobian.'''
    
    num_p = len(popt)
    res = fitted_expression(0, *popt)
    rel_dx = 0.0001
    dpopt = np.multiply(popt,rel_dx)
    jac = []
    for i in range(num_p):
                
        popt_u = np.zeros(num_p)
        popt_u[i] = popt_u[i] + dpopt[i]
        popt_u = popt + popt_u
        res_u = fitted_expression(0, *popt_u)
        dres_u = res - res_u
        dx_u = dpopt[i]
        dres_dx_u = np.divide(dres_u,dx_u)
        
        popt_d = np.zeros(num_p)
        popt_d[i] = popt_d[i] - dpopt[i]
        popt_d = popt + popt_d
        res_d = fitted_expression(0, *popt_d)
        dres_d = res - res_d
        dx_d = - dpopt[i]
        dres_dx_d = np.divide(dres_d,dx_d)
        
        dres_dx = np.add(dres_dx_u,dres_dx_d)/2      
        jac.append(dres_dx.tolist())
           
    JtxJ = np.dot(np.array(jac),np.array(jac).transpose())
    invJtxJ = np.linalg.pinv(JtxJ)
    res2 = (exp_amide - np.mean(exp_amide))**2
    sigma = np.sqrt(np.sum(res2)/2/len(exp_amide))
    
    pcov = sigma * invJtxJ
    serr = np.diag(pcov)
    
    return serr

def model_DFT_spectrum(theta_p, theta_ap, frac_ap, scaling):
    '''Modelled spectra are printed.'''
    
    mainpath = os.getcwd()
    p_path = os.path.join(mainpath, 'Spectra_Parallel_DFT')
    ap_path = os.path.join(mainpath, 'Spectra_AntiParallel_DFT')
       
    ang_p_u = 0
    ang_p_d = 0
    ang_ap_u = 0
    ang_ap_d = 0
    for ang in range(0,185,5):
        if theta_p < float(ang):
            ang_p_d = ang - 5
            ang_p_u = ang
            ang_p = [ang_p_d, ang_p_u]
            p_weights_d = abs(theta_p-ang_p_d)/5
            p_weights_u = abs(theta_p-ang_p_u)/5
            p_weight = [p_weights_d, p_weights_u]
            break 
        
    for ang in range(0,185,5):  
        if theta_ap < float(ang):
            ang_ap_d = ang - 5
            ang_ap_u = ang
            ang_ap = [ang_ap_d, ang_ap_u]
            ap_weights_d = abs(theta_ap-ang_ap_d)/5
            ap_weights_u = abs(theta_ap-ang_ap_u)/5
            ap_weight = [ap_weights_d, ap_weights_u]
            break
                
    os.chdir(p_path)
    allfiles = os.listdir()
    d_spec = []
    u_spec = []
    for file in allfiles:
        if '_'+str(ang_p[0])+'.' in file or '_0'+str(ang_p[0])+'.' in file:
            # print(file)
            d_spec = np.loadtxt(file).transpose()
        if '_'+str(ang_p[1])+'.' in file or '_0'+str(ang_p[1])+'.' in file:
            # print(file)
            u_spec = np.loadtxt(file).transpose()
        
    p_spec = p_weight[0]*d_spec + p_weight[1]*u_spec
        
    os.chdir(ap_path)
    allfiles = os.listdir()
    d_spec = []
    u_spec = []
    for file in allfiles:
        if '_'+str(ang_ap[0])+'.' in file or '_0'+str(ang_ap[0])+'.' in file:
            d_spec = np.loadtxt(file).transpose()
        if '_'+str(ang_ap[1])+'.' in file or '_0'+str(ang_ap[1])+'.' in file:
            u_spec = np.loadtxt(file).transpose()
            
    ap_spec = ap_weight[0]*d_spec + ap_weight[1]*u_spec    
        
    os.chdir(mainpath)
    p_spec[1]  = scaling*(1-frac_ap)*p_spec[1]
    ap_spec[1] = scaling*frac_ap*ap_spec[1]
 
    out_spec = p_spec[1] + ap_spec[1]
    
    waven = np.mean([p_spec[0], ap_spec[0]],axis=0)
    std_waven = np.std([p_spec[0], ap_spec[0]],axis=0)
    
    if all([x == 0 for x in np.round(std_waven,3)]):
        #print('All good! - similar wavenumber axes detected')
        pass
    else:
        print('WARNING - Wavenumber axes are different!')
        
    return waven, out_spec, p_spec[1], ap_spec[1]
    
def R2_RMSD(exp, fit):
    
    resid = np.subtract(exp,fit)
    RMSD  = np.sqrt(np.mean(resid**2,axis=0))
    
    SumSquaredRegression = np.sqrt(np.sum(resid**2))
    SumSquaredTotal      = np.sqrt(np.sum(exp**2))
    
    R2 = 1 - SumSquaredRegression/SumSquaredTotal
    
    # print(exp)
    # print(fit)
    # print(resid)
    # print(RMSD)
    # print(SumSquaredRegression)
    # print(SumSquaredTotal)
    # print(R2)

    return R2, RMSD

######################
# Main:

with open(filename+'_OUTPUT.txt', 'w') as file:
    file.write('Theta_para      Theta_antipara  frac_antipara   scaling                 1690            1630            1530 \n')
    print('Initialized file')

counter = 0

'''Here, the main part starts. Importantly, due to the periodicity of 180 deg
the fit will not be able to cross from 0 to 180 easily. Therefore, four separate fits are performed
starting from angles of 
- 45 and 45
- 45 and 135
- 135 and 135
- 135 and 45.
The results are projected back to a result in the rango of 0 and 90 orientation 
(we do not distinguish left and right-leaning solutions, as we couldn't from spectra)
are averaged if the fits converged above a certain R2. Otherwise a solution is discarded.'''

for amide in data:
    exp_amide = np.array(amide)
    exp_amide = exp_amide*OD_to_mOD
    exp_amide = np.concatenate((exp_amide, exp_amide), axis=0)
    
    # All four combination of left/left, left/right, etx
    
    popts = []
    ress  = []
    waven = []
    specs = []
    p_specs = []
    ap_specs = []
    info = ''
    
    # FIRST:
    p0 = [45, 45, 0.5, 0.0001]
    b_down = [0, 0, 0, 0]
    b_up = [180, 180, 1, np.inf]
    b0 = [b_down, b_up]
    popt, pcov = curve_fit(fitted_expression, np.zeros(len(exp_amide)), exp_amide, p0, bounds = b0, method = 'trf')
    res  = fitted_expression(np.zeros(len(exp_amide)),*popt)
    serr = stardard_errors(popt,exp_amide)
    
    R2, RMSD = R2_RMSD(exp_amide, res)
    
    if R2 >= R2_threshold:
        info = info + ' ' + '  45/45'
        popts.append(popt)
        ress.append(res)
        
        temp_specs = model_DFT_spectrum(*popt)
        waven.append(temp_specs[0])
        specs.append(temp_specs[1])
        p_specs.append(temp_specs[2])
        ap_specs.append(temp_specs[3])
    else:
        info = info + ' ' + '-------'

    if verbose == True:
        print('########################################################')
        print('NEXT TIMEPOINT:')
        print('########################################################')
        print('Results - 45/45:')
        
        print('Exp = ' + str(np.round(exp_amide[0:3],2)))
        print('Fit = ' + str(np.round(res[0:3],2)))
        print('R2, RMSD = ' + str(round(R2,3)) + ' ' + str(round(RMSD,3)))
        
        print('Theta_p  = ' + str(round(popt[0],1)) + ' +/- ' + str(round(serr[0],1)))
        print('Theta_ap = ' + str(round(popt[1],1)) + ' +/- ' + str(round(serr[1],1)))
        print('frac_ap  = ' + str(round(popt[2],2)) + ' +/- ' + str(round(serr[2],2)))
        print('scaling  = ' + str(round(popt[3],6)) + ' +/- ' + str(round(serr[3],6)))
        print('########')
    
    # SECOND:
    p0 = [45, 135, popt[-2], popt[-1]]
    b_down = [0, 0, 0, 0]
    b_up = [180, 180, 1, np.inf]
    b0 = [b_down, b_up]
    popt, pcov = curve_fit(fitted_expression, np.zeros(len(exp_amide)), exp_amide, p0, bounds = b0, method = 'trf')
    res  = fitted_expression(np.zeros(len(exp_amide)),*popt)
    serr = stardard_errors(popt,exp_amide)
    
    R2, RMSD = R2_RMSD(exp_amide, res)
    
    if R2 >= R2_threshold:
        info = info + ' ' + ' 45/135'
        popts.append(popt)
        ress.append(res)
        
        temp_specs = model_DFT_spectrum(*popt)
        waven.append(temp_specs[0])
        specs.append(temp_specs[1])
        p_specs.append(temp_specs[2])
        ap_specs.append(temp_specs[3])
    else:
        info = info + ' ' + '-------'
        
    if verbose == True:
        print('Results - 45/135:')
        
        print('Exp = ' + str(np.round(exp_amide[0:3],2)))
        print('Fit = ' + str(np.round(res[0:3],2)))
        print('R2, RMSD = ' + str(round(R2,3)) + ' ' + str(round(RMSD,3)))
        
        print('Theta_p  = ' + str(round(popt[0],1)) + ' +/- ' + str(round(serr[0],1)))
        print('Theta_ap = ' + str(round(popt[1],1)) + ' +/- ' + str(round(serr[1],1)))
        print('frac_ap  = ' + str(round(popt[2],2)) + ' +/- ' + str(round(serr[2],2)))
        print('scaling  = ' + str(round(popt[3],6)) + ' +/- ' + str(round(serr[3],6)))
        print('########')
    
    # THIRD:
    p0 = [135, 135, popt[-2], popt[-1]]
    b_down = [0, 0, 0, 0]
    b_up = [180, 180, 1, np.inf]
    b0 = [b_down, b_up]
    popt, pcov = curve_fit(fitted_expression, np.zeros(len(exp_amide)), exp_amide, p0, bounds = b0, method = 'trf')
    res  = fitted_expression(np.zeros(len(exp_amide)),*popt)
    serr = stardard_errors(popt,exp_amide)
    
    R2, RMSD = R2_RMSD(exp_amide, res)
    
    if R2 >= R2_threshold:
        info = info + ' ' + '135/135'
        popts.append(popt)
        ress.append(res)
        
        temp_specs = model_DFT_spectrum(*popt)
        waven.append(temp_specs[0])
        specs.append(temp_specs[1])
        p_specs.append(temp_specs[2])
        ap_specs.append(temp_specs[3])
    else:
        info = info + ' ' + '-------'
    
    if verbose == True:
        print('Results - 135/135:')
        
        print('Exp = ' + str(np.round(exp_amide[0:3],2)))
        print('Fit = ' + str(np.round(res[0:3],2)))
        print('R2, RMSD = ' + str(round(R2,3)) + ' ' + str(round(RMSD,3)))
        
        print('Theta_p  = ' + str(round(popt[0],1)) + ' +/- ' + str(round(serr[0],1)))
        print('Theta_ap = ' + str(round(popt[1],1)) + ' +/- ' + str(round(serr[1],1)))
        print('frac_ap  = ' + str(round(popt[2],2)) + ' +/- ' + str(round(serr[2],2)))
        print('scaling  = ' + str(round(popt[3],6)) + ' +/- ' + str(round(serr[3],6)))
        print('########')
    
    # FORTH:
    p0 = [135, 45, popt[-2], popt[-1]]
    b_down = [0, 0, 0, 0]
    b_up = [180, 180, 1, np.inf]
    b0 = [b_down, b_up]
    popt, pcov = curve_fit(fitted_expression, np.zeros(len(exp_amide)), exp_amide, p0, bounds = b0, method = 'trf')
    res  = fitted_expression(np.zeros(len(exp_amide)),*popt)
    serr = stardard_errors(popt,exp_amide)
    
    R2, RMSD = R2_RMSD(exp_amide, res)
    
    if R2 >= R2_threshold:
        info = info + ' ' + ' 135/45'
        popts.append(popt)
        ress.append(res)
        
        temp_specs = model_DFT_spectrum(*popt)
        waven.append(temp_specs[0])
        specs.append(temp_specs[1])
        p_specs.append(temp_specs[2])
        ap_specs.append(temp_specs[3])
    else:
        info = info + ' ' + '-------'
    
    if verbose == True:
        print('Results - 135/45:')
        
        print('Exp = ' + str(np.round(exp_amide[0:3],2)))
        print('Fit = ' + str(np.round(res[0:3],2)))
        print('R2, RMSD = ' + str(round(R2,3)) + ' ' + str(round(RMSD,3)))
        
        print('Theta_p  = ' + str(round(popt[0],1)) + ' +/- ' + str(round(serr[0],1)))
        print('Theta_ap = ' + str(round(popt[1],1)) + ' +/- ' + str(round(serr[1],1)))
        print('frac_ap  = ' + str(round(popt[2],2)) + ' +/- ' + str(round(serr[2],2)))
        print('scaling  = ' + str(round(popt[3],6)) + ' +/- ' + str(round(serr[3],6)))
        print('########')
    
    new_popts = []
    for i in range(len(popts)):
        new_popts_temp = []
        for j in range(len(popts[0])):
            if j < 2:
                if popts[i][j] > 90:
                    new_popts_temp.append(180 - popts[i][j])
                else:
                    new_popts_temp.append(popts[i][j])
            else:
                new_popts_temp.append(popts[i][j])
        new_popts.append(new_popts_temp)
        
    avg_popt = np.mean(new_popts,axis=0)
    std_popt = np.std(new_popts,axis=0)
    
    if verbose == True:
        print('Theta_p  = ' + str(round(avg_popt[0],1)) + ' +/- ' + str(round(std_popt[0],1)))
        print('Theta_ap = ' + str(round(avg_popt[1],1)) + ' +/- ' + str(round(std_popt[1],1)))
        print('frac_ap  = ' + str(round(avg_popt[2],2)) + ' +/- ' + str(round(std_popt[2],2)))
        print('scaling  = ' + str(round(avg_popt[3],4)) + ' +/- ' + str(round(std_popt[3],4)))  
    
    avg_ress = np.mean(ress,axis=0)
    std_ress = np.std(ress,axis=0)
    
    if verbose == True:
        print('1690 = ' + str(round(avg_ress[0],1)) + ' +/- ' + str(round(std_ress[0],1)))
        print('1630 = ' + str(round(avg_ress[1],1)) + ' +/- ' + str(round(std_ress[1],1)))
        print('1530 = ' + str(round(avg_ress[2],2)) + ' +/- ' + str(round(std_ress[2],2)))
    
    with open(filename+'_OUTPUT.txt', 'a') as file:

        string_1 = str(round(avg_popt[0],1)) + '  ' + str(round(std_popt[0],1))
        string_2 = str(round(avg_popt[1],1)) + '  ' + str(round(std_popt[1],1))
        string_3 = str(round(avg_popt[2],2)) + '  ' + str(round(std_popt[2],2))
        string_4 = str(round(avg_popt[3],6)) + '  ' + str(round(std_popt[3],6))
        
        string_5 = str(round(avg_ress[0],2)) + '  ' + str(round(std_ress[0],2))
        string_6 = str(round(avg_ress[1],2)) + '  ' + str(round(std_ress[1],2))
        string_7 = str(round(avg_ress[2],2)) + '  ' + str(round(std_ress[2],2))
        
        file.write(string_1.ljust(15) + ' ' + string_2.ljust(15) + ' ' + string_3.ljust(15) + ' ' + string_4.ljust(20) + '    '+\
                   string_5.ljust(15) + ' ' + string_6.ljust(15) + ' ' + string_7.ljust(15) + '    '+\
                   info + '\n')
    
    avg_waven = np.mean(waven,axis=0)
    std_waven = np.std(waven,axis=0)
    
    avg_spec = np.mean(specs,axis=0)
    std_spec = np.std(specs,axis=0)

    avg_p_spec = np.mean(p_specs,axis=0)
    std_p_spec = np.std(p_specs,axis=0)

    avg_ap_spec = np.mean(ap_specs,axis=0)
    std_ap_spec = np.std(ap_specs,axis=0)    
    
    if counter < 10:
        label = '00'+str(counter)
    elif counter < 100:
        label = '0'+str(counter)
    else:
        label = str(counter)
    
    np.savetxt('zzz_spec_'+label+'.txt', np.array([avg_waven,avg_spec,std_spec]).transpose())
    np.savetxt('zzz_spec_noSTD_'+label+'.txt', np.array([avg_waven,avg_spec]).transpose())

    np.savetxt('zzz_P_spec_'+label+'.txt', np.array([avg_waven,avg_p_spec,std_p_spec]).transpose())
    np.savetxt('zzz_P_spec_noSTD_'+label+'.txt', np.array([avg_waven,avg_p_spec]).transpose())

    np.savetxt('zzz_AP_spec_'+label+'.txt', np.array([avg_waven,avg_ap_spec,std_ap_spec]).transpose())
    np.savetxt('zzz_AP_spec_noSTD_'+label+'.txt', np.array([avg_waven,avg_ap_spec]).transpose())    
    
    counter = counter + 1
