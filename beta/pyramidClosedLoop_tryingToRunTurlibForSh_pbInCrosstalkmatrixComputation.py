#%% -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 13:42:39 2023

@author: flero
"""

import pickle

import time
import matplotlib.pyplot as plt
import numpy as np

from turlib import zernike_variance
from turlib import iterative_estimator

from OOPAO.Telescope import Telescope
from OOPAO.Source import Source
from OOPAO.Atmosphere import Atmosphere
from OOPAO.DeformableMirror import DeformableMirror
from OOPAO.Pyramid import Pyramid
from OOPAO.ShackHartmann import ShackHartmann
from OOPAO.Zernike import Zernike
from OOPAO.calibration.InteractionMatrix import InteractionMatrix
from OOPAO.calibration.InteractionMatrix import InteractionMatrixFromPhaseScreen

#%%%%%%%%%%%%%%%%%%%%%%%%%%%% Choose and import parameter file %%%%%%%%%%%%%%%%%

from parameterFilePapyrus_pyramidClosedLoop import initializeParameterFiles
#from parameterFile_pyramidClosedLoop import initializeParameterFiles

param = initializeParameterFiles()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TELESCOPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tel = Telescope(resolution=param['resolution'],
                diameter=param['diameter'],
                samplingTime=1/param['frequency'])

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

src = Source(optBand=param['opticalBand'],
             magnitude=param['magnitude'])
src*tel

#%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATMOSPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

atm = Atmosphere(telescope=tel,
                 r0=param['r0'],
                 L0=param['L0'],
                 windSpeed=param['windSpeed'],
                 fractionalR0=param['fractionnalR0'],
                 windDirection=param['windDirection'],
                 altitude=param['altitude'])

atm.initializeAtmosphere(tel)

atm.generateNewPhaseScreen(seed = 5) # get a new atmosphere

tel+atm

tel.computePSF(zeroPaddingFactor=param['zeroPaddingFactor'])

#%%%%%%%%%%%%%%%%%%%%% DM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dm = DeformableMirror(telescope=tel,
                      nSubap=param['nSubapDm'],
                      mechCoupling=0.35)

#%%%%%%%%%%%%%%%%%%%%%% PYWFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tel - atm

pywfs = Pyramid(nSubap=param['nSubap'],
                telescope=tel,
                modulation=param['modulationRadius'],
                lightRatio=param['lightRatio'])

#%%%%%%%%%%%%%%%%%%%%%% SH geometric %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shGeo = ShackHartmann(param['nSubap'], tel, param['lightRatio'], is_geometric=True)

#%%%%%%%%%%%%%%%%%%%%%% SH diffractive %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

shDiff = ShackHartmann(param['nSubap'], tel, param['lightRatio'], is_geometric=False)

#%%%%%%%%%%%% Zernike Modal Basis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# create Zernike Object
zernike = Zernike(tel,param['nModes'])
# compute polynomials for given telescope
zernike.computeZernike(tel)

# mode to command matrix to project Zernike Polynomials on DM
M2C_zernike = np.linalg.pinv(np.squeeze(dm.modes[tel.pupilLogical,:]))@zernike.modes

# compute Zernike variances on dm [rad^2]
varDm = np.zeros(param['nModes'])
for k in range(param['nModes']):
    dm.coefs = M2C_zernike[:,k]
    varDm[k] = np.var(dm.OPD*2*np.pi*src.wavelength)

#%%%%%%%%%%%%%%%%%%%%%% Calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if param['calibrate']:
    
    a = time.time()

    # for the loop
    calib_zernike_pywfs = InteractionMatrix(ngs=src,
                                      atm=atm,
                                      tel=tel,
                                      dm=dm,
                                      wfs=pywfs,
                                      M2C=M2C_zernike,
                                      stroke=param['stroke'],
                                      nMeasurements  = 1,
                                      noise='off')
    
    calib_zernike_shGeo = InteractionMatrix(ngs=src,
                                      atm=atm,
                                      tel=tel,
                                      dm=dm,
                                      wfs=shGeo,
                                      M2C=M2C_zernike,
                                      stroke=param['stroke'],
                                      nMeasurements  = 1,
                                      noise='off')
    
    calib_zernike_shDiff = InteractionMatrix(ngs=src,
                                      atm=atm,
                                      tel=tel,
                                      dm=dm,
                                      wfs=shDiff,
                                      M2C=M2C_zernike,
                                      stroke=param['stroke'],
                                      nMeasurements  = 1,
                                      noise='off')
    
    # for turlib
    zernikePlus = Zernike(tel,
                          J=param['nModesShown'])
    zernikePlus.computeZernike(tel)
   
#%%

    calib_zernike_turlib_shGeo = InteractionMatrixFromPhaseScreen(ngs=src,
                                          atm=atm,
                                          tel=tel,
                                          wfs=shGeo,
                                          phasScreens=zernikePlus.modesFullRes,
                                          stroke=param['stroke'])
    
#%%
 
    calib_zernike_turlib_shDiff = InteractionMatrixFromPhaseScreen(ngs=src,
                                          atm=atm,
                                          tel=tel,
                                          wfs=shDiff,
                                          phasScreens=zernikePlus.modesFullRes,
                                          stroke=param['stroke'])

#%%    
    calib_zernike_turlib_pywfs = InteractionMatrixFromPhaseScreen(ngs=src,
                                          atm=atm,
                                          tel=tel,
                                          wfs=pywfs,
                                          phasScreens=zernikePlus.modesFullRes,
                                          stroke=param['stroke'])
    
#%%    

    # save calib results
    with open(param['pathDataCalibration'] + 'calib_zernike_shGeo.pkl', 'wb') as outp:
        pickle.dump(calib_zernike_shGeo, outp, pickle.HIGHEST_PROTOCOL)
    
    with open(param['pathDataCalibration'] + 'calib_zernike_shDiff.pkl', 'wb') as outp:
        pickle.dump(calib_zernike_shDiff, outp, pickle.HIGHEST_PROTOCOL)
    
    with open(param['pathDataCalibration'] + 'calib_zernike_pywfs.pkl', 'wb') as outp:
        pickle.dump(calib_zernike_pywfs, outp, pickle.HIGHEST_PROTOCOL)
        
    with open(param['pathDataCalibration'] + 'calib_zernike_turlib_shGeo.pkl', 'wb') as outp:
        pickle.dump(calib_zernike_turlib_shGeo, outp, pickle.HIGHEST_PROTOCOL)
        
    with open(param['pathDataCalibration'] + 'calib_zernike_turlib_shDiff.pkl', 'wb') as outp:
        pickle.dump(calib_zernike_turlib_shDiff, outp, pickle.HIGHEST_PROTOCOL)
        
    with open(param['pathDataCalibration'] + 'calib_zernike_turlib_pywfs.pkl', 'wb') as outp:
        pickle.dump(calib_zernike_turlib_pywfs, outp, pickle.HIGHEST_PROTOCOL)

    # print time
    b = time.time()
    print('The calibration took ' + str(round(b-a)) + ' seconds.\n')

else: 
    
    #load previous calibration results
    with open(param['pathDataCalibration'] + 'calib_zernike_shGeo.pkl', 'rb') as inp:
        saved_calib_zernike_shGeo = pickle.load(inp)
        
    with open(param['pathDataCalibration'] + 'calib_zernike_shDiff.pkl', 'rb') as inp:
        saved_calib_zernike_shDiff = pickle.load(inp)
        
    with open(param['pathDataCalibration'] + 'calib_zernike_pywfs.pkl', 'rb') as inp:
        saved_calib_zernike_pywfs = pickle.load(inp)
    
    with open(param['pathDataCalibration'] + 'calib_zernike_turlib_shGeo.pkl', 'rb') as inp:
        saved_calib_zernike_turlib_shGeo = pickle.load(inp)
        
    with open(param['pathDataCalibration'] + 'calib_zernike_turlib_shDiff.pkl', 'rb') as inp:
        saved_calib_zernike_turlib_shDiff = pickle.load(inp)
        
    with open(param['pathDataCalibration'] + 'calib_zernike_turlib_pywfs.pkl', 'rb') as inp:
        saved_calib_zernike_turlib_pywfs = pickle.load(inp)

#%% compute cross-talk matrice (turlib input)

g_pywfs = calib_zernike_turlib_pywfs.D
gF_pywfs = g_pywfs[:,:param['nModes']] # modes2slopes for reconstructed modes
gR_pywfs = g_pywfs[:,param['nModes']:] # modes2slopes for higher order modes
gFplus_pywfs = np.linalg.pinv(np.transpose(gF_pywfs) @ gF_pywfs) @ np.transpose(gF_pywfs)

crosstalkMatrix_pywfs = gFplus_pywfs @ gR_pywfs

g_shGeo = calib_zernike_turlib_shGeo.D
gF_shGeo = g_shGeo[:,:param['nModes']] # modes2slopes for reconstructed modes
gR_shGeo = g_shGeo[:,param['nModes']:] # modes2slopes for higher order modes
gFplus_shGeo = np.linalg.pinv(np.transpose(gF_shGeo) @ gF_shGeo) @ np.transpose(gF_shGeo)

crosstalkMatrix_shGeo = gFplus_shGeo @ gR_shGeo

#%%%%%%%%%%%%%%%%%%%%%% Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# select and compute reconstructor

if param['wfs'] == 'shGeo':
    reconstructor = M2C_zernike @ calib_zernike_shGeo.M
elif param['wfs'] == 'shDiff':
    reconstructor = M2C_zernike @ calib_zernike_shDiff.M
elif param['wfs'] == 'pywfs':
    reconstructor = M2C_zernike @ calib_zernike_pywfs.M # choose calibration result from the one made only with controlled modes
    #reconstructor=M2C_zernike @ g.M[:param['nModes'],:] # choose calibration result from the one made for atmosphere reconstruction
else:
    raise ValueError('param[\'wfs\'] as been wrongly set')

# combine telescope with atmosphere
tel.resetOPD()
seed = int(np.random.rand()*1000)
atm.generateNewPhaseScreen(seed = seed) # get a new atmosphere each time
tel+atm

# initialize loop variables
dm.coefs=0
wfsSignal = np.zeros(pywfs.signal.shape)

# allow memory for saving loop data
bPywfs = np.zeros((param['nModes'],param['nIter']), dtype=float)
bShGeo = np.zeros((param['nModes'],param['nIter']), dtype=float)
bShDiff = np.zeros((param['nModes'],param['nIter']), dtype=float)
cDm = np.zeros((param['nModes'],param['nIter']), dtype=float)
phaseVar = np.zeros(param['nIter'], dtype=float)
resPhaseVar = np.zeros(param['nIter'], dtype=float)
psfLongExp = np.zeros(tel.PSF.shape, dtype=float)
strehl = np.zeros(param['nIter'], dtype=float)

if param['recordTurbulence'] == True:
    turbulence = np.zeros((np.sum(tel.pupil), param['nIter']))
    print('\nTurbulent phase screens will be recorded \n')

a=time.time()

for i in range(param['nIter']):
    
    # update turbulent phase screens
    atm.update()
    
    # save turbulence var
    phaseVar[i] = np.var(atm.OPD[np.where(tel.pupil==1)]*2*np.pi/src.wavelength) # var of the turbulent phase [rad]
    
    # save turbulent phase screen (optional)
    if param['recordTurbulence'] == True:
        turbulence[:,i] = atm.OPD[np.where(tel.pupil == True)]*2*np.pi/src.wavelength # [rad]
    
    # optical propagation involving dm shape of last iteration : 1 frame delay    
    src*tel*dm*shGeo*shDiff*pywfs
        
    # update wfs signal
    if param['wfs'] == 'pywfs':
        wfsSignal=pywfs.signal
    elif param['wfs'] == 'shGeo':
        wfsSignal=shGeo.signal
    else:
        raise ValueError('param[\'wfs\'] as been wrongly set')
    
    # update dm commands
    dm.coefs=dm.coefs-param['gain']*np.matmul(reconstructor,wfsSignal)
    
    # record resPhaseVar
    resPhaseVar[i] = np.var(tel.OPD[np.where(tel.pupil==1)]*2*np.pi/src.wavelength) # var of the residual phase [rad]
    
    # record strehl ratio
    strehl[i] = np.exp(-resPhaseVar[i])
    
    # record atmosphere modal decomposition [rad]
    bPywfs[:,i] = np.matmul(calib_zernike_pywfs.M, pywfs.signal)*2*np.pi/src.wavelength
    bShGeo[:,i] = np.matmul(calib_zernike_shGeo.M, shGeo.signal)*2*np.pi/src.wavelength
    bShDiff[:,i] = np.matmul(calib_zernike_shDiff.M, shDiff.signal)*2*np.pi/src.wavelength
    cDm[:,i] = np.matmul(np.linalg.pinv(M2C_zernike), dm.coefs)*2*np.pi/src.wavelength
    
    if i>20:
        #src*tel
        tel.computePSF(zeroPaddingFactor=param['zeroPaddingFactor'])
        psfLongExp = psfLongExp + tel.PSF
    
    # to ensure loop is running    
    if i%20 == 0: 
        print('iteration '+str(i)+'/' + str(param['nIter'])+'\n'
              +'turbulence (var) = ' + str(phaseVar[i])+' rad^2\n'
              +'wavefront error (var) = ' + str(resPhaseVar[i]) + ' rad^2\n'
              +'Strehl ratio = '+str(int(np.floor(strehl[i]*100)))+' %\n')    

b = time.time()

print('Closing on '+param['wfs']+'\n')
print('Simulating the closed loop took '+str(int(b-a)) +' seconds\n')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# compute atmosphere modal variances from telemetry data
varZernikePywfs = np.var(bPywfs + cDm, axis=1)
varZernikeShGeo = np.var(bShGeo + cDm, axis=1)
varZernikeShDiff = np.var(bShDiff + cDm, axis=1)

# get theoretical modal variances (Zernike)
varZernikeTh = zernike_variance(param['diameter'], (param['r0'], param['L0']), range(2,param['nModes']+2)) # in noll definition, piston is the mode n° 1

# compute empirical atmosphere modal variances
if param['recordTurbulence'] == True:
    varZernikeEmp = np.var(np.matmul(np.transpose(zernike.modes), turbulence)/np.sum(tel.pupil), axis=1)

#%% fit variances and extract (r0, L0) with turlib

a=time.time()

ai2pwfs=np.zeros(param['nModes']+1)
ai2pwfs[1:]=varZernikePywfs
outputPywfs = iterative_estimator(param['diameter'],
                              modes=np.arange(2,param['nModesFit']+1), 
                              #modes=zernike.modes[:,0:param['nModesFit']],
                              ai2=ai2pwfs, 
                              noise_estimate=np.zeros(param['nModes']+1), 
                              n_rec_modes=param['nModes'], 
                              m=param['nModesShown'],# crosstalkMatrix_pywfs.shape[0]*crosstalkMatrix_pywfs.shape[1], 
                              c_mat=crosstalkMatrix_pywfs,
                              hro=param['radialOrder'], # last radial order of the modes being fitted, e.g nModesFit
                              lro=1)

ai2shDiff=np.zeros(param['nModes']+1)
ai2shDiff[1:]=varZernikeShGeo
outputShDiff = iterative_estimator(param['diameter'],
                              modes=np.arange(2,param['nModesFit']+1), 
                              #modes=zernike.modes[:,0:param['nModesFit']],
                              ai2=ai2shDiff, 
                              noise_estimate=np.zeros(param['nModes']+1),
                              n_rec_modes=param['nModes'],
                              m=param['nModesShown'],# crosstalkMatrix_pywfs.shape[0]*crosstalkMatrix_pywfs.shape[1], 
                              c_mat=crosstalkMatrix_pywfs,
                              hro=param['radialOrder'], # last radial order of the modes being fitted, e.g nModesFit
                              lro=1)

ai2shGeo=np.zeros(param['nModes']+1)
ai2shGeo[1:]=varZernikeShGeo
outputShGeo = iterative_estimator(param['diameter'],
                              modes=np.arange(2,param['nModesFit']+1),
                              #modes=zernike.modes[:,0:param['nModesFit']],
                              ai2=ai2shGeo,
                              noise_estimate=np.zeros(param['nModes']+1),
                              n_rec_modes=param['nModes'],
                              m=param['nModesShown'],# crosstalkMatrix_pywfs.shape[0]*crosstalkMatrix_pywfs.shape[1], 
                              c_mat=crosstalkMatrix_pywfs,
                              hro=param['radialOrder'], # last radial order of the modes being fitted, e.g nModesFit
                              lro=1)

b=time.time()

#%%

print('r0: '+str(param['r0']*100) + ' cm')
print('L0: '+str(param['L0']) + ' m\n')
print('estimated r0 (pwfs): ' + str(outputPywfs[0]*100) + ' cm')
print('estimated L0(pwfs): ' + str(outputPywfs[1]) + ' m\n')
print('estimated r0 (shGeo): ' + str(outputShGeo[0]*100) + ' cm')
print('estimated L0(shGeo): ' + str(outputShGeo[1]) + ' m\n')
print('The DM was controlled with ' + str(param['nModes']) + ' modes' + 
      ', the reconstruction was done with ' + str(param['nModesShown'])+ ' modes' +
      ', the fitting of the variances was performed over ' + str(param['nModesFit']) + ' modes.' + '\n')
print('Turlib processing took ' + str(int(b-a)) + ' seconds.')

#%%

fig2, axs2 = plt.subplots(2)
axs2[0].plot(range(len(outputPywfs[2])), outputPywfs[2],'r', label='turliboutputPywfs')
axs2[1].plot(range(len(varZernikePywfs)),varZernikePywfs,'b',label='turlibInput')
axs2[0].set_yscale('log')
axs2[1].set_yscale('log')

#%%%%%%%%%%%%%%%%%%%%%%%%%%%% Save data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# we save the data along with the parameter file
data = param.copy()
data['varZernikePywfs'] = varZernikePywfs
data['varZernikeShGeo'] = varZernikeShGeo
data['varZernikeTh'] = varZernikeTh
if param['recordTurbulence'] == True:
    data['varZernikeEmp'] = varZernikeEmp

data['crosstalkMatrix_pywfs'] = crosstalkMatrix_pywfs
data['gF_pywfs'] = gF_pywfs
data['gR_pywfs'] = gR_pywfs
data['gFplus_pywfs'] = gFplus_pywfs

dataFileName = r'D:\FEUP\simuPyramidForTurblib\data\\'+str(param['nIter']/param['frequency'])+'s'\
             +'_loopFrequency-'+str(param['frequency'])+'Hz'\
             +'_r0-'+str(param['r0']*100)+'cm'\
             +'_band-'+param['opticalBand']\
             +'_nModes-'+str(param['nModes'])\
             +'_wfs-'+str(param['wfs'])+'.npy'
np.save(dataFileName, data)

read_data = np.load(dataFileName,allow_pickle='TRUE').item()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plots rad %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plt.close('all')

dirc = r'D:\FEUP\simuPyramidForTurblib\figure\\'

# create fig
fig1, axs1 = plt.subplots(nrows=1, ncols=3)
axs1[0].plot(np.arange(len(resPhaseVar))*1/param['frequency'] ,resPhaseVar, 'b', label='resPhaseVar')
axs1[0].plot(np.arange(len(resPhaseVar))*1/param['frequency'], phaseVar, 'r', label='phaseVar')
axs1[0].set_xlabel('time (s)')
axs1[0].set_ylabel('phase var [rad^2]')
axs1[0].legend()
axs1[0].set_title('Entrance and residual phase variance over time')
axs1[2].imshow(np.log(psfLongExp))
axs1[2].set_title('Long exposure Psf (log scale)')
axs1[1].plot(np.arange(len(strehl))*1/param['frequency'], 100*strehl)
axs1[1].set_ylabel('Strehl Ratio [%]')
axs1[1].set_title('Strehl Ratio over time, loop closed on '+param['wfs'])
fig1.colorbar(mappable=plt.imshow(np.log(psfLongExp)), ax=axs1[2])

# save fig
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
fig1.savefig(dirc+'loopOverviewPapyrus_texp-'+str(param['nIter']/param['frequency'])+'s'
             +'_loopFrequency-'+str(param['frequency'])+'Hz'
             +'_r0-'+str(param['r0']*100)+'cm'
             +'_band-'+param['opticalBand']
             +'_nModes-'+str(param['nModes'])
             +'_wfs-'+str(param['wfs'])
             +'.png',
             bbox_inches='tight')

# create fig
fig3, axs3 = plt.subplots(1,1)
axs3.plot(range(1,len(varZernikePywfs)+1), varZernikePywfs, 'b', label='varZernikePywfs')
axs3.plot(range(1,len(varZernikeShGeo)+1), varZernikeShGeo, 'g', label='varZernikeShGeo')
axs3.plot(range(1,len(varZernikeTh)+1), varZernikeEmp, 'm', label='varZernikeEmp')
axs3.plot(range(1,len(varZernikeTh)+1), varZernikeTh, 'r', label='varZernikeTh')
axs3.legend()
axs3.set_yscale('log')
axs3.set_xlabel('Noll index')
axs3.set_ylabel('Variance (log scale) [rad^2]')
axs3.set_title('Variance of the modal coeficients of the theoretical and reconstructed atmosphere\n'
               +'loop closed on '+param['wfs'])

# save fig
manager = plt.get_current_fig_manager()
manager.full_screen_toggle()
fig3.savefig(dirc+'modalVariancePapyrus_texp-'+str(param['nIter']/param['frequency'])+'s'
             +'_loopFrequency-'+str(param['frequency'])+'Hz'
             +'_r0-'+str(param['r0']*100)+'cm'
             +'_band-'+param['opticalBand']
             +'_nModes-'+str(param['nModes'])
             +'_wfs-'+str(param['wfs'])
             +'.png',
             bbox_inches='tight')


#%% close fig
if param['show'] == False:
    plt.close(fig1)
    plt.close(fig3)