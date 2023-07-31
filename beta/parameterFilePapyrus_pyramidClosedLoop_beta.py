# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 14:13:59 2023

@author: flero
"""

from turlib import n_modes_from_radial_order

def initializeParameterFiles():
    
    param = dict()
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SOURCE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['opticalBand'] = 'R' # star optical band (define both wavelength and zero point)
    param['magnitude'] = 0 # star magnitude
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATMOSPHERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['r0'] = 0.05 # value of r0 in the visibile in [m]
    param['L0'] = 30  # value of L0 in the visibile in [m]
    param['fractionnalR0'] = [0.45,0.1,0.1,0.25,0.1] # Cn2 profile
    param['windSpeed'] = [10,12,11,15,20] # wind speed of the different layers in [m.s-1]
    param['windDirection'] = [0,72,144,216,288] # wind direction of the different layers in [degrees]
    param['altitude'] = [0, 1000,5000,10000,12000] # altitude of the different layers in [m]
     
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TELESCOPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    param['diameter'] = 1.5 # M1 diameter [m]
    param['nSubap'] = 32 # number of subaps in the diameter of the pyramid pupil mask
    param['nPxPerSubap'] = 4 # ???? for pyramid
    param['resolution'] = param['nSubap']*param['nPxPerSubap'] # number of pixels in M1 diameter
    param['lenghtSubap'] = param['diameter']/param['nSubap']
     
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['nModes'] = 153 # number of modes of the modal callibration/command basis
    param['stroke'] = 1e-9 # standard deviation in [m] of each mode
    param['nSubapDm'] = 16 # nAct -1
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WFS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['wfs'] = 'pywfs' # choose who controls the loop : 'pywfs' or 'shGeo'
    param['modulationRadius'] = 3 # modulation radius in lambda / D
    param['lightRatio'] = 0.1 # light ratio for useful pixels selection
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
    param['frequency'] = 1000 # loop frequency in [Hz]
    param['nIter'] = 20000 # number of iterations
    param['gain'] = 0.5 # gain of the integrator
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOCAL PLANE SENSOR %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['zeroPaddingFactor'] = 2
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%% Simulation choices %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['calibrate'] = True # trigger calibration. If false, use previous calibration results.
    param['recordTurbulence'] = True # enable phase screen recording while closing the loop {bol}
    
    ##%%%%%%%%%%%%%%%%%%%%%%% Post Processing choices %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['nModesShown'] = 500 # number of modes being shown to the pywfs
    param['radialOrder'] = 8 # highest radial order included in the fit of the variances performed by turlib
    param['nModesFit'] = n_modes_from_radial_order(param['radialOrder']) # number of modes to be fitted by turlib, must be consistent with radial orders
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot preferences %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    param['nPsf'] = 9 # number of psf compute along the loop
    param['show'] = True # if True, figures are displayed    
    
    ##%%%%%%%%%%%%%%%%%%%%%%%%%%% Folder connections %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dirc = r'D:\FEUP\simuPyramidForTurblib\beta\\' # main directory
    param['pathDataCalibration'] = dirc + r'data\calibration\\'
    param['pathDataLoop'] = dirc + r'data\loop\\'
    
    return param