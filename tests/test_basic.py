def test_placeholder():
    assert 1 + 1 == 2


# def test_import():
#     from starrotator import starrotator 


def test_imports():
    #Test all dependencies
    from starrotator import StarRotator
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.special import sph_harm
    import requests
    import shutil
    import urllib.request as request
    from contextlib import closing
    import os.path
    import starrotator.lib.test as test
    import starrotator.lib.operations as ops
    import starrotator.lib.stellar_spectrum
    import starrotator.lib.vgrid as vgrid
    import starrotator.lib.integrate
    import starrotator.lib.util as util
    import astropy.io.fits as fits
    import astropy.units as u
    import scipy.interpolate as interp
    import astropy.constants as consts
    import numpy as np
    import copy
    import imp
    import sys
    from importlib.resources import files

"""
def test_demo_files():
    from importlib.resources import files
    import starrotator.lib.util as util
    data_path = files("starrotator.data").joinpath("demo_system.txt")
    obs_path = files("starrotator.data").joinpath("demo_observations.txt")
    paramdict = util.read_into_dictionary(data_path)
    util.check_integrity_input(paramdict)



def test_StarRotator():
    from starrotator import StarRotator
    import matplotlib.pyplot as plt
    from matplotlib import cm, colors
    from mpl_toolkits.mplot3d import Axes3D
    from scipy.special import sph_harm
    import requests
    import shutil
    import urllib.request as request
    from contextlib import closing
    import os.path
    import starrotator.lib.test as test
    import starrotator.lib.operations as ops
    import starrotator.lib.stellar_spectrum
    import starrotator.lib.vgrid as vgrid
    import starrotator.lib.integrate
    import astropy.io.fits as fits
    import astropy.units as u
    import scipy.interpolate as interp
    import astropy.constants as consts
    import numpy as np
    import copy
    import imp
    import sys
    from importlib.resources import files




    KELT2 = StarRotator(586.0,592.0,50.0)
    wl = KELT2.wl
    F_out = KELT2.stellar_spectrum
    spectra = KELT2.spectra
    residuals = KELT2.residual
    KELT2.convolve_spectral_resolution()

    assert KELT2.status == 'success'


    #Test that multiple-convolution is detected and blocked:
    try:
        KELT2.convolve_spectral_resolution()
        error_trigger=1
    except:
        pass
    if error_trigger==1:
        raise Exception("ERROR: Trying convolution twice in a row should be caught.")
    error_trigger=0

    #Testing another grid and another vartype.
    KELT3 = StarRotator(586,592.0,13)
    KELT3.convolve_spectral_resolution()

    assert KELT3.status == 'success'

    #Now test the whole thing with a dictionary as input. First test that the input is well
    #tested:
    in_dict = {'lala':1.0}
    try:
        KELT4 = StarRotator(586,592.0,13,input=in_dict)
        error_trigger=1
    except:
        pass
    if error_trigger==1:
        raise Exception("ERROR: Wrong input dictionary not caught.")
    error_trigger=0


    in_dict = {'veq':114000.0,
    'stelinc':90.0,
    'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
    'u1':0.93,'u2':-0.23,'R':115000.,'mus':0,'model':'PHOENIX','sma_Rs':3.153,
            'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
            'phases':[-0.02,-0.01,0.0,0.01,0.02]}
    KELT5 = StarRotator(586,592.0,13,input=in_dict)

    assert KELT5.status == 'success'

    


    pysme_error = 0
    try:
        from pysme.sme import SME_Structure as SME_Struct
        from pysme.abund import Abund
        from pysme.synthesize import synthesize_spectrum
        from pysme.solve import solve
        from pysme.linelist.vald import ValdFile
    except:
        print('WARNING: PYSME cannot be imported. PSME functionality will not be available.')
        n_warnings+=1
        pysme_error = 1


    if pysme_error == 0:
        in_dict = {'veq':114000.0,
        'stelinc':90.0,
        'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
        'u1':0.93,'u2':-0.23,'R':115000.,'mus':0,'model':'pySME','sma_Rs':3.153,
        'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
        'phases':[-0.02,-0.01,0.0,0.01,0.02]}#Without setting abund and grid_model keywords.
        try:
            KELT9 = StarRotator(586,592.0,13,input=in_dict)
            error_trigger=1
        except:
            pass
        if error_trigger==1:
            raise Exception("ERROR: Wrong input dictionary not caught.")
        error_trigger=0

        in_dict = {'veq':114000.0,
        'stelinc':90.0,
        'drr':0.0,'T':10000.0,'FeH':0.0,'logg':4.0,
        'u1':0.93,'u2':-0.23,'R':115000.,'mus':5,'model':'pySME','sma_Rs':3.153,
        'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
        'phases':[-0.02,-0.01,0.0,0.01,0.02],'grid_model':'atlas12.sav','abund':[],
        'linelist_path':'input/demo_linelist.dat'}
        KELT9 = StarRotator(586,592.0,13,input=in_dict)


        in_dict = {'veq':114000.0,
        'stelinc':90.0,
        'drr':0.0,'T':6000.0,'FeH':0.3,'logg':4.2,
        'u1':0.93,'u2':-0.23,'R':115000.,'mus':5,'model':'pySME','sma_Rs':3.153,
        'e':0.0,'omega':0.0,'inclination':86.79,'obliquity':-84.8,'RpRs':0.08228,'P':1.4811235,
        'phases':[-0.02,-0.01,0.0,0.01,0.02],'grid_model':'marcs2014.sav','abund':[],
        'linelist_path':'input/demo_linelist.dat'}
        KELT9 = StarRotator(586,592.0,13,input=in_dict)

    print('')
    print('')
    print('')
    print('Tests complete.')
    print(f'{n_warnings} warnings triggered.')

    """