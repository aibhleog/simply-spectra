import numpy as np
import pandas as pd
import scipy.integrate as integrate
from scipy.interpolate import interp1d
from calzetti import * # written by TAH

__author__ = 'Taylor Hutchison'
__email__ = 'aibhleog@tamu.edu'

# -- FILTERS ---------------------------- #
f105w = np.loadtxt('/home/aibhleog/Desktop/keckData/HST-WFC3_IR.F105W.dat')
f160w = np.loadtxt('/home/aibhleog/Desktop/keckData/HST-WFC3_IR.F105W.dat')

irac36 = np.loadtxt('/home/aibhleog/Desktop/keckData/information/Spitzer_3.6_IRAC.dat')
irac45 = np.loadtxt('/home/aibhleog/Desktop/keckData/information/Spitzer_4.5_IRAC.dat')
# --------------------------------------- #

def integrated_flux(nu,filt,flux):
    '''
    Purpose:
    Used to integrate a spectrum over a bandpass.

    Input:
    frequency, filter transmission curve, and spectrum to integrate
    in F_nu units of [erg/s/cm^2/Hz].

    Returns:
    Total integrated flux in units of [erg/s/cm^2/Hz].
    '''
    return integrate.trapz(filt*flux/nu,nu)/integrate.trapz(filt/nu,nu)


def reorder_spectrum(df):
    '''
    Purpose:
    Reorders the spectrum to go the opposite way.  Useful for the individual
    spectra because the models are in descending wavelength order originally.
    
    Input:
    dataframe of spectrum
    
    Returns:
    Reordered dataframe
    '''
    new_df = df.sort_index(axis=0,ascending=False)
    new_df.reset_index(inplace=True,drop=True)
    return new_df


def get_cloudy_table(foldername,age,z,zneb=False,high=False,imf=False,o3=False):
    '''
    Purpose:
    Gathers appropriate Cloudy output table (created by subroutine when
    creating all of the Cloudy modelling) given particular parameters.

    Input:
    Folder name which specifies the stellar population (SP) used, age of
    SP models used in Cloudy modelling, stellar metallicity, and nebular
    metallicity.
    ----->  the nebular metallicity defaults to matching stellar unless
            otherwise specified.
    ----->  if wanting to access the mixed S99 SPs, imf flag needs to be given
            float value of 2.0 or 1.7 (really these slopes are negative, just
            to be clear, but are positive for the purposes of this function).

    Returns:
    Dataframe (df) with all of the important lineflux information collected
    during the subroutines running the Cloudy modelling.
    '''
    home = '/home/aibhleog/Desktop/catalogs/cloudy/models/steidel/'
    if high == False:
        if o3 == False: filename = '%s/all-age.vary_ageUZ_wrtHbeta.txt'%(foldername)
        else: filename = '%s/all.vary_ageUZ_wrtHbeta_FEWER_LINES.txt'%(foldername)
    else:
        filename = '%s/higherU_vary_ageUZ_wrtHbeta.txt'%(foldername)
    df = pd.read_csv(home+filename,delimiter='\s+',header=0)
    if zneb == False: zneb = z
    if imf == 2.: df = df.loc[df.imf.values == 2.] # slices for imf when available
    if imf == 1.7: df = df.loc[df.imf.values == 1.7] # slices for imf when available
    if high == False:
        df = df.loc[df.age.values == age] # slicing for age
    df = df.loc[df.Z.values == z] # slicing for metallicity
    df = df.loc[df.Zneb.values == zneb] # slicing for nebular metallicity
    return df

def get_cloudy_table_agn(foldername,zneb):
    '''
    Purpose:
    Gathers appropriate Cloudy output table (created by subroutine when
    creating all of the Cloudy modelling) given particular parameters.

    Input:
    Folder name which specifies the AGN models and nebular
    metallicity.
    ----->  the nebular metallicity defaults to matching stellar unless
            otherwise specified.

    Returns:
    Dataframe (df) with all of the important lineflux information collected
    during the subroutines running the Cloudy modelling.
    '''
    home = '/home/aibhleog/Desktop/catalogs/cloudy/models/steidel/'
    filename = '%s/vary_ageUZ_wrtHbeta.txt'%(foldername)
    df = pd.read_csv(home+filename,delimiter='\s+',header=0)
    df = df.loc[df.Zneb.values == zneb] # slicing for nebular metallicity
    return df

def get_cloudy_line(line,input_df):
    '''
    Purpose:
    A quick and easy way to access particular lines in the Cloudy tables.

    Input:
    Line identifier name.

    Returns:
    Line flux of particular line, relative to Hbeta.
    '''
    lineidx = ['lya','heii','siiii1','siiii2','ciii1','ciii2','oii1','oii2',\
            'neiii','oiii1','oiii2','ha','nii','sii1','sii2']
    linename = ['Lya1216','HeII1640','[SiIII]1883','SiIII]1892','[CIII]1907','CIII]1909',\
            '[OII]3727','[OII]3729','[NeIII]3876','[OIII]4959','[OIII]5007',\
            'Ha6563','[NII]6585','[SII]6718','[SII]6732']
    index = lineidx.index(line)
    return input_df[linename[index]].values

# SP FOLDERS LOOPED INTO THE "ELSE" STATEMENT ON LINE 149:
# folder = ['binary_cont_100','binary_cont_300','single_cont_100','single_cont_300']
def get_cloudy_spec(foldername,mlimit,age,z,zneb=False,imf=1.7,ioni=False):
    '''
    Purpose:
    Gathers the Cloudy output spectra for the full range of ionization parameter (U)
    given particular parameters.

    Input:
    Name of folder where spectra are located, upper mass limit for IMF, age of model,
    stellar metallicity, and nebular metallicity.
    ----->  the nebular metallicity defaults to matching stellar unless
            otherwise specified.
    ----->  if wanting to access the mixed S99 SPs, imf flag needs to be given
            float value of 2.0 or 1.7 (really these slopes are negative, just
            to be clear, but are positive for the purposes of this function).

    Returns:
    Dataframe (df) of chosen spectra, covering the full range of U.
    '''
    home = '/home/aibhleog/Desktop/catalogs/cloudy/models/steidel/'
    df = pd.DataFrame({'u':[],'wavelength':[],'spectrum':[]})
    colnames = ['wavelength','spectrum']
    cols = ['u','wavelength','spectrum']
    if age != 7.477: age = int(age)

    if zneb == False: zneb = z
    ulow = -3.5
    if foldername == 'single_s99_kroupa' and z == 0.05 and zneb == 0.05: ulow = -3.3

    rangeit = np.concatenate((np.arange(ulow,-1.4,0.2),[-1.25,-1.]))
    for u in rangeit:
        if u == -1.25: u = round(u,2)
        elif u == -1.: u = int(u)
        else: u = round(u,1)
        if foldername == 'single_s99_miscIMF': 
            if imf == 2: imf = int(imf)
            filename = '%s/age%simf%sz%szneb%su%s_%s.con'%(foldername,age,imf,z,zneb,u,mlimit)
        elif foldername == 'single_s99_kroupa' and mlimit == 300: 
            filename = '%s/age%sz002zneb%su%s_%s.con'%(foldername,age,zneb,u,mlimit)
        else: filename = '%s/age%sz%szneb%su%s_%s.con'%(foldername,age,z,zneb,u,mlimit)
        fillerdf = pd.read_csv(home+filename,names=colnames,\
                                   delimiter='\s+',skiprows=1,usecols=[0,6])
        uarr = np.full(len(fillerdf),u)
        filldf = pd.DataFrame({'u':uarr,'wavelength':fillerdf.wavelength.values,\
                                  'spectrum':fillerdf.spectrum.values})
        df = df.append(filldf,ignore_index=True)
    df = df[cols]

    # returning the list of spectra or an individual spectrum
    if ioni == False: return df # if you want all U's 
    else:
        spec = df.loc[df.u.values == round(ioni,1)].copy()
        spec.reset_index(drop=True,inplace=True)
        return spec # spectrum for one particular ionization


def get_cloudy_spec_agn(foldername,zneb,manual=False,ioni=[1.]):
    '''
    Purpose:
    Gathers the Cloudy output spectra for the full range of ionization parameter (U)
    given particular parameters.

    Input:
    Name of folder where spectra are located and nebular metallicity.

    Returns:
    Dataframe (df) of chosen spectra, covering the full range of U.
    '''
    home = '/home/aibhleog/Desktop/catalogs/cloudy/models/steidel/'
    df = pd.DataFrame({'u':[],'wavelength':[],'spectrum':[]})
    colnames = ['wavelength','spectrum']
    cols = ['u','wavelength','spectrum']

    try:
        len(ioni)
    except TypeError: # if not an array or list
        ioni = np.asarray(ioni)

    if manual == True:
        for u in ioni:    
            filename = '%s/manualagn_age7zneb%su%s.con'%(foldername,zneb,u)
            fillerdf = pd.read_csv(home+filename,names=colnames,\
                                       delimiter='\s+',skiprows=1,usecols=[0,6])
            uarr = np.full(len(fillerdf),u)
            filldf = pd.DataFrame({'u':uarr,'wavelength':fillerdf.wavelength.values,\
                                      'spectrum':fillerdf.spectrum.values})
            df = df.append(filldf,ignore_index=True)
    else:
        ulow = -3.5
        rangeit = np.concatenate((np.arange(ulow,-1.4,0.2),[-1.25,-1.]))
        for u in rangeit:
            if u == -1.25: u = round(u,2)
            elif u == -1.: u = int(u)
            else: u = round(u,1)
            filename = '%s/agn_zneb%su%s.con'%(foldername,zneb,u)
            fillerdf = pd.read_csv(home+filename,names=colnames,\
                                       delimiter='\s+',skiprows=1,usecols=[0,6])
            uarr = np.full(len(fillerdf),u)
            filldf = pd.DataFrame({'u':uarr,'wavelength':fillerdf.wavelength.values,\
                                      'spectrum':fillerdf.spectrum.values})
            df = df.append(filldf,ignore_index=True)
    df = df[cols]
    return df

def get_cloudy_calzetti(ebv,input_df,manual=False,ioni=[1.]):
    sed = []
    try:
        len(ioni)
    except TypeError: # if not an array or list
        ioni = np.asarray(ioni)

    full_df = pd.DataFrame({'u':[],'wavelength':[],'spectrum':[]})
    colnames = ['wavelength','spectrum']
    cols = ['u','wavelength','spectrum']

    if manual == True: 
    # this just allows me to use this function when I manually run a Cloudy script
        for u in ioni:
            # testing
            df = input_df.loc[input_df.u.values == u].copy()
            below = df[df.wavelength.values <= 2.2]
            here = below[below.wavelength.values >= 0.12]

            wave = here.wavelength.values
            nu = 2.998e14/wave
            fnu = here.spectrum.values/nu

            calz = calzetti(wave,fnu,ebv)
            df['spectrum'][here.index.values] = calz
            full_df = full_df.append(df,ignore_index=True)
    else:
        for u in [-3.5,-3.3,-3.1,-2.9,-2.7,-2.5,-2.3,-2.1,-1.9,-1.7,-1.5,-1.25,-1.]:
            df = input_df.loc[input_df.u.values == u].copy()
            below = df[df.wavelength.values <= 2.2]
            here = below[below.wavelength.values >= 0.12]

            wave = here.wavelength.values
            nu = 2.998e14/wave
            fnu = here.spectrum.values/nu

            calz = calzetti(wave,fnu,ebv)
            df['spectrum'][here.index.values] = calz
            full_df = full_df.append(df,ignore_index=True)

    full_df = full_df[cols]
    return full_df # spectrum in units of vFv [erg/s/cm^2]

def get_cloudy_ew(line,input_df,manual=False,ioni=[1.]):
    '''
    Purpose:
    Calculates the equivalent width (EW) of a given line.

    Input:
    The line you wish to calculate the flux of (if listed, otherwise you'll
    need to identify the index points for your target line's emission), the
    dataframe (df) of spectra (varying by ionization parameter, U).
    ----->  important to note that the input df is specified to have a particular
            stellar and nebular metallicity by a different function.

    Returns:
    Calculated restframe EW of line in Angstroms.

    -------------------------------------------------------
    Other Lines:
        lo,hi,loc,hic,line = 2933,2940,2935,2938,'H-alpha'
        lo,hi,loc,hic,line = 2564,2568,2564,2567,'CIII]'
        lo,hi,loc,hic,line = 2427,2434,2429,2432,'Ly-alpha'
    -------------------------------------------------------
    '''
    line_ew = []
    if line == 'ciii':
        lo,hi,loc,hic = 2564,2568,2564,2567 # CIII Cloudy coordinates used
    elif line == 'lya':
        lo,hi,loc,hic = 2427,2434,2429,2432 # Ly-alpha Cloudy coordinates used
    elif line == 'nv':
        lo,hi,loc,hic = 2434,2441,2436,2439 # NV line (1243 line)
    else:
        print('Wrong line provided.',end='\n\n')

    try:
        len(ioni)
    except TypeError: # if not an array or list
        ioni = np.asarray(ioni)

    if manual == True: # this just allows me to use this function when I manually run a Cloudy script
        for u in ioni:
            df = input_df.loc[input_df.u.values == u]
            wave = df.wavelength.values*1e4
            nu = 2.998e18/wave
            fnu = df.spectrum.values/nu
            fnu = fnu[::-1] # values were backwards
            wave = wave[::-1] # values were backwards

            con = np.concatenate([fnu[lo:loc+1],fnu[hic-1:hi]])
            wavcon = np.concatenate([wave[lo:loc+1],wave[hic-1:hi]])
            f = interp1d(wavcon,con)
            con = f(wave[lo:hi])

            ew = np.trapz(fnu[lo:hi]/con-1,wave[lo:hi])
            #print('U: %s, W: %s'%(u,ew))
            line_ew.append(ew)
    else:
        ulow = -3.5
        rangeit = np.concatenate((np.arange(ulow,-1.4,0.2),[-1.25,-1.]))
        for u in rangeit:
            if u == -1.25: u = round(u,2)
            elif u == -1.: u = int(u)
            else: u = round(u,1)
            df = input_df.loc[input_df.u.values == u]
            wave = df.wavelength.values*1e4
            nu = 2.998e18/wave
            fnu = df.spectrum.values/nu
            fnu = fnu[::-1] # values were backwards
            wave = wave[::-1] # values were backwards

            con = np.concatenate([fnu[lo:loc+1],fnu[hic-1:hi]])
            wavcon = np.concatenate([wave[lo:loc+1],wave[hic-1:hi]])
            f = interp1d(wavcon,con)
            con = f(wave[lo:hi])

            ew = np.trapz(fnu[lo:hi]/con-1,wave[lo:hi])
            #print('U: %s, W: %s'%(u,ew))
            line_ew.append(ew)
    return line_ew

def get_cloudy_flux(line,input_df):
    '''
    Purpose:
    Currently unused but can be used to calculate the lineflux of particular
    lines for a given input dataframe (df).

    Input:
    The line you wish to calculate the flux of (if listed, otherwise you'll
    need to identify the index points for your target line's emission), 
    dataframe (df) of spectra (varying by ionization parameter, U).
    ----->  important to note that the input df is specified to have a particular
            stellar and nebular metallicity by a different function.

    Returns:
    The lineflux of the chosen line and the flux of the continuum at 1500A.
    '''
    lineflux,line1500 = [],[]
    mlo,mhi = 2492,2495 # M1500 flux
    if line == 'ciii':
        lo,hi = 2564,2567 # CIII Cloudy coordinates used
    elif line == 'lya':
        lo,hi = 2429,2432 # Ly-alpha Cloudy coordinates used
    else:
        print('Wrong line provided.',end='\n\n')

    for u in np.arange(-3.5,-1.4,0.2):
        df = input_df.loc[input_df.u.values == u]
        wave = df.wavelength.values*1e4
        nu = 2.998e18/wave
        fnu = df.spectrum.values/nu
        fnu = fnu[::-1] # values were backwards
        wave = wave[::-1] # values were backwards

        f1500 = np.trapz(fnu[mlo:mhi],wave[mlo:mhi])
        fline = np.trapz(fnu[lo:hi],wave[lo:hi])
        lineflux.append(fline)
        line1500.append(f1500)
    lineflux = np.asarray(lineflux)
    line1500 = np.asarray(line1500)
    return lineflux,line1500 # returns line flux and f_1500

def get_cloudy_norm_spec(input_df,norm=True,red=7.5027,manual=False,ioni=[1.]):
    '''
    Purpose:
    This function redshifts and normalizes the spectrum to match the GND_42912's
    F106W magnitude [AB] to calculate IRAC [3.6]-[4.5] color

    Input:
    Dataframe (df) of spectra (varying by ionization parameter, U), the redshift 
    (default GND_42912), and flags to allow for manual running of one or few
    spectra with different U's.
    ----->  important to note that the input df is specified to have a particular
            stellar and nebular metallicity by a different function.

    Returns:
    The magnitude(s) of IRAC CH1, [3.6], and IRAC CH2, [4.5] 
    '''
    arr_m36, arr_m45 = [],[]

    try:
        len(ioni)
    except TypeError: # if not an array or list
        ioni = np.asarray(ioni)

    if manual == True: # this just allows me to use this function when I manually run a Cloudy script
        for u in ioni:
            df = input_df.loc[input_df.u.values == u]
            wave = df.wavelength.values*1e4
            nu = 2.998e18/wave
            fnu = df.spectrum.values/nu
            fnu = fnu[::-1] # values were backwards
            wave = wave[::-1] # values were backwards
            znu = nu*(1+red)
            zwave = wave*(1+red)
            f = interp1d(zwave,fnu) #

            if norm == True:
                # F160W
                fnufilt = f(f160w[:,0])
                scale = np.power(10,-0.4*(25.48+48.6)) / integrated_flux(2.998e18/f160w[:,0],f160w[:,1],fnufilt)
            else : scale = 1.

            # IRAC [3.6]
            fnufilt = f(irac36[:,0])
            intflux = integrated_flux(2.998e18/irac36[:,0],irac36[:,1],fnufilt*scale)
            m36 = -2.5*np.log10(intflux)-48.6
            arr_m36.append(m36)
            # IRAC [4.5]
            fnufilt = f(irac45[:,0])
            intflux = integrated_flux(2.998e18/irac45[:,0],irac45[:,1],fnufilt*scale)
            m45 = -2.5*np.log10(intflux)-48.6
            arr_m45.append(m45)
    else:
        ulow = -3.5
        rangeit = np.concatenate((np.arange(ulow,-1.4,0.2),[-1.25,-1.]))
        for u in rangeit:
            if u == -1.25: u = round(u,2)
            elif u == -1.: u = int(u)
            else: u = round(u,1)
            df = input_df.loc[input_df.u.values == u]
            wave = df.wavelength.values*1e4
            nu = 2.998e18/wave
            fnu = df.spectrum.values/nu
            fnu = fnu[::-1] # values were backwards
            wave = wave[::-1] # values were backwards
            znu = nu*(1+red)
            zwave = wave*(1+red)
            f = interp1d(zwave,fnu)

            if norm == True:
                # F160W
                fnufilt = f(f160w[:,0])
                scale = np.power(10,-0.4*(25.48+48.6)) / integrated_flux(2.998e18/f160w[:,0],f160w[:,1],fnufilt)
            else : scale = 1.

            # IRAC [3.6]
            fnufilt = f(irac36[:,0])
            intflux = integrated_flux(2.998e18/irac36[:,0],irac36[:,1],fnufilt*scale)
            m36 = -2.5*np.log10(intflux)-48.6
            arr_m36.append(m36)
            # IRAC [4.5]
            fnufilt = f(irac45[:,0])
            intflux = integrated_flux(2.998e18/irac45[:,0],irac45[:,1],fnufilt*scale)
            m45 = -2.5*np.log10(intflux)-48.6
            arr_m45.append(m45)

    arr_m36, arr_m45 = np.asarray(arr_m36),np.asarray(arr_m45) # so that you can subtract them
    return arr_m36,arr_m45
