'''
functions needed for multiprocessing module ipymspectrum
'''
import chianti
def doFf(inpt):
    ''' multiprocessing helper for freefree'''
    ionS = inpt[0]
    temperature = inpt[1]
    wavelength = inpt[2]
    abund = inpt[3]
    cont = chianti.core.continuum(ionS, temperature, abund=abund)
    cont.freeFree(wavelength)
    return cont.FreeFree
    #
    # ----------------------------------------------
    #
def doFb(inpt):
    '''
    multiprocessing helper for freeBound
    '''
    ionS = inpt[0]
    temperature = inpt[1]
    wavelength = inpt[2]
    abund = inpt[3]
    cont = chianti.core.continuum(ionS, temperature, abund=abund)
    cont.freeBound(wavelength)
    return cont.FreeBound
    #
    # ----------------------------------------------
    #
def doIon(inpt):
    '''
    multiprocessing helper for ion, also does two-photon
    '''
 #     [ionS, temperature, eDensity, wavelength, filter, allLines, abund]
    ionS = inpt[0]
    temperature = inpt[1]
    density = inpt[2]
    wavelength = inpt[3]
    wvlRange = [wavelength.min(), wavelength.max()]
    filter = inpt[4]
    allLines = inpt[5]
    abund = inpt[6]
    thisIon = chianti.core.ion(ionS, temperature, density, abund=abund)
    thisIon.intensity(wvlRange = wvlRange, allLines = allLines)
    thisIon.spectrum(wavelength,  filter=filter)
#        outList = [ionS, thisIon.Spectrum]
    outList = [ionS, thisIon.Spectrum, thisIon.Intensity]
    if not thisIon.Dielectronic:
        if (thisIon.Z - thisIon.Ion) in [0, 1]:
            thisIon.twoPhoton(wavelength)
            outList.append(thisIon.TwoPhoton)
    return outList
    #
    # ----------------------------------------------
    #
