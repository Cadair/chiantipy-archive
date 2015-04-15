'''
functions needed for IPython multiprocessing module ipymspectrum
'''
import copy
import chianti
def doFf(inpt):
    ''' multiprocessing helper for freefree'''
    ionS = inpt[0]
    temperature = inpt[1]
    wavelength = inpt[2]
    abund = inpt[3]
    cont = chianti.core.continuum(ionS, temperature, abundance=abund)
    cont.freeFree(wavelength)
    return [ionS, copy.deepcopy(cont)]
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
    cont = chianti.core.continuum(ionS, temperature, abundance=abund)
    cont.freeBound(wavelength)
    return [ionS, copy.deepcopy(cont)]
    #
    # ----------------------------------------------
    #
def doIon1(inpt):
    '''
    multiprocessing helper for ion, also does two-photon
    '''
 #     [ionS, temperature, eDensity, wavelength, filter, allLines, abund]
    ionS = inpt[0]
    temperature = inpt[1]
    density = inpt[2]
    wavelength = inpt[3]
#    wvlRange = [wavelength.min(), wavelength.max()]
    filter = inpt[4]
    allLines = inpt[5]
    abund = inpt[6]
    em = input[7]
    thisIon = chianti.core.ion(ionS, temperature, density, abundance=abund)
#    thisIon.intensity(wvlRange = wvlRange, allLines = allLines, em=em)
    thisIon.spectrum(wavelength,  filter=filter, allLines=allLines,  em=em)
#        outList = [ionS, thisIon.Spectrum]
#    outList = [ionS, thisIon.Spectrum, thisIon.Intensity]
    outList = [ionS, thisIon]
#    outList = [ionS, ionS]
    if not thisIon.Dielectronic:
        if (thisIon.Z - thisIon.Ion) in [0, 1]:
            thisIon.twoPhoton(wavelength, em=em)
            outList.append(thisIon.TwoPhoton)
    return outList
    #
    # ----------------------------------------------
    #
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
    em = inpt[7]
    doContinuum = inpt[8]
    thisIon = chianti.core.ion(ionS, temperature, density, abundance=abund)
    thisIon.intensity(wvlRange = wvlRange, allLines = allLines, em=em)
    if 'errorMessage' not in thisIon.Intensity.keys():
        thisIon.spectrum(wavelength,  filter=filter, allLines=allLines,  em=em)
#        outList = [ionS, thisIon.Spectrum]
#    outList = [ionS, thisIon.Spectrum, thisIon.Intensity]
#    outList = [ionS, thisIon]
#    outList = [ionS, ionS]
    outList = [ionS, thisIon]
    if not thisIon.Dielectronic and doContinuum:
        if (thisIon.Z - thisIon.Ion) in [0, 1]:
            thisIon.twoPhoton(wavelength, em=em)
            outList.append(thisIon.TwoPhoton)
    return outList
    #
    # ----------------------------------------------
    #
