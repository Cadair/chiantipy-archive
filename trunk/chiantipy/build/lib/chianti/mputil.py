''' functions needed for multiprocessing module mspectrum'''
import chianti
def doFfQ(inQ, outQ):
    ''' multiprocessing helper for freefree'''
    for inputs in iter(inQ.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        abund = inputs[3]
        cont = chianti.core.continuum(ionS, temperature, abund=abund)
        cont.freeFree(wavelength)
        outQ.put(cont.FreeFree)
    return
    #
    # ----------------------------------------------
    #
def doFbQ(inQ, outQ):
    ''' multiprocessing helper for freeBound'''
    for inputs in iter(inQ.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        wavelength = inputs[2]
        abund = inputs[3]
        cont = chianti.core.continuum(ionS, temperature, abund=abund)
        cont.freeBound(wavelength)
        outQ.put(cont.FreeBound)
    return
    #
    # ----------------------------------------------
    #
def doIonQ(inQueue, outQueue):
    ''' multiprocessing helper for ion, also does two-photon'''
    for inputs in iter(inQueue.get, 'STOP'):
        ionS = inputs[0]
        temperature = inputs[1]
        density = inputs[2]
        wavelength = inputs[3]
        wvlRange = [wavelength.min(), wavelength.max()]
        filter = inputs[4]
        allLines = inputs[5]
        abund = inputs[6]
        thisIon = chianti.core.ion(ionS, temperature, density, abund=abund)
        thisIon.intensity(wvlRange = wvlRange, allLines = allLines)
        thisIon.spectrum(wavelength,  filter=filter)
#        outList = [ionS, thisIon.Spectrum]
        outList = [ionS, thisIon.Spectrum, thisIon.Intensity]
        if not thisIon.Dielectronic:
            if (thisIon.Z - thisIon.Ion) in [0, 1]:
                thisIon.twoPhoton(wavelength)
                outList.append(thisIon.TwoPhoton)
        outQueue.put(outList)
    return
    #
    # ----------------------------------------------
    #
