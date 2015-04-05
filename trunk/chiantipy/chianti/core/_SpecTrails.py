from datetime import datetime
from copy import copy
import numpy as np
import pylab as pl
import chianti.filters as chfilters
class _specTrails():
    '''
    a collection of methods for use in spectrum calculations
    '''
    def __init__(self):
        pass
        return
        #
        # ---------------------------------------------------------------------------
        #
    def convolve(self, filter=(chfilters.gaussianR, 1000.), label=0, verbose=0):
        '''
        the first application of spectrum calculates the line intensities within the specified wavelength range and for set of ions specified
        '''
        if type(label)!= type(0):
            if type(label) != str:
                print(' label must either be zero or a string')
                return
        #
        t1 = datetime.now()
        #
        nWvl = len(self.Wavelength)
        wavelength = self.Wavelength
        lineSpectrum = np.zeros((self.NTempDen, nWvl), 'float64').squeeze()
        for akey in self.IonInstances.iterkeys():
#            thisIon = self.IonInstances[akey]
            if not 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                if verbose:
                    print(' doing convolve on ion %s - %s'%(akey, self.IonInstances[akey].IonStr))
                self.IonInstances[akey].spectrum(wavelength, filter)
#                lineSpectrum = np.add(lineSpectrum, self.IonInstances[akey].Spectrum['intensity'])
                lineSpectrum += self.IonInstances[akey].Spectrum['intensity']
#                if self.NTempDen == 1:
#                    lineSpectrum += thisIon.Spectrum['intensity']
#                else:
#                    for iTempDen in range(self.NTempDen):
#                        lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
            else:
                if 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                    print(self.IonInstances[akey].Intensity['errorMessage'])
#
##        ijk = 0
##        for thisItem in self.IonInstances.iteritems():
##            ionS = thisItem[0]
##            thisIon = thisItem[1]
###        for thisIon in self.IonInstances.iteritems():
###            thisIon = self.IonInstances[akey]
##            if not 'errorMessage' in sorted(thisIon.Intensity.keys()):
##                if verbose:
##                    print(' doing convolve on ion %s - %s'%(ionS, thisIon.IonStr))
##                delattr(thisIon, 'Spectrum')
##                thisIon.spectrum(wavelength, filter)
##                ijk += 1
##                pl.figure(ijk)
##                pl.plot(wavelength, thisIon.Spectrum['intensity'][self.NTempDen/2])
##                pl.title(thisIon.IonStr)
###                lineSpectrum = np.add(lineSpectrum, self.IonInstances[akey].Spectrum['intensity'])
##                lineSpectrum += thisIon.Spectrum['intensity']
###                if self.NTempDen == 1:
###                    lineSpectrum += thisIon.Spectrum['intensity']
###                else:
###                    for iTempDen in range(self.NTempDen):
###                        lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
##            else:
##                if 'errorMessage' in sorted(thisIon.Intensity.keys()):
##                    print(thisIon.Intensity['errorMessage'])
##
##        ijk = 0
##        for thisIon in self.IonInstances:
##            ionS = thisIon.IonStr
###        for thisIon in self.IonInstances.iteritems():
###            thisIon = self.IonInstances[akey]
##            if not 'errorMessage' in sorted(thisIon.Intensity.keys()):
##                if verbose:
##                    print(' doing convolve on ion %s - %s'%(ionS, thisIon.IonStr))
###                delattr(thisIon, 'Spectrum')
##                thisIon.spectrum(wavelength, filter)
##                ijk += 1
##                pl.figure(ijk)
##                pl.plot(wavelength, thisIon.Spectrum['intensity'][self.NTempDen/2])
##                pl.title(thisIon.IonStr)
###                lineSpectrum = np.add(lineSpectrum, self.IonInstances[akey].Spectrum['intensity'])
##                lineSpectrum += thisIon.Spectrum['intensity']
###                if self.NTempDen == 1:
###                    lineSpectrum += thisIon.Spectrum['intensity']
###                else:
###                    for iTempDen in range(self.NTempDen):
###                        lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
##            else:
##                if 'errorMessage' in sorted(thisIon.Intensity.keys()):
##                    print(thisIon.Intensity['errorMessage'])
##
#        ijk = 0
#        for ion in range(len(self.IonInstances)):
#            ionS = self.IonInstances[ion].IonStr
##        for thisIon in self.IonInstances.iteritems():
##            thisIon = self.IonInstances[akey]
#            if not 'errorMessage' in sorted(self.IonInstances[ion].Intensity.keys()):
#                if verbose:
#                    print(' doing convolve on ion %s - %s'%(ionS, self.IonInstances[ion].IonStr))
##                delattr(thisIon, 'Spectrum')
#                self.IonInstances[ion].spectrum(wavelength, filter)
#                ijk += 1
#                pl.figure(ijk)
#                pl.plot(wavelength, self.IonInstances[ion].Spectrum['intensity'][self.NTempDen/2])
#                pl.title(self.IonInstances[ion].IonStr)
##                lineSpectrum = np.add(lineSpectrum, self.IonInstances[akey].Spectrum['intensity'])
#                lineSpectrum += self.IonInstances[ion].Spectrum['intensity']
##                if self.NTempDen == 1:
##                    lineSpectrum += thisIon.Spectrum['intensity']
##                else:
##                    for iTempDen in range(self.NTempDen):
##                        lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
#            else:
#                if 'errorMessage' in sorted(self.IonInstances[ion].Intensity.keys()):
#                    print(self.IonInstances[ion].Intensity['errorMessage'])
#        self.LineSpectrum = lineSpectrum

        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}

        total = self.FreeFree['intensity'] + self.FreeBound['intensity'] + self.TwoPhoton['intensity'] + self.LineSpectrum['intensity']
        self.Total = total
        #
        if self.NEm == 0:
            integrated = total
        else:
            if self.NEm == 1:
                integrated = total*self.Em
            else:
                integrated = np.zeros_like(wavelength)
                for iTempDen in range(self.NTempDen):
                    integrated += total[iTempDen]*self.Em[iTempDen]
        #
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %12.3e'%(dt.seconds))
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':self.Em, 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':self.Em, 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}
        return
