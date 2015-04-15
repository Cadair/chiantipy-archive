from datetime import datetime
import numpy as np
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
    def convolve(self, wavelength=0, filter=(chfilters.gaussianR, 1000.), label=0, verbose=0):
        '''
        the first application of spectrum calculates the line intensities within the specified wavelength range and for set of ions specified
        
        wavelength is will not be used if applied to 'spectrum' objects
        
        wavelength IS need for 'bunch' objects - in this case, the wavelength should not extend beyond the limits of the
        wvlRange used for the 'bunch' calculation
        
        '''
        if type(label)!= type(0):
            if type(label) != str:
                print(' label must either be zero or a string')
                return
        #
        t1 = datetime.now()
        #:
        if hasattr(self, 'Wavelength'):
            nWvl = len(self.Wavelength)
            wavelength = self.Wavelength
        elif type(wavelength) == int:
            print(' a wavelength array must be given')
            return
        else:
            nWvl = len(wavelength)
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
        #
        total = self.LineSpectrum['intensity']
        #
        # the following is required in order to be applied to both a 'spectrum' and a 'bunch' object
        #
        if hasattr(self, 'FreeFree'):
            total += self.FreeFree['intensity']
        if hasattr(self, 'FreeBound'):
            total += self.FreeBound['intensity']
        if hasattr(self, 'TwoPhoton'):
            total += self.TwoPhoton['intensity']
        self.Total = total
        #
        #
        if self.NTempDen == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
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
