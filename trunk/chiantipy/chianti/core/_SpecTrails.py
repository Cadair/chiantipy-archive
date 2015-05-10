from datetime import datetime
import numpy as np
import chianti.filters as chfilters
import chianti.util as util
import chianti.io as chio 
import chianti.data as chdata
import chianti.constants as const
#
defaults = chdata.Defaults
#
class _specTrails():
    '''
    a collection of methods for use in spectrum calculations
    '''
    def __init__(self, temperature, density):
        self.Temperature = temperature
        self.EDensity = density
        self.AbundanceName = defaults['abundfile']
        self.AbundAll = chdata.Abundance[self.AbundanceName]['abundance']
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
        if not hasattr(self, 'IonInstances'):
            print(' must set keepIons=1 in order to keep self.IonInstances')
            return
        #
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
            self.Wavelength = wavelength
            nWvl = len(wavelength)
        lineSpectrum = np.zeros((self.NTempDen, nWvl), 'float64').squeeze()
        for akey in self.IonInstances.iterkeys():
            if verbose:
                print( ' trying ion = %s'%(akey))
#            thisIon = self.IonInstances[akey]
            if not 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                if verbose:
                    print(' doing convolve on ion %s '%(akey))
                self.IonInstances[akey].spectrum(wavelength, filter)
#                lineSpectrum = np.add(lineSpectrum, self.IonInstances[akey].Spectrum['intensity'])
                if 'errorMessage' in sorted(self.IonInstances[akey].Spectrum.keys()):
                    print(self.IonInstances[akey].Spectrum['errorMessage'])
                else:
                    lineSpectrum += self.IonInstances[akey].Spectrum['intensity']
#                if self.NTempDen == 1:
#                    lineSpectrum += thisIon.Spectrum['intensity']
#                else:
#                    for iTempDen in range(self.NTempDen):
#                        lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
            else:
                if 'errorMessage' in sorted(self.IonInstances[akey].Intensity.keys()):
                    print(self.IonInstances[akey].Intensity['errorMessage'])

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
        dt = t2 - t1
        print(' elapsed seconds = %12.3e'%(dt.seconds))
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':self.Em,  'Abundance':self.AbundanceName}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':self.Em,  'Abundance':self.AbundanceName}}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'Abundance':self.AbundanceName}
        return
        #
        # ---------------------------------------------------------------------------
        #
    def ionGate(self, elementList = 0, ionList = 0, minAbund=0, doContinuum=1, verbose = 0):
        '''
        creates a list of ions for free-free, free-bound, and line intensity calculations
        '''
        #
        masterlist = chdata.MasterList
        abundAll = self.AbundAll
        #
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        if minAbund:
            if minAbund < minAbundAll:
                minAbund = minAbundAll
        ionInfo = chio.masterListInfo()
        #
        if hasattr(self, 'Wavelength'):
            wvlRange = [self.Wavelength.min(), self.Wavelength.max()]
        elif hasattr(self, 'WvlRange'):
            wvlRange = self.WvlRange
        else:
            print(' need a wavelength range in ionGate ')
        #
        temperature = self.Temperature
        #
        # use the ionList but make sure the ions are in the database
        self.Todo = {}
        #
        if minAbund:
            if doContinuum:
                for iz in range(1, 31):
                    abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
                    if abundance >= minAbund:
                        if verbose:
                            print(' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance))
                        #
                        for ionstage in range(1, iz+2):
                            ionS = util.zion2name(iz, ionstage)
                            masterListTest = ionS in masterlist
                            masterListInfoTest = ionS in sorted(ionInfo.keys())
                            if masterListTest or masterListInfoTest:
                                wvlTestMin = wvlRange[0] <= ionInfo[ionS]['wmax']
                                wvlTestMax = wvlRange[1] >= ionInfo[ionS]['wmin']
                                ioneqTest = (temperature.max() >= ionInfo[ionS]['tmin']) and (temperature.min() <= ionInfo[ionS]['tmax'])
                            # construct similar test for the dielectronic files
                            ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                            masterListTestD = ionSd in masterlist
                            masterListInfoTestD = ionSd in sorted(ionInfo.keys())
                            if masterListTestD or masterListInfoTestD:
                                wvlTestMinD = wvlRange[0] <= ionInfo[ionSd]['wmax']
                                wvlTestMaxD = wvlRange[1] >= ionInfo[ionSd]['wmin']
                                ioneqTestD = (temperature.max() >= ionInfo[ionSd]['tmin']) and (temperature.min() <=ionInfo[ionSd]['tmax'])
                            ionstageTest = ionstage > 1
                            if ionstageTest and ioneqTest and doContinuum:
                                # ionS is the target ion, cannot be the neutral for the continuum
                                if verbose:
                                    print(' setting up continuum calculation for %s  '%(ionS))
                                if ionS in sorted(self.Todo.keys()):
                                    self.Todo[ionS] += '_ff'
                                else:
                                    self.Todo[ionS] = 'ff'
                                    if iz +1 != ionstage:
                                        self.Todo[ionS] += '_fb'
                                if verbose:
                                    print(' for ion %s do : %s'%(ionS, self.Todo[ionS]))
        #
        if elementList:
            for i,  one in enumerate(elementList):
                elementList[i] = one.lower()
            for one in masterlist:
                stuff = util.convertName(one)
                bare = stuff['Z'] == stuff['Ion']
                if stuff['Element'] in  elementList:
                    self.Todo[one] = 'line'
                    if doContinuum and not stuff['Dielectronic']:
                        self.Todo[one]+= '_ff'
                        if not bare:
                            self.Todo[one] += '_fb'
        if ionList: 
            for one in ionList:
                stuff = util.convertName(one)
                bare = stuff['Z'] == stuff['Ion']
                if masterlist.count(one):
                    self.Todo[one] = 'line'
                    if doContinuum and not stuff['Dielectronic']:
                        self.Todo[one]+= '_ff'
                        if not bare:
                            self.Todo[one] += '_fb'
                else:
                    if verbose:
                        pstring = ' %s not in CHIANTI database'%(one)
                        print(pstring)
        #
        #
        #
        if minAbund:
            for iz in range(1, 31):
                abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
                if abundance >= minAbund:
                    if verbose:
                        print(' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance))
                    #
                    for ionstage in range(1, iz+2):
                        ionS = util.zion2name(iz, ionstage)
                        masterListTest = ionS in masterlist
                        masterListInfoTest = ionS in sorted(ionInfo.keys())
                        if masterListTest or masterListInfoTest:
                            wvlTestMin = wvlRange[0] <= ionInfo[ionS]['wmax']
                            wvlTestMax = wvlRange[1] >= ionInfo[ionS]['wmin']
                            ioneqTest = (temperature.max() >= ionInfo[ionS]['tmin']) and (temperature.min() <= ionInfo[ionS]['tmax'])
                        # construct similar test for the dielectronic files
                        ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                        masterListTestD = ionSd in masterlist
                        masterListInfoTestD = ionSd in sorted(ionInfo.keys())
                        if masterListTestD or masterListInfoTestD:
                            wvlTestMinD = wvlRange[0] <= ionInfo[ionSd]['wmax']
                            wvlTestMaxD = wvlRange[1] >= ionInfo[ionSd]['wmin']
                            ioneqTestD = (temperature.max() >= ionInfo[ionSd]['tmin']) and (temperature.min() <=ionInfo[ionSd]['tmax'])
                            #
                        if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                            #if verbose:
                                #print(' setting up spectrum calculation for  %s'%(ionS))
                            if ionS in sorted(self.Todo.keys()):
                                self.Todo[ionS] += '_line'
                            else:
                                self.Todo[ionS] = 'line'
                        # get dielectronic lines
                            if verbose:
                                print(' for ion %s do : %s'%(ionS, self.Todo[ionS]))
                        if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                            #if verbose:
                                #print(' setting up  spectrum calculation for  %s '%(ionSd))
                            if ionSd in sorted(self.Todo.keys()):
                                self.Todo[ionSd] += '_line'
                            else:
                                self.Todo[ionSd] = 'line'
                            if verbose:
                                print(' for ion %s do : %s'%(ionSd, self.Todo[ionSd]))
        return


