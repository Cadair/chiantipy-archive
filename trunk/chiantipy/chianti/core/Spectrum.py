import copy
from datetime import datetime
import numpy as np
import chianti
import chianti.data as chdata
import chianti.constants as const
import chianti.filters as chfilters
import chianti.util as util
import chianti.io as chio
#
import chianti.Gui as chGui
#
from ._IonTrails import _ionTrails

defaults = chdata.Defaults

class spectrum(_ionTrails):
    '''
    Calculate the emission spectrum as a function of temperature and density.

    includes elemental abundances and ionization equilibria

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 1000.  Other filters,
    such as gaussian, box and lorentz, are available in chianti.filters.  When using the box filter,
    the width should equal the wavelength interval to keep the units of the continuum and line
    spectrum the same.

    A selection of elements can be make with elementList a list containing the names of elements
    that are desired to be included, e.g., ['fe','ni']

    A selection of ions can be make with ionList containing the names of
    the desired lines in Chianti notation, i.e. C VI = c_6

    Both elementList and ionList can not be specified at the same time

    a minimum abundance can be specified so that the calculation can be speeded up by excluding
    elements with a low abundance. With solar photospheric abundances -

    setting minAbund = 1.e-4 will include H, He, C, O, Ne
    setting minAbund = 2.e-5 adds  N, Mg, Si, S, Fe
    setting minAbund = 1.e-6 adds  Na, Al, Ar, Ca, Ni

    Setting doContinuum =0 will skip the continuum calculation.

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density
    '''
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), elementList = 0, ionList = 0, minAbund=0, doContinuum=1, em = None, abundanceName=0, verbose=0, allLines=1):
        #
        t1 = datetime.now()
        # creates Intensity dict from first ion calculated
        setupIntensity = 0
        #
        masterlist = chdata.MasterList
        # use the ionList but make sure the ions are in the database
        if elementList:
            for i,  one in enumerate(elementList):
                elementList[i] = one.lower()
            alist = []
            for one in masterlist:
                stuff = util.convertName(one)
                if stuff['Element'] in  elementList:
                    alist.append(one)
            masterlist = alist
        elif ionList:
            alist=[]
            for one in ionList:
                if masterlist.count(one):
                    alist.append(one)
                else:
                    if verbose:
                        pstring = ' %s not in CHIANTI database'%(one)
                        print(pstring)
            masterlist = alist
        self.Defaults=defaults
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.EDensity = np.asarray(eDensity, 'float64')
        nDen = self.EDensity.size
        nTempDen = max([nTemp, nDen])
        if type(em) != type(None):
            if type(em) == type(None):
                if nTempDen > 1:
                    em = np.ones_like(self.Temperature)*em
                    nEm = nTempDen
                else:
                    nEm = 1
            else:
                em = np.asarray(em, 'float64')
                nEm = em.size
                if nEm != nTempDen:
                    print(' the emission measure array must be the same size as the temperature/density array')
                    return

        #self.AbundanceName = defaults['abundfile']
        #self.AbundanceAll = chdata.AbundanceAll
        #
        if abundanceName:
            if abundanceName in list(chdata.Abundance.keys()):
                self.AbundanceName = abundanceName
            else:
                abundChoices = list(chdata.Abundance.keys())
#                for one in wvl[topLines]:
#                    wvlChoices.append('%12.3f'%(one))
                abundChoice = chGui.gui.selectorDialog(abundChoices,label='Select Abundance name')
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                abundanceName = self.AbundanceName
                print((' Abundance chosen:  %s '%(self.AbundanceName)))
        else:
            self.AbundanceName = self.Defaults['abundfile']
        #
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        #
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        if minAbund < minAbundAll:
            minAbund = minAbundAll
        self.minAbund = minAbund
        ionInfo = chio.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        wvlRange = [wavelength.min(), wavelength.max()]
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
#        self.Intensity = {'ionS':[], 'lvl1':[], 'lvl2':[], 'wvl':np.ndarray, 'pretty1':np.ndarray, 'pretty2':np.ndarray, 'intensity':np.zeros((nTempDen, 0),'float64'), 'obs':np.ndarray }
        ionsCalculated = []
        #
        for iz in range(31):
            abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
            if abundance >= minAbund:
                print((' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance)))
                #
                for ionstage in range(1, iz+2):
                    ionS = util.zion2name(iz, ionstage)
#                   print ' ionS = ', ionS
                    masterListTest = ionS in masterlist
                    masterListInfoTest = ionS in list(ionInfo.keys())
                    if masterListTest or masterListInfoTest:
                        wvlTestMin = self.Wavelength.min() <= ionInfo[ionS]['wmax']
                        wvlTestMax = self.Wavelength.max() >= ionInfo[ionS]['wmin']
                        ioneqTest = (self.Temperature.max() >= ionInfo[ionS]['tmin']) and (self.Temperature.min() <= ionInfo[ionS]['tmax'])
                    # construct similar test for the dielectronic files
                    ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                    masterListTestD = ionSd in masterlist
                    masterListInfoTestD = ionSd in list(ionInfo.keys())
                    if masterListTestD or masterListInfoTestD:
                        wvlTestMinD = self.Wavelength.min() <= ionInfo[ionSd]['wmax']
                        wvlTestMaxD = self.Wavelength.max() >= ionInfo[ionSd]['wmin']
                        ioneqTestD = (self.Temperature.max() >= ionInfo[ionSd]['tmin']) and (self.Temperature.min() <=ionInfo[ionSd]['tmax'])
                    ionstageTest = ionstage > 1
                    if ionstageTest and ioneqTest and doContinuum:
                        # ionS is the target ion, cannot be the neutral for the continuum
                        if verbose:
                            print(' calculating continuum for :  %s'%(ionS))
                        cont = chianti.core.continuum(ionS, temperature, abundanceName=self.AbundanceName)
                        cont.freeFree(wavelength)
    #                   print dir(thisIon)
    #                   print ' wvl = ', thisIon.FreeFree['wvl']
                        if nTempDen ==1:
                            freeFree += cont.FreeFree['rate']
                        else:
                            for iTempDen in range(nTempDen):
                                freeFree[iTempDen] += cont.FreeFree['rate'][iTempDen]
                    #
                        cont.freeBound(wavelength)
                        if 'errorMessage' not in list(cont.FreeBound.keys()):
                            #  an fblvl file exists for this ions
                            if nTempDen == 1:
                                freeBound += cont.FreeBound['rate']
                            else:
                                for iTempDen in range(nTempDen):
                                    freeBound[iTempDen] += cont.FreeBound['rate'][iTempDen]
                    if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                        if verbose:
                            print(' calculating spectrum for  :  %s'%(ionS))
                        #
                        thisIon = chianti.core.ion(ionS, temperature, eDensity, abundanceName=self.AbundanceName)
                        ionsCalculated.append(ionS)
#                       print ' dir = ', dir(thisIon)
#                        thisIon.emiss(wvlRange = wvlRange, allLines=allLines)
                        thisIon.intensity(wvlRange = wvlRange, allLines=allLines)
#                        print(' intensity shape %5i %5i '%(thisIon.Intensity['intensity'].shape[0], thisIon.Intensity['intensity'].shape[1]))
                        # check that there are lines in this wavelength range
                        if 'errorMessage' not in  list(thisIon.Intensity.keys()):
                            thisIon.spectrum(wavelength, filter=filter)
                            if setupIntensity:
                                for akey in self.Intensity:
                                    self.Intensity[akey] = np.hstack((copy.copy(self.Intensity[akey]), thisIon.Intensity[akey]))
                            else:
                                setupIntensity = 1
                                print(' creating Intensity dict from ion %s'%(ionS))
                                self.Intensity  = thisIon.Intensity
#                           intensity = thisIon.Intensity['intensity']
                            if nTempDen == 1:
                                lineSpectrum += thisIon.Spectrum['intensity']
                            else:
                                for iTempDen in range(nTempDen):
                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
                        else:
                            print((' error with ion = %s'%(ionS)))
                        # get 2 photon emission for H and He sequences
                        if (iz - ionstage) in [0, 1]:
                            thisIon.twoPhoton(wavelength)
                            twoPhoton += thisIon.TwoPhoton['rate']
                    # get dielectronic lines
                    if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                        if verbose:
                            print(' calculating spectrum for  :  ', ionSd)
                        #
                        thisIon = chianti.core.ion(ionSd, temperature, eDensity, abundanceName=self.AbundanceName)
                        ionsCalculated.append(ionSd)
#                       print ' dir = ', dir(thisIon)
#                       have to do all lines for the dielectronic satellites
#                        thisIon.emiss(allLines=1)
                        thisIon.intensity(wvlRange = wvlRange, allLines=allLines)
                        # check that there are lines in this wavelength range - probably not redundant
                        if 'errorMessage' not in  list(thisIon.Intensity.keys()):
                            thisIon.spectrum(wavelength, filter=filter)
                            if setupIntensity:
                                for akey in self.Intensity:
                                    self.Intensity[akey] = np.hstack((self.Intensity[akey], thisIon.Intensity[akey]))
                            else:
                                setupIntensity = 1
                                self.Intensity = thisIon.Intensity
                            if nTempDen == 1:
                                lineSpectrum += thisIon.Spectrum['intensity']
                            else:
                                for iTempDen in range(nTempDen):
                                    lineSpectrum[iTempDen] += thisIon.Spectrum['intensity'][iTempDen]
        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
        #
        total = freeFree + freeBound + lineSpectrum + twoPhoton
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = ', dt.seconds)
        if type(em) != type(None):
            if nEm == 1:
                integrated = total*em
            else:
                integrated = np.zeros_like(wavelength)
                for iTempDen in range(nTempDen):
                    integrated += total[iTempDen]*em[iTempDen]
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'ions':ionsCalculated, 'Abundance':self.AbundanceName}
        else:
            self.Spectrum ={'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'ions':ionsCalculated, 'Abundance':self.AbundanceName}
#    #
#    # ---------------------------------------------------------------------------
#    #
#    def intensityList(self, index=-1,  wvlRange=None, wvlRanges=None,   top=10, relative=0, outFile=0 ):
#        '''
#        List the line intensities
#
#        wvlRange, a 2 element tuple, list or array determines the wavelength range
#
#        Top specifies to plot only the top strongest lines, default = 10
#
#        normalize = 1 specifies whether to normalize to strongest line, default = 0
#        rewrite of emissList
#        this has been directly copied from ion -- not the right way to do this
#        '''
#        #
#        #
#        #
#        if not hasattr(self, 'Intensity'):
#            try:
#                self.intensity()
#            except:
#                print(' intensities not calculated and emiss() is unable to calculate them')
#                print(' perhaps the temperature and/or eDensity are not set')
#                return
#        #
#        # everything in self.Intensity should be a numpy array
#        #
#        intens = copy.copy(self.Intensity)
#        intensity = intens['intensity']
#        ionS = intens['ionS']
#        wvl = intens['wvl']
#        lvl1 = intens['lvl1']
#        lvl2 = intens['lvl2']
#        pretty1 = intens['pretty1']
#        pretty2 = intens['pretty2']
#        obs = intens['obs']
#        avalue = intens['avalue']
#        #
#        temperature = self.Temperature
#        eDensity = self.EDensity
#        #
#            #
#        ndens = eDensity.size
#        ntemp = temperature.size
#        #
#        if ndens == 1 and ntemp == 1:
#            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity)
#            tstr = ' -  T = %10.2e (K)' %(temperature)
#            print(dstr+tstr)
#        elif ndens == 1 and ntemp > 1:
#            if index < 0:
#                index = ntemp/2
#            print('using index = %5i specifying temperature =  %10.2e'%(index, temperature[index]))
#
#            self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
#            intensity=intensity[index]
#        elif ndens > 1 and ntemp == 1:
#            if index < 0:
#                index = ndens/2
#            print('using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index]))
#            self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
#            intensity=intensity[index]
#        elif ndens > 1 and ntemp > 1:
#            if index < 0:
#                index = ntemp/2
#            print('using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index]))
#            self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
#            intensity=intensity[index]
#        #
#        if wvlRange:
#            wvlIndex=util.between(wvl,wvlRange)
#        elif wvlRanges:
#            wvlIndex = []
#            for awvlRange in wvlRanges:
#                wvlIndex.extend(util.between(wvl,awvlRange))
#        else:
#            wvlIndex = list(range(wvl.size))
#        #
#        #
#        #  get lines in the specified wavelength range
#        #
#        intensity = intensity[wvlIndex]
#        ionS = ionS[wvlIndex]
#        wvl = wvl[wvlIndex]
#        lvl1 = lvl1[wvlIndex]
#        lvl2 = lvl2[wvlIndex]
#        avalue = avalue[wvlIndex]
#        pretty1 = pretty1[wvlIndex]
#        pretty2 = pretty2[wvlIndex]
#        obs = obs[wvlIndex]
#        #
#        self.Error = 0
#        if wvl.size == 0:
#            print('No lines in this wavelength interval')
#            self.Error = 1
#            self.Message = 'No lines in this wavelength interval'
#            return
#        #
#        elif top == 0:
#            top = wvl.size
#        elif top > wvl.size:
#            top = wvl.size
##
#        #
#        # sort by intensity
#        #
#        isrt = np.argsort(intensity)
#        #
#        ionS = ionS[isrt[-top:]]
#        wvl = wvl[isrt[-top:]]
#        lvl1 = lvl1[isrt[-top:]]
#        lvl2 = lvl2[isrt[-top:]]
#        obs = obs[isrt[-top:]]
#        intensity = intensity[isrt[-top:]]
#        avalue = avalue[isrt[-top:]]
#        pretty1 = pretty1[isrt[-top:]]
#        pretty2 = pretty2[isrt[-top:]]
#        #
#    # must follow setting top
#        #
#        if relative:
#            intensity = intensity/intensity[:top].max()
#        #
#        #
#        idx = np.argsort(wvl)
#        fmt = '%5s %5i %5i %25s - %25s %12.4f %12.3e %12.2e %1s'
#        print('   ')
#        print(' ------------------------------------------')
#        print('   ')
#        print(' Ion   lvl1  lvl2         lower                     upper                   Wvl(A)   Intensity       Obs')
#        for kdx in idx:
#            print(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], intensity[kdx], avalue[kdx], obs[kdx]))
#        print('   ')
#        print(' ------------------------------------------')
#        print('   ')
#        #
#        self.Intensity['wvlTop'] = wvl[idx]
#        self.Intensity['intensityTop'] = intensity[idx]
#        if outFile:
#            fmt = '%5s %5i %5i %25s - %25s %12.4f %12.3e %1s \n'
#            outpt = open(outFile, 'w')
#            outpt.write('Ion lvl1  lvl2         lower                       upper                   Wvl(A)   Intensity       Obs \n')
#            for kdx in idx:
#                outpt.write(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], intensity[kdx], avalue[kdx], obs[kdx]))
#            outpt.close()
#        return
