from datetime import datetime
import copy
#
#
import numpy as np
import pylab as pl
import chianti.data as chdata
import chianti.constants as const
import chianti.filters as chfilters
import chianti.util as util
import chianti.io as chio
import chianti.Gui as chgui

#
from ._IonTrails import _ionTrails
from ._SpecTrails import _specTrails

try:
    from IPython import parallel
#    from chianti import mputil
    import chianti.ipymputil as mputil
except:
    print(' your version of IPython does not support multiprocessing \n you will not be able to use ipymspectrum')
#
defaults = chdata.Defaults
#
#
#class ipymspectrum(_ionTrails, _specTrails):
#    '''
#    this is the multiprocessing version of spectrum for using inside an IPython Qtconsole or notebook.
#
#    be for creating an instance, it is necessary to type something like the following into a console
#    > ipcluster start --profile=notebook  --n=3
#    this is the way to invoke things under the IPython 2.1 notation, a bit different in 2.0
#
#    Calculate the emission spectrum as a function of temperature and density.
#
#    temperature and density can be arrays but, unless the size of either is one (1),
#    the two must have the same size
#
#    the returned spectrum will be convolved with a filter of the specified width on the
#    specified wavelength array
#
#    the default filter is gaussianR with a resolving power of 100.  Other filters,
#    such as gaussian, box and lorentz, are available in chianti.filters.  When using the box filter,
#    the width should equal the wavelength interval to keep the units of the continuum and line
#    spectrum the same.
#
#    A selection of elements can be make with elementList a list containing the names of elements
#    that are desired to be included, e.g., ['fe','ni']
#
#    A selection of ions can be make with ionList containing the names of
#    the desired lines in Chianti notation, i.e. C VI = c_6
#
#    Both elementList and ionList can not be specified at the same time
#
#    a minimum abundance can be specified so that the calculation can be speeded up by excluding
#    elements with a low abundance. With solar photospheric abundances -
#
#    setting minAbund = 1.e-4 will include H, He, C, O, Ne
#    setting minAbund = 2.e-5 adds  N, Mg, Si, S, Fe
#    setting minAbund = 1.e-6 adds  Na, Al, Ar, Ca, Ni
#
#    Setting em will multiply the spectrum at each temperature by the value of em.
#
#    em [for emission measure], can be a float or an array of the same length as the
#    temperature/density.
#    allLines = 1 will include lines with either theoretical or observed wavelengths.  allLines=0 will
#    include only those lines with observed wavelengths
#
#    proc = the number of processors to use
#    timeout - a small but non-zero value seems to be necessary
#    '''
#    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), label=0, elementList = 0, ionList = 0, minAbund=0., doContinuum=1, allLines = 1, em = 1.,  proc=3,  abundanceName=0, verbose = 0,  timeout=0.1):
#        #
#        t1 = datetime.now()
#        #
#        rcff = parallel.Client()
#        rcfb = parallel.Client()
#        rcion = parallel.Client()
#        all_engines = rcff[:]
#        all_engines = rcfb[:]
#        all_engines = rcion[:]
#        lbvff = rcff.load_balanced_view()
#        lbvfb = rcfb.load_balanced_view()
#        lbvion = rcion.load_balanced_view()
#        #
#        # creates Intensity dict from first ion calculated
#        #
#        setupIntensity = 0
#        #
#        masterlist = chdata.MasterList
#        # use the ionList but make sure the ions are in the database
#        if elementList:
#            for i,  one in enumerate(elementList):
#                elementList[i] = one.lower()
#            alist = []
#            for one in masterlist:
#                stuff = util.convertName(one)
#                if stuff['Element'] in  elementList:
#                    alist.append(one)
#            masterlist = alist
#        elif ionList:
#            alist=[]
#            for one in ionList:
#                if masterlist.count(one):
#                    alist.append(one)
#                else:
#                    if verbose:
#                        pstring = ' %s not in CHIANTI database'%(one)
#                        print(pstring)
#            masterlist = alist
#        self.Defaults = defaults
#        self.Temperature = np.asarray(temperature, 'float64')
#        nTemp = self.Temperature.size
#        self.EDensity = np.asarray(eDensity, 'float64')
#        nDen = self.EDensity.size
#        nTempDen = max([nTemp, nDen])
#        self.NTempDen = nTempDen
##        em = np.asarray(em, 'float64')
##        self.Em = em
#        self.AllLines = allLines
#        #
#        if not abundanceName:
#            self.AbundanceName = self.Defaults['abundfile']
#        else:
#            if abundanceName in chdata.Abundance:
#                self.AbundanceName = abundanceName
#            else:
#                abundChoices = list(chdata.Abundance.keys())
##                for one in wvl[topLines]:
##                    wvlChoices.append('%12.3f'%(one))
#                abundChoice = chgui.gui.selectorDialog(abundChoices,label='Select Abundance name')
#                abundChoice_idx = abundChoice.selectedIndex
#                self.AbundanceName = abundChoices[abundChoice_idx[0]]
#                abundanceName = self.AbundanceName
#                print(' Abundance chosen:  %s '%(self.AbundanceName))
#        #
#        #
#        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
#        #
#        nonzed = abundAll > 0.
#        minAbundAll = abundAll[nonzed].min()
#        if minAbund < minAbundAll:
#            minAbund = minAbundAll
#        ionInfo = chio.masterListInfo()
#        wavelength = np.asarray(wavelength)
#        nWvl = wavelength.size
#        self.Wavelength = wavelength
##        wvlRange = [wavelength.min(), wavelength.max()]
#        #
#        #
#        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
#        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
#        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
#        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
#        #
#         #
#        ffInpt = []
#        fbInpt = []
#        ionInpt = []
#        #
#        self.IonsCalculated = []
#        self.IonInstances = {}
#        self.Finished = []
#        #
#        self.Todo = []
#        for iz in range(31):
#            abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
#            if abundance >= minAbund:
#                if verbose:
#                    print(' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance))
#                #
#                for ionstage in range(1, iz+2):
#                    ionS = util.zion2name(iz, ionstage)
#                    masterListTest = ionS in masterlist
#                    masterListInfoTest = ionS in sorted(ionInfo.keys())
#                    if masterListTest or masterListInfoTest:
#                        wvlTestMin = self.Wavelength.min() <= ionInfo[ionS]['wmax']
#                        wvlTestMax = self.Wavelength.max() >= ionInfo[ionS]['wmin']
#                        ioneqTest = (self.Temperature.max() >= ionInfo[ionS]['tmin']) and (self.Temperature.min() <= ionInfo[ionS]['tmax'])
#                    # construct similar test for the dielectronic files
#                    ionSd = util.zion2name(iz, ionstage, dielectronic=1)
#                    masterListTestD = ionSd in masterlist
#                    masterListInfoTestD = ionSd in sorted(ionInfo.keys())
#                    if masterListTestD or masterListInfoTestD:
#                        wvlTestMinD = self.Wavelength.min() <= ionInfo[ionSd]['wmax']
#                        wvlTestMaxD = self.Wavelength.max() >= ionInfo[ionSd]['wmin']
#                        ioneqTestD = (self.Temperature.max() >= ionInfo[ionSd]['tmin']) and (self.Temperature.min() <=ionInfo[ionSd]['tmax'])
#                    ionstageTest = ionstage > 1
#                    if ionstageTest and ioneqTest and doContinuum:
#                        # ionS is the target ion, cannot be the neutral for the continuum
#                        if verbose:
#                            print(' setting up continuum calculation for %s  '%(ionS))
#                        ffinpt.append([ionS, temperature, wavelength, abundance])
#                        fbinpt.append([ionS, temperature, wavelength, abundance])
#                        #
#                    if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
#                        if verbose:
#                            print(' setting up spectrum calculation for  %s'%(ionS))
#                        ionInpt.append([ionS, temperature, eDensity, wavelength, filter, allLines, abundance, em])
#                        self.Todo.append(ionS)
#                        self.IonsCalculated.append(ionS)
#                    # get dielectronic lines
#                    if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
#                        if verbose:
#                            print(' setting up  spectrum calculation for  %s '%(ionSd))
##                        dielWorkerQ.put((ionSd, temperature, density, wavelength, filter))
#                        # set allLines fo dielectronic
#                        ionInpt.append([ionSd, temperature, eDensity, wavelength, filter, allLines, abundance, em])
#                        self.Todo.append(ionSd)
#                        self.IonsCalculated.append(ionS)
#        #
#        if doContinuum:
#            for anInpt in ffInpt:
#                lvbff.apply(mputil.doFf, anInpt)
#            lvbff.wait()
#            lvbff.get_result()
#            if verbose:
#                print(' got ff result')
#            #
#            for iff in range(ffInpt):
#                thisFreeFree = list(lvbff.results.values()[iff])
#                if nTempDen ==1:
#                    freeFree += thisFreeFree['rate']
#                else:
#                    for iTempDen in range(nTempDen):
#                        freeFree[iTempDen] += thisFreeFree['rate'][iTempDen]
#        #
#            for anInpt in fbInpt:
#                lvbfb.apply(mputil.doFb, anInpt)
#            lvbfb.wait()
#            lvbfb.get_result()
#            if verbose:
#                print(' got fb result')
#            #
#            for ifb in range(fbInpt):
#                thisFreeBound = list(lvbfb.results.values()[ifb])
#                if nTempDen ==1:
#                    freeBound += thisFreeBoound['rate']
#                else:
#                    for iTempDen in range(nTempDen):
#                        freeBound[iTempDen] += thisFreeBound['rate'][iTempDen]
#        #
#        #
#        for anInpt in ionInpt:
#            lbvion.apply(mputil.doIon, anInpt)
#        lbvion.wait()
#        lbvion.get_result()
#        if verbose:
#            print(' got ion results')
#        #
#        #
#        #
##        for out in lbvion.results.values():
#        for ijk in range(len(list(lbvion.results.values()))):
#            out = list(lbvion.results.values())[ijk]
#            ions = out[0]
#            if verbose:
#                print('type of result %5i = %s %s'%(ijk, type(out[0]), type(out[1])))
#                if type(out[0]) == str:
#                    print('out[0] = %s'%(out[0]))
#                if type(out[1]) == str:
#                    print('out[1] = %s'%(out[1]))
##                print(' collecting calculation for %s'%(ions))
#            thisIon = out[1]
##            thisSpectrum = thisIon.Spectrum
##            thisIntensity = thisIon.Intensity
#            if not 'errorMessage' in sorted(thisIntensity.keys()):
#                self.Finished.append(ions)
#                self.IonInstances[ions] = copy.deepcopy(thisIon)
##                self.IonInstances.append(copy.deepcopy(thisIon))
#                if setupIntensity:
#                    for akey in sorted(self.Intensity.keys()):
#                        self.Intensity[akey] = np.hstack((self.Intensity[akey], thisIntensity[akey]))
#                else:
#                    setupIntensity = 1
#                    self.Intensity  = thisIntensity
#                #
#                lineSpectrum += thisSpectrum['intensity']
##                if nTempDen == 1:
##                    lineSpectrum += thisSpectrum['intensity']
##                else:
##                    for iTempDen in range(nTempDen):
##                        lineSpectrum[iTempDen] += thisSpectrum['intensity'][iTempDen]
#               # check for two-photon emission
#                if len(out) == 3:
#                    tp = out[2]
#                    if nTempDen == 1:
#                        twoPhoton += tp['rate']
#                    else:
#                        for iTempDen in range(nTempDen):
#                            twoPhoton[iTempDen] += tp['rate'][iTempDen]
#            else:
#                if 'errorMessage' in sorted(thisIntensity.keys()):
#                    print(thisIntensity['errorMessage'])
#        #
#        #
#        #
#        #
#        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
#        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
#        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
#        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
#        #
#        total = freeFree + freeBound + lineSpectrum + twoPhoton
#        #
#        t2 = datetime.now()
#        dt=t2-t1
#        print(' elapsed seconds = %12.3e'%(dt.seconds))
#        rcff.purge_results('all')
#        rcfb.purge_results('all')
#        rcion.purge_results('all')
#        #
#        if nEm == 0:
#            integrated = total
#        else:
#            if nEm == 1:
#                integrated = total*em
#            else:
#                integrated = np.zeros_like(wavelength)
#                for iTempDen in range(nTempDen):
#                    integrated += total[iTempDen]*em[iTempDen]
#        #
#        if type(label) == type(''):
#            if hasattr(self, 'Spectrum'):
#                print(' hasattr = true')
#                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}
#            else:
#                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}}
#        else:
#            self.Spectrum = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}
#
    #
    # -------------------------------------------------------------------------
    #
#
#
class ipymspectrum(_ionTrails, _specTrails):
    '''
    this is the multiprocessing version of spectrum for using inside an IPython Qtconsole or notebook.

    be for creating an instance, it is necessary to type something like the following into a console
    > ipcluster start --profile=notebook  --n=3
    this is the way to invoke things under the IPython 2.1 notation, a bit different in 2.0

    Calculate the emission spectrum as a function of temperature and density.

    temperature and density can be arrays but, unless the size of either is one (1),
    the two must have the same size

    the returned spectrum will be convolved with a filter of the specified width on the
    specified wavelength array

    the default filter is gaussianR with a resolving power of 100.  Other filters,
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

    Setting em will multiply the spectrum at each temperature by the value of em.

    em [for emission measure], can be a float or an array of the same length as the
    temperature/density.
    allLines = 1 will include lines with either theoretical or observed wavelengths.  allLines=0 will
    include only those lines with observed wavelengths

    proc = the number of processors to use
    timeout - a small but non-zero value seems to be necessary
    '''
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), label=0, elementList = 0, ionList = 0, minAbund=0., doContinuum=1, allLines = 1, em = 1.,  proc=3,  abundanceName=0, verbose = 0,  timeout=0.1):
        #
        t1 = datetime.now()
        #
        rcff = parallel.Client()
        rcfb = parallel.Client()
        rcion = parallel.Client()
        all_engines = rcff[:]
        all_engines = rcfb[:]
        all_engines = rcion[:]
        lbvff = rcff.load_balanced_view()
        lbvfb = rcfb.load_balanced_view()
        lbvion = rcion.load_balanced_view()
        #
        #
        # creates Intensity dict from first ion calculated
        #
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
        self.Defaults = defaults
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.EDensity = np.asarray(eDensity, 'float64')
        nDen = self.EDensity.size
        nTempDen = max([nTemp, nDen])
        self.NTempDen = nTempDen
        em = np.asarray(em, 'float64')
        if em.size != nTempDen:
            if em.size == 1:
                em = np.ones(nTempDen, 'float64')*em[0]
        self.Em = em
        self.AllLines = allLines
        #
        if not abundanceName:
            self.AbundanceName = self.Defaults['abundfile']
        else:
            if abundanceName in chdata.Abundance:
                self.AbundanceName = abundanceName
            else:
                abundChoices = list(chdata.Abundance.keys())
#                for one in wvl[topLines]:
#                    wvlChoices.append('%12.3f'%(one))
                abundChoice = chgui.gui.selectorDialog(abundChoices,label='Select Abundance name')
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                abundanceName = self.AbundanceName
                print(' Abundance chosen:  %s '%(self.AbundanceName))
        #
        #
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        #
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        if minAbund < minAbundAll:
            minAbund = minAbundAll
        ionInfo = chio.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
#        wvlRange = [wavelength.min(), wavelength.max()]
        #
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
         #
        ffInpt = []
        fbInpt = []
        ionInpt = []
        #
        self.IonsCalculated = []
        self.IonInstances = {}
        self.Finished = []
        #
        self.Todo = []
        for iz in range(31):
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
                        wvlTestMin = self.Wavelength.min() <= ionInfo[ionS]['wmax']
                        wvlTestMax = self.Wavelength.max() >= ionInfo[ionS]['wmin']
                        ioneqTest = (self.Temperature.max() >= ionInfo[ionS]['tmin']) and (self.Temperature.min() <= ionInfo[ionS]['tmax'])
                    # construct similar test for the dielectronic files
                    ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                    masterListTestD = ionSd in masterlist
                    masterListInfoTestD = ionSd in sorted(ionInfo.keys())
                    if masterListTestD or masterListInfoTestD:
                        wvlTestMinD = self.Wavelength.min() <= ionInfo[ionSd]['wmax']
                        wvlTestMaxD = self.Wavelength.max() >= ionInfo[ionSd]['wmin']
                        ioneqTestD = (self.Temperature.max() >= ionInfo[ionSd]['tmin']) and (self.Temperature.min() <=ionInfo[ionSd]['tmax'])
                    ionstageTest = ionstage > 1
                    if ionstageTest and ioneqTest and doContinuum:
                        # ionS is the target ion, cannot be the neutral for the continuum
                        if verbose:
                            print(' setting up continuum calculation for %s  '%(ionS))
                        ffInpt.append([ionS, temperature, wavelength, abundance])
                        fbInpt.append([ionS, temperature, wavelength, abundance])
                        #
                    if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                        if verbose:
                            print(' setting up spectrum calculation for  %s'%(ionS))
                        ionInpt.append([ionS, temperature, eDensity, wavelength, filter, allLines, abundance, em, doContinuum])
                        self.Todo.append(ionS)
                        self.IonsCalculated.append(ionS)
                    # get dielectronic lines
                    if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                        if verbose:
                            print(' setting up  spectrum calculation for  %s '%(ionSd))
#                        dielWorkerQ.put((ionSd, temperature, density, wavelength, filter))
                        # set allLines fo dielectronic
                        ionInpt.append([ionSd, temperature, eDensity, wavelength, filter, allLines, abundance, em, doContinuum])
                        self.Todo.append(ionSd)
                        self.IonsCalculated.append(ionS)
        #
        if doContinuum:
            for anInpt in ffInpt:
                lbvff.apply(mputil.doFf, anInpt)
            lbvff.wait()
            lbvff.get_result()
            if verbose:
                print(' got ff result')
            #
#            for iff in range(len(ffInpt)):
            for iff in range(len(list(lbvff.results.values()))):
                aFreeFree = list(lbvff.results.values()[iff])
                ionS = aFreeFree[0]
                thisFfCont = aFreeFree[1]
                if 'errorMessage' not in sorted(thisFfCont.FreeFree.keys()):
                    if nTempDen == 1:
                        freeFree += thisFfCont.FreeFree['rate']*em[0]
                    else:
                        for iTempDen in range(nTempDen):
                            freeFree[iTempDen] += thisFfCont.FreeFree['rate'][iTempDen]*em[iTempDen]
                else:
                    print(thisFfCont.FreeFree['errorMessage'])
        #
            for anInpt in fbInpt:
                lbvfb.apply(mputil.doFb, anInpt)
            lbvfb.wait()
            lbvfb.get_result()
            if verbose:
                print(' got fb result')
            #
#            for ifb in range(fbInpt):
            for ifb in range(len(list(lbvfb.results.values()))):
                aFreeBound = list(lbvfb.results.values()[ifb])
                ionS = aFreeBound[0]
                thisFbCont = aFreeBound[1]
                print(' fb ion = %s'%(ionS))
                for akey in sorted(thisFbCont.FreeBound.keys()):
                    print('FreeBound key = %s'%(akey))
                if 'errorMessage' not in sorted(thisFbCont.FreeBound.keys()):
                    if nTempDen == 1:
                        freeBound += thisFbCont.FreeBound['rate']*em[0]
                    else:
                        for iTempDen in range(nTempDen):
                            freeBound[iTempDen] += thisFbCont.FreeBound['rate'][iTempDen]*em[iTempDen]
                else:
                    print(thisFbCont.FreeBound['errorMessage'])
        #
        #
        for anInpt in ionInpt:
            lbvion.apply(mputil.doIon, anInpt)
        lbvion.wait()
        lbvion.get_result()
        if verbose:
            print(' got ion results')
        #
        #
        #
#        for out in lbvion.results.values():
        for ijk in range(len(list(lbvion.results.values()))):
            out = list(lbvion.results.values())[ijk]
            ions = out[0]
            print(' collecting calculation for %s'%(ions))
#            if verbose:
#                print('type of result %5i = %s %s'%(ijk, type(out[0]), type(out[1])))
#                if type(out[0]) == str:
#                    print('out[0] = %s'%(out[0]))
#                if type(out[1]) == str:
#                    print('out[1] = %s'%(out[1]))
#                print(' collecting calculation for %s'%(ions))
            thisIon = out[1]
#            self.Finished.append(ions)
#            thisSpectrum = thisIon.Spectrum
#            thisIntensity = thisIon.Intensity
            if not 'errorMessage' in sorted(thisIon.Intensity.keys()):
                self.Finished.append(ions)
                self.IonInstances[ions] = copy.deepcopy(thisIon)
                thisIntensity = thisIon.Intensity
##                self.IonInstances.append(copy.deepcopy(thisIon))
                if setupIntensity:
                    for akey in sorted(self.Intensity.keys()):
                        self.Intensity[akey] = np.hstack((self.Intensity[akey], thisIntensity[akey]))
                else:
                    setupIntensity = 1
                    self.Intensity  = thisIntensity
                #
                lineSpectrum += thisIon.Spectrum['intensity']
#                if nTempDen == 1:
#                    lineSpectrum += thisSpectrum['intensity']
#                else:
#                    for iTempDen in range(nTempDen):
#                        lineSpectrum[iTempDen] += thisSpectrum['intensity'][iTempDen]
               # check for two-photon emission
                if len(out) == 3:
                    tp = out[2]
                    if nTempDen == 1:
                        twoPhoton += tp['rate']
                    else:
                        for iTempDen in range(nTempDen):
                            twoPhoton[iTempDen] += tp['rate'][iTempDen]
            else:
                if 'errorMessage' in sorted(thisIon.Intensity.keys()):
                    print(thisIon.Intensity['errorMessage'])
        #
        #
        #
        #
        self.FreeFree = {'wavelength':wavelength, 'intensity':freeFree.squeeze()}
        self.FreeBound = {'wavelength':wavelength, 'intensity':freeBound.squeeze()}
        self.LineSpectrum = {'wavelength':wavelength, 'intensity':lineSpectrum.squeeze()}
        self.TwoPhoton = {'wavelength':wavelength, 'intensity':twoPhoton.squeeze()}
        #
        total = freeFree + freeBound + lineSpectrum + twoPhoton
        #
        t2 = datetime.now()
        dt=t2-t1
        print(' elapsed seconds = %12.3e'%(dt.seconds))
        rcff.purge_results('all')
        rcfb.purge_results('all')
        rcion.purge_results('all')
        #
        if nTempDen == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                print(' hasattr = true')
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}}
        else:
            self.Spectrum = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'ions':self.IonsCalculated, 'Abundance':self.AbundanceName}

    #
    # -------------------------------------------------------------------------
    #
    def lineSpectrumPlot(self, saveFile=0, plotContinuum=0, linLog = 'lin'):
        '''
        to plot the spectrum as a function of wavelength
        '''
        # must follow setting top
        #
        pl.figure()
        ylabel = 'Intensity'
        #
        xlabel = 'Wavelength ('+self.Defaults['wavelength'] +')'
        #
#        ymin = 10.**(np.log10(emiss.min()).round(0))
        #
        pl.ion()
        #
        pl.plot(self.LineSpectrum['wavelength'], self.LineSpectrum['intensity'])
        pl.xlabel(xlabel)
        pl.ylabel(ylabel)
        if saveFile:
            pl.savefig(saveFile)
    #
    # -------------------------------------------------------------------------
    #
