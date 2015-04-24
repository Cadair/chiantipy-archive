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
    #
    # -------------------------------------------------------------------------
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
        rcAll = parallel.Client()
        all_engines = rcAll[:]
        lbvAll = rcAll.load_balanced_view()
        #
        #
        # creates Intensity dict from first ion calculated
        #
        setupIntensity = 0
        #
        masterlist = chdata.MasterList
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
                em = np.ones(nTempDen, 'float64')*em
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
        self.AbundAll = abundAll
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
        allInpt = []
        #
        self.IonInstances = {}
        # ionGate creates the self.Todo list
        #
        self.ionGate(elementList = elementList, ionList = ionList, minAbund=minAbund, doContinuum=doContinuum, verbose = verbose)
        #
        for akey in sorted(self.Todo.keys()):
            zStuff = util.convertName(akey)
            Z = zStuff['Z']
            abundance = chdata.Abundance[self.AbundanceName]['abundance'][Z - 1]
            if verbose:
                print(' %5i %5s abundance = %10.2e '%(Z, const.El[Z-1],  abundance))
            if verbose:
                print(' doing ion %s for the following processes %s'%(akey, self.Todo[akey]))
            if 'ff' in self.Todo[akey]:
#                if verbose:
#                    print(' doing ff')
                allInpt.append([akey, 'ff', temperature, wavelength, abundance])
            if 'fb' in self.Todo[akey]:
#                if verbose:
#                    print(' doing fb')
                allInpt.append([akey, 'fb', temperature, wavelength, abundance])
            if 'line' in self.Todo[akey]:
#                if verbose:
#                    print(' doing line')
                allInpt.append([akey, 'line', temperature, eDensity, wavelength, filter, allLines, abundance, em, doContinuum])
        #
        for anInpt in allInpt:
            lbvAll.apply(mputil.doAll, anInpt)
        lbvAll.wait()
        lbvAll.get_result()
        if verbose:
            print(' got all ff, fb, line results')
        #
        for ijk in range(len(list(lbvAll.results.values()))):
            out = list(lbvAll.results.values())[ijk]
            ionS = out[0]
            print(' collecting calculation for %s'%(ionS))
            calcType = out[1]
            print(' processing %s results'%(calcType))
            #
            if calcType == 'ff':
                thisFf = out[2]
                if 'errorMessage' not in sorted(thisFf.keys()):
                    if nTempDen == 1:
                        freeFree += thisFf['rate']*em[0]
                    else:
                        for iTempDen in range(nTempDen):
                            freeFree[iTempDen] += thisFf['rate'][iTempDen]*em[iTempDen]
                elif type(thisFf) == str:
                    print(' error in FfCont %s'%(thisFf))
                else:
                    print(thisFf['errorMessage'])
            #
            elif calcType == 'fb':
                thisFbCont = out[2]
                print(' fb ion = %s'%(ionS))
                if hasattr(thisFbCont, 'FreeBound'):
                    if 'errorMessage' not in sorted(thisFbCont.FreeBound.keys()):
                        if nTempDen == 1:
                            freeBound += thisFbCont.FreeBound['rate']*em[0]
                        else:
                            for iTempDen in range(nTempDen):
                                freeBound[iTempDen] += thisFbCont.FreeBound['rate'][iTempDen]*em[iTempDen]
                    else:
                        print(thisFbCont.FreeBound['errorMessage'])
            #
            elif calcType == 'line':
                thisIon = out[2]
                if not 'errorMessage' in sorted(thisIon.Intensity.keys()):
                    self.IonInstances[ionS] = thisIon
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
                    if len(out) == 4:
                        tp = out[3]
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
        rcAll.purge_results('all')
        #
        if nTempDen == 1:
            integrated = total
        else:
            integrated = total.sum(axis=0)
        #
        if type(label) == type(''):
            if hasattr(self, 'Spectrum'):
                print(' hasattr = true')
                self.Spectrum[label] = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em,  'Abundance':self.AbundanceName}
            else:
                self.Spectrum = {label:{'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'Abundance':self.AbundanceName}}
        else:
            self.Spectrum = {'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'Abundance':self.AbundanceName}

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
