from datetime import datetime
import copy
#

#
import numpy as np
import chianti.data as chdata
import chianti.constants as const
import chianti.filters as chfilters
import chianti.util as util
#
#chInteractive = chdata.chInteractive
#if chInteractive:
#    import pylab as pl
#else:
#    import matplotlib
#    matplotlib.use('Agg')
#    import matplotlib.pyplot as pl

try:
    import multiprocessing as mp
    from chianti import mputil
    import chianti.mputil as mputil
except:
    print ' your version of Python does not support multiprocessing \n you will not be able to use mspectrum'
#
defaults = chdata.Defaults
#
# the following is necessary to make chiantipy non interactive for the web
#try:
#    chInteractive = int(os.environ['CHIANTIPY_INTERACTIVE'])
#except:
#    chInteractive = 1
#
class mspectrum:
    ''' this is the multiprocessing version of spectrum
    set proc to the desired number of processors, default=3

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
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), elementList = 0, ionList = 0, minAbund=0., abundanceName=0,  doContinuum=1, allLines = 1, em = None,  proc=3, verbose = 0,  timeout=0.1):
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
                        print('')
            masterlist = alist
        self.Defaults = defaults
        self.Temperature = np.asarray(temperature, 'float64')
        nTemp = self.Temperature.size
        self.EDensity = np.asarray(eDensity, 'float64')
        nDen = self.EDensity.size
        nTempDen = max([nTemp, nDen])
        if type(em) != type(None):
            if isinstance(em, float):
                if nTempDen > 1:
                    em = np.ones_like(self.Temperature)*em
                    nEm = nTempDen
                else:
                    nEm = 1
            else:
                em = np.asarray(em, 'float64')
                nEm = em.size
                if nEm != nTempDen:
                    print ' the emission measure array must be the same size as the temperature/density array'
                    return
            self.Em = em
        self.AllLines = allLines
        #
        if not abundName:
            self.AbundanceName = self.Defaults['abundfile']
        else:
            if abundName in chdata.Abundance.keys():
                self.AbundanceName = abundName
            else:
                abundChoices = chdata.Abundance.keys()
#                for one in wvl[topLines]:
#                    wvlChoices.append('%12.3f'%(one))
                abundChoice = gui.selectorDialog(abundChoices,label='Select Abundance name')
                abundChoice_idx = abundChoice.selectedIndex
                self.AbundanceName = abundChoices[abundChoice_idx[0]]
                abund = self.AbundanceName
                print(' Abundance chosen:  %s '%(self.AbundanceName))
        #
        abundAll = chdata.Abundance[self.AbundanceName]['abundance']
        #
        nonzed = abundAll > 0.
        minAbundAll = abundAll[nonzed].min()
        if minAbund < minAbundAll:
            minAbund = minAbundAll
        ionInfo = util.masterListInfo()
        wavelength = np.asarray(wavelength)
        nWvl = wavelength.size
        self.Wavelength = wavelength
        wvlRange = [wavelength.min(), wavelength.max()]
        #
        proc = min([proc, mp.cpu_count()])
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
        #  free-free multiprocessing setup
        ffWorkerQ = mp.Queue()
        ffDoneQ = mp.Queue()
        #
        #  free-bound multiprocessing setup
        #
        fbWorkerQ = mp.Queue()
        fbDoneQ = mp.Queue()
        #
        #  ion multiprocessing setup
        ionWorkerQ = mp.Queue()
        ionDoneQ = mp.Queue()
        #
        ionsCalculated = []
        #
        self.Todo = []
        for iz in range(31):
            abundance = chdata.Abundance[self.AbundanceName]['abundance'][iz-1]
            if abundance >= minAbund:
                if verbose:
                    print ' %5i %5s abundance = %10.2e '%(iz, const.El[iz-1],  abundance)
                #
                for ionstage in range(1, iz+2):
                    ionS = util.zion2name(iz, ionstage)
                    masterListTest = ionS in masterlist
                    masterListInfoTest = ionS in ionInfo.keys()
                    if masterListTest or masterListInfoTest:
                        wvlTestMin = self.Wavelength.min() <= ionInfo[ionS]['wmax']
                        wvlTestMax = self.Wavelength.max() >= ionInfo[ionS]['wmin']
                        ioneqTest = (self.Temperature.max() >= ionInfo[ionS]['tmin']) and (self.Temperature.min() <= ionInfo[ionS]['tmax'])
                    # construct similar test for the dielectronic files
                    ionSd = util.zion2name(iz, ionstage, dielectronic=1)
                    masterListTestD = ionSd in masterlist
                    masterListInfoTestD = ionSd in ionInfo.keys()
                    if masterListTestD or masterListInfoTestD:
                        wvlTestMinD = self.Wavelength.min() <= ionInfo[ionSd]['wmax']
                        wvlTestMaxD = self.Wavelength.max() >= ionInfo[ionSd]['wmin']
                        ioneqTestD = (self.Temperature.max() >= ionInfo[ionSd]['tmin']) and (self.Temperature.min() <=ionInfo[ionSd]['tmax'])
                    ionstageTest = ionstage > 1
                    if ionstageTest and ioneqTest and doContinuum:
                        # ionS is the target ion, cannot be the neutral for the continuum
                        if verbose:
                            print ' setting up continuum calculation for :  ',  ionS
                        ffWorkerQ.put((ionS, temperature, wavelength, abund))
                        fbWorkerQ.put((ionS, temperature, wavelength, abund))
#                        fbInputs.append([ionS, temperature, wavelength])
                        #
                    if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                        if verbose:
                            print ' setting up spectrum calculation for  :  ', ionS
                        ionWorkerQ.put((ionS, temperature, eDensity, wavelength, filter, allLines, abund))
                        self.Todo.append(ionS)
                        ionsCalculated.append(ionS)
                    # get dielectronic lines
                    if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                        if verbose:
                            print ' setting up  spectrum calculation for  :  ', ionSd
#                        dielWorkerQ.put((ionSd, temperature, density, wavelength, filter))
                        # set allLines fo dielectronic
                        ionWorkerQ.put((ionSd, temperature, eDensity, wavelength, filter, 1, abund))
                        self.Todo.append(ionSd)
                        ionsCalculated.append(ionS)
        #
        ffWorkerQSize = ffWorkerQ.qsize()
        fbWorkerQSize = fbWorkerQ.qsize()
        ionWorkerQSize = ionWorkerQ.qsize()
        if doContinuum:
            ffProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFfQ, args=(ffWorkerQ, ffDoneQ))
                p.start()
                ffProcesses.append(p)
    #       timeout is not necessary
            for p in ffProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
#            for i in range(proc):
#                ffProcesses.append('STOP')
            #
            for iff in range(ffWorkerQSize):
                thisFreeFree =ffDoneQ.get()
                if nTempDen ==1:
                    freeFree += thisFreeFree['rate']
                else:
                    for iTempDen in range(nTempDen):
                        freeFree[iTempDen] += thisFreeFree['rate'][iTempDen]
            for p in ffProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
            fbProcesses = []
            for i in range(proc):
                p = mp.Process(target=mputil.doFbQ, args=(fbWorkerQ, fbDoneQ))
                p.start()
                fbProcesses.append(p)
    #       timeout is not necessary
            for p in fbProcesses:
                if p.is_alive():
                    p.join(timeout=timeout)
#            for i in range(proc):
#                fbProcesses.append('STOP')
            #
            for ifb in range(fbWorkerQSize):
                thisFreeBound = fbDoneQ.get()
                if thisFreeBound.has_key('rate'):
                    if nTempDen ==1:
                        freeBound += thisFreeBound['rate']
                    else:
                        for iTempDen in range(nTempDen):
                            freeBound[iTempDen] += thisFreeBound['rate'][iTempDen]
            for p in fbProcesses:
                if not isinstance(p, str):
                    p.terminate()
        #
        ionProcesses = []
        if ionWorkerQSize < proc:
            nproc = ionWorkerQSize
        for i in range(proc):
            p = mp.Process(target=mputil.doIonQ, args=(ionWorkerQ, ionDoneQ))
            p.start()
            ionProcesses.append(p)
#            ionWorkerQ.put('STOP')
#       timeout is not necessary
        for p in ionProcesses:
#            print' process is alive:  ', p.is_alive()
            if p.is_alive():
#                p.join()
                p.join(timeout=timeout)
#        for i in range(proc):
#            ionProcesses.append('STOP')
        self.Finished = []
        #
        for ijk in range(ionWorkerQSize):
            out = ionDoneQ.get()
            ions = out[0]
            if verbose:
                print(' collecting calculation for %s'%(ions))
            aspectrum = out[1]
            if not 'errorMessage' in out[2].keys():
                self.Finished.append(ions)
                try:
                    if setupIntensity:
                        for akey in self.Intensity.keys():
                            self.Intensity[akey] = np.hstack((copy.copy(self.Intensity[akey]), out[2][akey]))
                    else:
                        setupIntensity = 1
                        self.Intensity  = out[2]
                    #
                    if nTempDen == 1:
                        lineSpectrum += aspectrum['intensity']
                    else:
                        for iTempDen in range(nTempDen):
                            lineSpectrum[iTempDen] += aspectrum['intensity'][iTempDen]
                   # check for two-photon emission
                    if len(out) == 4:
                        tp = out[3]
                        if nTempDen == 1:
                            twoPhoton += tp['rate']
                        else:
                            for iTempDen in range(nTempDen):
                                twoPhoton[iTempDen] += tp['rate'][iTempDen]
                except:
                    print('  error with ion %s'%(ions))
                    for akey in out[2]:
                        print(akey)
        #
        for p in ionProcesses:
            if not isinstance(p, str):
                p.terminate()
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
        print ' elapsed seconds = ', dt.seconds
        #
        if type(em) != type(None):
            if nEm == 1:
                integrated = total*em
            else:
                integrated = np.zeros_like(wavelength)
                for iTempDen in range(nTempDen):
                    integrated += total[iTempDen]*em[iTempDen]
            self.Spectrum ={'temperature':temperature, 'eDensity':eDensity, 'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'integrated':integrated, 'em':em, 'minAbund':minAbund, 'masterlist':masterlist, 'ions':ionsCalculated, 'Abundance':self.AbundanceName}
        else:
            self.Spectrum ={'temperature':temperature, 'eDensity':eDensity, 'wavelength':wavelength, 'intensity':total.squeeze(), 'filter':filter[0].__name__,   'width':filter[1], 'minAbund':minAbund, 'masterlist':masterlist, 'ions':ionsCalculated, 'abundance':self.AbundanceName}
    #
    # ---------------------------------------------------------------------------
    #
    def intensityList(self, index=-1,  wvlRange=None, wvlRanges=None,   top=10, relative=0, outFile=0 ):
        '''
        List the line intensities

        wvlRange, a 2 element tuple, list or array determines the wavelength range

        Top specifies to plot only the top strongest lines, default = 10

        normalize = 1 specifies whether to normalize to strongest line, default = 0
        rewrite of emissList
        this has been directly copied from ion -- not the right way to do this
        '''
        #
        #
        #
        if not hasattr(self, 'Intensity'):
            print ' intensities not calculated and emiss() is unable to calculate them'
            print ' perhaps the temperature and/or eDensity are not set'
            return
        #
        # everything in self.Intensity should be a numpy array
        #
        intens = copy.copy(self.Intensity)
        intensity = intens['intensity']
        ionS = intens['ionS']
        wvl = intens['wvl']
        lvl1 = intens['lvl1']
        lvl2 = intens['lvl2']
        pretty1 = intens['pretty1']
        pretty2 = intens['pretty2']
        obs = intens['obs']
        avalue = intens['avalue']
        #
        temperature = self.Temperature
        eDensity = self.EDensity
        #
            #
        ndens = eDensity.size
        ntemp = temperature.size
        #
        if ndens == 1 and ntemp == 1:
            dstr = ' -  Density = %10.2e (cm$^{-3}$)' %(eDensity)
            tstr = ' -  T = %10.2e (K)' %(temperature)
        elif ndens == 1 and ntemp > 1:
            if index < 0:
                index = ntemp/2
            print 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
            self.Message = 'using index = %5i specifying temperature =  %10.2e'%(index, temperature[index])
            intensity=intensity[index]
        elif ndens > 1 and ntemp == 1:
            if index < 0:
                index = ndens/2
            print 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
            self.Message = 'using index =%5i specifying eDensity = %10.2e'%(index, eDensity[index])
            intensity=intensity[index]
        elif ndens > 1 and ntemp > 1:
            if index < 0:
                index = ntemp/2
            print 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
            self.Message = 'using index = %5i specifying temperature = %10.2e, eDensity =  %10.2e'%(index, temperature[index], eDensity[index])
            intensity=intensity[index]
        #
        if wvlRange:
            wvlIndex=util.between(wvl,wvlRange)
        elif wvlRanges:
            wvlIndex = []
            for awvlRange in wvlRanges:
                wvlIndex.extend(util.between(wvl,awvlRange))
        else:
            wvlIndex = range(wvl.size)
        #
        #
        #  get lines in the specified wavelength range
        #
        intensity = intensity[wvlIndex]
        ionS = ionS[wvlIndex]
        wvl = wvl[wvlIndex]
        lvl1 = lvl1[wvlIndex]
        lvl2 = lvl2[wvlIndex]
        avalue = avalue[wvlIndex]
        pretty1 = pretty1[wvlIndex]
        pretty2 = pretty2[wvlIndex]
        obs = obs[wvlIndex]
        #
        self.Error = 0
        if wvl.size == 0:
            print 'No lines in this wavelength interval'
            self.Error = 1
            self.Message = 'No lines in this wavelength interval'
            return
        #
        elif top == 0:
            top = wvl.size
        elif top > wvl.size:
            top = wvl.size
#
        #
        # sort by intensity
        #
        isrt = np.argsort(intensity)
        #
        ionS = ionS[isrt[-top:]]
        wvl = wvl[isrt[-top:]]
        lvl1 = lvl1[isrt[-top:]]
        lvl2 = lvl2[isrt[-top:]]
        obs = obs[isrt[-top:]]
        intensity = intensity[isrt[-top:]]
        avalue = avalue[isrt[-top:]]
        pretty1 = pretty1[isrt[-top:]]
        pretty2 = pretty2[isrt[-top:]]
        #
    # must follow setting top
        #
        if relative:
            intensity = intensity/intensity[:top].max()
        #
        #
        idx = np.argsort(wvl)
        fmt = '%5s %5i %5i %25s - %25s %12.4f %12.3e %12.2e %1s'
        print '   '
        print ' ------------------------------------------'
        print '   '
        print ' Ion   lvl1  lvl2         lower                     upper                   Wvl(A)   Intensity       Obs'
        for kdx in idx:
            print(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], intensity[kdx], avalue[kdx], obs[kdx]))
        print '   '
        print ' ------------------------------------------'
        print '   '
        #
        self.Intensity['wvlTop'] = wvl[idx]
        self.Intensity['intensityTop'] = intensity[idx]
        if outFile:
            fmt = '%5s %5i %5i %25s - %25s %12.4f %12.3e %1s \n'
            outpt = open(outFile, 'w')
            outpt.write('Ion lvl1  lvl2         lower                       upper                   Wvl(A)   Intensity       Obs \n')
            for kdx in idx:
                outpt.write(fmt%(ionS[kdx], lvl1[kdx], lvl2[kdx], pretty1[kdx], pretty2[kdx], wvl[kdx], intensity[kdx], avalue[kdx], obs[kdx]))
            outpt.close()
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
