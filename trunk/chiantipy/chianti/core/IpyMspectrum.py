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
#

try:
    from IPython import parallel
#    from chianti import mputil
    import chianti.ipymputil as mputil
except:
    print ' your version of Python does not support multiprocessing \n you will not be able to use mspectrum'
#
defaults = chdata.Defaults
#
#
class ipymspectrum:
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
    def __init__(self, temperature, eDensity, wavelength, filter=(chfilters.gaussianR, 1000.), elementList = 0, ionList = 0, minAbund=0., doContinuum=1, allLines = 1, em = None,  proc=3,  abund=0, verbose = 0,  timeout=0.1):
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
        if not abund:
            self.AbundanceName = self.Defaults['abundfile']
        else:
            if abund in chdata.Abundance.keys():
                self.AbundanceName = abund
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
        #
        freeFree = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        freeBound = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        twoPhoton = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        lineSpectrum = np.zeros((nTempDen, nWvl), 'float64').squeeze()
        #
        #  free-free multiprocessing setup
#        ffWorkerQ = mp.Queue()
#        ffDoneQ = mp.Queue()
        #
        #  free-bound multiprocessing setup
        #
#        fbWorkerQ = mp.Queue()
#        fbDoneQ = mp.Queue()
        #
        #  ion multiprocessing setup
#        ionWorkerQ = mp.Queue()
#        ionDoneQ = mp.Queue()
        #
        ffInpt = []
        fbInpt = []
        ionInpt = []
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
                        ffinpt.append([ionS, temperature, wavelength, abund])
                        fbinpt.append([ionS, temperature, wavelength, abund])
                        #
                    if masterListTest and wvlTestMin and wvlTestMax and ioneqTest:
                        if verbose:
                            print ' setting up spectrum calculation for  :  ', ionS
                        ionInpt.append([ionS, temperature, eDensity, wavelength, filter, allLines, abund])
                        self.Todo.append(ionS)
                        ionsCalculated.append(ionS)
                    # get dielectronic lines
                    if masterListTestD and wvlTestMinD and wvlTestMaxD and ioneqTestD:
                        if verbose:
                            print ' setting up  spectrum calculation for  :  ', ionSd
#                        dielWorkerQ.put((ionSd, temperature, density, wavelength, filter))
                        # set allLines fo dielectronic
                        ionInpt.append([ionSd, temperature, eDensity, wavelength, filter, 1, abund])
                        self.Todo.append(ionSd)
                        ionsCalculated.append(ionS)
        #
        if doContinuum:
            for anInpt in ffInpt:
                lvbff.apply(mputil.doFf, anInpt)
            lvbff.wait()
            lvbff.get_result()
            if verbose:
                print' got ff result'
            #
            for iff in range(ffInpt):
                thisFreeFree = lvbff.results.values()[iff]
                if nTempDen ==1:
                    freeFree += thisFreeFree['rate']
                else:
                    for iTempDen in range(nTempDen):
                        freeFree[iTempDen] += thisFreeFree['rate'][iTempDen]
        #
            for anInpt in fbInpt:
                lvbfb.apply(mputil.doFb, anInpt)
            lvbfb.wait()
            lvbfb.get_result()
            if verbose:
                print' got fb result'
            #
            for ifb in range(fbInpt):
                thisFreeBound = lvbfb.results.values()[ifb]
                if nTempDen ==1:
                    freeBound += thisFreeBoound['rate']
                else:
                    for iTempDen in range(nTempDen):
                        freeBound[iTempDen] += thisFreeBound['rate'][iTempDen]
        #
        #
        for anInpt in ionInpt:
            lbvion.apply(mputil.doIon, anInpt)
        lbvion.wait()
        lbvion.get_result()
        if verbose:
            print' got ion result'
        #
        #
        #
#        for out in lbvion.results.values():
        for ijk in range(len(lbvion.results.values())):
            out = lbvion.results.values()[ijk]
            ions = out[0]
            if verbose:
                print(' collecting calculation for %s'%(ions))
            aspectrum = out[1]
            if not 'errorMessage' in out[2].keys():
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
        rcff.purge_results('all')
        rcfb.purge_results('all')
        rcion.purge_results('all')
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
        # -------------------------------------------------------------------------------------
        #
    def intensityRatio(self,wvlRange=None, wvlRanges=None,top=10):
        """
        Plot the ratio of 2 lines or sums of lines.
        Shown as a function of density and/or temperature.
        For a single wavelength range, set wvlRange = [wMin, wMax]
        For multiple wavelength ranges, set wvlRanges = [[wMin1,wMax1],[wMin2,wMax2], ...]
        A plot of relative emissivities is shown and then a dialog appears for the user to
        choose a set of lines.
        """
        #
        #        self.Emiss={"temperature":temperature,"density":density,"wvl":wvl,"emiss":em,
        #        "plotLabels":plotLabels}
        #
        if hasattr(self, 'Intensity'):
            doIntensity = False
            intens = self.Intensity
        else:
            doIntensity = True
        #
        #
        if doIntensity:
            # new values of temperature or eDensity
            self.intensity()
            intens = self.Intensity
        #
        #
        fontsize=14
        #
        temperature = self.Temperature
        eDensity = self.EDensity
        intensity = intens['intensity']
        ionS = intens['ionS']
        wvl = intens["wvl"]
#        plotLabels = intens["plotLabels"]
#        xLabel = plotLabels["xLabel"]
#        yLabel = plotLabels["yLabel"]
        #
        # find which lines are in the wavelength range if it is set
        #
        #
        if wvlRange:
            igvl=util.between(wvl,wvlRange)
        elif wvlRanges:
            igvl = []
            for awvlRange in wvlRanges:
                igvl.extend(util.between(wvl,awvlRange))
        else:
            igvl=range(len(wvl))
        #
        nlines=len(igvl)
        #
#        print ' nlines = ',nlines
#        print ' iglv = ',igvl
        igvl=np.take(igvl,wvl[igvl].argsort())
        # find the top most intense lines
        #
        if top > nlines:
            top=nlines
            #
        maxIntens = np.zeros(nlines,'Float64')
        for iline in range(nlines):
            maxIntens[iline] = intensity[:, igvl[iline]].max()
        for iline in range(nlines):
            if maxIntens[iline]==maxIntens.max():
                maxAll=intensity[:, igvl[iline]]
        line=range(nlines)
        igvlsort=np.take(igvl,np.argsort(maxIntens))
#        print 'igvlsort = ', igvlsort
        topLines=igvlsort[-top:]
#        print ' topLines = ', topLines
        maxWvl='%5.3f' % wvl[topLines[-1]]
        maxline=topLines[-1]
        #
        topLines=topLines[wvl[topLines].argsort()]
        #
        #
        # need to make sure there are no negative values before plotting
        good = intensity > 0.
        intensMin = intensity[good].min()
        bad = intensity <= 0.
        intensity[bad] = intensMin
        #
        #
        ntemp=self.Temperature.size
        #
        ndens=self.EDensity.size
        #
        ylabel='Emissivity relative to '+maxWvl
#        title=self.Spectroscopic
        title = ' No title'
        #
        #
        if ndens==1 and ntemp==1:
            print ' only a single temperature and eDensity'
            return
        elif ndens == 1:
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            outTemperature=self.Temperature
            outDensity=np.zeros(ntemp,'Float64')
            outDensity.fill(self.EDensity)
            desc_str=' at  Density = %10.2e (cm)$^{-3}$' % self.EDensity
        elif ntemp == 1:
            xvalues=self.EDensity
            outTemperature=np.zeros(ndens,'Float64')
            outTemperature.fill(self.Temperature)
            outDensity=self.EDensity
            xlabel=r'$\rm{Electron Density (cm)^{-3}}$'
            desc_str=' at Temp = %10.2e (K)' % self.Temperature
        else:
            outTemperature=self.Temperature
            outDensity=self.EDensity
            xlabel='Temperature (K)'
            xvalues=self.Temperature
            desc_str=' for variable Density'
        #
        # put all actual plotting here
        #
        pl.ion()
#        if chInteractive:
#            pl.ion()
#        else:
#            pl.ioff()
        #
        #  maxAll is an array
        ymax = np.max(intensity[:, topLines[0]]/maxAll)
        ymin = ymax
        pl.figure()
        ax = pl.subplot(111)
        nxvalues=len(xvalues)
        for iline in range(top):
            tline=topLines[iline]
            pl.loglog(xvalues,intensity[:, tline]/maxAll)
            if np.min(intensity[:, tline]/maxAll) < ymin:
                ymin = np.min(intensity[:, tline]/maxAll)
            if np.max(intensity[:, tline]/maxAll) > ymax:
                ymax = np.max(intensity[:, tline]/maxAll)
            skip=2
            start=divmod(iline,nxvalues)[1]
            for ixvalue in range(start,nxvalues,nxvalues/skip):
                pl.text(xvalues[ixvalue], intensity[ixvalue, tline]/maxAll[ixvalue], str(wvl[tline]))
        pl.xlim(xvalues.min(),xvalues.max())
#        pl.ylim(ymin, ymax)
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel(ylabel,fontsize=fontsize)
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,title, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(eDensity,intensity[:, topLines[top-1]]/maxAll, visible=False)
            ax2.xaxis.tick_top()
            pl.ylim(ymin/1.2, 1.2*ymax)
        else:
            pl.ylim(ymin/1.2, 1.2*ymax)
            pl.title(title+desc_str,fontsize=fontsize)
        pl.draw()
        #  need time to let matplotlib finish plotting
        time.sleep(0.5)
        #
        # get line selection
        #
        selectTags = []
        for itop in topLines:
            selectTags.append(ionS[itop]+ ' '+ str(wvl[itop]))
        #
        numden = gui.choice2Dialog(wvl[topLines])
        #
        # num_idx and den_idx are tuples
        #
        num_idx=numden.numIndex
        if len(num_idx) == 0:
            print ' no numerator lines were selected'
            return
        #
        den_idx=numden.denIndex
        if len(den_idx) == 0:
            print ' no denominator lines were selected'
            return
        #
        numIntens=np.zeros(len(xvalues),'Float64')
        for aline in num_idx:
            numIntens += intensity[:, topLines[aline]]
        #
        denIntens = np.zeros(len(xvalues),'Float64')
        for aline in den_idx:
            denIntens += intensity[:, topLines[aline]]
        #
        # plot the desired ratio
        #  maxAll is an array
        pl.figure()
        ax = pl.subplot(111)
        pl.loglog(xvalues,numIntens/denIntens)
        pl.xlim(xvalues.min(),xvalues.max())
        pl.xlabel(xlabel,fontsize=fontsize)
        pl.ylabel('Ratio ('+self.Defaults['flux']+')',fontsize=fontsize)
        desc = title + ':'
        for aline in num_idx:
            desc += ' ' + str(wvl[topLines[aline]])
        desc +=' / '
        for aline in den_idx:
            desc += ' ' + str(wvl[topLines[aline]])
        if ndens == ntemp and ntemp > 1:
            pl.text(0.07, 0.5,desc, horizontalalignment='left', verticalalignment='center', fontsize=fontsize,  transform = ax.transAxes)
            #
            ax2 = pl.twiny()
            xlabelDen=r'Electron Density (cm$^{-3}$)'
            pl.xlabel(xlabelDen, fontsize=fontsize)
            pl.loglog(eDensity,numIntens/denIntens, visible=False)
            ax2.xaxis.tick_top()
        else:
#            pl.ylim(ymin, ymax)
            pl.title(desc,fontsize=fontsize)
#       desc=title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str
#        pl.title(desc, fontsize=fontsize)
#       pl.title(title+' '+str(wvl[num_line])+' / '+str(wvl[den_line])+' '+desc_str,fontsize=fontsize)
#        pl.draw()
#        pl.ioff()
#        pl.show()
        #
        intensityRatioFileName=self.IonStr
        for aline in num_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='_2'
        for aline in den_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='.rat'
        self.IntensityRatio={'ratio':numIntens/denIntens,'desc':desc,
                'temperature':outTemperature,'eDensity':outDensity,'filename':intensityRatioFileName, 'numIdx':num_idx, 'denIdx':den_idx}
        #
        # -------------------------------------------------------------------------------------
        #
    def intensityRatioSave(self,outFile=''):
        '''Save the intensity ratio to a file.

        The intensity ratio as a function to temperature and eDensity is saved to an asciii file.

        Descriptive information is included at the top of the file.'''
        if outFile == '':
            outfile=self.IntensityRatio['filename']
#            if chInteractive:
            print ' saving ratio to filename = ',outfile
        if hasattr(self, 'IntensityRatio'):
            temperature=self.IntensityRatio['temperature']
            eDensity=self.IntensityRatio['eDensity']
            ratio=self.IntensityRatio['ratio']
            out=open(outFile,'w')
            nvalues=len(ratio)
            #
            #  need to add 7 lines to maintain IDL like files
            #
            out.write(outFile+'\n')    #1
            out.write(self.IntensityRatio['desc']+'\n') #2
            out.write(' created with ChiantiPy version '+ chdata.__version__ +'\n')   #3
            out.write(' columns are temperature, eDensity, ratio'+'\n')  #5
            tunit = 'K'
            out.write(' temperature in '+tunit+', electron eDensity in cm^(-3)'+'\n')  #6
            out.write(' ratio given in '+self.Defaults['flux']+'\n')   #4
            out.write(' '+'\n') #7
            for ivalue in range(nvalues):
                s='%12.3e %12.3e  %12.3e ' % (temperature[ivalue],eDensity[ivalue],ratio[ivalue])
                out.write(s+os.linesep)
            out.close()
        else:
#            if chInteractive:
            print ' in .intensityRatioSave(), no IntensityRatio is found'
        intensityRatioFileName=self.IonStr
        for aline in num_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='_2'
        for aline in den_idx:
            intensityRatioFileName+= '_%3i'%(wvl[topLines[aline]])
        intensityRatioFileName+='.rat'
        self.IntensityRatio={'ratio':numEmiss/denEmiss,'desc':desc,
                'temperature':outTemperature,'eDensity':outDensity,'filename':intensityRatioFileName, 'numIdx':num_idx, 'denIdx':den_idx}
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
