import chianti.util as util
class _archivalTrails():
    '''
    a set of methods for using archival data sets
    '''
        #
        # -------------------------------------------------------------------------
        #
    def setupScups(self, dir=0, verbose=0):
        '''
        if ion is initiated with setup=0, this allows the setup to be done at a later point
        perhaps, more importantly,  by setting dir to a directory cotaining the necessary files
        for a ChiantiPy ion, it allows one to setup an ion with files not in the current
        Chianti directory
        '''
        #
        # read in all data if in masterlist
        #  if not, there should still be ionization and recombination rates
        #
        MasterList = chdata.MasterList
        #
        if self.IonStr in MasterList:
            if dir:
                fileName = os.path.join(dir, self.IonStr)
                self.Elvlc = io.elvlcRead('',filename=fileName+'.elvlc')
                self.Wgfa = io.wgfaRead('',filename=fileName+'.wgfa', elvlcname=fileName+'.elvlc')
                self.Nwgfa=len(self.Wgfa['lvl1'])
                nlvlWgfa = max(self.Wgfa['lvl2'])
                nlvlList =[nlvlWgfa]
    #                print 'fileName = ', fileName
                scupsfile = fileName + '.scups'
                # read the splups file
                if os.path.isfile(scupsfile):
                    # happens the case of fe_3 and prob. a few others
                    self.Scups = io.scupsRead('', filename=scupsfile)
                    self.Nscups=len(self.Scups['lvl1'])
                    nlvlScups = max(self.Scups['lvl2'])
                    nlvlList.append(nlvlScups)
                else:
                    self.Nscups = 0
                    nlvlScups = 0
            else:
                fileName = util.ion2filename(self.IonStr)
                self.Elvlc = io.elvlcRead(self.IonStr)
                self.Wgfa = io.wgfaRead(self.IonStr)
                self.Nwgfa=len(self.Wgfa['lvl1'])
                nlvlWgfa = max(self.Wgfa['lvl2'])
                nlvlList =[nlvlWgfa]
    #                print 'fileName = ', fileName
                scupsfile = fileName + '.scups'
                if os.path.isfile(scupsfile):
                    # happens the case of fe_3 and prob. a few others
                    self.Scups = io.scupsRead(self.IonStr)
                    self.Nscups=len(self.Scups['lvl1'])
                    nlvlScups = max(self.Scups['lvl2'])
                    nlvlList.append(nlvlScups)
                else:
                    self.Nscups = 0
                    nlvlScups = 0
##                self.Nlvls = nlvlElvlc
            #
            file = fileName +'.cilvl'
            if os.path.isfile(file):
                self.Cilvl = io.cireclvlRead('',filename = fileName, cilvl=1)
                self.Ncilvl=len(self.Cilvl['lvl1'])
                nlvlCilvl = max(self.Cilvl['lvl2'])
                nlvlList.append(nlvlCilvl)
            else:
                self.Ncilvl = 0
            #  .reclvl file may not exist
            reclvlfile = fileName +'.reclvl'
            if os.path.isfile(reclvlfile):
                self.Reclvl = io.cireclvlRead('',filename=fileName, reclvl=1)
                self.Nreclvl = len(self.Reclvl['lvl1'])
                nlvlReclvl = max(self.Reclvl['lvl2'])
                nlvlList.append(nlvlReclvl)
            else:
                self.Nreclvl = 0
            #  .dielsplups files not longer exist
#            dielsplupsfile = fileName +'.splups'
#            if self.Dielectronic and os.path.isfile(dielsplupsfile):
#                self.DielSplups = io.splupsRead('', filename=dielsplupsfile, diel=1)
#                self.Ndielsplups=len(self.DielSplups["lvl1"])
#                nlvlDielSplups = max(self.DielSplups['lvl2'])
#                nlvlList.append(nlvlDielSplups)
#            else:
#                self.Ndielsplups = 0
            #
            #  psplups file may not exist
            psplupsfile = fileName +'.psplups'
            if os.path.isfile(psplupsfile):
                self.Psplups = io.splupsRead('', filename=psplupsfile,  prot=True)
                self.Npsplups=len(self.Psplups["lvl1"])
            else:
                self.Npsplups = 0
            #
            drparamsFile = fileName +'.drparams'
            if os.path.isfile(drparamsFile):
                self.DrParams = io.drRead(self.IonStr)
            #
            rrparamsFile = fileName +'.rrparams'
            if os.path.isfile(rrparamsFile):
                self.RrParams = io.rrRead(self.IonStr)

            #  not needed for ion, only phion
#                photoxfile = util.ion2filename(self.IonStr)+'.photox'
#                if os.path.isfile(photoxfile):
#                    self.Photox = util.photoxRead(self.IonStr)
            #
            # need to determine the number of levels that can be populated
            nlvlElvlc = len(self.Elvlc['lvl'])
#                print ' nlvlElvlc = ', nlvlElvlc
#                print ' other nlvls = ',  nlvlList
#                nlvlWgfa = max(self.Wgfa['lvl2'])
            #  elvlc file can have more levels than the rate level files
            self.Nlvls = min([nlvlElvlc, max(nlvlList)])
        else:
            try:
                self.Elvlc = io.elvlcRead(self.IonStr, verbose=verbose)
            except:
                print(' the ion %s is not in the CHIANTI masterlist '%(self.IonStr))
                print(' elvlc file NOT available for %s'%(self.IonStr))
                return
        #
        # -------------------------------------------------------------------------------------
        #
    def populateScups(self, popCorrect=1, verbose=0, **kwargs):
        """
        Calculate level populations for specified ion.
        possible keyword arguments include temperature, eDensity, pDensity, radTemperature and rStar
        """
        #
        #
        for one in kwargs.keys():
            if one not in chdata.keywordArgs:
                print(' following keyword is not understood - %20s '%(one))
        #
        nlvls=self.Nlvls
        nwgfa=self.Nwgfa
        nscups=self.Nscups
        npsplups=self.Npsplups
        #
        if 'temperature' in kwargs.keys():
            self.Temperature = np.asarray(kwargs['temperature'])
            temperature = self.Temperature
        elif hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
                print(' no temperature values have been set')
                return
        #
        if 'eDensity' in kwargs.keys():
            self.EDensity = np.asarray(kwargs['eDensity'])
            eDensity = self.EDensity
        elif hasattr(self, 'EDensity'):
            eDensity = self.EDensity
        else:
            print(' no eDensity values have been set')
            return
        #
        if 'pDensity' in kwargs.keys():
            if kwargs['pDensity'] == 'default':
                self.p2eRatio()
                protonDensity = self.ProtonDensityRatio*self.EDensity
            else:
                try:
                    self.PDensity = np.asarray(kwargs['pDensity'])
                except:
                    print(' could not interpret value for keyword pDensity')
                    print(' should be either "default" or a number or array')
                    return
        else:
            if hasattr(self, 'PDensity'):
                protonDensity = self.PDensity
            else:
                self.p2eRatio()
                self.PDensity = self.ProtonDensityRatio*self.EDensity
                protonDensity = self.PDensity
                print(' proton density not specified, set to \"default\" ')
        #
        if 'radTemperature' in kwargs.keys() and 'rStar' in kwargs.keys():
            self.RadTemperature = np.asarray(kwargs['radTemperature'])
            radTemperature = np.array(self.RadTemperature)
            self.RStar = np.asarray(kwargs['rStar'])
            rStar = np.asarray(self.RStar)
        elif hasattr(self, 'RadTemperature') and hasattr(self, 'RStar'):
            radTemperature = self.RadTemperature
            rStar = self.RStar
        #
        #
        # the Dielectronic test should eventually go away
        rec = 0
        ci = 0
        if popCorrect and (not self.Dielectronic):
            if self.Ncilvl:
                ci = 1
                cilvl = self.Cilvl
                if hasattr(self, 'CilvlRate'):
                    cilvlRate = self.CilvlRate
                else:
                    self.cireclvlDescale('cilvl')
                    cilvlRate = self.CilvlRate
                self.recombRate()
                #
                lowers = util.zion2name(self.Z, self.Ion-1)
                # get the lower ionization stage
                lower = ion(lowers, temperature=self.Temperature, eDensity = self.EDensity)
                lower.ionizRate()
                # need to get multiplicity of lower ionization stage
                lowMult = lower.Elvlc['mult']
#            else:
#                ci = 0
#            try:
            rec = 0
            if self.Nreclvl:
                rec = 1
                reclvl = self.Reclvl
                if hasattr(self, 'ReclvlRate'):
                    reclvlRate = self.ReclvlRate
                else:
                    self.cireclvlDescale('reclvl')
                    reclvlRate = self.ReclvlRate
#            if self.Ndielsplups:
#                rec = 1
#                if hasattr(self, 'DielUpsilon'):
#                    dielexRate = self.DielUpsilon['exRate']
#                else:
#                    print(' doing upsilonDescale')
#                    self.upsilonDescale(diel=1)
#                    dielexRate = self.DielUpsilon['exRate']

#            except:
#                self.Reclvl = io.cireclvlRead(self.IonStr,'reclvl' )
#                reclvl = self.Reclvl
#                self.reclvlDescale()
#                if type(self.Reclvl) != type(None):
#                    rec = 1
#                    self.ionizRate()
#                    #  get the higher ionization stage
#                    highers = util.zion2name(self.Z, self.Ion+1)
##                   print ' highers = ', highers
#                    higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
#                    higher.recombRate()
#                else:
#                    rec = 0
            #
#        elif self.Dielectronic:
##            self.Ndielsplups and self.Dielectronic:
#            rec = 0
#            if hasattr(self, 'DielUpsilon'):
#                dielexRate = self.DielUpsilon['exRate']
#            else:
##                print(' doing upsilonDescale')
#                self.upsilonDescaleSplups(diel=1)
#                dielexRate = self.DielUpsilon['exRate']
         #
        if rec:
#            if self.Ndielsplups:
#                self.upsilonDescale(diel=1)
#                dielexRate = self.DielUpsilon['exRate']
            # get ionization rate of this iion
            self.ionizRate()
            #  get the higher ionization stage
            highers = util.zion2name(self.Z, self.Ion+1)
            higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
            higher.recombRate()
#        print ' nlvls, ci, rec = ', nlvls, ci, rec
        #
        rad=np.zeros((nlvls+ci+rec,nlvls+ci+rec),"float64")  #  the populating matrix for radiative transitions
        #
        #
        for iwgfa in range(nwgfa):
            l1 = self.Wgfa["lvl1"][iwgfa]-1
            l2 = self.Wgfa["lvl2"][iwgfa]-1
            rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]
            rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]
            # photo-excitation and stimulated emission
            if self.RadTemperature:
                if not self.RStar:
                    dilute = 0.5
                else:
                    dilute = util.dilute(self.RStar)
                # next - don't include autoionization lines
                if abs(self.Wgfa['wvl'][iwgfa]) > 0.:
                    de = const.invCm2Erg*(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                    dekt = de/(const.boltzmann*self.RadTemperature)
                    # photoexcitation
                    phexFactor = dilute*(float(self.Elvlc['mult'][l2])/float(self.Elvlc['mult'][l1]))/(np.exp(dekt) -1.)
                    rad[l2+ci,l1+ci] += self.Wgfa["avalue"][iwgfa]*phexFactor
                    rad[l1+ci,l1+ci] -= self.Wgfa["avalue"][iwgfa]*phexFactor
                    # stimulated emission
                    stemFactor = dilute/(np.exp(-dekt) -1.)
                    rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]*stemFactor
                    rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]*stemFactor

        #
        #
        if self.Nscups:
            self.upsilonDescale(diel=self.Dielectronic)
            ups = self.Upsilon['upsilon']
            exRate = self.Upsilon['exRate']
            dexRate = self.Upsilon['dexRate']                
        #
        if npsplups:
            self.upsilonDescaleSplups(prot=1)
#            pups = self.PUpsilon['upsilon']
            pexRate = self.PUpsilon['exRate']
            pdexRate = self.PUpsilon['dexRate']
            
        #
        temp=temperature
        ntemp=temp.size
        #
        cc=const.collision*self.EDensity
        ndens=cc.size
        if npsplups:
            cp=const.collision*protonDensity
        if ntemp > 1 and ndens >1 and ntemp != ndens:
            print(' unless temperature or eDensity are single values')
            print(' the number of temperatures values must match the ')
            print(' the number of eDensity values')
            return
        #
        # get corrections for recombination and excitation
        #
        nscups = self.Nscups
        #
        #  first, for ntemp=ndens=1
        if ndens==1 and ntemp==1:
            popmat=np.copy(rad)
            if not self.Dielectronic:
                for iscups in range(0,nscups):
                    l1=self.Scups["lvl1"][iscups]-1
                    l2=self.Scups["lvl2"][iscups]-1
                    #
                    popmat[l1+ci,l2+ci] += self.EDensity*dexRate[iscups]
                    popmat[l2+ci,l1+ci] += self.EDensity*exRate[iscups]
                    popmat[l1+ci,l1+ci] -= self.EDensity*exRate[iscups]
                    popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[iscups]
                #
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                 #
                popmat[l1+ci,l2+ci] += self.PDensity*pdexRate[isplups]
                popmat[l2+ci,l1+ci] += self.PDensity*pexRate[isplups]
                popmat[l1+ci,l1+ci] -= self.PDensity*pexRate[isplups]
                popmat[l2+ci,l2+ci] -= self.PDensity*pdexRate[isplups]
           # now include ionization rate from
            if ci:
#                print ' ci = ', ci
                #
                # the ciRate can be computed for all temperatures
                #
                ciTot = 0.
                for itrans in range(len(cilvl['lvl1'])):
                    lvl1 = cilvl['lvl1'][itrans]-1
                    lvl2 = cilvl['lvl2'][itrans]-1
#                    de = cilvl['de'][itrans]
#                    ekt = (de*1.57888e+5)/temperature
#                    mult = lowMult[lvl1-1]
                    # this is kind of double booking the ionization rate components
                    popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans]
                    popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans]
                    ciTot += self.EDensity*self.CilvlRate['rate'][itrans]
                #
                popmat[1, 0] += (self.EDensity*lower.IonizRate['rate'] - ciTot)
                popmat[0, 0] -= (self.EDensity*lower.IonizRate['rate'] - ciTot)
                popmat[0, 1] += self.EDensity*self.RecombRate['rate']
                popmat[1, 1] -= self.EDensity*self.RecombRate['rate']
            if rec:
                #
#                print ' rec, dielTot  = ', rec,  dielTot
                #
                for itrans in range(self.Nreclvl):
#                    lvl1 = reclvl['lvl1'][itrans]-1
                    lvl2 = reclvl['lvl2'][itrans]-1
                    popmat[lvl2+ci, -1] += self.EDensity*reclvlRate['rate'][itrans]
                    popmat[-1, -1] -= self.EDensity*reclvlRate['rate'][itrans]
                if self.Nreclvl:
                    recTot = reclvlRate['rate'].sum(axis=0)
                else:
                    recTot = 0.

                #
                popmat[-1,  ci] += self.EDensity*self.IonizRate['rate']
                popmat[ci, ci] -= self.EDensity*self.IonizRate['rate']
                # next 2 line take care of overbooking
                popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate']- recTot)
                popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate']- recTot)
                #
            if self.Dielectronic:
                dielTot = 0.
                for iscups in range(0,nscups):
                    l1=self.Scups["lvl1"][iscups]-1
                    l2=self.Scups["lvl2"][iscups]-1
                    #
                    popmat[l2+ci,-1] += self.EDensity*exRate[iscups]
                    popmat[-1, -1] -= self.EDensity*exRate[iscups]
                # for dielectronic ions, l1 = ground level of the ion itself
#                    branch = np.zeros(self.Ndielsplups, 'float64')
#                    for isplups in range(0,self.Ndielsplups):
#                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
#                        l2 = self.DielSplups["lvl2"][isplups]-1
#                        auto = rad[l1+ci, l2+ci]
#                        avalue = rad[:, l2+ci]
#                        good = avalue > 0.
#                        avalueTot = avalue[good].sum()
#                        branch[isplups] = (avalueTot-auto)/avalueTot
##                        print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
#                        self.DielUpsilon['branch'] =  branch
#                    #
#                    dielTot = 0.
#                    print ' Ndielsplups > 0 '
#                for isplups in range(0,self.Ndielsplups):
#                    l1 = self.DielSplups["lvl1"][isplups]-1
#                    l2 = self.DielSplups["lvl2"][isplups]-1
#                     #
##                        print ' l1, l2, dielexRate = ', l1, l2, dielexRate[isplups]
#                    popmat[l2+ci,l1+ci] += self.EDensity*dielexRate[isplups]
#                    popmat[l1+ci,l1+ci] -= self.EDensity*dielexRate[isplups]
                    #
                dielTot += self.EDensity*exRate[isplups]*branch[isplups]
            else:
                dielTot = 0.

#                print ' higher, rec , dieltot = ',  self.EDensity*higher.RecombRate['rate'], self.EDensity*reclvlRate['rate'].sum(axis=0),  dielTot
            # normalize to unity
#            print(' rec =  %5i  ci = %5i'%(rec, ci))
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[nlvls+ci+rec-1] = 0.
            if self.Dielectronic:
                norm[nlvls-1] = 0.
            popmat[nlvls+ci+rec-1]=norm
#            popmata = np.copy(popmat)
#            popmata[nlvls+ci+rec-1]=norm
            #popmata[nlvls+ci+rec-1]=norm
            b=np.zeros(nlvls+ci+rec,'float64')
            b[nlvls+ci+rec-1]=1.
#            print ' norm = ', norm
#            print 'popmat, last line',  popmat[-1]
#            print ' b = ', b
#            popmat[nlvls/2]=norm
#            b=np.zeros(nlvls+ci+rec,'float64')
#            b[nlvls/2]=1.
#            if rec:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:ci+nlvls+rec-1]
#            else:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:]
#            fullpop = np.linalg.solve(popmat,b)
            try:
                fullpop=np.linalg.solve(popmat,b)
                pop = fullpop[ci:ci+nlvls]
            except np.linalg.LinAlgError:
                pop = np.zeros(nlvls, 'float64')
#                print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature)
            #
            # ----------------------------------------------------------------------------------
        #   next, in case of a single eDensity value
#            pop = np.linalg.solve(popmat,b)
        elif ndens == 1:
            pop=np.zeros((ntemp, nlvls),"float64")
#            pop=np.zeros((ntemp,ci + nlvls + rec),"float64")
            for itemp in range(0,ntemp):
                popmat=np.copy(rad)
                for iscups in range(0,nscups):
                    l1=self.Scups["lvl1"][iscups]-1
                    l2=self.Scups["lvl2"][iscups]-1
                    popmat[l1+ci,l2+ci] += self.EDensity*dexRate[iscups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity*exRate[iscups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity*exRate[iscups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[iscups, itemp]
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans]-1
                        lvl2 = cilvl['lvl2'][itrans]-1
#                        de = cilvl['de'][itrans]
#                        ekt = (de*1.57888e+5)/temperature
                        #mult = lowMult[lvl1-1]
                        # this is kind of double booking the ionization rate components
                        popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
#                        popmat[lvl2, lvl1-1] += self.EDensity*cirate[itemp]
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity*cirate[itemp]
                    popmat[1, 0] += (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity*self.RecombRate['rate'][itemp]
                if rec:
                #
                    if self.Nreclvl:
                        recTot = self.ReclvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
                    popmat[-1,  ci] += self.EDensity*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate'][itemp]- recTot)
                    popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate'][itemp]- recTot)
#                    popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate'][itemp]- self.ReclvlRate['rate'][:, itemp].sum()) - dielTot
#                    popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate'][itemp]- self.ReclvlRate['rate'][:, itemp].sum()) + dielTot
#                    popmat[ci, -1] += self.EDensity*higher.RecombRate['rate'][itemp]
#                    popmat[-1, -1] -= self.EDensity*higher.RecombRate['rate'][itemp]
                    #
#                    for itrans in range(len(reclvl['lvl1'])):
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                    #
#                if self.Dielectronic:
#                    branch = np.zeros(self.Ndielsplups, 'float64')
#                    for isplups in range(0,self.Ndielsplups):
#                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
#                        l2 = self.DielSplups["lvl2"][isplups]-1
#                        auto = rad[l1+ci, l2+ci]
#                        avalue = rad[:, l2+ci]
#                        good = avalue > 0.
#                        avalueTot = avalue[good].sum()
#                        branch[isplups] = (avalueTot-auto)/avalueTot
##                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
#                        self.DielUpsilon['branch'] =  branch
#                    #
#                    dielTot = 0.
#                        print ' Ndielsplups > 0 '
#                    for isplups in range(0,self.Ndielsplups):
#                        l1 = self.DielSplups["lvl1"][isplups]-1
#                        l2 = self.DielSplups["lvl2"][isplups]-1
#                         #
#                        popmat[l2+ci,l1+ci] += self.EDensity*dielexRate[isplups, itemp]
#                        popmat[l1+ci,l1+ci] -= self.EDensity*dielexRate[isplups, itemp]
                        #
#                        dielTot += self.EDensity*dielexRate[isplups, itemp]*branch[isplups]
#                else:
#                    dielTot = 0.
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                if self.Dielectronic:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
            #
        elif ntemp == 1:
#            pop=np.zeros((ndens,nlvls),"float64")
            pop=np.zeros((ndens,nlvls),"float64")
            for idens in range(0,ndens):
                popmat=np.copy(rad)
                for isplups in range(0,nscups):
                    l1=self.Scups["lvl1"][isplups]-1
                    l2=self.Scups["lvl2"][isplups]-1
#                    if self.Dielectronic:
#                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
#                    else:
#                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
                #
                    popmat[l1+ci,l2+ci] += self.EDensity[idens]*dexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.EDensity[idens]*exRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.EDensity[idens]*exRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.EDensity[idens]*dexRate[isplups]
                #
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
#                    # for proton excitation, the levels are all below the ionization potential
#                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
                 #
                    popmat[l1+ci,l2+ci] += self.PDensity[idens]*pdexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.PDensity[idens]*pexRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.PDensity[idens]*pexRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.PDensity[idens]*pdexRate[isplups]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
#                        de = cilvl['de'][itrans]
#                        ekt = (de*1.57888e+5)/temperature
                        # this is kind of double booking the ionization rate components
#                        popmat[lvl2, lvl1-1] += self.EDensity[idens]*cirate
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity[idens]*cirate
                        popmat[lvl2+ci, lvl1] += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        popmat[lvl1, lvl1] -= self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        ciTot += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                    popmat[1, 0] += (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 0] -= (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 1] += self.EDensity[idens]*self.RecombRate['rate']
                    popmat[1, 1] -= self.EDensity[idens]*self.RecombRate['rate']
                if rec:
#                    print ' rec = ', rec
                    if self.Nreclvl:
#                        print ' ReclvlRate.shape = ', self.ReclvlRate['rate'].shape
                        recTot = self.ReclvlRate['rate'].sum()
                    else:
                        recTot = 0.
                    #
                    popmat[-1,  ci] += self.EDensity[idens]*self.IonizRate['rate']
                    popmat[ci, ci] -= self.EDensity[idens]*self.IonizRate['rate']
                    popmat[ci, -1] += self.EDensity[idens]*(higher.RecombRate['rate'] - recTot)
                    popmat[-1, -1] -= self.EDensity[idens]*(higher.RecombRate['rate'] - recTot)
#                    popmat[ci, -1] += self.EDensity[idens]*higher.RecombRate['rate']
#                    popmat[-1, -1] -= self.EDensity[idens]*higher.RecombRate['rate']
                    #
#                    for itrans in range(len(reclvl['lvl1'])):
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
                        popmat[-1, -1] -= self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
                    #
#                if self.Dielectronic:
#                    branch = np.zeros(self.Ndielsplups, 'float64')
#                    for isplups in range(0,self.Ndielsplups):
#                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
#                        l2 = self.DielSplups["lvl2"][isplups]-1
#                        auto = rad[l1+ci, l2+ci]
#                        avalue = rad[:, l2+ci]
#                        good = avalue > 0.
#                        avalueTot = avalue[good].sum()
#                        branch[isplups] = (avalueTot-auto)/avalueTot
##                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
#                        self.DielUpsilon['branch'] =  branch
#                    #
#                    dielTot = 0.
##                        print ' Ndielsplups > 0 '
#                    for isplups in range(0,self.Ndielsplups):
#                        l1 = self.DielSplups["lvl1"][isplups]-1
#                        l2 = self.DielSplups["lvl2"][isplups]-1
#                         #
#                        popmat[l2+ci,l1+ci] += self.EDensity[idens]*dielexRate[isplups]
#                        popmat[l1+ci,l1+ci] -= self.EDensity[idens]*dielexRate[isplups]
                        #
#                        dielTot += self.EDensity[idens]*dielexRate[isplups]*branch[isplups]
#                else:
#                    dielTot = 0.

# normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                if self.Dielectronic:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[idens] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[idens] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at eDensity = ', ('%8.2e')%(eDensity[idens])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[idens] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[idens] = thispop[ci:]
#                pop[idens] = thispop[ci:ci+nlvls]
                #
        elif ntemp>1  and ntemp==ndens:
            pop=np.zeros((ntemp,nlvls),"float64")
#            pop=np.zeros((ntemp,ci+nlvls+rec),"float64")
            for itemp in range(0,ntemp):
                temp=self.Temperature[itemp]
                popmat=np.copy(rad)
                for isplups in range(0,nscups):
                    l1=self.Scups["lvl1"][isplups]-1
                    l2=self.Scups["lvl2"][isplups]-1
#                    if self.Dielectronic:
#                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
#                    else:
#                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
                    #
                    popmat[l1+ci,l2+ci] += self.EDensity[itemp]*dexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity[itemp]*dexRate[isplups, itemp]
                # proton rates
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
#                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                   print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
                        # this is kind of double booking the ionization rate components
#                        popmat[lvl2, lvl1-1] += self.EDensity[itemp]*cirate[itemp]
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity[itemp]*cirate[itemp]
                        popmat[lvl2+ci, lvl1] += self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity[itemp]*self.CilvlRAte['rate'][itrans, itemp]
                    popmat[1, 0] += (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                if rec:
                #
                    if self.Nreclvl:
                        recTot = self.ReclvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
#                   print ' rec = ', rec
                    popmat[-1,  ci] += self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.EDensity[itemp]*(higher.RecombRate['rate'][itemp] - recTot)
                    popmat[-1, -1] -= self.EDensity[itemp]*(higher.RecombRate['rate'][itemp] - recTot)
#                    popmat[ci, -1] += self.EDensity[itemp]*higher.RecombRate['rate'][itemp]
#                    popmat[-1, -1] -= self.EDensity[itemp]*higher.RecombRate['rate'][itemp]
                    #
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                # normalize to unity
#                if self.Dielectronic:
#                    branch = np.zeros(self.Ndielsplups, 'float64')
#                    for isplups in range(0,self.Ndielsplups):
#                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
#                        l2 = self.DielSplups["lvl2"][isplups]-1
#                        auto = rad[l1+ci, l2+ci]
#                        avalue = rad[:, l2+ci]
#                        good = avalue > 0.
#                        avalueTot = avalue[good].sum()
#                        branch[isplups] = (avalueTot-auto)/avalueTot
##                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
#                        self.DielUpsilon['branch'] =  branch
#                    #
#                    dielTot = 0.
#                    for isplups in range(0,self.Ndielsplups):
#                        l1 = self.DielSplups["lvl1"][isplups]-1
#                        l2 = self.DielSplups["lvl2"][isplups]-1
#                         #
#                        popmat[l2+ci,l1+ci] += self.EDensity[itemp]*dielexRate[isplups, itemp]
#                        popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*dielexRate[isplups, itemp]
                        #
#                        dielTot += self.EDensity[itemp]*dielexRate[isplups, itemp]*branch[isplups]
#                else:
#                    dielTot = 0.

                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                if self.Dielectronic:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[itemp] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[itemp] = thispop[ci:]
#                pop[itemp] = thispop[ci:ci+nlvls]
            #
        pop=np.where(pop >0., pop,0.)
        self.Population={"temperature":temperature,"eDensity":eDensity,"population":pop, "protonDensity":protonDensity, "ci":ci, "rec":rec}
        #
        return
        #
        # -------------------------------------------------------------------------------------
        #
    def populateSplups(self, popCorrect=1, verbose=0, **kwargs):
        """
        Calculate level populations for specified ion.  This is a new version that will enable the calculation
        of dielectronic satellite lines without resorting to the dielectronic ions, such as c_5d
        possible keyword arguments include temperature, eDensity, pDensity, radTemperature and rStar
        """
        #
        #
        for one in kwargs.keys():
            if one not in chdata.keywordArgs:
                print ' following keyword is not understood - ',  one
        #
        nlvls=self.Nlvls
        nwgfa=self.Nwgfa
        nsplups=self.Nsplups
        npsplups=self.Npsplups
        #
        if kwargs.has_key('temperature'):
            self.Temperature = np.asarray(kwargs['temperature'])
            temperature = self.Temperature
        elif hasattr(self, 'Temperature'):
            temperature=self.Temperature
        else:
                print ' no temperature values have been set'
                return
        #
        if kwargs.has_key('eDensity'):
            self.EDensity = np.asarray(kwargs['eDensity'])
            eDensity = self.EDensity
        elif hasattr(self, 'EDensity'):
            eDensity = self.EDensity
        else:
            print ' no eDensity values have been set'
            return
        #
        if kwargs.has_key('pDensity'):
            if kwargs['pDensity'] == 'default':
                self.p2eRatio()
                protonDensity = self.ProtonDensityRatio*self.EDensity
            else:
                try:
                    self.PDensity = np.asarray(kwargs['pDensity'])
                except:
                    print ' could not interpret value for keyword pDensity'
                    print ' should be either "default" or a number or array'
                    return
        else:
            if hasattr(self, 'PDensity'):
                protonDensity = self.PDensity
            else:
                self.p2eRatio()
                self.PDensity = self.ProtonDensityRatio*self.EDensity
                protonDensity = self.PDensity
                print ' proton density not specified, set to "default"'
        #
        if 'radTemperature' in kwargs.keys() and 'rStar' in kwargs.keys():
            self.RadTemperature = np.asarray(kwargs['radTemperature'])
            radTemperature = np.array(self.RadTemperature)
            self.RStar = np.asarray(kwargs['rStar'])
            rStar = np.asarray(self.RStar)
        elif hasattr(self, 'RadTemperature') and hasattr(self, 'RStar'):
            radTemperature = self.RadTemperature
            rStar = self.RStar
        #
        #
        # the Dielectronic test should eventually go away
        if popCorrect and (not self.Dielectronic):
            if self.Ncilvl:
                ci = 1
                cilvl = self.Cilvl
                if hasattr(self, 'CilvlRate'):
                    cilvlRate = self.CilvlRate
                else:
                    self.cireclvlDescale('cilvl')
                    cilvlRate = self.CilvlRate
                self.recombRate()
                #
                lowers = util.zion2name(self.Z, self.Ion-1)
                # get the lower ionization stage
                lower = ion(lowers, temperature=self.Temperature, eDensity = self.EDensity)
                lower.ionizRate()
                # need to get multiplicity of lower ionization stage
                lowMult = lower.Elvlc['mult']
            else:
                ci = 0
#            try:
            if self.Nreclvl:
                rec = 1
                reclvl = self.Reclvl
                if hasattr(self, 'ReclvlRate'):
                    reclvlRate = self.ReclvlRate
                else:
#                    print ' doing reclvlDescale in populate'
                    self.cireclvlDescale('reclvl')
                    reclvlRate = self.ReclvlRate
            elif self.Ndielsplups:
                self.upsilonDescale(diel=1)
                dielexRate = self.DielUpsilon['exRate']
                rec = 1
            else:
                rec = 0
#            except:
#                self.Reclvl = util.cireclvlRead(self.IonStr,'reclvl' )
#                reclvl = self.Reclvl
#                self.reclvlDescale()
#                if type(self.Reclvl) != type(None):
#                    rec = 1
#                    self.ionizRate()
#                    #  get the higher ionization stage
#                    highers = util.zion2name(self.Z, self.Ion+1)
##                   print ' highers = ', highers
#                    higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
#                    higher.recombRate()
#                else:
#                    rec = 0
            #
        else:
            ci = 0
            rec = 0
        #
        if rec:
            if self.Ndielsplups:
                self.upsilonDescale(diel=1)
                dielexRate = self.DielUpsilon['exRate']
            # get ionization rate of this iion
            self.ionizRate()
            #  get the higher ionization stage
            highers = util.zion2name(self.Z, self.Ion+1)
            higher = ion(highers, temperature=self.Temperature, eDensity=self.EDensity)
            higher.recombRate()
#        print ' nlvls, ci, rec = ', nlvls, ci, rec
        #
        rad=np.zeros((nlvls+ci+rec,nlvls+ci+rec),"float64")  #  the populating matrix for radiative transitions
        #
        #
        for iwgfa in range(nwgfa):
            l1 = self.Wgfa["lvl1"][iwgfa]-1
            l2 = self.Wgfa["lvl2"][iwgfa]-1
            rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]
            rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]
            # photo-excitation and stimulated emission
            if self.RadTemperature:
                if not self.RStar:
                    dilute = 0.5
                else:
                    dilute = util.dilute(self.RStar)
                # next - don't include autoionization lines
                if abs(self.Wgfa['wvl'][iwgfa]) > 0.:
                    de = const.invCm2Erg*(self.Elvlc['ecm'][l2] - self.Elvlc['ecm'][l1])
                    dekt = de/(const.boltzmann*self.RadTemperature)
                    # photoexcitation
                    phexFactor = dilute*(float(self.Elvlc['mult'][l2])/float(self.Elvlc['mult'][l1]))/(np.exp(dekt) -1.)
                    rad[l2+ci,l1+ci] += self.Wgfa["avalue"][iwgfa]*phexFactor
                    rad[l1+ci,l1+ci] -= self.Wgfa["avalue"][iwgfa]*phexFactor
                    # stimulated emission
                    stemFactor = dilute/(np.exp(-dekt) -1.)
                    rad[l1+ci,l2+ci] += self.Wgfa["avalue"][iwgfa]*stemFactor
                    rad[l2+ci,l2+ci] -= self.Wgfa["avalue"][iwgfa]*stemFactor

        #
        #
        if self.Nsplups:
            self.upsilonDescale()
            ups = self.Upsilon['upsilon']
            exRate = self.Upsilon['exRate']
            dexRate = self.Upsilon['dexRate']
        #
        if npsplups:
            self.upsilonDescale(prot=1)
#            pups = self.PUpsilon['upsilon']
            pexRate = self.PUpsilon['exRate']
            pdexRate = self.PUpsilon['dexRate']
        #
        temp=temperature
        ntemp=temp.size
        #
        cc=const.collision*self.EDensity
        ndens=cc.size
        if npsplups:
            cp=const.collision*protonDensity
        if ntemp > 1 and ndens >1 and ntemp != ndens:
            print ' unless temperature or eDensity are single values'
            print ' the number of temperatures values must match the '
            print ' the number of eDensity values'
            return
        #
        # get corrections for recombination and excitation
        #
        #
        #  first, for ntemp=ndens=1
        if ndens==1 and ntemp==1:
            popmat=np.copy(rad)
            for isplups in range(0,nsplups):
                l1=self.Splups["lvl1"][isplups]-1
                l2=self.Splups["lvl2"][isplups]-1
                #
                popmat[l1+ci,l2+ci] += self.EDensity*dexRate[isplups]
                popmat[l2+ci,l1+ci] += self.EDensity*exRate[isplups]
                popmat[l1+ci,l1+ci] -= self.EDensity*exRate[isplups]
                popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[isplups]
                #
            for isplups in range(0,npsplups):
                l1=self.Psplups["lvl1"][isplups]-1
                l2=self.Psplups["lvl2"][isplups]-1
                 #
                popmat[l1+ci,l2+ci] += self.PDensity*pdexRate[isplups]
                popmat[l2+ci,l1+ci] += self.PDensity*pexRate[isplups]
                popmat[l1+ci,l1+ci] -= self.PDensity*pexRate[isplups]
                popmat[l2+ci,l2+ci] -= self.PDensity*pdexRate[isplups]
           # now include ionization rate from
            if ci:
#                print ' ci = ', ci
                #
                # the ciRate can be computed for all temperatures
                #
                ciTot = 0.
                for itrans in range(len(cilvl['lvl1'])):
                    lvl1 = cilvl['lvl1'][itrans]-1
                    lvl2 = cilvl['lvl2'][itrans]-1
#                    de = cilvl['de'][itrans]
#                    ekt = (de*1.57888e+5)/temperature
#                    mult = lowMult[lvl1-1]
                    # this is kind of double booking the ionization rate components
                    popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans]
                    popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans]
                    ciTot += self.EDensity*self.CilvlRate['rate'][itrans]
                #
                popmat[1, 0] += (self.EDensity*lower.IonizRate['rate'] - ciTot)
                popmat[0, 0] -= (self.EDensity*lower.IonizRate['rate'] - ciTot)
                popmat[0, 1] += self.EDensity*self.RecombRate['rate']
                popmat[1, 1] -= self.EDensity*self.RecombRate['rate']
            if rec:
                #
                if self.Ndielsplups:
                    branch = np.zeros(self.Ndielsplups, 'float64')
                    for isplups in range(0,self.Ndielsplups):
                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                        l2 = self.DielSplups["lvl2"][isplups]-1
                        auto = rad[l1+ci, l2+ci]
                        avalue = rad[:, l2+ci]
                        good = avalue > 0.
                        avalueTot = avalue[good].sum()
                        branch[isplups] = (avalueTot-auto)/avalueTot
#                        print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                        self.DielUpsilon['branch'] =  branch
                    #
                    dielTot = 0.
#                    print ' Ndielsplups > 0 '
                    for isplups in range(0,self.Ndielsplups):
                        l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                        l2 = self.DielSplups["lvl2"][isplups]-1
                         #
#                        print ' l1, l2, dielexRate = ', l1, l2, dielexRate[isplups]
                        popmat[l2+ci,l1+ci] += self.EDensity*dielexRate[isplups]
                        popmat[l1+ci,l1+ci] -= self.EDensity*dielexRate[isplups]
                        #
                        dielTot += self.EDensity*dielexRate[isplups]*branch[isplups]
                else:
                    dielTot = 0.
                #
#                print ' rec, dielTot  = ', rec,  dielTot
                #
                for itrans in range(self.Nreclvl):
#                    lvl1 = reclvl['lvl1'][itrans]-1
                    lvl2 = reclvl['lvl2'][itrans]-1
                    popmat[lvl2+ci, -1] += self.EDensity*reclvlRate['rate'][itrans]
                    popmat[-1, -1] -= self.EDensity*reclvlRate['rate'][itrans]
                if self.Nreclvl:
                    reclvlRateTot = reclvlRate['rate'].sum(axis=0)
                else:
                    reclvlRateTot = 0.

                #
                popmat[-1,  ci] += self.EDensity*self.IonizRate['rate']
                popmat[ci, ci] -= self.EDensity*self.IonizRate['rate']
                # next 2 line take care of overbooking
                popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate']- reclvlRateTot - dielTot)
                popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate']- reclvlRateTot - dielTot)
#                print ' higher, rec , dieltot = ',  self.EDensity*higher.RecombRate['rate'], self.EDensity*reclvlRate['rate'].sum(axis=0),  dielTot
            # normalize to unity
            norm=np.ones(nlvls+ci+rec,'float64')
            if ci:
                norm[0] = 0.
            if rec:
                norm[nlvls+ci+rec-1] = 0.
            popmat[nlvls+ci+rec-1]=norm
            b=np.zeros(nlvls+ci+rec,'float64')
            b[nlvls+ci+rec-1]=1.
#            print ' norm = ', norm
#            print 'popmat, last line',  popmat[-1]
#            print ' b = ', b
#            popmat[nlvls/2]=norm
#            b=np.zeros(nlvls+ci+rec,'float64')
#            b[nlvls/2]=1.
#            if rec:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:ci+nlvls+rec-1]
#            else:
#                fullpop = np.linalg.solve(popmat,b)
#                pop = fullpop[ci:]
#            fullpop = np.linalg.solve(popmat,b)
            try:
                fullpop=np.linalg.solve(popmat,b)
                pop = fullpop[ci:ci+nlvls]
            except np.linalg.LinAlgError:
                pop = np.zeros(nlvls, 'float64')
#                print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature)
            #
        #   next, in case of a single eDensity value
#            pop = np.linalg.solve(popmat,b)
        elif ndens == 1:
            pop=np.zeros((ntemp, nlvls),"float64")
#            pop=np.zeros((ntemp,ci + nlvls + rec),"float64")
            for itemp in range(0,ntemp):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
                    popmat[l1+ci,l2+ci] += self.EDensity*dexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity*exRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity*exRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity*dexRate[isplups, itemp]
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans]-1
                        lvl2 = cilvl['lvl2'][itrans]-1
#                        de = cilvl['de'][itrans]
#                        ekt = (de*1.57888e+5)/temperature
                        mult = lowMult[lvl1-1]
                        # this is kind of double booking the ionization rate components
                        popmat[lvl2+ci, lvl1] += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity*self.CilvlRate['rate'][itrans, itemp]
#                        popmat[lvl2, lvl1-1] += self.EDensity*cirate[itemp]
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity*cirate[itemp]
                    popmat[1, 0] += (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity*self.RecombRate['rate'][itemp]
                if rec:
                #
                    if self.Ndielsplups:
                        branch = np.zeros(self.Ndielsplups, 'float64')
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                            auto = rad[l1+ci, l2+ci]
                            avalue = rad[:, l2+ci]
                            good = avalue > 0.
                            avalueTot = avalue[good].sum()
                            branch[isplups] = (avalueTot-auto)/avalueTot
#                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                            self.DielUpsilon['branch'] =  branch
                        #
                        dielTot = 0.
#                        print ' Ndielsplups > 0 '
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                             #
                            popmat[l2+ci,l1+ci] += self.EDensity*dielexRate[isplups, itemp]
                            popmat[l1+ci,l1+ci] -= self.EDensity*dielexRate[isplups, itemp]
                            #
                            dielTot += self.EDensity*dielexRate[isplups, itemp]*branch[isplups]
                    else:
                        dielTot = 0.
                    if self.Nreclvl:
                        recTot = self.ReclvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
                    popmat[-1,  ci] += self.EDensity*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate'][itemp]- recTot - dielTot)
                    popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate'][itemp]- recTot - dielTot)
#                    popmat[ci, -1] += self.EDensity*(higher.RecombRate['rate'][itemp]- self.ReclvlRate['rate'][:, itemp].sum()) - dielTot
#                    popmat[-1, -1] -= self.EDensity*(higher.RecombRate['rate'][itemp]- self.ReclvlRate['rate'][:, itemp].sum()) + dielTot
#                    popmat[ci, -1] += self.EDensity*higher.RecombRate['rate'][itemp]
#                    popmat[-1, -1] -= self.EDensity*higher.RecombRate['rate'][itemp]
                    #
#                    for itrans in range(len(reclvl['lvl1'])):
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity*self.ReclvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
            #
        elif ntemp == 1:
#            pop=np.zeros((ndens,nlvls),"float64")
            pop=np.zeros((ndens,nlvls),"float64")
            for idens in range(0,ndens):
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
#                    if self.Dielectronic:
#                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
#                    else:
#                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cc[idens]*ups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cc[idens]*ups[isplups]/(fmult2*np.sqrt(temp))
                #
                    popmat[l1+ci,l2+ci] += self.EDensity[idens]*dexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.EDensity[idens]*exRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.EDensity[idens]*exRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.EDensity[idens]*dexRate[isplups]
                #
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
#                    # for proton excitation, the levels are all below the ionization potential
#                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cp[idens]*pups[isplups]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cp[idens]*pups[isplups]/(fmult2*np.sqrt(temp))
                 #
                    popmat[l1+ci,l2+ci] += self.PDensity[idens]*pdexRate[isplups]
                    popmat[l2+ci,l1+ci] += self.PDensity[idens]*pexRate[isplups]
                    popmat[l1+ci,l1+ci] -= self.PDensity[idens]*pexRate[isplups]
                    popmat[l2+ci,l2+ci] -= self.PDensity[idens]*pdexRate[isplups]
                # now include ionization rate from
                if ci:
#                    print ' ci = ', ci
                    #
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
#                        de = cilvl['de'][itrans]
#                        ekt = (de*1.57888e+5)/temperature
                        # this is kind of double booking the ionization rate components
#                        popmat[lvl2, lvl1-1] += self.EDensity[idens]*cirate
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity[idens]*cirate
                        popmat[lvl2+ci, lvl1] += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        popmat[lvl1, lvl1] -= self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                        ciTot += self.EDensity[idens]*self.CilvlRate['rate'][itrans]
                    popmat[1, 0] += (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 0] -= (self.EDensity[idens]*lower.IonizRate['rate'] -ciTot)
                    popmat[0, 1] += self.EDensity[idens]*self.RecombRate['rate']
                    popmat[1, 1] -= self.EDensity[idens]*self.RecombRate['rate']
                if rec:
#                    print ' rec = ', rec
                    if self.Ndielsplups:
                        branch = np.zeros(self.Ndielsplups, 'float64')
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                            auto = rad[l1+ci, l2+ci]
                            avalue = rad[:, l2+ci]
                            good = avalue > 0.
                            avalueTot = avalue[good].sum()
                            branch[isplups] = (avalueTot-auto)/avalueTot
#                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                            self.DielUpsilon['branch'] =  branch
                        #
                        dielTot = 0.
#                        print ' Ndielsplups > 0 '
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                             #
                            popmat[l2+ci,l1+ci] += self.EDensity[idens]*dielexRate[isplups]
                            popmat[l1+ci,l1+ci] -= self.EDensity[idens]*dielexRate[isplups]
                            #
                            dielTot += self.EDensity[idens]*dielexRate[isplups]*branch[isplups]
                    else:
                        dielTot = 0.
                    if self.Nreclvl:
#                        print ' ReclvlRate.shape = ', self.ReclvlRate['rate'].shape
                        recTot = self.ReclvlRate['rate'].sum()
                    else:
                        recTot = 0.
                    #
                    popmat[-1,  ci] += self.EDensity[idens]*self.IonizRate['rate']
                    popmat[ci, ci] -= self.EDensity[idens]*self.IonizRate['rate']
                    popmat[ci, -1] += self.EDensity[idens]*(higher.RecombRate['rate']
                        - recTot - dielTot)
                    popmat[-1, -1] -= self.EDensity[idens]*(higher.RecombRate['rate']
                        - recTot - dielTot)
#                    popmat[ci, -1] += self.EDensity[idens]*higher.RecombRate['rate']
#                    popmat[-1, -1] -= self.EDensity[idens]*higher.RecombRate['rate']
                    #
#                    for itrans in range(len(reclvl['lvl1'])):
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
                        popmat[-1, -1] -= self.EDensity[idens]*self.ReclvlRate['rate'][itrans]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[idens] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[idens] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at eDensity = ', ('%8.2e')%(eDensity[idens])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[idens] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[idens] = thispop[ci:]
#                pop[idens] = thispop[ci:ci+nlvls]
                #
        elif ntemp>1  and ntemp==ndens:
            pop=np.zeros((ntemp,nlvls),"float64")
#            pop=np.zeros((ntemp,ci+nlvls+rec),"float64")
            for itemp in range(0,ntemp):
                temp=self.Temperature[itemp]
                popmat=np.copy(rad)
                for isplups in range(0,nsplups):
                    l1=self.Splups["lvl1"][isplups]-1
                    l2=self.Splups["lvl2"][isplups]-1
#                    if self.Dielectronic:
#                        de=np.abs((self.Elvlc["eryd"][l2]-self.Ip/const.ryd2Ev)-self.Elvlc["eryd"][l1])
#                    else:
#                        de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cc[itemp]*ups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cc[itemp]*ups[isplups,itemp]/(fmult2*np.sqrt(temp))
                    #
                    popmat[l1+ci,l2+ci] += self.EDensity[itemp]*dexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*exRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.EDensity[itemp]*dexRate[isplups, itemp]
                # proton rates
                for isplups in range(0,npsplups):
                    l1=self.Psplups["lvl1"][isplups]-1
                    l2=self.Psplups["lvl2"][isplups]-1
                    # for proton excitation, the levels are all below the ionization potential
#                    de=np.abs(self.Elvlc["eryd"][l2]-self.Elvlc["eryd"][l1])
#                    ekt=(de*1.57888e+5)/temp
#                    fmult1=float(self.Elvlc["mult"][l1])
#                    fmult2=float(self.Elvlc["mult"][l2])
#                    popmat[l1+ci,l2+ci]+=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
#                    popmat[l2+ci,l1+ci]+=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l1+ci,l1+ci]-=cp[itemp]*pups[isplups,itemp]*np.exp(-ekt)/(fmult1*np.sqrt(temp))
#                    popmat[l2+ci,l2+ci]-=cp[itemp]*pups[isplups,itemp]/(fmult2*np.sqrt(temp))
                     #
                    popmat[l1+ci,l2+ci] += self.PDensity[itemp]*pdexRate[isplups, itemp]
                    popmat[l2+ci,l1+ci] += self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l1+ci,l1+ci] -= self.PDensity[itemp]*pexRate[isplups, itemp]
                    popmat[l2+ci,l2+ci] -= self.PDensity[itemp]*pdexRate[isplups, itemp]
                # now include ionization rate from
                if ci:
#                   print ' ci = ', ci
                    #
                    # the ciRate can be computed for all temperatures
                    #
                    ciTot = 0.
                    for itrans in range(len(cilvl['lvl1'])):
                        lvl1 = cilvl['lvl1'][itrans] -1
                        lvl2 = cilvl['lvl2'][itrans] -1
                        # this is kind of double booking the ionization rate components
#                        popmat[lvl2, lvl1-1] += self.EDensity[itemp]*cirate[itemp]
#                        popmat[lvl1-1, lvl1-1] -= self.EDensity[itemp]*cirate[itemp]
                        popmat[lvl2+ci, lvl1] += self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        popmat[lvl1, lvl1] -= self.EDensity[itemp]*self.CilvlRate['rate'][itrans, itemp]
                        ciTot += self.EDensity[itemp]*self.CilvlRAte['rate'][itrans, itemp]
                    popmat[1, 0] += (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 0] -= (self.EDensity[itemp]*lower.IonizRate['rate'][itemp] - ciTot)
                    popmat[0, 1] += self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                    popmat[1, 1] -= self.EDensity[itemp]*self.RecombRate['rate'][itemp]
                if rec:
                #
                    if self.Ndielsplups:
                        branch = np.zeros(self.Ndielsplups, 'float64')
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                            auto = rad[l1+ci, l2+ci]
                            avalue = rad[:, l2+ci]
                            good = avalue > 0.
                            avalueTot = avalue[good].sum()
                            branch[isplups] = (avalueTot-auto)/avalueTot
#                            print ' l1 %4i l2 %4i auto %10.2e  avalue %10.2e tot %10.3f'%( l1,  l2,  auto,  avalueTot,  branch[isplups])
                            self.DielUpsilon['branch'] =  branch
                        #
                        dielTot = 0.
                        for isplups in range(0,self.Ndielsplups):
                            l1 = self.DielSplups["lvl1"][isplups]-1 + nlvls
                            l2 = self.DielSplups["lvl2"][isplups]-1
                             #
                            popmat[l2+ci,l1+ci] += self.EDensity[itemp]*dielexRate[isplups, itemp]
                            popmat[l1+ci,l1+ci] -= self.EDensity[itemp]*dielexRate[isplups, itemp]
                            #
                            dielTot += self.EDensity[itemp]*dielexRate[isplups, itemp]*branch[isplups]
                    else:
                        dielTot = 0.
                    if self.Nreclvl:
                        recTot = self.ReclvlRate['rate'][:, itemp].sum()
                    else:
                        recTot = 0.
                #
#                   print ' rec = ', rec
                    popmat[-1,  ci] += self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, ci] -= self.EDensity[itemp]*self.IonizRate['rate'][itemp]
                    popmat[ci, -1] += self.EDensity[itemp]*(higher.RecombRate['rate'][itemp]
                        - recTot - dielTot)
                    popmat[-1, -1] -= self.EDensity[itemp]*(higher.RecombRate['rate'][itemp]
                        - recTot - dielTot)
#                    popmat[ci, -1] += self.EDensity[itemp]*higher.RecombRate['rate'][itemp]
#                    popmat[-1, -1] -= self.EDensity[itemp]*higher.RecombRate['rate'][itemp]
                    #
                    for itrans in range(self.Nreclvl):
                        lvl1 = reclvl['lvl1'][itrans]-1
                        lvl2 = reclvl['lvl2'][itrans]-1
                        popmat[lvl2+ci, -1] += self.EDensity[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                        popmat[-1, -1] -= self.EDensity[itemp]*self.ReclvlRate['rate'][itrans, itemp]
                # normalize to unity
                norm=np.ones(nlvls+ci+rec,'float64')
                if ci:
                    norm[0] = 0.
                if rec:
                    norm[-1] = 0.
                popmat[nlvls+ci+rec-1]=norm
                b=np.zeros(nlvls+ci+rec,'float64')
                b[nlvls+ci+rec-1]=1.
                try:
                    thispop=np.linalg.solve(popmat,b)
                    pop[itemp] = thispop[ci:ci+nlvls]
                except np.linalg.LinAlgError:
                    pop[itemp] = np.zeros(nlvls, 'float64')
#                    print ' error in matrix inversion, setting populations to zero at T = ', ('%8.2e')%(temperature[itemp])
#                thispop=np.linalg.solve(popmat,b)
#                if rec:
#                    pop[itemp] = thispop[ci:ci+nlvls+rec-1]
#                else:
#                    pop[itemp] = thispop[ci:]
#                pop[itemp] = thispop[ci:ci+nlvls]
            #
        pop=np.where(pop >0., pop,0.)
        self.Population={"temperature":temperature,"eDensity":eDensity,"population":pop, "protonDensity":protonDensity, "ci":ci, "rec":rec}
        #
        return
        #
        # -------------------------------------------------------------------------
        #
    def setupSplups(self, dir=0, verbose=0):
        '''
        if ion is initiated with setup=0, this allows the setup to be done at a later point
        perhaps, more importantly,  by setting dir to a directory cotaining the necessary files
        for a ChiantiPy ion, it allows one to setup an ion with files not in the current
        Chianti directory
        '''
        #
        # read in all data if in masterlist
        #  if not, there should still be ionization and recombination rates
        #
        MasterList = chdata.MasterList
        #
        if dir:
            fileName = os.path.join(dir, self.IonStr)
        else:
            fileName = util.ion2filename(self.IonStr)
        if self.IonStr in MasterList:
            self.Elvlc = util.elvlcRead('', filename=fileName+'.elvlc',  verbose=verbose)
            self.Wgfa = util.wgfaRead('', filename=fileName+'.wgfa')
            self.Nwgfa=len(self.Wgfa['lvl1'])
            nlvlWgfa = max(self.Wgfa['lvl2'])
            nlvlList =[nlvlWgfa]
#                print 'fileName = ', fileName
            splupsfile = fileName + '.splups'
            if os.path.isfile(splupsfile):
                # happens the case of fe_3 and prob. a few others
                self.Splups = util.splupsRead('', filename=fileName+'.splups')
                self.Nsplups=len(self.Splups['lvl1'])
                nlvlSplups = max(self.Splups['lvl2'])
                nlvlList.append(nlvlSplups)
            else:
                self.Nsplups = 0
                nlvlSplups = 0
##                self.Nlvls = nlvlElvlc
            #
            file = fileName +'.cilvl'
            if os.path.isfile(file):
                self.Cilvl = util.cireclvlRead('',filename = fileName, cilvl=1)
                self.Ncilvl=len(self.Cilvl['lvl1'])
                nlvlCilvl = max(self.Cilvl['lvl2'])
                nlvlList.append(nlvlCilvl)
            else:
                self.Ncilvl = 0
            #  .reclvl file may not exist
            reclvlfile = fileName +'.reclvl'
            if os.path.isfile(reclvlfile):
                self.Reclvl = util.cireclvlRead('',filename=fileName, reclvl=1)
                self.Nreclvl = len(self.Reclvl['lvl1'])
                nlvlReclvl = max(self.Reclvl['lvl2'])
                nlvlList.append(nlvlReclvl)
            else:
                self.Nreclvl = 0
            #  .dielsplups file may not exist
            dielsplupsfile = fileName +'.splups'
            if os.path.isfile(dielsplupsfile):
                self.DielSplups = util.splupsRead('', filename=dielsplupsfile, diel=1)
                self.Ndielsplups=len(self.DielSplups["lvl1"])
                nlvlDielSplups = max(self.DielSplups['lvl2'])
                nlvlList.append(nlvlDielSplups)
            else:
                self.Ndielsplups = 0
            #
            #  psplups file may not exist
            psplupsfile = fileName +'.psplups'
            if os.path.isfile(psplupsfile):
                self.Psplups = util.splupsRead('', filename=psplupsfile,  prot=True)
                self.Npsplups=len(self.Psplups["lvl1"])
            else:
                self.Npsplups = 0
            #
            drparamsFile = fileName +'.drparams'
            if os.path.isfile(drparamsFile):
                self.DrParams = util.drRead(self.IonStr)
            #
            rrparamsFile = fileName +'.rrparams'
            if os.path.isfile(rrparamsFile):
                self.RrParams = util.rrRead(self.IonStr)

            #  not needed for ion, only phion
#                photoxfile = util.ion2filename(self.IonStr)+'.photox'
#                if os.path.isfile(photoxfile):
#                    self.Photox = util.photoxRead(self.IonStr)
            #
            # need to determine the number of levels that can be populated
            nlvlElvlc = len(self.Elvlc['lvl'])
#                print ' nlvlElvlc = ', nlvlElvlc
#                print ' other nlvls = ',  nlvlList
#                nlvlWgfa = max(self.Wgfa['lvl2'])
            #  elvlc file can have more levels than the rate level files
            self.Nlvls = min([nlvlElvlc, max(nlvlList)])
        else:
            print ' the ion ' + self.IonStr + ' is not in the CHIANTI masterlist '
            try:
                self.Elvlc = util.elvlcRead(self.IonStr, verbose=verbose)
                print ' elvlc file available '
            except:
                print ' elvlc file NOT available '
