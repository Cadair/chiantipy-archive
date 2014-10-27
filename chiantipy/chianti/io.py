'''
a collection of reading and writing functions
'''
import os, fnmatch
try:
    # for Python 3 import
    import configparser
except ImportError:
    # for Python 2 import
    import ConfigParser as configparser
#from ConfigParser import *
import numpy as np
import chianti.util as util
import chianti.constants as const
from .FortranFormat import *
    #
    # -------------------------------------------------------------------------------------
    #
def abundanceRead(abundancename=''):
    """ read an abundanc file and returns the abundance values relative to hydrogen"""
    abundance=np.zeros((50),'Float64')
    xuvtop=os.environ["XUVTOP"]
    if abundancename!='':
        # a specific abundance file name has been specified
        fname=os.path.join(xuvtop,'abundance',abundancename+'.abund')
    else:
        # the user will select an abundance file
        abundir=os.path.join(xuvtop,'abundance')
        abundlabel = 'ChiantiPy - Select an abundance file'
        fname = chianti.gui.chpicker(abundir, filter='*.abund', label=abundlabel)
        if fname == None:
            print((' no abundance file selected'))
            return 0
        else:
            abundancefilename=os.path.basename(fname)
            abundancename,ext=os.path.splitext(abundancefilename)
#    else:
#        # the default abundance file will be used
#        abundancename=self.Defaults['abundfile']
#        fname=os.path.join(xuvtop,'abundance',abundancename+'.abund')
    input=open(fname,'r')
    s1=input.readlines()
    input.close()
    nlines=0
    idx=-1
    while idx <= 0:
        aline=s1[nlines][0:5]
        idx=aline.find('-1')
        nlines+=1
    nlines-=1
    for line in range(nlines):
        z,ab,element=s1[line].split()
        abundance[int(z)-1]=float(ab)
    gz=np.nonzero(abundance)
    abs=10.**(abundance[gz]-abundance[0])
    abundance.put(gz,abs)
    abundanceRef=s1[nlines+1:-2]
    return {'abundancename':abundancename,'abundance':abundance,'abundanceRef':abundanceRef}
    #
    # -------------------------------------------------------------------------------------
    #
def cireclvlRead(ions, filename=0, cilvl=0, reclvl=0, rrlvl=0):
    '''
    to read Chianti cilvl and reclvl files and return data
    must specify type as either cilvl, reclvl or rrlvl
    '''
    if filename:
        fname = filename
    else:
        fname = ion2filename(ions)
    if cilvl:
        paramname=fname+'.cilvl'
    elif reclvl:
        paramname = fname + '.reclvl'
    elif rrlvl:
        paramname = fname + '.rrlvl'
    else:
        print('either "cilvl", "reclvl" ir "rrlvl" must be specified')
        return {}
    if os.path.exists(paramname):
        input=open(paramname,'r')
        lines = input.readlines()
        input.close()
    else:
        print(('file does not exist:  ', paramname))
        return {'error':'file does not exist: ' + paramname}
    #
    iline = 0
    idx = -1
    while idx < 0:
        aline=lines[iline][0:5]
        idx=aline.find('-1')
        iline += 1
    ndata = iline - 1
    ntrans = ndata/2
    #
#    nref = 0
#    idx = -1
#    while idx < 0:
#        aline=lines[iline][0:5]
#        idx=aline.find('-1')
#        iline += 1
#        nref += 1
#    nref -= 1
    #
    # need to find the maximum number of temperatures, not all lines are the same
    #
    ntemp = np.zeros(ntrans, 'int32')
    iline = 0
    for jline in range(0, ndata, 2):
        dummy = lines[jline].replace(os.linesep, '').split()
        ntemp[iline] = len(dummy[4:])
        iline += 1
    maxNtemp = ntemp.max()
#   print ' maxNtemp = ', maxNtemp
    temp = np.zeros((ntrans,maxNtemp), 'float64')
    iline = 0
    for jline in range(0, ndata, 2):
        recdat = lines[jline].replace(os.linesep, '').split()
        shortT = np.asarray(recdat[4:], 'float64')
        # the result of the next statement is to continue to replicate t
        t = np.resize(shortT, maxNtemp)
        if rrlvl:
            temp[iline] = t
        else:
            temp[iline] = 10.**t
        iline += 1
    #
    lvl1 = np.zeros(ntrans, 'int64')
    lvl2 = np.zeros(ntrans, 'int64')
    ci = np.zeros((ntrans, maxNtemp), 'float64')
    #
    idat = 0
    for jline in range(1, ndata, 2):
        cidat = lines[jline].replace(os.linesep, '').split()
        shortCi = np.asarray(cidat[4:], 'float64')
        lvl1[idat] = int(cidat[2])
        lvl2[idat] = int(cidat[3])
        ci[idat] = np.resize(shortCi, maxNtemp)
        idat += 1
    return {'temperature':temp, 'ntemp':ntemp,'lvl1':lvl1, 'lvl2':lvl2, 'rate':ci,'ref':lines[ndata+1:], 'ionS':ions}
    #
    # -------------------------------------------------------------------------------------
    #
def defaultsRead(verbose=0):
    '''
    possibleDefaults = {'wavelength':['angstrom', 'kev', 'nm']}
    symbolDefaults = {'wavelength':['A', 'keV', 'nm']}
    '''
    initDefaults={'abundfile': 'sun_photospheric_1998_grevesse','ioneqfile': 'chianti', 'wavelength': 'angstrom', 'flux': 'energy','gui':False}
    rcfile=os.path.join(os.environ['HOME'],'.chianti/chiantirc')
    if os.path.isfile(rcfile):
        print((' reading chiantirc file'))
        config = configparser.RawConfigParser(initDefaults)
        config.read(rcfile)
        defaults = {}
        for anitem in config.items('chianti'):
            defaults[anitem[0]] = anitem[1]
        if defaults['gui'].lower() in ('t', 'y', 'yes', 'on', 'true', '1', 1, True):
            defaults['gui'] = True
        elif defaults['gui'].lower() in ('f', 'n', 'no', 'off', 'false', '0', 0, False):
            defaults['gui'] = False
    else:
        defaults = initDefaults
        if verbose:
            print((' chiantirc file (/HOME/.chianti/chiantirc) does not exist'))
            print((' using the following defaults'))
            for akey in list(defaults.keys()):
                print((' %s = %s'%(akey, defaults[akey])))
    return defaults
    #
    #-----------------------------------------------------------
    #
def diRead(ions, filename=0):
    """
    read chianti direct ionization .params files and return
        {"info":info,"btf":btf,"ev1":ev1,"xsplom":xsplom,"ysplom":ysplom,"ref":hdr}
        info={"iz":iz,"ion":ion,"nspl":nspl,"neaev":neaev}
    cannot read dilvlparams files
    """
    #
    if filename:
        paramname = filename
    else:
        zion = util.convertName(ions)
        if zion['Z'] < zion['Ion']:
            print(' this is a bare nucleus that has no ionization rate')
            return
        #
        fname = util.ion2filename(ions)
        paramname=fname+'.diparams'
    #
    input=open(paramname,'r')
    #  need to read first line and see how many elements
    line1=input.readline()
    indices=line1.split()
    iz=int(indices[0])
    ion=int(indices[1])
    nspl=indices[2]
    nfac=int(indices[3])
    neaev=int(indices[4])
    nspl=int(nspl)
    format=FortranFormat(str(nspl+1)+'E10.2')
    #
    ev1=np.zeros(nfac,'Float64')
    btf=np.zeros(nfac,'Float64')
    xsplom=np.zeros([nfac, nspl],'Float64')
    ysplom=np.zeros([nfac, nspl],'Float64')
    #
    for ifac in range(nfac):
        line=input.readline()
        paramdat=FortranLine(line,format)
        btf[ifac]=paramdat[0]
        xsplom[ifac]=paramdat[1:]
        line=input.readline()
        paramdat=FortranLine(line,format)
        ev1[ifac]=paramdat[0]
        ysplom[ifac]=paramdat[1:]
    if neaev:
        line=input.readline()
        eacoef=line.split()
#            print ' eaev = ', type(eacoef), eacoef
        eaev=[float(avalue) for avalue in eacoef]
#            print ' eaev = ', type(eaev), eaev
#            print ' eaev = ', type(eaev), eaev
#            if len(eaev) == 1:
#                eaev=float(eaev[0])
#                eaev=np.asarray(eaev, 'float32')
#            else:
#                eaev=np.asarray(eaev, 'float32')
    else:
        eaev=0.
    hdr=input.readlines()
    input.close()
    info={"iz":iz,"ion":ion,"nspl":nspl,"neaev":neaev, 'nfac':nfac}
    if neaev:
        info['eaev'] = eaev
    DiParams={"info":info,"btf":btf,"ev1":ev1,"xsplom":xsplom,"ysplom":ysplom, 'eaev':eaev,"ref":hdr}
    return DiParams
    #
    # -------------------------------------------------------------------------------------
    #
def drRead(ions):
    """
    read chianti dielectronic recombination .drparams files and return
        {'rrtype','params','ref'}
        """
    #
    #
    fname = util.ion2filename(ions)
    paramname=fname+'.drparams'
    if os.path.isfile(paramname):
        input=open(paramname,'r')
        #  need to read first line and see how many elements
        lines=input.readlines()
        input.close()
        drtype=int(lines[0])
        ref=lines[4:-1]
        #
        if drtype == 1:
            # a Badnell type
            fmt=FortranFormat('2i5,8e12.4')
            eparams=np.asarray(FortranLine(lines[1],fmt)[2:], 'float64')
            cparams=np.asarray(FortranLine(lines[2],fmt)[2:], 'float64')
            DrParams={'drtype':drtype, 'eparams':eparams,'cparams':cparams,  'ref':ref}
        elif drtype == 2:
            # shull type
            fmt=FortranFormat('2i5,4e12.4')
            params=np.asarray(FortranLine(lines[1],fmt)[2:], 'float64')
            DrParams={'drtype':drtype, 'params':params, 'ref':ref}
        else:
            DrParams = None
            print((' for ion %5s unknown DR type = %5i' %(ions, drtype)))
    else:
        DrParams=None
    return DrParams
    #
    # -------------------------------------------------------------------------------------
    #
def eaRead(ions, filename=0):
    '''
    read a chianti excitation-autoionization file and return the EA ionization rate data
    derived from splupsRead
    {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups,"bsplups":bsplups,"ref":ref}
    '''
    if filename:
        splupsname = filename
    else:
        zion = util.convertName(ions)
        if zion['Z'] < zion['Ion']:
            print(' this is a bare nucleus that has no ionization rate')
            return
        #
        fname = util.ion2filename(ions)
        splupsname = fname+'.easplups'
    if not os.path.exists(splupsname):
        print((' could not find file:  ', splupsname))
        self.Splups={"lvl1":-1}
        return {"lvl1":-1}
    # there is splups/psplups data
    else:
        input=open(splupsname,'r')
        s1=input.readlines()
        dum=input.close()
        nsplups=0
        ndata=2
        while ndata > 1:
            s1a=s1[nsplups][:]
            s2=s1a.split()
            ndata=len(s2)
            nsplups=nsplups+1
        nsplups=nsplups-1
        lvl1=[0]*nsplups
        lvl2=[0]*nsplups
        ttype=[0]*nsplups
        gf=[0.]*nsplups
        de=[0.]*nsplups
        cups=[0.]*nsplups
        nspl=[0]*nsplups
        splups=np.zeros((nsplups,9),'Float64')
        splupsFormat1='(6x,3i3,8e10.0)'
        splupsFormat2='(6x,3i3,12e10.0)'
        #
        for i in range(0,nsplups):
            try:
                inpt=FortranLine(s1[i],splupsFormat1)
            except:
                inpt=FortranLine(s1[i],splupsFormat2)
            lvl1[i]=inpt[0]
            lvl2[i]=inpt[1]
            ttype[i]=inpt[2]
            gf[i]=inpt[3]
            de[i]=inpt[4]
            cups[i]=inpt[5]
            if len(inpt)  > 13:
                nspl[i]=9
                splups[i].put(list(range(9)),inpt[6:])
            else:
                nspl[i]=5
                splups[i].put(list(range(5)),inpt[6:])
        #
        ref=[]
        for i in range(nsplups+1,len(s1)-1):
            s1a=s1[i][:-1]
            ref.append(s1a.strip())
#        self.EaParams={"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
#                ,"nspl":nspl,"splups":splups,"ref":ref}
        return {"lvl1":lvl1,"lvl2":lvl2,"ttype":ttype,"gf":gf,"de":de,"cups":cups
                ,"nspl":nspl,"splups":splups,"ref":ref}
    #
    # -----------------------------------------------------------------------
    #
def elvlcRead(ions, filename=0, getExtended=0, verbose=0,  useTh=1):
    """
    reads the new format elvlc files
    read a chianti energy level file that has 6 energy columns
    and returns
    {"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
    ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,
    "ecmx":ecmx,"erydx":erydx,"ref":ref,"pretty":pretty, 'ionS':ions}
    if a energy value for ecm or eryd is zero(=unknown), the theoretical values
    (ecmth and erydth) are inserted if useTh is true
    Python 2.7 and 3 compliant
    """
    #
    #
    '%7i%30s%5s%5i%5s%5.1f%15.3f%15.3f \n'
    #
    fstring='i7,a30,a5,i5,a5,f5.1,2f15.3'
    elvlcFormat  = FortranFormat(fstring)
    #
    #
    if filename:
        elvlname = filename
        bname = os.path.basename(filename)
        ions = bname.split('.')[0]
    else:
        fname = util.ion2filename(ions)
        elvlname=fname+'.elvlc'
    if not os.path.isfile(elvlname):
        print((' elvlc file does not exist:  %s'%(elvlname)))
        return {'status':0}
    status = 1
    input=open(elvlname,'r')
    s1=input.readlines()
    input.close()
    nlvls=0
    ndata=2
    while ndata > 1:
        s1a=s1[nlvls][:-1]
        s2=s1a.split()
        ndata=len(s2)
        nlvls=nlvls+1
    nlvls-=1
    if verbose:
        print((' nlvls = %i'%(nlvls)))
    lvl=[0]*nlvls
    conf = [0]*nlvls
    term=[0]*nlvls
    label = [0]*nlvls
    spin=[0]*nlvls
    spd=[0]*nlvls
    l = ['']*nlvls
    j = [0.]*nlvls
    mult = [0.]*nlvls
    ecm=[0]*nlvls
    ecmth=[0]*nlvls
    pretty=[0]*nlvls
    if getExtended:
        extended = [' ']*nlvls
    for i in range(0,nlvls):
        if verbose:
            print((s1[i][0:115]))
        inpt = FortranLine(s1[i][0:115],elvlcFormat)
        lvl[i]=inpt[0]
        term[i]=inpt[1].strip()
        label[i] = inpt[2]
        spin[i]=inpt[3]
        spd[i]=inpt[4].strip()
        l[i] = const.Spd.index(spd[i])
        j[i]=inpt[5]
        mult[i] = 2.*inpt[5] + 1.
        ecm[i]=inpt[6]
        ecmth[i]=inpt[7]
        if ecm[i] == 0.:
            if useTh:
                ecm[i] = ecmth[i]
        stuff = term[i].strip() + ' %1i%1s%3.1f'%( spin[i], spd[i], j[i])
        pretty[i] = stuff.strip()
        if getExtended:
            cnt = s1[i].count(',')
            if cnt > 0:
                idx = s1[i].index(',')
                extended[i] = s1[i][idx+1:]
    eryd = [x*const.invCm2ryd for x in ecm]
    erydth = [x*const.invCm2ryd for x in ecmth]
    ref=[]
    # this should skip the last '-1' in the file
    for i in range(nlvls+1,len(s1) -1):
        s1a=s1[i][:-1]
        ref.append(s1a.strip())
#    self.const.Elvlc={"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
#            ,"mult":mult,"ecm":ecm,"eryd":eryd,"ecmth":ecmth,"erydth":erydth,"ref":ref}
    info = {"lvl":lvl,"conf":conf, "term":term,'label':label, "spin":spin, "spd":spd, "l":l, "j":j,
             'mult':mult, "ecm":ecm, 'eryd':eryd,'erydth':erydth, "ecmth":ecmth, "ref":ref,
             "pretty":pretty, 'status':status, 'filename':elvlname}
    if getExtended:
        info['extended'] = extended
    return info
    #
    # -------------------------------------------------------------------------------------
    #
def elvlcWrite(info, outfile=0, addLvl=0, includeRyd=0):
    '''
    for files created after elvlc format change in November 2012
    creates a .elvlc in the current directory
    info is a dictionary that must contain the following keys
    ionS, the Chianti style name of the ion such as c_4
    term, a string showing the configuration
    spin, an integer of the spin of the state in LS coupling
    l, an integer of the angular momentum quantum number
    spd, an string for the alphabetic symbol of the angular momemtum, S, P, D, etc
    j, a floating point number, the total angular momentum
    ecm, the observed energy in inverse cm, if unknown, the value is 0.
    eryd, the observed energy in Rydbergs, if unknown, the value is 0.
    ecmth, the calculated energy from the scattering calculation, in inverse cm
    erydth, the calculated energy from the scattering calculation in Rydbergs
    ref, the references in the literature to the data in the input info

    the output filename will be ionS+'.elvlc' unless outfile is specified
    addLvl is to add a constant value to the index of all levels
    setting includeRyd will also write the Rydberg energies in the extended area, demarked by a comma
    '''
    if outfile:
        elvlcName = outfile
    else:
        try:
            gname = info['ionS']
        except:
            print(' ''ionS'' not included in input dict')
            return
        elvlcName = gname + '.elvlc'
    print((' elvlc file name = ', elvlcName))
    #
#    if not info.has_key('ecmx'):
#        info['ecmx'] = np.zeros_like(info['ecm'])
#    if not info.has_key('erydx'):
#        info['erydx'] = np.zeros_like(info['eryd'])
    if 'label' not in info:
        nlvl = len(info['ecm'])
        info['label'] = [' ']*nlvl
    if 'eryd' not in info:
        info['eryd'] = [x*const.invCm2ryd for x in info['ecm']]
    if 'erydth 'not in info:
        info['erydth'] = [x*const.invCm2ryd for x in info['ecmth']]
   #
    out = open(elvlcName, 'w')
    for i,  aterm in enumerate(info['term']):
        thisTerm = aterm.ljust(29)
        thisLabel = info['label'][i].ljust(4)
#        print, ' len of thisTerm = ', len(thisTerm)
        if includeRyd:
            pstring = '%7i%30s%5s%5i%5s%5.1f%15.3f%15.3f , %15.8f , %15.8f \n'%(i+1+addLvl, thisTerm, thisLabel, info['spin'][i], info['spd'][i],info['j'][i],  info['ecm'][i], info['ecmth'][i], info['eryd'][i], info['erydth'][i])
        else:
            pstring = '%7i%30s%5s%5i%5s%5.1f%15.3f%15.3f \n'%(i+1+addLvl, thisTerm, thisLabel, info['spin'][i], info['spd'][i],info['j'][i],  info['ecm'][i], info['ecmth'][i])
        out.write(pstring)
    out.write(' -1\n')
    out.write('%filename:  ' + os.path.split(elvlcName)[1] + '\n')
#    info['ref'].append(' produced as a part of the \'CHIANTI\' atomic database for astrophysical spectroscopy')
#    today = date.today()
#    info['ref'].append(' K. Dere (GMU) - ' + today.strftime('%Y %B %d'))
#    for one in info['ref']:
#        out.write(one+'\n')
    for aref in info['ref']:
        out.write(aref + '\n')
    out.write(' -1\n')
    out.close()
    return
    #
    # ----------------------------------------------------------------------------------------
    #
def fblvlRead(filename, verbose=0):
    """
    read a chianti energy level file and returns
    {"lvl":lvl,"conf":conf,"term":term,"spin":spin,"l":l,"spd":spd,"j":j
    ,"mult":mult,"ecm":ecm,"eryd":eryd,"ref":ref}
    """
#        #  ,format='(i5,a20,2i5,a3,i5,2f20.3)'
    fstring='i5,a20,2i5,a3,i5,2f20.3'
    elvlcFormat=FortranFormat(fstring)
    #
    if os.path.exists(filename):
        input=open(filename,'r')
        s1=input.readlines()
        input.close()
        nlvls=0
        ndata=2
        while ndata > 1:
            s1a=s1[nlvls][:-1]
            s2=s1a.split()
            ndata=len(s2)
            nlvls=nlvls+1
        nlvls-=1
        if verbose:
            print((' nlvls = %5i'%(nlvls)))
        lvl=[0]*nlvls
        conf=[0]*nlvls
        pqn=[0]*nlvls
        l=[0]*nlvls
        spd=[0]*nlvls
        mult=[0]*nlvls
        ecm=[0]*nlvls
        ecmth=[0]*nlvls
        for i in range(0,nlvls):
            if verbose:
                print((s1[i]))
            inpt=FortranLine(s1[i],elvlcFormat)
            lvl[i]=inpt[0]
            conf[i]=inpt[1].strip()
            pqn[i]=inpt[2]
            l[i]=inpt[3]
            spd[i]=inpt[4].strip()
            mult[i]=inpt[5]
            if inpt[6] == 0.:
                ecm[i]=inpt[7]
            else:
                ecm[i]=inpt[6]
                ecmth[i]=inpt[7]
        ref=[]
        for i in range(nlvls+1,len(s1)-1):
            s1a=s1[i][:-1]
            ref.append(s1a.strip())
        return {"lvl":lvl,"conf":conf,'pqn':pqn,"l":l,"spd":spd,"mult":mult,
            "ecm":ecm,'ecmth':ecmth, 'ref':ref}
    else:
        return {'errorMessage':' fblvl file does not exist'}
    #
    # -------------------------------------------------------------------------------------
    #
def wgfaRead(ions, filename=0, elvlcname=-1, total=0, verbose=0):
    """
    reads chianti wgfa file and returns
    {"lvl1":lvl1,"lvl2":lvl2,"wvl":wvl,"gf":gf,"avalue":avalue,"ref":ref}
    if elvlcname is specified, the lsj term labels are returned as 'pretty1' and 'pretty2'
    """
    #
    if filename:
        wgfaname = filename
        if elvlcname < 0:
            elvlcnamee = 0
            elvlc = 0
        elif not elvlcname:
            elvlcname = os.path.splitext(wgfaname)[0] + '.elvlc'
            if os.path.isfile(elvlcname):
                elvlc = elvlcRead('', elvlcname)
            else:
                elvlc = 0
        else:
            elvlc = elvlcRead('',elvlcname)

    else:
        fname = util.ion2filename(ions)
        wgfaname=fname+'.wgfa'
        elvlcname = fname + '.elvlc'
        if os.path.isfile(elvlcname):
            elvlc = elvlcRead('', elvlcname)
        else:
            elvlc = 0
    if verbose:
        if elvlc:
            print(' have elvlc data')
        else:
            print(' do not have elvlc data')
    #
    input=open(wgfaname,'r')
    s1=input.readlines()
    input.close()
    nwvl=0
    ndata=2
    while ndata > 1:
        s1a=s1[nwvl]
        s2=s1a.split()
        ndata=len(s2)
        nwvl += 1
    nwvl -= 1
    if verbose:
        print((' nwvl = %10i ndata = %4i'%(nwvl, ndata)))
    lvl1=[0]*nwvl
    lvl2=[0]*nwvl
    wvl=[0.]*nwvl
    gf=[0.]*nwvl
    avalue=[0.]*nwvl
    if elvlc:
        pretty1 = ['']*nwvl
        pretty2 = ['']*nwvl
    #
    if verbose:
        print((' nwvl = %10i'%(nwvl)))
    #
    wgfaFormat='(2i5,f15.3,2e15.3)'
    for ivl in range(nwvl):
        inpt=FortranLine(s1[ivl],wgfaFormat)
        lvl1[ivl]=inpt[0]
        lvl2[ivl]=inpt[1]
        wvl[ivl]=inpt[2]
        gf[ivl]=inpt[3]
        avalue[ivl]=inpt[4]
        if elvlc:
            pretty1[ivl] = elvlc['pretty'][inpt[0] - 1]
            pretty2[ivl] = elvlc['pretty'][inpt[1] - 1]

    ref=[]
    # should skip the last '-1' in the file
    for i in range(nwvl+1,len(s1) -1):
        s1a=s1[i][:-1]
        ref.append(s1a.strip())
    Wgfa={"lvl1":lvl1,"lvl2":lvl2,"wvl":wvl,"gf":gf,"avalue":avalue,"ref":ref, 'ionS':ions, 'filename':wgfaname}
    if total:
        avalueLvl = [0.]*max(lvl2)
        for iwvl in range(nwvl):
            avalueLvl[lvl2[iwvl] -1] += avalue[iwvl]
        Wgfa['avalueLvl'] = avalueLvl

    if elvlc:
        Wgfa['pretty1'] = pretty1
        Wgfa['pretty2'] = pretty2
    #
    return Wgfa
    #
    # --------------------------------------
    #
def scupsRead(ions, filename=0, verbose=0):
    '''
    to read the new format ~ version 8 scups file containing the Burgess and Tully scaled temperature and upsilons.
    Python 2.7/3 compliant
    '''
    #
    if filename:
        scupsFileName = filename
        bname = os.path.basename(scupsFileName)
        ions = bname.split('.')[0]
    else:
        fname = util.ion2filename(ions)
        scupsFileName = fname+'.scups'
    if not os.path.isfile(scupsFileName):
        print((' elvlc file does not exist:  %s'%(scupsFileName)))
        return {'status':0}
    status = 1
    #
    if os.path.isfile(scupsFileName):
        inpt = open(scupsFileName)
        lines = inpt.readlines()
        inpt.close()
    else:
        print(('file does not exist: '+str(scupsFileName)))
        return {'errorMessage':'file does not exist' +str(scupsFileName)}
        return
    ll = lines[1].split()
    temp = np.asarray(ll[3:], 'float64')
    minusOne = 0
    counter = 0
    while not minusOne:
        if '-1' in lines[counter][:4]:
            minusOne = 1
        else:
            counter += 1
    ntrans = (counter)/3
    #print(' counter %i4 ntrans %i4'%(counter, ntrans))
    lvl1 = []
    lvl2 = []
    de = []
    gf = []
    lim = []
    ttype = []
    cups = []
    ntemp = []
    btemp = []
    bscups = []
    counter = 0
    for itrans in range(ntrans):
        if verbose:
            print((lines[counter]))
            print((lines[counter+1]))
            print((lines[counter+2]))
        ll1 = lines[counter].split()
        lvl1.append(int(ll1[0]))
        lvl2.append(int(ll1[1]))
        de.append(float(ll1[2]))
        gf.append(float(ll1[3]))
        lim.append(float(ll1[4]))
        ntemp.append(int(ll1[5]))
        ttype.append(int(ll1[6]))
        cups.append(float(ll1[7]))
        ll2 = lines[counter+1].split()
        ll3 = lines[counter+2].split()
#        print ' ll2 = ', ll2
#        print ' gf = ', ll2[2]
        btemp.append(np.asarray(ll2, 'float64'))
        bscups.append(np.asarray(ll3, 'float64'))
        counter += 3
    counter += 1
    ref = []
    for aline in lines[counter:]:
        ref.append(aline.strip('\n'))
    return {'lvl1':lvl1, 'lvl2':lvl2, 'de':de, 'gf':gf, 'lim':lim, 'ttype':ttype,'cups':cups,'ntemp':ntemp, 'btemp':btemp, 'bscups':bscups, 'ntrans':ntrans, 'ref':ref}
