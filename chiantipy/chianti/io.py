'''
a collection of reading and writing functions
'''
import os, fnmatch
import numpy as np
import chianti.util as util
    #
    # -------------------------------------------------------------------------------------
    #
def scupsRead(ions, filename=0, verbose=0):
    '''
    to read the new format ~ version 8 scups file containing the Burgess and Tully scaled temperature and upsilons.
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
        print('file does not exist: '+str(scupsFileName))
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
            print(lines[counter])
            print(lines[counter+1])
            print(lines[counter+2])
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
