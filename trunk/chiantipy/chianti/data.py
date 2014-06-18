'''
module for collecting various top-level Chianti data
for the keyword arguments below
temperature = temperature in Kelvin
eDensity is the electron density per cubic cm
hDensity is the hydrogen density per cubic cm
pDensity is the proton density per cubic cm
radTemperature is the radiation temperature of central source
rStar is the distance of the plasma from the source in units of the sources radius
distance is the distance from the central source
'''
import os, fnmatch
from  . import util

###
xuvtop = os.environ['XUVTOP']
#chInteractive=1
Defaults = util.defaultsRead()
Ip = util.ipRead()
MasterList = util.masterListRead()
#AbundanceAll = util.abundanceRead(abundancename = Defaults['abundfile'])
IoneqAll = util.ioneqRead(ioneqname = Defaults['ioneqfile'])
# gets the ChianitPy version
# gets the version of the CHIANTI database
ChiantiVersion = util.versionRead()
keywordArgs = ['temperature','eDensity','hDensity', 'pDensity','radTemperature', 'rStar', 'distance']
#
abunddir = os.path.join(xuvtop,'abundance')
filelist = os.listdir(abunddir)
#
abundList = []
for one in filelist:
    fname = os.path.join(abunddir,one)
    if os.path.isfile(fname):
        abundList.append(os.path.splitext(one)[0])
#for one in abundList:
#    print(one)
Abundance = {abundList[0]:util.abundanceRead(abundancename = abundList[0])}
for one in abundList[1:]:
    Abundance[one] = util.abundanceRead(abundancename = one)

print(' importing data')

