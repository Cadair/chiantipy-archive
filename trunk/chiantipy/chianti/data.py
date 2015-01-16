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
import os

from . import io

###
xuvtop = os.environ['XUVTOP']
#chInteractive=1
Defaults = io.defaultsRead()
Ip = io.ipRead()
MasterList = io.masterListRead()
IoneqAll = io.ioneqRead(ioneqname = Defaults['ioneqfile'])
# gets the ChianitPy version
# gets the version of the CHIANTI database
ChiantiVersion = io.versionRead()
keywordArgs = ['temperature','eDensity','hDensity', 'pDensity','radTemperature', 'rStar', 'distance']
#
AbundanceDefault = io.abundanceRead(abundancename = Defaults['abundfile'])
abunddir = os.path.join(xuvtop,'abundance')
filelist = os.listdir(abunddir)
#
AbundanceList = []
for one in filelist:
    fname = os.path.join(abunddir,one)
    if os.path.isfile(fname):
        AbundanceList.append(os.path.splitext(one)[0])
#for one in abundList:
#    print(one)
Abundance = {AbundanceList[0]:io.abundanceRead(abundancename = AbundanceList[0])}
for one in AbundanceList[1:]:
    Abundance[one] = io.abundanceRead(abundancename = one)

