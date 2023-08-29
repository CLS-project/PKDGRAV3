from sys import argv,exit
from MASTER import MSR # Once MSR is imported, simulation mode is no longer automatically entered
from math import pi

import numpy as np

# Parse parameters. Name and grid size are required; box size is optional
if len(argv)<=2:
    print("Usage: {} <filename> <dTau> [nMinMembers]".format(argv[0]))
    exit(1)
name=argv[1]
dTau = float(argv[2])
nMinMembers = 10 if len(argv)<4 else int(argv[3])

msr = MSR()
msr.setParameters(bFindGroups = True,bMemGlobalGid = True)

# Load the file and setup the tree, then measure the power
time = msr.Load(name)
msr.DomainDecomp()
msr.BuildTree()

msr.Fof(dTau,nMinMembers)
msr.Reorder()

# Compare group assignments to some standard reference.
msr.WriteArray(f"{name}.fof",113)  # 113 = GLOBAL Group ids

#b = msr.GetArray(field=17,time=1) # 17 = GLOBAL Group id (oGlobalGid)

