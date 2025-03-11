from sys import argv,exit
import PKDGRAV as msr
from PKDGRAV import OUT_TYPE
import numpy as np

# Parse parameters. Name and grid size are required; box size is optional
if len(argv)<=2:
    print("Usage: {} <filename> <dTau> [nMinMembers]".format(argv[0]))
    exit(1)
name=argv[1]
dTau = float(argv[2])
nMinMembers = 10 if len(argv)<4 else int(argv[3])

# Load the file and setup the tree, then measure the power
time = msr.load(name,bFindGroups = True,bMemGlobalGid = True)
msr.domain_decompose()
msr.build_tree()

msr.fof(dTau,nMinMembers)
msr.reorder()

# Compare group assignments to some standard reference.
msr.write_array(f"{name}.fof",OUT_TYPE.OUT_GLOBALGID_ARRAY)

#b = msr.GetArray(field=msr.FIELD_GLOBAL_GID,time=time)
