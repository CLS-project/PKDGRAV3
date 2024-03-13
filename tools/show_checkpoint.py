import PKDGRAV
from sys import argv,exit
from dill import load

if len(argv) != 2:
    print(f"Usage: {argv[0]} <checkpoint.pkl>")
    exit(1)

with open(argv[1],'rb') as fp:
    version = load(fp)
    if version != 1:
        print("Invalid checkpoint version")
        exit(1)
    print( f'species = {load(fp)}')
    print( f'classes = {load(fp)}')
    print( f'step    = {load(fp)}')
    print( f'steps   = {load(fp)}')
    print( f'time    = {load(fp)}')
    print( f'delta   = {load(fp)}')
    E = load(fp)
    U = load(fp)
    Told = load(fp)
    a = load(fp)
    s = load(fp)
    print(f"{'Parameters (explicitly set only):':<20}")
    for k in vars(a):
    	if (getattr(s,k,False)):
    		print("  ",f"{k:<20}",getattr(a,k))