#!/usr/bin/env python3

import numpy as np

def main():
    print("This tool lets you calculate the values for")
    print()
    print("dKpcUnit (length unit in kpc)")
    print("and")
    print("dMsolUnit (mass unit in solar masses)")
    print()
    print("used in pkdgrav3 to define the unit system together with G = 1.")
    print()
    print("You need to choose 2 values from the following:")
    print()
    print("1: Mass unit and Length unit")
    print("2: Mass unit and Velocity unit")
    print("3: Length unit and Velocity unit")
    print("")
    print("Choice: ",end="")
    choice = int(input())
    if (choice < 1 or choice > 3):
        print("Invalid choice! Exiting!")
        return

# we do everything in cgs
    Gcgs = 6.67430e-8
# Mass units
    g = 1.0
    kg = 1000 * g
    ME = 5.97237e24 * kg
    MJ = 1.8986e27 * kg
    MS = 1.9884e30 * kg
# Length units
    cm = 1.0
    m = 100.0 * cm
    RE = 6.3710084e6 * m
    AU = 1.495978707e11 * m
    ly = 9.4607e15 * m
    pc = 3.0857e16 * m
    kpc = 1000 * pc
    Mpc = 1000 * kpc
# Velocity units
    cmps = 1.0
    mps = 100.0 * cmps
    kmps = 1000.0 * mps
    
    massUnitString = "Units for mass: g, kg, ME (Earth), MJ (Jupiter), MS (Sun)"
    lengthUnitString = "Units for length: cm, m, RE (Earth), AU, ly, pc, kpc, Mpc"
    velocityUnitString = "Units for velocity: cmps, mps, kmps"

    print()
    print("Enter the values as [value] * [unit] pairs.")
    print("[value] can be any valid python or numpy (np.xxx) expression, for example 1.0 / 62.5 or np.sqrt(8*np.pi)")
    print("A list of the possible units is provided before each prompt.")
    print()

    if choice == 1:
        print(massUnitString)
        print("Provide the mass unit: ",end="")
        massinput = input()
        massvalue = eval(massinput)
        print(lengthUnitString)
        print("Provide the length unit: ",end="")
        lengthinput = input()
        lengthvalue = eval(lengthinput)
        lengthunit = lengthvalue / kpc
        massunit = massvalue / MS
    elif choice == 2:
        print(massUnitString)
        print("Provide the mass unit: ",end="")
        massinput = input()
        massvalue = eval(massinput)
        print(velocityUnitString)
        print("Provide the velocity unit: ",end="")
        velocityinput = input()
        velocityvalue = eval(velocityinput)
        massunit = massvalue / MS
        lengthunit = (Gcgs * MS * massunit) / (kpc * velocityvalue**2)
    elif choice == 3:
        print(lengthUnitString)
        print("Provide the length unit: ",end="")
        lengthinput = input()
        lengthvalue = eval(lengthinput)
        print(velocityUnitString)
        print("Provide the velocity unit: ",end="")
        velocityinput = input()
        velocityvalue = eval(velocityinput)
        lengthunit = lengthvalue / kpc
        massunit = kpc * lengthunit * velocityvalue**2.0 / (Gcgs * MS)
    else:
        print("Unexpected choice: {}! Exiting!".format(choice))
        return
    
    print("Result (this can be copied direclty into the parameter file):")
    print("dKpcUnit = {}".format(lengthunit))
    print("dMsolUnit = {}".format(massunit))

if __name__ == "__main__":
    main()