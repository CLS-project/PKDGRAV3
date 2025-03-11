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
    print("You need to choose 2 values from the following combinations:")
    print()
    print("1: Mass unit and Length unit")
    print("2: Mass unit and Velocity unit")
    print("3: Length unit and Velocity unit")
    print("4: Mass unit and Time unit")
    print("5: Length unit and Time unit")
    print("6: Velocity unit and Time unit")
    print("")
    print("Choice: ",end="")
    inputvalue = 0
    try:
        inputvalue = input()
        choice = int(inputvalue)
    except:
        print("Input {} is not a valid choice! Exiting".format(inputvalue))
        return
    if (choice < 1 or choice > 6):
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
    kpc = 1000.0 * pc
    Mpc = 1e6 * pc
# Velocity units
    cmps = 1.0
    mps = 100.0 * cmps
    kmps = 1000.0 * mps
# Time units
    s = 1.0
    min = 60.0 * s
    h = 60.0 * min
    d = 24.0 * h
    y = 365.25 * d
    ky = 1000.0 * y
    My = 1e6 * y
    Gy = 1e9 * y
    
    massUnitString = "Units for mass: g, kg, ME (Earth), MJ (Jupiter), MS (Sun)"
    lengthUnitString = "Units for length: cm, m, RE (Earth), AU, ly, pc, kpc, Mpc"
    velocityUnitString = "Units for velocity: cmps, mps, kmps"
    timeUnitString = "Units for time: s, min, h, d, y, ky, My, Gy"

    print()
    print("Enter the values as [value] * [unit] pairs.")
    print("[value] can be any valid python or numpy (np.xxx) expression, for example 1.0 / 62.5 or np.sqrt(8*np.pi) or 1.7256e16")
    print("A list of the possible units is provided before each prompt.")
    print("Units are actually just multipliers to the base units in cgs.")
    print("g for mass, cm for length, cmps for velocity and s for time correspond to 1.0.")
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
    elif choice == 4:
        print(massUnitString)
        print("Provide the mass unit: ",end="")
        massinput = input()
        massvalue = eval(massinput)
        print(timeUnitString)
        print("Provide the time unit: ",end="")
        timeinput = input()
        timevalue = eval(timeinput)
        massunit = massvalue / MS
        lengthunit = ((Gcgs * MS * massunit * timevalue**2.0) / kpc**3.0)**(1.0/3.0)
    elif choice == 5:
        print(lengthUnitString)
        print("Provide the length unit: ",end="")
        lengthinput = input()
        lengthvalue = eval(lengthinput)
        print(timeUnitString)
        print("Provide the time unit: ",end="")
        timeinput = input()
        timevalue = eval(timeinput)
        lengthunit = lengthvalue / kpc
        massunit = kpc**3.0 * lengthunit**3.0 / (Gcgs * MS * timevalue**2.0)
    elif choice == 6:
        print(velocityUnitString)
        print("Provide the velocity unit: ",end="")
        velocityinput = input()
        velocityvalue = eval(velocityinput)
        print(timeUnitString)
        print("Provide the time unit: ",end="")
        timeinput = input()
        timevalue = eval(timeinput)
        massunit = timevalue * velocityvalue**3.0 / (Gcgs * MS)
        lengthunit = timevalue * velocityvalue / kpc
    else:
        print("Unexpected choice: {}! Exiting!".format(choice))
        return
    
    print("Result (this can be copied direclty into the parameter file):")
    print("dKpcUnit = {}".format(lengthunit))
    print("dMsolUnit = {}".format(massunit))

if __name__ == "__main__":
    main()