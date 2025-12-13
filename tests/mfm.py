from __future__ import division
from ddt import ddt
import numpy as np
import unittest
import xmlrunner
import PKDGRAV as msr
import h5py
import os
from pathlib import Path
from PKDGRAV import HYDRO_MODEL, PKD_FIELD


# Need to think how to save this outside of here
def save(data_dict, filename, time=0., bMemSoft=False):
    """!
    Save data into a hdf5 snapshot prepared for PKDGRAV3

    This function takes the labels from the data_dict structure, such that
    it must already represent the tree structure of the snapshot. For example:

    PartType0
        |- Masses
        |- Coordinates
        |- Velocities
        |- Temperature

    PartType1
        |- Masses
        |- Coordinates
        |- Velocities

    NOTE: Even in the case of having the same mass for all particles,
    the 'Masses' array must be present to make it easier for PKDGRAV (and us).
    NOTE2: For this function to work when running simulations with gravity,
    dSoft must be given to the parameters files.

    @param data_dic Dictionary with the data and structure of the file to save.
    @param filename Name of the file which will be created.
    @param time Time of the snapshot. If using comoving coordinates, the
        scale-factor is to be given instead.
    @param bMemSoft If set to True, the softening of the particles is directly
        stored in the PartType group as a dataset. If set to False, it stores
        the softening in classes. In this case the maximum number of unique
        softening values has to be less than 256. This parameter must also be
        set accordingly in the parameter file of PKDGRAV3
    """
    snap = h5py.File(filename, 'w')

    gHeader = snap.create_group('Header')
    gHeader.attrs['Time'] = time

    dTypes = {'default': "<f4", 'Velocites': "<f4", 'Coordinates': "<f8",
              'ParticlesIDs': "<u8"}

    class_type = {'names': ['class', 'mass', 'soft', 'start'],
                  'formats': ['u1', '<f8', '<f8', '<u8'],
                  'offsets': [0, 16, 24, 8], 'itemsize': 40}

    classes = np.zeros(256, dtype=class_type)
    iType0 = 0
    iStart = 0

    for partType in data_dict.keys():
        gPartType = snap.create_group(partType)

        # We can 'safely' assume that at least 'Coordinates' is present
        npart = data_dict[partType]['Coordinates'].shape[0]

        # If not provided, make sure that the particleIDs are stored to assert
        # that unique IDs are generated
        if 'ParticleIDs' not in data_dict[partType].keys():
            data_dict[partType]['ParticleIDs'] = np.arange(npart) + iStart

        for var in data_dict[partType].keys():
            if var in dTypes:
                var_type = dTypes[var]
            else:
                var_type = dTypes['default']

            gPartType.create_dataset(var, data=data_dict[partType][var],
                                     dtype=var_type)

        # We use the class to store information about the softening of the
        # particles in the case that they are not provided by the param file.
        # This means that a 'Softening' field must be present.
        # In order to avoid too many classes, there has to be less than
        # 256 different softening values
        if 'Softening' in data_dict[partType].keys() and not bMemSoft:
            print("Storing softening as classes")
            softs, inverse = np.unique(data_dict[partType]['Softening'],
                                       return_inverse=True)
            if softs.shape[0] + iType0 > 255:
                print("Too many softening values!")
                exit()

            for iType in range(softs.shape[0]):
                classes[iType0 + iType]['class'] = iType0 + iType
                classes[iType0 + iType]['soft'] = softs[iType]
                classes[iType0 + iType]['mass'] = 0.  # Should be provided as 'Masses' field
                classes[iType0 + iType]['start'] = 0  # Not used

            inverse += iType0
            gPartType.create_dataset('class', data=inverse.astype("|u1"))

            iType0 += softs.shape[0]
        else:
            iType = iType0

            gPartType.create_dataset('class', data=np.full(npart, iType, dtype="|u1"))

            classes[iType]['class'] = iType
            classes[iType]['soft'] = 0.  # Should be provided by dSoft in parameter file
            classes[iType]['mass'] = 0.  # Should be provided as 'Masses' field
            classes[iType]['start'] = 0  # Not used

            iType0 += 1

        iStart += npart

    for partType in data_dict.keys():
        snap[partType].create_dataset("classes", data=classes[:iType0])

    snap.close()


@ddt
class TestMFM(unittest.TestCase):
    n = 32

    def setUp(self):

        # Size of the box
        L = 1.

        # The particles are equally distributed along the three axis, from -0.5 up to 0.5
        N = self.n**3

        x = np.linspace(-0.5, 0.5, num=self.n, endpoint=False)
        mesh = np.meshgrid(x, x, x)

        x = mesh[0].flatten()
        y = mesh[1].flatten()
        z = mesh[2].flatten()
        del mesh

        gamma = 5. / 3.  # gas adiabatic constant
        rho0 = 1.
        mass0 = rho0 * L * L * L / N

        # Durier & Dalla Vecchia
        rho0 = 1.
        E0 = 1e-5  # Total
        E = 1. / mass0
        Uint0 = E0 * (gamma - 1.) / N / mass0

        # Position of the blast
        offsetX = 0.
        offsetY = 0.
        offsetZ = 0.

        # We set the particle mass to have rho0
        mass = np.ones_like(x) * mass0

        Uint = np.ones(N) * Uint0

        r = np.sqrt((x - offsetX)**2 + (y - offsetY)**2 + (z - offsetZ)**2)
        pos = np.stack((x, y, z), axis=1)
        vel = np.zeros_like(pos)

        N_center = 7
        ord_r = np.argsort(r)
        mask = ord_r[:N_center]

        h = 0.015625 * 0.01

        fac = 1. / np.pi
        h1 = 1. / h
        q = r * h1
        fac = fac * h1 * h1 * h1
        tmp2 = 2. - q
        val = np.zeros_like(q)
        mask = q < 2.
        val[mask] = 0.25 * tmp2[mask] * tmp2[mask] * tmp2[mask]
        mask = q < 1.
        val[mask] = 1 - 1.5 * q[mask] * q[mask] * (1 - 0.5 * q[mask])
        w = val * fac
        print("Kicking %d particles" % np.sum(w > 0))

        # Thermal energy injection
        norm = w / np.sum(w)
        print("Background temperature: ", np.mean(Uint * (gamma - 1.)))
        E_perP = E * norm
        Uint += E_perP
        print("Central Temperatures: ", Uint[mask] * (gamma - 1.))

        print("Thermal energy: ", np.sum(Uint))

        gas_dict = {'Coordinates': pos,
                    'Masses': mass.flatten(),
                    'Velocities': vel,
                    'InternalEnergy': Uint.flatten(),
                    'SmoothingLength': np.full(N, np.power(3. / (4. * np.pi * N), 1. / 3.))}

        save({'PartType0': gas_dict}, 'sedov.hdf5')

    def test_MFMSedov(self):
        msr.simulate(
            hydro_model=HYDRO_MODEL.MFM,
            bDoGravity=0,
            bMemUnordered=1,
            bNewKDK=1,
            achInFile='sedov.hdf5',
            nSmooth=64,
            dNeighborsStd=0.64,
            bGlobalDt=0,
            bMemMass=1,
            iStartStep=0,
            nSteps=20,
            dDelta=5e-3,
            dEtaCourant=0.1,
            dCFLacc=0.05,
            bHDF5=1,
            bVDetails=1,
            bParaRead=0,
            bParaWrite=0,
            bComove=0,
            bPeriodic=1,
            bEwald=0,
        )
        mass = msr.get_array(field=PKD_FIELD.FIELD_MASS)
        vel = msr.get_array(field=PKD_FIELD.FIELD_VELOCITY)
        pos = msr.get_array(field=PKD_FIELD.FIELD_POSITION)

        linMom = np.sqrt(((vel * mass[:, np.newaxis]).sum(axis=0)**2).sum())
        cross_rv = np.cross(pos, vel)
        L_particles = mass[:, np.newaxis] * cross_rv
        L_total = np.sqrt(np.sum(L_particles.sum(axis=0)**2))
        com = np.sqrt(((pos * mass[:, np.newaxis]).sum(axis=0)**2).sum())

        shock_pos = pos[np.argmax(vel[:, 0]), :][0]

        print(f"""
        Linear mometum: {linMom}
        Angular momentum: {L_total}
        Center of mass: {com}
        Shock position: {shock_pos}
        """)

        assert linMom < 5e-6
        assert L_total < 5e-6
        assert com < 1. / self.n
        try:
            # This would be the result of the full scheme (2nd order)
            assert np.isclose(shock_pos, 0.456, atol=0.01)
        except:
            # This would be the result of the "lean" scheme (1st order)
            assert np.isclose(shock_pos, 0.431, atol=0.01)

    def tearDown(self):
        os.remove('sedov.hdf5')


if __name__ == '__main__':
    print('Running test')
    unittest.main(verbosity=2, testRunner=xmlrunner.XMLTestRunner(
        output='test-reports'))
