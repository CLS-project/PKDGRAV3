# distutils: language = c++
# cython: always_allow_keywords=True
import cython

@cython.cclass
cdef class classic_theta_switch:
    before20: cython.double
    between20and2: cython.double
    after2: cython.double
    def __init__(self,before20=0.4,between20and2=0.55,after2=0.7):
        self.before20 = before20
        self.between20and2 = between20and2
        self.after2 = after2
    def __call__(self,a,**kwargs):
        z = 1.0 / a - 1.0
        if z>=20.0:
            return self.before20
        elif z>=2.0:
            return self.between20and2
        else:
            return self.after2
