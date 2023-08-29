# distutils: language = c++
# cython: always_allow_keywords=True
import cython

@cython.cclass
class classic_theta_switch:
    """
    This is the classic theta switch. Use this by setting
    the parameter `dTheta` to an instance of this class.

    :param number before20: theta to use before redshift 20 
    :param number between20and2: theta to use between redshift 20 and 2
    :param number after2: theta to use after redshift 2
    :return: theta
    :rtype: number
    """
    before20: cython.double
    between20and2: cython.double
    after2: cython.double
    def __init__(self,before20=0.4,between20and2=0.55,after2=0.7):
        self.before20 = before20
        self.between20and2 = between20and2
        self.after2 = after2
    def __call__(self,a,**kwargs):
        z = 1.0 / a - 1.0
        if z>20.0:
            return self.before20
        elif z>2.0:
            return self.between20and2
        else:
            return self.after2
    def __repr__(self):
        return f"classic_theta_switch({self.before20},{self.between20and2},{self.after2})"

@cython.cclass
class classic_replicas_switch:
    """
    This is the classic replicas switch. Use this by setting
    the parameter `nReplicas` to an instance of this class.

    :param number theta_switch: theta threshold to switch between 1 and 2 replicas
    :return: number of replicas
    :rtype: integer
    """
    theta_switch: cython.double
    def __init__(self,theta_switch=0.52):
        self.theta_switch = theta_switch
    def __call__(self,theta,**kwargs):
        return 1 if theta >= self.theta_switch else 2
    def __repr__(self):
        return f"classic_replicas_switch({self.theta_switch})"

