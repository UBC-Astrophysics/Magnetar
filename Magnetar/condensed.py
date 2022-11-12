    
from numpy import sin, cos, arccos, radians, exp, expm1, log, sqrt, minimum, maximum, where, pi, logspace, array, zeros
from Magnetar.utils import atmosphere
from scipy.integrate import simps
from Magnetar.simple_atmospheres import bbfunk


import ctypes
_sum = ctypes.CDLL('Magnetar/emissivity.so')
_sum.emissivity.argtypes = (ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double), ctypes.c_int,
                            ctypes.c_double,ctypes.c_double, ctypes.c_double,ctypes.c_double, ctypes.c_double,ctypes.c_int, ctypes.c_int,
                            ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double))
_sum.emissivity.restype = None
 
class condensed_surface(atmosphere):
    def __init__(self,effective_temperature,mag_strength,mag_inclination,density,fixed_ions=True):
        self.effective_temperature=effective_temperature
        self.mag_strength=mag_strength
        self.mag_inclination=mag_inclination
        self.fixed_ions=fixed_ions
        self.dens=density
        self.Z=26.
        self.A=56.
        self.surface_temperature=self.effective_temperature
        self.adjust_surface_temperature()
    def __str__(self):
        outstring='''#
# class condensed_surface
#
# effective_temperature %12g keV
# surface_temperature   %12g keV
# mag_strength          %12g Gauss
# mag_inclination       %12g radians
# fixed_ions            %12s
# density               %12g g/cc
# Z                     %12g
# A                     %12g
#       
''' % (self.effective_temperature, self.surface_temperature, self.mag_strength, self.mag_inclination, self.fixed_ions, self.dens, self.Z, self.A )
        return outstring+atmosphere.__str__(self)
    #
    # This computes the emissivity from a condensed surface using the approximate treatment
    # by Potekhin et al. (2012, A&A, 546, A121) https://arxiv.org/pdf/1208.6582.pdf
    #
    

    def _emissivity_xo(dataarray,mag_strength,mag_inclination,dens,fixed_ions,Z,A,mixing=True):
        thetak=radians(dataarray[-3])
        phik=radians(dataarray[-2])
        ene=dataarray[-1]
        lx=len(thetak)
        array_type=ctypes.c_double * lx
        eXout=array_type(*zeros(lx))
        eOout=array_type(*zeros(lx))
             
        _sum.emissivity(array_type(*thetak),array_type(*phik),array_type(*ene),ctypes.c_int(lx),
                                  ctypes.c_double(mag_strength),ctypes.c_double(radians(mag_inclination)),ctypes.c_double(dens),
                                  ctypes.c_double(Z),ctypes.c_double(A),ctypes.c_int(1 if fixed_ions else 0),ctypes.c_int(1 if mixing else 0),eXout,eOout)
        return array(eXout),array(eOout)


    def emissivity_xo(self,dataarray):
        #a,b = condensed_surface._emissivity_xo(array(dataarray),self.mag_strength,self.mag_inclination,self.dens,self.fixed_ions,self.Z,self.A)
        #print(a,b)
        #return a,b
        return condensed_surface._emissivity_xo(array(dataarray),self.mag_strength,self.mag_inclination,self.dens,self.fixed_ions,self.Z,self.A)
    def calcIQ(self, dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=bbfunk(array(dataarray[-1]),self.surface_temperature)
        return (eo+ex)*bbintens,(eo-ex)*bbintens
    def xintensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=bbfunk(array(dataarray[-1]),self.surface_temperature)
        return ex*bbintens
    def ointensity(self,dataarray):
        ex,eo = self.emissivity_xo(dataarray)
        bbintens=bbfunk(array(dataarray[-1]),self.surface_temperature)
        return eo*bbintens
