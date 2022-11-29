try:
    from numba import jit
except:
    def jit(**kwargs):
        def jitd(func):
            def helper(*args, **kwargs):
                return func(*args, **kwargs)
            return helper
        return jitd
    
from Magnetar.utils import atmosphere
import numpy as np

@jit(nopython=True,parallel=True)
def bbfunk( ee, tt):  # per mode
    return 208452.792 * ee**3 / np.expm1(ee / tt) / 2


def bb_atmo_f(temp):
    return atmosphere() * (lambda ee: bbfunk(ee,temp))


#
# blackbody atmosphere (convenience class with parallel structure to condensed_surface)
#
class bb_atmo(atmosphere):
    def __init__(self,teff,mag_strength,mag_inclination,*args,ofraction=0.5,**kwargs):
        self.teff  = teff
        self.ofraction = np.clip(ofraction,0,1)
        self.xtemp = (1-self.ofraction)**0.25
        self.otemp = self.ofraction**0.25
        rat=teff*(2.0/(self.xtemp**4+self.otemp**4))**0.25
        self.xtemp*=rat
        self.otemp*=rat
        self.mag_inclination = mag_inclination
        self.mag_strength=mag_strength
    def __str__(self):
        outstring='''#
# class bb_atmo
#
# effective_temperature %12g keV
# O fraction            %12g
# X temperature         %12g keV
# O temperature         %12g keV
# mag_strength          %12g Gauss [not used]
# mag_inclination       %12g radians [not used]
#\n''' % (self.teff, self.ofraction, self.xtemp, self.otemp, self.mag_strength, self.mag_inclination)
        return outstring+atmosphere.__str__(self)

    def xintensity(self, dataarray):
        if self.xtemp>0:
            return bbfunk(dataarray[-1], self.xtemp)
        else:
            return 0*dataarray[-1]

    def ointensity(self, dataarray):
        if self.otemp>0:
            return bbfunk(dataarray[-1], self.otemp)
        else:
            return 0*dataarray[-1]

#
# pure X blackbody atmosphere (convenience function with parallel structure to condensed_surface)
#
def bb_atmo_purex(teff,mag_strength,mag_inclination,*args,**kwargs):
    return bb_atmo(teff,mag_strength,mag_inclination,ofraction=0)

#
# Thompson_Kostenko_Magnetosphere (convenience class function with parallel structure to condensed_surface)
#
class Thompson_Kostenko_Magnetosphere(atmosphere):
    def __init__(self,teff,mag_strength,mag_inclination,normalization=1,alpha=1):
        self.mag_inclination=mag_inclination
        self.normalization=normalization
        self.alpha=alpha
    def __str__(self):
        outstring='''#
# class modified_bb_atmo
#
# normalization %12g at 10 keV
# alpha         %12g
''' % (normalization,alpha)
        return outstring+atmosphere.__str__(self)
    @jit(nopython=True,parallel=True)
    def _ointensity(dataarray,mag_inclination,normalization,alpha):
        ee=dataarray[-1]/10.0
        coskb2=(np.cos(np.radians(mag_inclination)
                           ) * np.cos(np.radians(dataarray[-3])) +
                    np.sin(np.radians(mag_inclination)
                           ) * np.sin(np.radians(dataarray[-3])
                                  ) * np.cos(np.radians(dataarray[-2])))**2
        return normalization*ee**alpha*(1.0-coskb2)
    def ointensity(self, dataarray):
        return Thompson_Kostenko_Magnetosphere._ointensity(np.array(dataarray),self.mag_inclination,self.normalization,self.alpha)
    def xintensity(self, dataarray):
        return 0.0 # np.zeros((len(dataarray[-1])))
      
        

#
# modified blackbody atmosphere (convenience class with parallel structure to condensed_surface)
#
# freq_power quantifies how the opacity depends on photon energy: kappa goes as (1/energy)**freq_power
# freq_power is 2 (E<T) to 3 (E>t) for free-free opacity and zero for scattering
#
# sigma_power quantifies how the temperature depends on column density: T goes as sigma**sigma_power
# sigma_power is (alpha+1)/(4+alpha-beta) if opacity is proportional to density**alpha temperature**beta
#
#    unmagnetized free-free : alpha=1, beta=-3.5 -> 0.23529411764705882 (4/17)
#    magnetized free-free   : alpha=1, beta=-1.5 -> 0.3076923077 (4/13)         
#    electron-scattering    : alpha=0, beta=0    -> 0.25 (4/16)
#
# Based on the power-law atmospheres for neutron stars in
#
# https://ui.adsabs.harvard.edu/abs/1998MNRAS.300..599H for the magnetized case and
#
# https://ui.adsabs.harvard.edu/abs/1984ApJ...287..244H for unmagnetized case
#
#

class modified_bb_atmo(atmosphere):
    def __init__(self,teff,mag_strength,mag_inclination,*args,freq_power=2,sigma_power=4./13.,kb_suppress=True,limb_darkening=True,**kwargs):
        self.effective_temperature  = teff
        self.freq_power = freq_power
        self.sigma_power = sigma_power
        self.surface_temperature = teff
        self.mag_inclination = mag_inclination
        self.mag_strength=mag_strength
        self.ecyc = mag_strength/4.4e13*511.
        self.kb_suppress=kb_suppress
        self.limb_darkening=limb_darkening
        self.adjust_surface_temperature()
        
    def __str__(self):
        outstring='''#
# class modified_bb_atmo
#
# effective_temperature %12g keV
# surface_temperature   %12g keV [tau=1 for upgoing O-mode with E=surface_temperature, neglecting k dot b suppression]
# mag_strength          %12g Gauss
# mag_inclination       %12g degrees
# cyclotron energy      %12g keV
# k dot b supression    %12s 
# limb darkening        %12s 
# freq_power            %12g [cross-section goes as 1/freq**freqpower for O-mode]
# sigma_power           %12g [temperature goes as column-density**sigma_power]
#
# ./atm -o XXX_%g_%g_%g -B %g -T %g -b %g -m 11 -p 5 -a 1.5 -M 2 -D 5 
#
''' % (self.effective_temperature, self.surface_temperature, self.mag_strength, self.mag_inclination, self.ecyc, 
       self.kb_suppress, self.limb_darkening, self.freq_power, self.sigma_power,
       np.log10(self.mag_strength),np.log10(self.effective_temperature)+7.06462726,self.mag_inclination,
       np.log10(self.mag_strength),np.log10(self.effective_temperature)+7.06462726,self.mag_inclination)
        return outstring+atmosphere.__str__(self)
    
    # Auxiliary functions defined for numba
    @jit(nopython=True,parallel=True)
    def _xintensity(dataarray,surface_temperature,freq_power,sigma_power,ecyc,limb_darkening):
        ee=dataarray[-1]
        sigmax=np.abs(surface_temperature/ee)**freq_power
        sigmax*=np.where(ee<ecyc,(ee/ecyc)**2,1)
        if limb_darkening:
            tempx=surface_temperature*np.abs(np.cos(np.radians(dataarray[-3]))/sigmax)**sigma_power
        else:
            tempx=surface_temperature*np.abs(1/sigmax)**sigma_power
        return bbfunk(ee,tempx)

    @jit(nopython=True,parallel=True)
    def _ointensity(dataarray,surface_temperature,freq_power,sigma_power,ecyc,limb_darkening,kb_suppress,mag_inclination):
        ee=dataarray[-1]
        sigmao=np.abs(surface_temperature/ee)**freq_power
        if kb_suppress:
            coskb2=(np.cos(np.radians(mag_inclination)
                       ) * np.cos(np.radians(dataarray[-3])) +
                np.sin(np.radians(mag_inclination)
                       ) * np.sin(np.radians(dataarray[-3])
                                  ) * np.cos(np.radians(dataarray[-2])))**2
            sigmao*=np.where(ee<ecyc,(1-coskb2)+coskb2*(ee/ecyc)**2,1)
        if limb_darkening:
             tempo=surface_temperature*np.abs(np.cos(np.radians(dataarray[-3]))/sigmao)**sigma_power
        else:
             tempo=surface_temperature*np.abs(1/sigmao)**sigma_power
        return bbfunk(ee,tempo)


    
    def xintensity(self, dataarray):
        return modified_bb_atmo._xintensity(np.array(dataarray),self.surface_temperature,self.freq_power,self.sigma_power,self.ecyc,self.limb_darkening)
 
    def ointensity(self, dataarray):
        return modified_bb_atmo._ointensity(np.array(dataarray),self.surface_temperature,self.freq_power,self.sigma_power,
                                            self.ecyc,self.limb_darkening,self.kb_suppress,self.mag_inclination)


_absstring='''0.20738 68.644
0.21805 68.561
0.22926 68.789
0.24106 68.881
0.25345 69.235
0.26649 69.087
0.28100 69.568
0.28929 69.087
0.29058 72.665
0.29726 74.831
0.31255 74.557
0.32863 74.831
0.34553 74.934
0.36330 75.210
0.38199 75.279
0.40162 76.639
0.42226 78.779
0.44398 78.779
0.46681 79.087
0.49083 79.178
0.51607 79.178
0.53402 79.679
0.53612 104.13
0.53627 92.947
0.53760 85.408
0.54625 117.60
0.57520 119.82
0.60477 121.62
0.63585 123.34
0.66854 125.05
0.70291 126.61
0.72357 140.37
0.77686 143.31
0.81680 145.30
0.85291 148.61
0.88636 168.76
0.93191 172.04
0.97981 174.75
1.0302 177.54
1.0831 180.25
1.1388 182.79
1.1974 185.03
1.2589 187.07
1.3012 190.39
1.3417 201.38
1.4107 203.66
1.4832 205.91
1.5594 208.58
1.6396 211.27
1.7239 213.17
1.8076 215.95
1.8881 232.55
1.9852 234.69
2.0873 236.86
2.1946 239.15
2.3075 241.08
2.4206 242.98
2.5273 261.97
2.6573 264.93
2.7939 267.19
2.9376 269.47
3.0886 271.52
3.2473 278.65
3.4143 281.55
3.5898 283.43
3.7745 285.12
3.9685 288.41
4.1724 294.98
4.3870 296.54
4.6126 298.11
4.8498 299.48
5.0993 301.13
5.3615 302.31
5.6373 303.35
5.9272 305.73
6.2320 307.98
6.5524 310.32
6.8894 312.47
7.0669 480.08
7.0698 590.78
7.0745 447.52
7.0754 423.34
7.0763 400.46
7.0772 378.83
7.0782 358.36
7.0791 338.99
7.0800 320.68
7.0873 544.92
7.0880 522.03
7.0906 616.51
7.2445 648.30
7.6031 653.56
7.9940 661.25
8.4044 687.76
8.8366 695.34
9.2909 703.21
9.7463 710.11
'''

_absarr=np.fromstring(_absstring,sep=' ').reshape((97,2)).T

def absorption(ee,nH):
    ssabs=np.interp(ee,_absarr[0],_absarr[1])/(ee)**3*1e-24
    return np.exp(-nH*ssabs)    
