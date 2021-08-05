#!/usr/bin/env python3
'''
    Zodiacal light calculator
       Based on Tsumura-san's table
        Ver. 1.0  2020/11/02   H. Akitaya
        Ver. 1.1  2021/07/06   H. Akitaya; comments implemented
          Unit: lambda*I(lambda) [nW/m^2/sr]
        Ver. 1.2  2021/07/28   H. Akitaya; output with astropy.units
        Ver. 1.3  2021/07/29   H. Akitaya; integration in wavelength range; photon flux
        Ver. 1.4  2021/07/30   H. Akitaya; interpolation using scipy
        Ver. 1.5  2021/08/05   H. Akitaya; change wavelength range for o and j-bands.

   Sample usage:
       from zodiacallight import ZodiacalLight
       zd = ZodiacalLight()
       zd(1.2, 0)
       714.1626404659994
'''

import os
import sys
import astropy.units as u
import astropy.constants as c
import numpy as np
from scipy import integrate
from scipy import interpolate

# Table of zodiacal light wavelength dependency
ZL_SPEC_FN = os.path.dirname(__file__) + '/' + 'ZLspectrum.txt'  

# Table of zodiacal light latitude dependency
ZL_ECLIPLAT_FN = os.path.dirname(__file__) + '/' + 'ZL_ecliptic_profile.txt'

# normalization flux (flux at ecliptic latitude at 0 deg)
ZL_ECLIPLAT_NORM = 691.82538

# Zodiacal light photon number flux at zlat=0 deg.
# Unit: photons/s/str/m^2
# Calculated in 2021-07-30 using:
# zd = ZodiacalLight()
# wl1 = [0.5, 1.0, 1.5, 2.0]
# dwl = 0.5
# for wl in wl1:
#    phflux = zd.get_zl_photonflux_wlinteg(wl, wl+dwl, 0.0, unit=True)
#    print('{:.2f} {:.2f} {:.7e}'.format(wl, wl+dwl, phflux[0]))

ZL_PHOTONFLUX_BAND = {'o': 1.9220403980285938e+12,
                      'j': 2.143305343449887e+12,
                      'h': 1.1003014e+12,
                      'k': 5.6836829e+11,
                      }

class ZodiacalLight(object):
    ''' Class for Zodiacal Light Calculation.
    '''
    
    def __init__(self, data_read=True):
        self.zl_spec = []
        self._zl_spec_wl = []
        self._zl_spec_flux = []
        self.zl_ecliplat = []
        self._zl_ecliplat_lat = []
        self._zl_ecliplat_scale = []
        self._f_wscipy_get_zl_wavelength_at = None
        self._f_wscipy_get_zl_normalized_ecliplat_at = None
        
        if data_read is True:
            self.read_data_all()

    def __call__(self, wl, lat, unit=False, photon=False, w_scipy=True):
        if unit is False:
            return self.get_zl_at(wl, lat, w_scipy=w_scipy)
        else:
            return self.get_zl_at(wl, lat, w_scipy=w_scipy) *(1e-9*u.W/(u.m**2)/u.sr)

    def __str__(self):
        return 'ZodiacalLight(wl[um], latitude[deg]) [nW/m^2/str]'

    def read_data_zl_spec(self, fn=ZL_SPEC_FN):
        ''' Read zodiacal light table (wavelength dependency).
        '''
        with open(fn, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                items = line.strip().split()
                if len(items) != 3:
                    continue
                self.zl_spec.append([float(items[0]), float(items[1]),
                                     float(items[2])])
            _tmp_array = np.array(self.zl_spec)
            self._zl_spec_wl = _tmp_array[:,0]
            self._zl_spec_flux = _tmp_array[:,1]
            self._f_wscipy_get_zl_wavelength_at = \
                interpolate.interp1d(self._zl_spec_wl,
                                     self._zl_spec_flux,
                                     kind='cubic')

    def read_data_zl_ecliplat(self, fn=ZL_ECLIPLAT_FN):
        ''' Read zodiacal light table (latitude dependency).
        '''
        with open(fn, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                items = line.strip().split()
                if len(items) != 2:
                    continue
                self.zl_ecliplat.append([float(items[0]), float(items[1])])
            _tmp_array = np.array(self.zl_ecliplat)
            self._zl_ecliplat_lat = _tmp_array[:,0]
            self._zl_ecliplat_scale = _tmp_array[:,1]
            self._f_wscipy_get_zl_normalized_ecliplat_at = \
                interpolate.interp1d(self._zl_ecliplat_lat,
                                     self._zl_ecliplat_scale,
                                     kind='cubic')

    def read_data_all(self):
        ''' Read zodiacal light tables.
        '''
        self.read_data_zl_spec()
        self.read_data_zl_ecliplat()

    def get_zl_wavelength_at(self, wl, w_scipy=True):
        ''' Calculate scale of zodiacal light at wavelength wl.
        '''
        if w_scipy is True:  # Use scipy interpolation
            return self._f_wscipy_get_zl_wavelength_at(wl)
        else:
            return self._function_from_table(self.zl_spec, wl)

    def get_zl_normalized_ecliplat_at(self, ecliplat, w_scipy=True):
        ''' Calculate scale of zodiacal light at latitude ecliplat.
        '''
        if w_scipy is True:
            return self._f_wscipy_get_zl_normalized_ecliplat_at(ecliplat) / ZL_ECLIPLAT_NORM
        else:
            return self._function_from_table(self.zl_ecliplat,
                                                ecliplat) / ZL_ECLIPLAT_NORM

    def get_zl_at(self, wl, ecliplat, unit=False, w_scipy=True):
        ''' Calculate zodiacal light flux at wavelength wl and latitude ecliplat.
        '''
        try:
            f =  self.get_zl_wavelength_at(wl, w_scipy=w_scipy) * \
                self.get_zl_normalized_ecliplat_at(ecliplat, w_scipy=w_scipy)
        except TypeError:
            sys.stderr.write('Zodiacal light calcuration error\n')
            raise(TypeError)
        if unit is False:
            return f
        else:
            return 1e-9 * f * (u.W/u.m**2/u.sr)  # [W/m^2/str]

    def get_zl_photonflux_at(self, wl, ecliplat, unit=False):
        ''' Calculate zodiacal light photon flux at wavelength wl 
        and latitude ecliplat.
        '''
        f = self.get_zl_at(wl, ecliplat, unit=True)
        pe = ZodiacalLight.get_photon_energy(wl, unit=True)
        pf = (f/pe).si
        if unit is False:
            return pf.value
        else:
            return pf

    def get_zl_photonflux_wlinteg(self, wl1, wl2, ecliplat, unit=False):
        ''' Calculate zodiacal light photon flux integrated between 
        wavelength wl1 and wl2, atlatitude ecliplat.
        '''
        pf_integ = integrate.quad(lambda x: self.get_zl_photonflux_at(
            x, ecliplat, unit=False)/(x*1e-6)*1e-6, wl1, wl2)
        return pf_integ

    def get_zl_photonflux_band(self, band, ecliplat, unit=False, w_scipy=True):
        if not band in ZL_PHOTONFLUX_BAND:
            return None
        phflux = ZL_PHOTONFLUX_BAND[band] * self.get_zl_normalized_ecliplat_at(ecliplat, w_scipy=w_scipy)
        if unit is False:
            return phflux
        else:
            return phflux * u.photon / u.m**2 / u.s / u.sr
            

    def get_effective_wl_for_photonflux(self, wl1, wl2, unit=False):
        ''' Calculate effective wavelength within the flat efficiency.
        '''
        ecliplat = 0.0  # Not dependent on the ecliptic latitude in this calc.
        pf_integ, _ = integrate.quad(lambda x: self.get_zl_photonflux_at(
            x, ecliplat, unit=False)/x, wl1, wl2)
        pf_wl_integ, _ = integrate.quad(lambda x: self.get_zl_photonflux_at(
            x, ecliplat, unit=False)*1e-6, wl1, wl2)
        wl_eff = pf_wl_integ / pf_integ
        if unit is False:
            return wl_eff
        else:
            return wl_eff * u.m

    def _function_from_table(self, table_data, x):
        ''' Linear interporation for discrete values in a table.
        '''
        if len(table_data) < 2:
            raise ValueError
        xf_before = table_data[0][0]
        yf_before = table_data[0][1]
        for elms in table_data[1:]:
            if len(elms) < 2:
                raise ValueError
            xf = elms[0]
            yf = elms[1]
            if (xf_before <= x) and (x < xf):
                y = ZodiacalLight.function_interpolate(xf_before, xf, yf_before, yf, x)
                return y
            xf_before = xf
            yf_before = yf
        return None

    @staticmethod
    def function_interpolate(xf_before, xf, yf_before, yf, x):
        ''' Calculate linear interporated value between two vectors.
        '''
        y = (yf-yf_before)/(xf-xf_before)*(x-xf_before) + yf_before
        return y

    @staticmethod
    def get_photon_energy(wl, unit=True):
        ''' Get phtoton energy fron wavelength wl (m).
        unit: True -> with astropy.units (J), False: w/o units.
        '''
        photon_energy = (c.h * c.c / (wl * u.um)).to(u.J)
        if unit is True:
            return photon_energy
        else:
            return photon_energy.value

if __name__ == '__main__':
    if len(sys.argv) < 3:
        sys.exit(1)
    wl = float(sys.argv[1])  # wavelength (um)
    ecliplat = float(sys.argv[2])  # latitude (deg)
    zl = ZodiacalLight()
    zl.read_data_all()
    flux = zl.get_zl_at(1.25, 40.0)
    # print('{:13.4f}'.format(flux))

    lat = -90.0
    while (lat <= 90.0):
        flux = zl.get_zl_at(wl, lat)
        if flux is not None:
            print('{:13.4f} {:13.4f}'.format(lat, flux))
        lat += (5-1e-8)


'''
./zodiacallight.py 1.25 0 > zod_1.25.xy
./zodiacallight.py 0.75 0 > zod_0.75.xy
./zodiacallight.py 1.75 0 > zod_1.75.xy
./zodiacallight.py 2.25 0 > zod_2.25.xy
'''
