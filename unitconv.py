#!/usr/bin/env python3

import sys
import math

from astropy import constants as const
from astropy import units as u

class UnitConv(object):
    def __init__(self):
        self.flux = 0.0
        pass

    def set_lambda_i_zl(self, flux):
        self.flux = flux * 1e-9 * u.W / (u.m*u.m) / u.sr

    def get_photons(self, dtel, wl, dwl, expt, pixscale):
        stel = (dtel/2.0)**2 * math.pi  # telescope area
        spix = (pixscale**2).to(u.sr)  # pixel scale (arcsec)
        phe = const.h * const.c / wl  # photon energy
        return self.flux * stel * expt / wl * dwl * spix / phe


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit(1)
    flux = float(sys.argv[1])

    dtel = 0.3 * u.m
    wl = 1.25 * 1e-6 * u.m
    dwl = 0.5 * 1e-6 * u.m
    expt = 120.0 * u.s
    uc = UnitConv()
    uc.set_lambda_i_zl(flux)
    pixscale = 1.0 * u.arcsec
    print(uc.flux)
    a = uc.get_photons(dtel, wl, dwl, expt, pixscale)
    print(a.to(u.dimensionless_unscaled))
    
    
