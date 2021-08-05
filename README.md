# Zodical Light Calculator

2021-08-25 H. Akitaya (PERC, CIT)

# Install from pip

```
$ git clone https://github.com/h-akitaya/zodicallight.git
$ cd zodicallight
$ pip install -e .
```

# Usage


* Calculate energy flux density at wavelength wl (um) and ecliptical latitude lat (deg).
```python:test.py
from zodicallight import ZodicalLight
zd = ZodicalLight()

wl = 1.0  # Wavelength (um)
lat = 20.0  # Ecliptical Latitude (deg)
print(zd(wl, lat))  # nW/m^2/str
617.2478965691489
```

* Calculate photon flux density at wavelength 'optical-band (o)' and ecliptical latitude lat (deg).
```python:test02.py
from zodicallight import ZodicalLight
zd = ZodicalLight()

band = 'o'  # optical band ('o'(0.5-0.9 um), 'j(0.9-1.5 um)', 
            # 'h'(1.5-2.0 um), 'k'(2.0-2.5 um))
lat = 0.0  # Ecliptical Latitude (deg)

print(zd.get_zl_photonflux_band(band, lat, unit=True))
1922040398028.5942 ph / (m2 s sr)

print(zd.get_zl_photonflux_band(band, lat, unit=False))
1922040398028.5942
```

* Calculate electron number for each pixel.
```python:test03.py
import numpy as np
import astropy.units as u

from zodicallight import ZodicalLight

zd = ZodicalLight()

band = 'o'  # band name
lat = 0.0  # Ecliptical Latitude (deg)

# Telescope and observation parameters
tel_d = 0.3*u.m  # Telescope diameter.
tel_s = np.pi*(tel_d/2)**2  # Telescope effective area.
tel_e = 0.6  # Telescope and instrument efficiency.
pix_fov = (2.0*u.arcsec)**2  # Pixel scale.
expt = 120*u.s  # Exposure time.

ph = zd.get_zl_photonflux_band(band, lat, unit=True) * tel_s * tel_e * pix_fov * expt

print(ph.decompose())
919.6808082096942 ph
```
