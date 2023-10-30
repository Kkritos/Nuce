"""
 Copyright (C) 2023  Konstantinos Kritos <kkritos1@jhu.edu>

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import numpy as np
from scipy.special import jv
from scipy.stats import poisson
import scipy.integrate as integrate
from astropy.cosmology import Planck18
from tqdm import tqdm

# Newton constant:
G_Newt = 1/232

# speed of light:
c_light = 3e5

# Salpeter time:
t_Sal = 50.0

# solar radius
r_sun = 2.3e-8

# high-mass IMF index:
alphaIMF = -2.3

# stellar evolution parameter:
nu = 0.07

# relaxation efficiency:
zeta = 0.10

# stellar evolution timescale:
t_sev = 2.0

# potential shape parameter:
kappa = 0.4

# hardening rate parameter:
beta = 4/7

# sound speed:
c_s = 10.0

# TDE efficiency:
f_TDE = 0.50

# Bondi parameter:
lamb = 1.0

# star ejection parameter:
xi_e = 0.0074

# critical Spitzer parameter:
S_BH_star_critical = 0.16

# end of file