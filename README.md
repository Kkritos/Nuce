# Nuce

The repository contains a suite of functions required for NUclear star Cluster Evolution (NUCE).

### Requirements

The following Python packages are required,

- $\tt numpy$
- $\tt scipy$
- $\tt astropy$
- $\tt tqdm$

### Running
Import the ``simulate_cluster`` function from the other_functions.py file into your Python code and run

  > cluster = simulate_cluster(N0=N0, rh0=rh0, Mgas0=Mgas0, mBH0=mBH0, MH0=MH0, tmax=tmax, dt=dt, tge=tge, seed=seed)

Below we give a description of the input parameters:
| Variable | Description |
|:--- |:--- |
| N0 | initial number of stars |
| rh0 | initial half-mass radius (pc) |
| Mgas0 | initial gas amount ($M_\odot$) |
| mBH0 | black hole mass ($M_\odot$) |
| MH0 | initial runaway black hole mass ($M_\odot$) |
| tmax | maximum simulation time (Myr) |
| dt | simulation time step (Myr) |
| tge | gas expulsion timescale (Myr) |
| seed | seed number |

The output ``cluster`` above is a dictionary with the following keys:
| Key | Description |
|:--- |:--- |
| 't' | time (Myr) |
| 'M_cl' | total cluster mass ($M_\odot$) |
| 'M_H' | runaway black hole mass ($M_\odot$) |
| 'S_H' | runaway black hole dimensionless spin parameter |
| 'N_star' | number of stars |
| 'N_BH' | number of stellar black holes |
| 'm_star' | stellar mass ($M_\odot$) |
| 'm_BH' | stellar black hole mass ($M_\odot$) |
| 'M_gas' | residual gas mass ($M_\odot$) |
| 'r_h' | half-mass radius (pc) |
| 'm_r' | runaway star mass ($M_\odot$) |
| 't_rh' | half-mass relaxation time (Myr) |
| 'psi' | multimas relaxation factor |
| 'S_BH_star' | Spitzer parameter |
| 'R_HBH' | three-body binary formation rate (Myr$^{-1}$) |
| 'R_TDE' | tidal disruption rate (Myr$^{-1}$) |
| 'R_col' | stellar collision rate (Myr$^{-1}$) |
| 'v_esc' | escape velocity (km/s) |
| 'a_ejs' | ejecter radius (pc) |
| 'a_ejb' | self-ejecter radius (pc) |
| 'a_GW' | gravitational-wave radius (pc) |
| 'a0' | initial binary semimajor axis (pc) |
| 'xi_BH_star' | temperature ratio between black holes and stars |
| 'v_k' | gravitational-wave kick velocity (km/s) |
| 'k_HBH' | number of binaries formed in current time step |
| 'k_in_me_re' | number of in-cluster mergers retained |
| 'k_in_me_ej' | number of in-cluster mergers ejected |
| 'k_out_me' | number of ejected mergers (in current time step) |
| 'k_cap' | number of captures (in current time step) |
| 'dt' | time step (Myr) |

See also the example.ipynb notebook in this repository.

### Contact

email: kkritos1@jhu.edu
