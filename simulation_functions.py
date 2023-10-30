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

from global_constants import *
from merger_remnant import *

def eccentricity_function(e):
    """
    Auxiliary function of eccentricity, from Peters & Mathews (1963).
    
    Inputs:
    @in e: eccentricity [0,1)
    """
    
    F = (1 - e**2)**(-7 / 2) * (1 + 73 / 24 * e**2 + 37 / 96 * e**4)
    
    return F

def relaxation_time(M_cl, r_h, m_avg, lnLambda, psi):
    """
    Half-mass relaxation timescale.
    
    Inputs:
    @in M_cl: cluster mass
    @in r_h: half-mass radius
    @in m_avg: average cluster mass
    @in lnLambda: Coulomb logarithm
    @in psi: multimass relaxation factor
    """
    
    t_rh = 0.138 * np.sqrt(M_cl * r_h**3 / G_Newt) / m_avg / lnLambda / psi
    
    return t_rh

def IMF_kroupa(m):
    """
    Kroupa (2002) initial mass function.
    Returns the number dN of stars in mass bin [m, m + dm] (1/Msun)
    
    Inputs:
    @in m : stellar mass array (Msun)
    """
    
    # mass boundaries (in solar masses):
    m1 = 0.08
    m2 = 0.50
    m3 = 1.00
    
    # spectral indices (broken power law; central values):
    a0 = -0.3
    a1 = -1.3
    a2 = -2.3
    a3 = alphaIMF
    
    # normalization constants:
    c1 = m1**a0 / m1**a1
    c2 = c1 * m2**a1 / m2**a2
    c3 = c2 * m3**a2 / m3**a3
    
    out = np.zeros(m.size)
    
    for i in range(0, m.size):
        
        if  (m[i] <= m1):
            
            out[i] = m[i]**a0
            
        elif(m[i] <= m2 and m[i] > m1):
            
            out[i] = c1 * m[i]**a1
            
        elif(m[i] <= m3 and m[i] > m2):
            
            out[i] = c2 * m[i]**a2
            
        elif(m[i] >= m3):
            
            out[i] = c3 * m[i]**a3
            
    return out

def stellar_number_fraction(m_a, m_b, m_min, m_max):
    """
    Fraction of stars with masses in the range [m_a, m_b],
    assuming a Kroupa IMF in [m_min, m_max].
    
    Inputs:
    @in m_a: lower limit
    @in m_b: upper limit
    @in m_min: mimimum star mass
    @in m_max: maximum star mass
    """
    
    f_ab = integrate.quad(lambda x: IMF_kroupa(np.array([x])), m_a, m_b)[0] / \
           integrate.quad(lambda x: IMF_kroupa(np.array([x])), m_min, m_max)[0]
    
    return f_ab

def mean_stellar_mass(m_a, m_b, m_min, m_max):
    """
    Mean stellar mass in range [m_a, m_b],
    assuming a Kroupa IMF in [m_min, m_max].
    
    Inputs:
    @in m_a: lower limit
    @in m_b: upper limit
    @in m_min: mimimum star mass
    @in m_max: maximum star mass
    """
    
    m_mean = integrate.quad(lambda x: IMF_kroupa(np.array([x])) * x, m_a, m_b)[0] / \
             integrate.quad(lambda x: IMF_kroupa(np.array([x])), m_min, m_max)[0]
    
    return m_mean

def T_coal(m1, m2, a0, e0):
    """
    GW coalescence timescale (I. Mandel 2021 fit to Peters timescale),
    including 1st order pN correction effects [Zwick et al., MNRAS 495, 2321 (2020)]
    
    Inputs:
    @in m1: primary   mass component
    @in m2: secondary mass component
    @in a0: initial semimajor axis
    @in e0: initial eccentricity in range [0,1)
    """
    
    # coalescence timescale for circular orbit:
    Tc = 5*c_light**5*a0**4/(256*G_Newt**3*m1*m2*(m1+m2))
    
    # 1st order pN correction:
    S = 8**(1-np.sqrt(1-e0)) * np.exp( 5*G_Newt*(m1+m2)/c_light**2/a0/(1-e0) )
    
    return Tc*(1+0.27*e0**10+0.33*e0**20+0.2*e0**1000)*(1-e0**2)**(7/2) * S

def forb(e, f0, e0):
    """
    Orbital frequency.
    
    Inputs:
    @in e: eccentricity
    @in f0: initial orbital frequency
    @in e0: initial eccentricity
    """
    
    term = (1-e0**2)/(1-e**2) * (e/e0)**(12/19) * ((1+121/304*e**2)/(1+121/304*e0**2))**(870/2299)
    forb = f0 * term**(-3/2)
    return forb

def ecce(f, f0, e0):
    """
    Eccentricity at a specific frequency given initial conditions.
    
    # Inputs:
    @in f: orbital frequency
    @in f0: initial orbital frequency
    @in e0: initial eccentricity
    """
    
    ecce = e0
    
    while f > forb(ecce, f0, e0):
        
        if ecce > 0.9:
            
            de = (1-ecce)/1000
            
        else:
            
            de = ecce/1000
            
        ecce = ecce - de
        
    return ecce

def g(n, e):
    """
    From Peters & Mathews (1963).
    
    Inputs:
    @in n: harmonic number
    @in e: eccentricity
    """
    
    if e==0:
        return 1 if n==2 else 0
    else:
        x = n*e
        g1 = jv(n-2, x) - 2*e*jv(n-1, x) + 2/n*jv(n, x) + 2*e*jv(n+1, x) - jv(n+2, x)
        g2 = jv(n-2, x) - 2*jv(n, x) + jv(n+2, x)
        g = n**4/32 * (g1**2 + (1-e**2)*g2**2 + 4/3/n**2*jv(n, x)**2)
        return g

def Phi(f, f0, e0, flso):
    """
    From Bonetti & Sesana (2020).
    
    Inputs:
    @in f: orbital frequency
    @in f0: initial orbital frequency
    @in e0: initial eccentricity
    @in flso: last stable orbit frequency
    """
    
    nmax = n_max(e0)
    s = 0
    for n in range(1, nmax):
        
        if f/n > f0 and f/n < flso:
            
            e = ecce(f/n, f0, e0)
            s += g(n, e) / n**(2/3) / eccentricity_function(e)
    Phi = 2**(2/3) * s
    
    return Phi

def n_max(e):
    """
    Maximum harmonic number (for satisfactory accuracy).
    
    Inputs:
    @in e: eccentricity
    """
    
    n_max = round(5*(1+e)**(1/2)/(1-e)**(3/2))
    
    return n_max

def hc(f, m1, m2, z, f0, e0):
    """
    Characteristic GW strain.
    
    Inputs:
    @in f: GW frequency
    @in m1: primary mass
    @in m2: secondary mass
    @in z: merger redshift
    @in f0: initial orbital frequency
    @in e0: initial eccentricity
    """
    
    Mt = m1+m2
    Mc = (m1*m2)**(3/5)/Mt**(1/5)
    Mcz = Mc*(1+z)
    dL = Planck18.luminosity_distance(z).value*1e6
    flso = c_light**3 / 2/np.pi /G_Newt / Mt / 6**(3/2)
    const = 2*G_Newt**(5/3)*np.pi**(2/3)*Mcz**(5/3)*f**(-1/3)/3/c_light**3/np.pi**2/dL**2
    hc2 = const * Phi(f*(1+z), f0, e0, flso)
    hc = np.sqrt(hc2)
    return hc

def dYdt(Y, t, extras):
    """
    Inputs:
    @in Y: current solution vector
    @in t: current time
    @in extras: extra variables
    """
    
    # get extras:
    Y0, dt, seed, discrete_evolution, tge, t_max = extras
    
    # random state:
    np.random.seed(seed)
    
    # unpack variables:
    M_H, S_H, N_star, N_BH, m_star, m_BH, M_gas, r_h, m_r = Y
    
    # unpack initial conditions:
    M_H0, S_H0, N_star0, N_BH0, m_star0, m_BH0, M_gas0, r_h0, m_r0 = Y0
    
    # star mass:
    M_star = N_star * m_star
    
    # stellar BH mass:
    M_BH = N_BH * m_BH
    
    if M_gas < 0: M_gas = 0
    
    # cluster mass:
    M_cl = M_star + M_BH + M_H + M_gas
    
    # cluster energy:
    E_cl = - kappa / 2 * G_Newt * M_cl**2 / r_h
    
    # 3D stellar velocity dispersion:
    v_star = np.sqrt(kappa * G_Newt * M_cl / r_h)
    
    # escape velocity:
    v_esc = 2 * v_star
    
    # Coulomb logarithms:
    lnLambda, lnLambda_BH = 10.0, 1.0
    
    # Spitzer parameter:
    S_BH_star = (N_BH / N_star) * (m_BH / m_star)**(5 / 2)
    
    # multimass factor:
    psi = 1 + S_BH_star
    
    # relaxation time:
    t_rh = relaxation_time(M_cl, r_h, m_star, lnLambda, psi)
    
    # core collapse time:
    t_cc = 0.2 * relaxation_time(N_star0 * m_star0 + M_gas0, r_h0, m_star0, lnLambda, 1)
    
    # heat released by a stellar collision:
    Q_col = G_Newt * m_r * m_star / 2 / r_sun
    
    # ejecter radius:
    a_ejs = beta * G_Newt * M_H / (2 + M_H / m_BH) / v_esc**2
    
    # self-ejecter radius:
    a_ejb = a_ejs / (1 + M_H / m_BH)**2
    
    C = - 2 * zeta * E_cl / G_Newt / m_BH / M_H / t_rh
    D = 64 / 5 * G_Newt**3 * M_H * m_BH * (M_H + m_BH) / c_light**5
    
    # thermal eccentricity:
    e = np.sqrt(np.random.rand())
    
    # GW radius:
    a_GW = (D / C * eccentricity_function(e))**(1 / 5)
    
    # equipartition parameter:
    xi_BH_star = (m_BH / m_star)**(3 / 5) * (M_BH / M_star)**(2 / 5) * (lnLambda_BH/lnLambda)**(-2 / 5) if S_BH_star > S_BH_star_critical else 1
    
    # initial binary hardness:
    eta0 = (1-np.random.rand())**(-2/7) #get_hardness_sample(seed)
    
    # initial radius:
    a0 = 4 * G_Newt * M_H * m_BH / m_star / xi_BH_star / eta0 / v_esc**2 if N_BH>0 else 0
    
    # number of singles ejected:
    N_BH_ejs = np.log(a_ejs / np.max([a_ejb, a_GW])) / np.log(1 + beta / (1 + M_H / m_BH)) \
    if a_ejs > np.max([a_ejs, a_GW]) else 0
    
    # heat released by a HBH:
    Q_HBH = G_Newt * M_H * m_BH / 2 / a_GW if a_GW > a_ejs \
    else G_Newt * M_H * m_BH / 2 / (M_H + m_BH) * (m_BH / np.max([a_GW, a_ejb]) + M_H / a_ejs)
    
    # tidal radius:
    r_tidal = r_sun * (M_H / m_star)**(1 / 3)
    
    # heat released by a TDE:
    Q_TDE = G_Newt * M_H * m_star / 2 / r_tidal
    
    # power demand:
    dEdt = - zeta * E_cl / t_rh
    
    # collision rate:
    R_col = dEdt / Q_col if np.max([t_cc, t]) < t_sev else 0.0
    
    # HBH rate:
    R_HBH = dEdt / Q_HBH if N_BH  > 0 and t > np.max([t_cc, t_sev]) else 0.0
    
    # TDE rate:
    R_TDE = dEdt / Q_TDE if N_BH == 0 and t > np.max([t_cc, t_sev]) else 0.0
    
    # gas density:
    rho_gas = 3 * M_gas / (4 * np.pi * (r_h / 1.3)**3)
    
    # ----------------------------------------------------------------------------------------------
    # M_H continuous evolution:
    
    dM_Hdt__acc = np.min([4 * np.pi * lamb * (G_Newt * M_H)**2 * rho_gas / c_s**3, M_H / t_Sal])
    dM_Hdt__TDE = f_TDE * m_star * R_TDE
    
    dM_Hdt = dM_Hdt__acc + dM_Hdt__TDE
    
    # ----------------------------------------------------------------------------------------------
    # S_H continuous evolution:
    
    dS_Hdt = 0.0
    
    # ----------------------------------------------------------------------------------------------
    # N_star continuous evolution:
    
    dN_stardt__rel = - xi_e * N_star / t_rh
    dN_stardt__TDE = - R_TDE
    dN_stardt__col = - R_col
    
    dN_stardt = dN_stardt__rel + dN_stardt__TDE + dN_stardt__col
    
    # ----------------------------------------------------------------------------------------------
    # N_BH continuous evolution:
    
    dN_BHdt = - N_BH_ejs * R_HBH
    
    if N_BH < 0: N_BH = 0
    
    # ----------------------------------------------------------------------------------------------
    # m_star continuous evolution:
    
    dm_stardt = - nu * m_star / t if t > t_sev else 0.0
    
    # ----------------------------------------------------------------------------------------------
    # m_BH continuous evolution:
    
    dm_BHdt = 0.0
    
    # ----------------------------------------------------------------------------------------------
    # M_gas continuous evolution:
    
    dM_gasdt__fee = - M_gas / tge
    dM_gasdt__acc = - dM_Hdt__acc
    
    dM_gasdt = dM_gasdt__fee + dM_gasdt__acc
    
    # ----------------------------------------------------------------------------------------------
    # r_h continuous evolution:
    
    dM_cldt = m_star*dN_stardt + N_star*dm_stardt + m_BH*dN_BHdt + dM_gasdt + dM_Hdt
    
    dM_cldt__sev = N_star * dm_stardt
    dM_cldt__fee = dM_gasdt__fee
    
    dr_hdt__rel = zeta * r_h / t_rh + 2 * r_h * dM_cldt / M_cl
    dr_hdt__adi = - (dM_cldt__sev + dM_cldt__fee) * r_h / M_cl
    
    dr_hdt = dr_hdt__rel + dr_hdt__adi if t > t_cc else dr_hdt__adi
    
    # ----------------------------------------------------------------------------------------------
    # m_r continuous evolution:
    
    dm_rdt = m_star * R_col
    
    if t > t_sev:
        m_r = 0.0
        
    # ----------------------------------------------------------------------------------------------
    # M_H & S_H discrete evolution:
    
    v_k=0
    k_HBH=0
    k_in_me_re=0
    k_in_me_ej=0
    k_out_me=0
    k_cap=0
    
    if discrete_evolution:
        
        if N_BH > 0 and R_HBH > 0:
            
            # number of HBH binaries formed:
            k_HBH = poisson.rvs(mu=dt * R_HBH)
            
            for k in range(k_HBH):
                
                if a_GW < a_ejb: # binary ejected
                    
                    M_H = m_BH
                    S_H = 0.0
                    
                    N_BH = N_BH - 2
                    
                    if T_coal(M_H, m_BH, a_ejb, e) < t_max - t:
                        # binary merges in field within the available time
                        k_out_me += 1
                        
                else: # binary merges in cluster
                    
                    v_k = remnant_kick(m_BH/ M_H, S_H, 0.0) # remnant kick
                    
                    N_BH = N_BH - 1
                    
                    if a_GW > a0: # GW capture
                        
                        k_cap += 1
                    
                    if v_k > v_esc: # merger remnant ejected
                        
                        M_H = m_BH
                        S_H = 0.0
                        
                        N_BH = N_BH - 1
                        
                        k_in_me_ej += 1
                        
                    else: # merger remnant retained
                        
                        M_H = remnant_mass(M_H, m_BH, S_H, 0.0) # remnant mass
                        S_H = remnant_spin(m_BH/ M_H, S_H, 0.0) # remnant spin
                        
                        k_in_me_re += 1
                        
    # ----------------------------------------------------------------------------------------------
    
    dYdt = np.array([dM_Hdt, dS_Hdt, dN_stardt, dN_BHdt, dm_stardt, dm_BHdt, dM_gasdt, dr_hdt, dm_rdt])
    
    other = {'t_rh':t_rh, 'psi':psi, 'S_BH_star':S_BH_star, 'R_HBH':R_HBH, 'R_TDE':R_TDE, 'R_col':R_col, 
             'M_cl':M_cl, 'v_esc':v_esc, 'a_ejs':a_ejs, 'a_ejb':a_ejb, 'a_GW':a_GW, 'a0':a0, 'xi_BH_star':xi_BH_star, 
             'v_k':v_k, 'k_HBH':k_HBH, 'k_in_me_re':k_in_me_re, 'k_in_me_ej':k_in_me_ej, 'k_out_me':k_out_me, 'k_cap':k_cap}
    
    Y = [M_H, S_H, N_star, N_BH, m_star, m_BH, M_gas, r_h, m_r]
    
    return dYdt, other, Y

def solve_cluster(Y0, t_max, dt, every, tge, seed):
    """
    Inputs:
    @in Y0: initial condition
    @in t_max: maximum simulation time (Myr)
    @in dt: time step (Myr)
    @in every: save data every _every_
    @in tge: gas expulsion timescale (Myr)
    @in seed: random seed number
    """
    
    np.random.seed(seed)
    
    # initialize outputs:
    time = []
    M_H = []
    S_H = []
    N_star = []
    N_BH = []
    m_star = []
    m_BH = []
    M_gas = []
    r_h = []
    m_r = []
    
    t_rh = []
    psi = []
    S_BH_star = []
    R_HBH = []
    R_TDE = []
    R_col = []
    M_cl = []
    v_esc = []
    a_ejs = []
    a_ejb = []
    a_GW = []
    a0 = []
    xi_BH_star = []
    v_k = []
    k_HBH = []
    k_in_me_re = []
    k_in_me_ej = []
    k_out_me = []
    k_cap = []
    
    timestep = []
    
    # auxiliary integer:
    i_BH = 0
    
    # set initial condition:
    Y = Y0
    
    # set initial time:
    t = 0.0
    
    N_steps = int(t_max / dt)
    
    # fourth order Runge-Kutta:
    #while t < t_max:
    for i in tqdm(range(N_steps)):
        
        seed = np.random.randint(99999999)
        
        extras = [Y0, dt, seed, False, tge, t_max]
        
        # save data:
        if i % every == 0:
            
            time.append(t)
            M_H.append(Y[0])
            S_H.append(Y[1])
            N_star.append(Y[2])
            N_BH.append(Y[3])
            m_star.append(Y[4])
            m_BH.append(Y[5])
            M_gas.append(Y[6])
            r_h.append(Y[7])
            m_r.append(Y[8])
            
            t_rh.append(dYdt(Y, t, extras)[1]['t_rh'])
            psi.append(dYdt(Y, t, extras)[1]['psi'])
            S_BH_star.append(dYdt(Y, t, extras)[1]['S_BH_star'])
            R_HBH.append(dYdt(Y, t, extras)[1]['R_HBH'])
            R_TDE.append(dYdt(Y, t, extras)[1]['R_TDE'])
            R_col.append(dYdt(Y, t, extras)[1]['R_col'])
            M_cl.append(dYdt(Y, t, extras)[1]['M_cl'])
            v_esc.append(dYdt(Y, t, extras)[1]['v_esc'])
            a_ejs.append(dYdt(Y, t, extras)[1]['a_ejs'])
            a_ejb.append(dYdt(Y, t, extras)[1]['a_ejb'])
            a_GW.append(dYdt(Y, t, extras)[1]['a_GW'])
            a0.append(dYdt(Y, t, extras)[1]['a0'])
            xi_BH_star.append(dYdt(Y, t, extras)[1]['xi_BH_star'])
            
            timestep.append(dt)
            
        # RK4 slopes:
        k1 = dYdt(Y, t, extras)[0]
        k2 = dYdt(Y + dt / 2 * k1, t + dt / 2, extras)[0]
        k3 = dYdt(Y + dt / 2 * k2, t + dt / 2, extras)[0]
        k4 = dYdt(Y + dt * k3, t + dt, extras)[0]
        
        # continuous evolution:
        Y = Y + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6
        
        # discrete evolution:
        discrete_data = dYdt(Y, t, [Y0, dt, seed, True, tge, t_max])
        Y = discrete_data[2]
        
        if i % every == 0:
            
            v_k.append(discrete_data[1]['v_k'])
            k_HBH.append(discrete_data[1]['k_HBH'])
            k_in_me_re.append(discrete_data[1]['k_in_me_re'])
            k_in_me_ej.append(discrete_data[1]['k_in_me_ej'])
            k_out_me.append(discrete_data[1]['k_out_me'])
            k_cap.append(discrete_data[1]['k_cap'])
            
        # update time:
        t = t + dt
        
        # update auxiliary integer:
        i += 1
        
    time = np.array(time)
    M_H = np.array(M_H)
    S_H = np.array(S_H)
    N_star = np.array(N_star)
    N_BH = np.array(N_BH)
    m_star = np.array(m_star)
    m_BH = np.array(m_BH)
    M_gas = np.array(M_gas)
    r_h = np.array(r_h)
    m_r = np.array(m_r)
    
    t_rh = np.array(t_rh)
    psi = np.array(psi)
    S_BH_star = np.array(S_BH_star)
    R_HBH = np.array(R_HBH)
    R_TDE = np.array(R_TDE)
    R_col = np.array(R_col)
    M_cl = np.array(M_cl)
    v_esc = np.array(v_esc)
    a_ejs = np.array(a_ejs)
    a_ejb = np.array(a_ejb)
    a_GW = np.array(a_GW)
    a0 = np.array(a0)
    xi_BH_star = np.array(xi_BH_star)
    v_k = np.array(v_k)
    k_HBH = np.array(k_HBH)
    k_in_me_re = np.array(k_in_me_re)
    k_in_me_ej = np.array(k_in_me_ej)
    k_out_me = np.array(k_out_me)
    k_cap = np.array(k_cap)
    
    timestep = np.array(timestep)
    
    data = {'t': time, 'M_cl':M_cl, 'M_H': M_H, 'S_H':S_H, 'N_star':N_star, 'N_BH':N_BH, 'm_star':m_star, 'm_BH':m_BH, 'M_gas':M_gas, 'r_h':r_h, 'm_r':m_r, 't_rh':t_rh, 'psi':psi, 'S_BH_star':S_BH_star, 'R_HBH':R_HBH, 'R_TDE':R_TDE, 'R_col':R_col, 'v_esc':v_esc, 'a_ejs':a_ejs, 'a_ejb':a_ejb, 'a_GW':a_GW, 'a0':a0, 'xi_BH_star':xi_BH_star, 'v_k':v_k, 'k_HBH':k_HBH, 'k_in_me_re':k_in_me_re, 'k_in_me_ej':k_in_me_ej, 'k_out_me':k_out_me, 'k_cap':k_cap, 'dt':timestep}
    
    return data

def simulate_cluster(N0, rh0, Mgas0=0.0, mBH0=10.0, MH0=10.0, tmax=14000.0, dt=0.01, tge=3.0, seed=348529):
    
    M_H0 = MH0 # seed BH mass
    S_H0 = 0.0 # seed BH spin
    N_star0 = N0 # initial star number
    m_r0 = 150.0 # seed runaway star mass
    N_BH0 = stellar_number_fraction(20, m_r0, 0.08, m_r0) * N_star0 # initial BH number
    m_star0 = mean_stellar_mass(0.08, m_r0, 0.08, m_r0) # initial mean star mass
    m_BH0 = mBH0 # initial BH mass
    M_gas0 = Mgas0 # initial gas mass
    r_h0 = rh0 # initial half-mass radius
    
    # initial condition:
    Y0 = np.array([M_H0, S_H0, N_star0, N_BH0, m_star0, m_BH0, M_gas0, r_h0, m_r0])
    
    # maximum time:
    t_max = tmax
    
    # every:
    every = 1
    
    cluster = solve_cluster(Y0, t_max, dt, every, tge, seed)
    
    return cluster

# end of file