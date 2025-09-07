## Load libraries
import numpy as np
import pynbody as pnb




## Simulation and observation scaling relations

## log(M500*E(z)/Msun) - log(Tspec [keV]) scaling relation for the 512/50 Simba-C box
# tspecs = np.linspace(-0.6,0.48, 20)
# y=1.935*tspecs+13.32
# y1=1.935*tspecs+13.62
# y2=1.935*tspecs+13.02

h_sun09 = 0.73
h_simba = 0.68

M3_sun09 = pnb.array.SimArray(1.27e14 * h_simba**(-1), units='Msol')  # Msun
alpha_sun09 = pnb.array.SimArray(1.67, units='1')
r3_sun09 = pnb.array.SimArray(0.602 * h_simba**(-1), units='Mpc')  # Mpc

def M500_Tspec_simbac_N1024L100_renier(Tspec):
    ## tspec,corr -> log(Tspec [keV])
    ## Returns: M500/Msun -> log(M500/Msun)
    
    logTspec = np.log10(Tspec.in_units('keV'))
    return pnb.array.SimArray(10**(1.935*logTspec + 13.32), units='Msol')

def Tspec_M500_simbac_N1024L100_renier(M500):
    ## M500/Msun -> log(M500/Msun)
    ## Returns: tspec,corr -> log(Tspec [keV])
    
    logM500 = np.log10(M500.in_units('Msol'))
    return pnb.array.SimArray(10**((logM500 - 13.32)/1.935), units='keV')

def M500_Tspec_simba_N1024L100_aviv(Tspec):
    ## tspec,corr -> log(Tspec [keV])
    ## Returns: M500/Msun -> log(M500/Msun)
    
    logTspec = np.log10(Tspec.in_units('keV'))
    return pnb.array.SimArray(10**(2.109*logTspec + 13.61), units='Msol')

def Tspec_M500_simba_N1024L100_aviv(M500):
    ## M500/Msun -> log(M500/Msun)
    ## Returns: tspec,corr -> log(Tspec [keV])
    
    logM500 = np.log10(M500.in_units('Msol'))
    return pnb.array.SimArray(10**((logM500 - 13.61)/2.109), units='keV')

def M500_Tspec_simbac_N1024L100_aviv(Tspec):
    ## tspec,corr -> log(Tspec [keV])
    ## Returns: M500/Msun -> log(M500/Msun)
    
    logTspec = np.log10(Tspec.in_units('keV'))
    return pnb.array.SimArray(10**(2.134*logTspec + 13.63), units='Msol')

def Tspec_M500_simbac_N1024L100_aviv(M500):
    ## M500/Msun -> log(M500/Msun)
    ## Returns: tspec,corr -> log(Tspec [keV])
    
    logM500 = np.log10(M500.in_units('Msol'))
    return pnb.array.SimArray(10**((logM500 - 13.63)/2.134), units='keV')


def M500_T500_sun09(T500, h):
    ## T500 in keV
    ## M500 in Msun
    
#     M3_sun09 = pnb.array.SimArray(1.27e14 * h**(-1), units='Msol')  # Msun
# #     M3_sun09 = 1.27e14 * h**(-1)  # Msun
#     alpha_sun09 = 1.67
#     r3_sun09 = pnb.array.SimArray(0.602 * h**(-1), units='Mpc')  # Mpc
    
#     print(T500.units)
#     print(kB*T500)
#     print((kB*T500).units)
#     print((kB*T500).in_units('keV'))
    
#     return pnb.array.SimArray(M3_sun09 * pnb.array.SimArray((T500.in_units('keV') * pnb.array.SimArray(3., units='keV')**(-1))**alpha_sun09, units='1'), units='Msol')
    return pnb.array.SimArray(M3_sun09 * (T500.in_units('keV') * pnb.array.SimArray(3., units='keV')**(-1))**alpha_sun09, units='Msol')

def R500_T500_sun09(T500, h):
    ## T500 in keV
    ## R500 in Mpc
    
#     M3_sun09 = pnb.array.SimArray(1.27e14 * h**(-1), units='Msol')  # Msun
#     alpha_sun09 = 1.67
#     r3_sun09 = pnb.array.SimArray(0.602 * h**(-1), units='Mpc')  # Mpc
    
    return r3_sun09 * ((T500).in_units('keV')/pnb.array.SimArray(3., units='keV'))**(alpha_sun09/3.)

def T500_R500_sun09(R500, h):
    ## R500 in Mpc
    ## T500 in K
    
#     M3_sun09 = pnb.array.SimArray(1.27e14 * h**(-1), units='Msol')  # Msun
#     alpha_sun09 = 1.67
#     r3_sun09 = pnb.array.SimArray(0.602 * h**(-1), units='Mpc')  # Mpc
    
#     return (pnb.array.SimArray(3., units='keV') * (R500.in_units('Mpc')/r3_sun09)**(3./alpha_sun09) / kB).in_units('K')
    return pnb.array.SimArray(pnb.array.SimArray(3., units='keV') * (R500.in_units('Mpc')/r3_sun09)**(3./alpha_sun09), units='keV')