## Load libraries
import pynbody as pnb
from pynbody.array import SimArray
import numpy as np
import periodictable as pt
import astropy.constants as c
import astropy.units as u



## Constants

ABOUT_ZERO = 1e-9
SMALL_OVERLAP = 1e-4

# Avogadro's constant
NA_no_units = c.N_A.value  # particles/mol
NA = pnb.array.SimArray(NA_no_units, units='')

# Boltzmann's constant
kB_no_units = c.k_B.to('keV/K').value
# kB_with_units = pnb.units.Unit('keV K**-1') * kB
kB = pnb.array.SimArray(kB_no_units, units='keV K**-1')
# print(kB.units)
# print((kB*pnb.array.SimArray(2, units='K')).units)

# Hydrogen mass in Msun (want proton mass though)
mH = pnb.array.SimArray(pt.H.mass/NA, units='g').in_units('Msol')

# Proton mass in Msun
mp_no_units = c.m_p.to('Msun').value
mp = pnb.array.SimArray(mp_no_units, units='Msol')

# mean molecular weight per free electron (from Oppenheimer+21 <-- McCarthy+08)
mu_e = pnb.array.SimArray(1.14, units='1')

# mean molecular weight (from Oppenheimer+21 --> McCarthy+08)
mu = pnb.array.SimArray(0.59, units='1')
# Above two values are appropriate for gas with metallicity Z=0.3*Zsol

#G = u.G  # Newton's Gravitational constant (set value?)
# G_no_units = c.G.to('m**3 kg**-1 s**-2').value
G_no_units = c.G.to('keV Msun**-2 kpc').value
# print(G_no_units)
# G = pnb.array.SimArray(G_no_units, units='m**3 kg**-1 s**-2').in_units('keV Msol**-2 kpc')
G = pnb.array.SimArray(G_no_units, units='keV Msol**-2 kpc')
# print(G)
# print(G.units)



## Should be using R200-R500 and R2500-R500 scaling relations from Simba or Simba-C for these conversions??

# R500_OVER_R200 = 0.65  # Approximately for NFW profiles, from Reiprich+13
# R200_OVER_R500 = 1/R500_OVER_R200  # Approximately for NFW profiles, from Reiprich+13
# R2500_OVER_R200 = 0.3  # Porciani/Basu (Uni Bonn)
# R200_OVER_R2500 = 1/R2500_OVER_R200
# R2500_OVER_R500 = R2500_OVER_R200 * R200_OVER_R500
# R500_OVER_R2500 = 1/R2500_OVER_R500
# # R500_OVER_R200 = 0.7  # Porciani/Basu (Uni Bonn)
# # R2500_OVER_R200 = 0.4  # Liu+2023


## Median values of Rxxx and Rxxx ratios from simba and simba-c for all groups (N_lgal>=3)) M500>=1e13 Msun, and in 3 mass bins
## From scaling_relations_with_xigrm.ipynb
median_R_values = {'simba': {'ratios': {'all': {'R200/R2500': SimArray(3.46687471, '1.00e+00'),
                              'R200/R500': SimArray(1.52865472, '1.00e+00'),
                              'R2500/R200': SimArray(0.28844423, '1.00e+00'),
                              'R2500/R500': SimArray(0.44297859, '1.00e+00'),
                              'R500/R200': SimArray(0.65416997, '1.00e+00'),
                              'R500/R2500': SimArray(2.25744545, '1.00e+00')},
                      'bin_0': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                'R200/R500': SimArray(np.nan, '1.00e+00'),
                                'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                'R500/R200': SimArray(np.nan, '1.00e+00'),
                                'R500/R2500': SimArray(np.nan, '1.00e+00')},
                      'bin_1': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                'R200/R500': SimArray(np.nan, '1.00e+00'),
                                'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                'R500/R200': SimArray(np.nan, '1.00e+00'),
                                'R500/R2500': SimArray(np.nan, '1.00e+00')},
                      'bin_2': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                'R200/R500': SimArray(np.nan, '1.00e+00'),
                                'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                'R500/R200': SimArray(np.nan, '1.00e+00'),
                                'R500/R2500': SimArray(np.nan, '1.00e+00')},
                      'bin_3': {'R200/R2500': SimArray(3.43868275, '1.00e+00'),
                                'R200/R500': SimArray(1.52638365, '1.00e+00'),
                                'R2500/R200': SimArray(0.29080903, '1.00e+00'),
                                'R2500/R500': SimArray(0.44607821, '1.00e+00'),
                                'R500/R200': SimArray(0.65514328, '1.00e+00'),
                                'R500/R2500': SimArray(2.24175933, '1.00e+00')},
                      'bin_4': {'R200/R2500': SimArray(3.54047536, '1.00e+00'),
                                'R200/R500': SimArray(1.53815798, '1.00e+00'),
                                'R2500/R200': SimArray(0.28244804, '1.00e+00'),
                                'R2500/R500': SimArray(0.43196709, '1.00e+00'),
                                'R500/R200': SimArray(0.65013529, '1.00e+00'),
                                'R500/R2500': SimArray(2.31508285, '1.00e+00')},
                      'bin_5': {'R200/R2500': SimArray(3.47913028, '1.00e+00'),
                                'R200/R500': SimArray(1.53201475, '1.00e+00'),
                                'R2500/R200': SimArray(0.28742816, '1.00e+00'),
                                'R2500/R500': SimArray(0.44100594, '1.00e+00'),
                                'R500/R200': SimArray(0.65273523, '1.00e+00'),
                                'R500/R2500': SimArray(2.26754313, '1.00e+00')},
                      'bin_6': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                'R200/R500': SimArray(np.nan, '1.00e+00'),
                                'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                'R500/R200': SimArray(np.nan, '1.00e+00'),
                                'R500/R2500': SimArray(np.nan, '1.00e+00')}},
           'value_ratios': {'all': {'R200/R2500': SimArray(3.50762883, '1.00e+00'),
                                    'R200/R500': SimArray(1.56146942, '1.00e+00'),
                                    'R2500/R200': SimArray(0.28509288, '1.00e+00'),
                                    'R2500/R500': SimArray(0.44516381, '1.00e+00'),
                                    'R500/R200': SimArray(0.6404224, '1.00e+00'),
                                    'R500/R2500': SimArray(2.24636409, '1.00e+00')},
                            'bin_0': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                      'R200/R500': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                      'R500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R500/R2500': SimArray(np.nan, '1.00e+00')},
                            'bin_1': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                      'R200/R500': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                      'R500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R500/R2500': SimArray(np.nan, '1.00e+00')},
                            'bin_2': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                      'R200/R500': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                      'R500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R500/R2500': SimArray(np.nan, '1.00e+00')},
                            'bin_3': {'R200/R2500': SimArray(3.45619257, '1.00e+00'),
                                      'R200/R500': SimArray(1.54616799, '1.00e+00'),
                                      'R2500/R200': SimArray(0.28933573, '1.00e+00'),
                                      'R2500/R500': SimArray(0.44736165, '1.00e+00'),
                                      'R500/R200': SimArray(0.64676025, '1.00e+00'),
                                      'R500/R2500': SimArray(2.23532797, '1.00e+00')},
                            'bin_4': {'R200/R2500': SimArray(3.57719874, '1.00e+00'),
                                      'R200/R500': SimArray(1.55362052, '1.00e+00'),
                                      'R2500/R200': SimArray(0.27954835, '1.00e+00'),
                                      'R2500/R500': SimArray(0.43431205, '1.00e+00'),
                                      'R500/R200': SimArray(0.64365782, '1.00e+00'),
                                      'R500/R2500': SimArray(2.30249195, '1.00e+00')},
                            'bin_5': {'R200/R2500': SimArray(3.52682605, '1.00e+00'),
                                      'R200/R500': SimArray(1.55535125, '1.00e+00'),
                                      'R2500/R200': SimArray(0.28354106, '1.00e+00'),
                                      'R2500/R500': SimArray(0.44100594, '1.00e+00'),
                                      'R500/R200': SimArray(0.64294159, '1.00e+00'),
                                      'R500/R2500': SimArray(2.26754313, '1.00e+00')},
                            'bin_6': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                      'R200/R500': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                      'R500/R200': SimArray(np.nan, '1.00e+00'),
                                      'R500/R2500': SimArray(np.nan, '1.00e+00')}},
           'values': {'all': {'R200': SimArray(645.14886419, 'kpc'),
                              'R2500': SimArray(183.92734668, 'kpc'),
                              'R500': SimArray(413.16778666, 'kpc')},
                      'bin_0': {'R200': SimArray(np.nan, 'kpc'),
                                'R2500': SimArray(np.nan, 'kpc'),
                                'R500': SimArray(np.nan, 'kpc')},
                      'bin_1': {'R200': SimArray(np.nan, 'kpc'),
                                'R2500': SimArray(np.nan, 'kpc'),
                                'R500': SimArray(np.nan, 'kpc')},
                      'bin_2': {'R200': SimArray(np.nan, 'kpc'),
                                'R2500': SimArray(np.nan, 'kpc'),
                                'R500': SimArray(np.nan, 'kpc')},
                      'bin_3': {'R200': SimArray(597.48015273, 'kpc'),
                                'R2500': SimArray(172.87235634, 'kpc'),
                                'R500': SimArray(386.42641419, 'kpc')},
                      'bin_4': {'R200': SimArray(914.88333837, 'kpc'),
                                'R2500': SimArray(255.75412665, 'kpc'),
                                'R500': SimArray(588.87181799, 'kpc')},
                      'bin_5': {'R200': SimArray(1510.74257879, 'kpc'),
                                'R2500': SimArray(428.35755402, 'kpc'),
                                'R500': SimArray(971.31923031, 'kpc')},
                      'bin_6': {'R200': SimArray(np.nan, 'kpc'),
                                'R2500': SimArray(np.nan, 'kpc'),
                                'R500': SimArray(np.nan, 'kpc')}}},
 'simba-c': {'ratios': {'all': {'R200/R2500': SimArray(3.44463263, '1.00e+00'),
                                'R200/R500': SimArray(1.52932679, '1.00e+00'),
                                'R2500/R200': SimArray(0.29030672, '1.00e+00'),
                                'R2500/R500': SimArray(0.4459453, '1.00e+00'),
                                'R500/R200': SimArray(0.6538825, '1.00e+00'),
                                'R500/R2500': SimArray(2.24242776, '1.00e+00')},
                        'bin_0': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                  'R200/R500': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                  'R500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R500/R2500': SimArray(np.nan, '1.00e+00')},
                        'bin_1': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                  'R200/R500': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                  'R500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R500/R2500': SimArray(np.nan, '1.00e+00')},
                        'bin_2': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                  'R200/R500': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                  'R500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R500/R2500': SimArray(np.nan, '1.00e+00')},
                        'bin_3': {'R200/R2500': SimArray(3.41042255, '1.00e+00'),
                                  'R200/R500': SimArray(1.52810375, '1.00e+00'),
                                  'R2500/R200': SimArray(0.2932188, '1.00e+00'),
                                  'R2500/R500': SimArray(0.44946462, '1.00e+00'),
                                  'R500/R200': SimArray(0.65440583, '1.00e+00'),
                                  'R500/R2500': SimArray(2.2248692, '1.00e+00')},
                        'bin_4': {'R200/R2500': SimArray(3.5487451, '1.00e+00'),
                                  'R200/R500': SimArray(1.53773192, '1.00e+00'),
                                  'R2500/R200': SimArray(0.28178975, '1.00e+00'),
                                  'R2500/R500': SimArray(0.43484928, '1.00e+00'),
                                  'R500/R200': SimArray(0.65030841, '1.00e+00'),
                                  'R500/R2500': SimArray(2.29964735, '1.00e+00')},
                        'bin_5': {'R200/R2500': SimArray(3.50339804, '1.00e+00'),
                                  'R200/R500': SimArray(1.52508255, '1.00e+00'),
                                  'R2500/R200': SimArray(0.28543716, '1.00e+00'),
                                  'R2500/R500': SimArray(0.43581811, '1.00e+00'),
                                  'R500/R200': SimArray(0.65570221, '1.00e+00'),
                                  'R500/R2500': SimArray(2.29453519, '1.00e+00')},
                        'bin_6': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                  'R200/R500': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                  'R500/R200': SimArray(np.nan, '1.00e+00'),
                                  'R500/R2500': SimArray(np.nan, '1.00e+00')}},
             'value_ratios': {'all': {'R200/R2500': SimArray(3.51515619, '1.00e+00'),
                                      'R200/R500': SimArray(1.54830985, '1.00e+00'),
                                      'R2500/R200': SimArray(0.28448238, '1.00e+00'),
                                      'R2500/R500': SimArray(0.44046687, '1.00e+00'),
                                      'R500/R200': SimArray(0.64586556, '1.00e+00'),
                                      'R500/R2500': SimArray(2.27031831, '1.00e+00')},
                              'bin_0': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                        'R200/R500': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                        'R500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R500/R2500': SimArray(np.nan, '1.00e+00')},
                              'bin_1': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                        'R200/R500': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                        'R500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R500/R2500': SimArray(np.nan, '1.00e+00')},
                              'bin_2': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                        'R200/R500': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                        'R500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R500/R2500': SimArray(np.nan, '1.00e+00')},
                              'bin_3': {'R200/R2500': SimArray(3.41226968, '1.00e+00'),
                                        'R200/R500': SimArray(1.5304933, '1.00e+00'),
                                        'R2500/R200': SimArray(0.29306007, '1.00e+00'),
                                        'R2500/R500': SimArray(0.44852648, '1.00e+00'),
                                        'R500/R200': SimArray(0.65338411, '1.00e+00'),
                                        'R500/R2500': SimArray(2.22952278, '1.00e+00')},
                              'bin_4': {'R200/R2500': SimArray(3.49385092, '1.00e+00'),
                                        'R200/R500': SimArray(1.54419255, '1.00e+00'),
                                        'R2500/R200': SimArray(0.28621713, '1.00e+00'),
                                        'R2500/R500': SimArray(0.44197437, '1.00e+00'),
                                        'R500/R200': SimArray(0.64758764, '1.00e+00'),
                                        'R500/R2500': SimArray(2.26257466, '1.00e+00')},
                              'bin_5': {'R200/R2500': SimArray(3.45486384, '1.00e+00'),
                                        'R200/R500': SimArray(1.54999876, '1.00e+00'),
                                        'R2500/R200': SimArray(0.28944701, '1.00e+00'),
                                        'R2500/R500': SimArray(0.4486425, '1.00e+00'),
                                        'R500/R200': SimArray(0.64516181, '1.00e+00'),
                                        'R500/R2500': SimArray(2.2289462, '1.00e+00')},
                              'bin_6': {'R200/R2500': SimArray(np.nan, '1.00e+00'),
                                        'R200/R500': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R2500/R500': SimArray(np.nan, '1.00e+00'),
                                        'R500/R200': SimArray(np.nan, '1.00e+00'),
                                        'R500/R2500': SimArray(np.nan, '1.00e+00')}},
             'values': {'all': {'R200': SimArray(641.1774632, 'kpc'),
                                'R2500': SimArray(182.40369082, 'kpc'),
                                'R500': SimArray(414.11443838, 'kpc')},
                        'bin_0': {'R200': SimArray(np.nan, 'kpc'),
                                  'R2500': SimArray(np.nan, 'kpc'),
                                  'R500': SimArray(np.nan, 'kpc')},
                        'bin_1': {'R200': SimArray(np.nan, 'kpc'),
                                  'R2500': SimArray(np.nan, 'kpc'),
                                  'R500': SimArray(np.nan, 'kpc')},
                        'bin_2': {'R200': SimArray(np.nan, 'kpc'),
                                  'R2500': SimArray(np.nan, 'kpc'),
                                  'R500': SimArray(np.nan, 'kpc')},
                        'bin_3': {'R200': SimArray(588.42783886, 'kpc'),
                                  'R2500': SimArray(172.4447052, 'kpc'),
                                  'R500': SimArray(384.46939811, 'kpc')},
                        'bin_4': {'R200': SimArray(886.03334618, 'kpc'),
                                  'R2500': SimArray(253.59792552, 'kpc'),
                                  'R500': SimArray(573.7842392, 'kpc')},
                        'bin_5': {'R200': SimArray(1502.24274674, 'kpc'),
                                  'R2500': SimArray(434.81966803, 'kpc'),
                                  'R500': SimArray(969.18964485, 'kpc')},
                        'bin_6': {'R200': SimArray(np.nan, 'kpc'),
                                  'R2500': SimArray(np.nan, 'kpc'),
                                  'R500': SimArray(np.nan, 'kpc')}}}}

# Rxxx_values = median_R_values['simba-c']

## Conversion ratios for Rxxx values
Rxxx_conversions = median_R_values#['simba-c']


## Median values of Rxxx and Rxxx ratios from simba and simba-c for all groups (N_lgal>=3)) M500>=1e13 Msun
# median_R_values = {
#     'simba': {'ratios': {'R200/R2500': pnb.array.SimArray(3.46687471, '1.00e+00'),
#                       'R200/R500': pnb.array.SimArray(1.52865472, '1.00e+00'),
#                       'R2500/R200': pnb.array.SimArray(0.28844423, '1.00e+00'),
#                       'R2500/R500': pnb.array.SimArray(0.44297859, '1.00e+00'),
#                       'R500/R200': pnb.array.SimArray(0.65416997, '1.00e+00'),
#                       'R500/R2500': pnb.array.SimArray(2.25744545, '1.00e+00')},
#            'values': {'R200': pnb.array.SimArray(645.14886419, 'kpc'),
#                       'R2500': pnb.array.SimArray(183.92734668, 'kpc'),
#                       'R500': pnb.array.SimArray(413.16778666, 'kpc'),
#                       'ratios': {'R200/R2500': pnb.array.SimArray(3.50762883, '1.00e+00'),
#                                  'R200/R500': pnb.array.SimArray(1.56146942, '1.00e+00'),
#                                  'R2500/R200': pnb.array.SimArray(0.28509288, '1.00e+00'),
#                                  'R2500/R500': pnb.array.SimArray(0.44516381, '1.00e+00'),
#                                  'R500/R200': pnb.array.SimArray(0.6404224, '1.00e+00'),
#                                  'R500/R2500': pnb.array.SimArray(2.24636409, '1.00e+00')}}},
#  'simba-c': {'ratios': {'R200/R2500': pnb.array.SimArray(3.44463263, '1.00e+00'),
#                         'R200/R500': pnb.array.SimArray(1.52932679, '1.00e+00'),
#                         'R2500/R200': pnb.array.SimArray(0.29030672, '1.00e+00'),
#                         'R2500/R500': pnb.array.SimArray(0.4459453, '1.00e+00'),
#                         'R500/R200': pnb.array.SimArray(0.6538825, '1.00e+00'),
#                         'R500/R2500': pnb.array.SimArray(2.24242776, '1.00e+00')},
#              'values': {'R200': pnb.array.SimArray(641.1774632, 'kpc'),
#                         'R2500': pnb.array.SimArray(182.40369082, 'kpc'),
#                         'R500': pnb.array.SimArray(414.11443838, 'kpc'),
#                         'ratios': {'R200/R2500': pnb.array.SimArray(3.51515619, '1.00e+00'),
#                                    'R200/R500': pnb.array.SimArray(1.54830985, '1.00e+00'),
#                                    'R2500/R200': pnb.array.SimArray(0.28448238, '1.00e+00'),
#                                    'R2500/R500': pnb.array.SimArray(0.44046687, '1.00e+00'),
#                                    'R500/R200': pnb.array.SimArray(0.64586556, '1.00e+00'),
#                                    'R500/R2500': pnb.array.SimArray(2.27031831, '1.00e+00')}}}}

# Rxxx_values = median_R_values['simba-c']['values']

# ## Conversion ratios for Rxxx values
# Rxxx_conversions = median_R_values['simba-c']['values']['ratios']