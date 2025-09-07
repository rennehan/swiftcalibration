USE_VERSION = 'v0.2b'
import os
# os.environ['PYNBODY_RAMS_LOAD'] = '0'
import sys
import yt
import unyt
from unyt import unyt_quantity, unyt_array, km, s, Mpc, kg, m, Msun, kpc
from unyt import Gyr, Mpc, km, s
from unyt.physical_constants import G as G_unyt, kboltz, mp as mp_unyt  # G, k_B, m_p

import caesar

import argparse
import numpy as np
# import pynbody as pnb
# import pynbody.units as u
# import pynbody.array as pna

import gen_sim_data as gen

from velociraptor.observations.objects import ObservationalData
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates entropy profile calculations from caesar files for the given simulation in a given directory.
#### Call it using python gen_entr.py --modeldir=path/to/model --simdir=SIMDIR --sim=SIM --targethalo=caesar halo ID
##########################################################################################################################


parser = argparse.ArgumentParser()
parser.add_argument('--modeldir')
parser.add_argument('--simdir')
parser.add_argument('--sim', type=str)
parser.add_argument('--targethalo', type=int)
args = parser.parse_args()

# caesar file needs to be formatted as snapshot_caesar_151.hdf5 in the sim directory
MODELDIR = args.modeldir
SIMDIR = args.simdir
SIM = args.sim
HALO = args.targethalo

SNAPLIST = [151]

print(SNAPLIST)


def gen_observable(SNAPLIST, SNAPDIR, SIM):
    print(SNAPDIR)
    
    for j in SNAPLIST:
        print(j)
        SNAP = os.path.join(SNAPDIR, '%s_%03d.hdf5' % (SIM, j))
        print(SNAP)
        if not os.path.exists(SNAP):
            print(SNAP, "does not exist")
            continue

        CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%03d.hdf5' % (SIM, j))

        snap = yt.load(SNAP)#, unit_system='cgs')
        data = snap.all_data()
        obj = caesar.load(CAESARFILE)
        # s = pnb.load(SNAP)

        output_directory = SNAPDIR

        h0 = obj.simulation.hubble_constant
        omega_matter = snap.cosmology.omega_matter
        cosmology = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=h0*100. * u.km / u.s / u.Mpc, Om0=omega_matter)

        # ================================================================================================================================================

        # Get data, and some functions

        def gas_entropy(field, data):
            return data["gas", "kT"] * data["gas", "El_number_density"] ** (-2.0/3.0)

        def gas_entropy2(field, data):
            return data["gas", "kT"] * data["PartType0", "El_number_density"] ** (-2.0/3.0)

        snap.add_field(
            ("gas", "gas_entropy"),
            units="keV*cm**2",
            function=gas_entropy,
            sampling_type="cell",
        )

        snap.add_field(
            ("PartType0", "gas_entropy"),
            units="keV*cm**2",
            function=gas_entropy2,
            sampling_type="particle",
            force_override=True
        )

        gas_entropy = data['gas', 'entropy']

        redshift = obj.simulation.redshift
        halo_data = obj.halos

        halo_center = halo_data[HALO].minpotpos
        halo_center = halo_center.in_units('Mpc')
        print(halo_center)

        igrm = data.include_above(('PartType0', 'temperature'), 5e5, units='K')
        igrm = igrm.include_below(('PartType0', 'H_nuclei_density'), 0.13, units='cm**-3')
        igrm = igrm.include_equal(('PartType0', 'DelayTime'), 0)

        # Radial Binning
        radial_bins = unyt.unyt_array(np.append([0], np.logspace(-4, np.log10(5), 45)), 'Mpc')

        halo_center.in_units('Mpc')

        igrm_pos_x = igrm['gas', 'particle_position_x']
        igrm_pos_y = igrm['gas', 'particle_position_y']
        igrm_pos_z = igrm['gas', 'particle_position_z']

        igrm_halo_center = ((igrm_pos_x-halo_center[0])**2 + (igrm_pos_y-halo_center[1])**2 + (igrm_pos_z-halo_center[2])**2)**0.5

        igrm_zoom_entropies = igrm['gas', 'entropy']
        igrm_zoom_masses = igrm['gas', 'mass']

        entropy_profile_median = []
        entropy_profile_mean = []
        entropy_profile_mass_weighted = []
        entropy_profile_mass_weighted_mean_std = []
        for bin_lo, bin_hi in zip(radial_bins[:-1], radial_bins[1:]):
            radial_filter = np.logical_and(igrm_halo_center.in_units('Mpc').value>=bin_lo.in_units('Mpc').value,
                                        igrm_halo_center.in_units('Mpc').value<=bin_hi.in_units('Mpc').value)
            entropies = igrm_zoom_entropies[radial_filter]
            masses = igrm_zoom_masses[radial_filter]
            
            entropy_profile_median.append(np.nanmedian(entropies))
            entropy_profile_mean.append(np.nanmean(entropies))
            if sum(masses)>0:
                K_mass = np.average(entropies, weights=masses)
                entropy_profile_mass_weighted.append(K_mass)
            else:
                entropy_profile_mass_weighted.append(np.nan)
            
            if np.sum(masses) > 0 and len(entropies) > 0:
                mean = np.average(entropies, weights=masses)
                variance = np.average((entropies - mean)**2, weights=masses)
                std_dev = np.sqrt(variance)
                entropy_profile_mass_weighted_mean_std.append(std_dev)
            else:
                entropy_profile_mass_weighted_mean_std.append(np.nan)



        from unyt import Msun, kpc, cm, keV, unyt_quantity
        from unyt.physical_constants import mp, kboltz, G  # ✅ use unyt's version


        # Constants
        mu = 0.59        # Mean molecular weight
        mu_e = 1.14      # Mean molecular weight per free electron

        # Load CAESAR catalog
        h = obj.simulation.hubble_constant
        omega_b = obj.simulation.omega_baryon
        omega_m = obj.simulation.omega_matter

        # Update f_b if desired
        f_b = omega_b / omega_m

        # Critical density
        H0_si = (100 * h * 1e3) / (3.086e22)  # H0 in 1/s
        rho_crit_si = 3 * H0_si**2 / (8 * np.pi * 6.6743e-11)  # kg/m^3
        rho_crit = unyt_quantity(rho_crit_si, "kg/m**3").to("Msun/kpc**3")

        # Choose your halo
        halo = obj.halos[HALO]  # Replace with desired index or search method

        M500 = unyt_quantity(halo.virial_quantities['m500c'].in_units("Msun"), "Msun")
        R500 = unyt_quantity(halo.virial_quantities['r500c'].in_units("m"), "m") # IN METERS FOR TEMP CALC
        ne_500 = (500 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")
        T500 = (G * M500 * mu * mp / R500).to("keV") # removing 1/kb here
        K500 = (T500 * ne_500**(-2/3)).to("keV*cm**2") # removing kb here
        R500 = halo_data[HALO].virial_quantities['r500c'].in_units('Mpc') # redefine r500 and r200 with correct units before binning
        
        radial_bin_centres = (radial_bins[:-1] + radial_bins[1:])/2.

        normalized_r = radial_bin_centres/R500.value
        normalized_k_mean = np.array(entropy_profile_mean)/K500.value
        normalized_k_std = np.array(entropy_profile_mass_weighted_mean_std) / K500.value

        log_r = np.log10(normalized_r)
        log_k = np.log10(normalized_k_mean)
        log_k_error = np.log10(1 + normalized_k_std / normalized_k_mean)

        # Define a common binning across all runs:
        common_log_r_bins = np.linspace(np.min(log_r), np.max(log_r), 30)

        # Interpolate your entropy profile onto these bins:
        interp_log_k = np.interp(common_log_r_bins, log_r, log_k, left=np.nan, right=np.nan)


        # ================================================================================================================================================

        save_data = {
            'log_data': {
                'x': common_log_r_bins,
                'xerr': np.zeros(len(log_r)),
                'y': interp_log_k,
                'yerr':  log_k_error,
                'x_label': r'log($R / R_{500}$)',
                'y_label': r'$\log(K_{\mathrm{mass}} / K_{500})$',
                'name': 'Entropy Profile'
            }
        }


        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        for key, val in save_data.items():
            output_filename = f"{SIM}_{j:03d}_entropy_profile_common.hdf5"
            comment = "Halo3224 - BAL Variations"
            citation = ""
            bibcode = ""
            name = val['name']
            plot_as = "line"
            redshift_lower = redshift
            redshift_upper = redshift

            x = val['x'] * unyt.dimensionless
            xerr = val['xerr'] * unyt.dimensionless
            y = val['y'] * unyt.dimensionless
            yerr = val['yerr'] * unyt.dimensionless

            processed = ObservationalData()
            processed.associate_x(x, scatter=xerr, comoving=True, description=val['x_label'])
            processed.associate_y(y, scatter=yerr, comoving=True, description=val['y_label'])
            processed.associate_citation(citation, bibcode)
            processed.associate_name(name)
            processed.associate_comment(comment)
            processed.associate_redshift(redshift, redshift_lower, redshift_upper)
            processed.associate_plot_as(plot_as)
            processed.associate_cosmology(cosmology)

            output_path = os.path.join(output_directory, output_filename)

            if os.path.exists(output_path):
                os.remove(output_path)

            processed.write(filename=output_path)

        print()

    print("\nDONE\n")

gen_observable(SNAPLIST, os.path.join(MODELDIR, SIMDIR), SIM)