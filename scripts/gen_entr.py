USE_VERSION = 'v0.2b'
import os
import sys
import yt
import unyt
import numpy as np
import caesar

from unyt import unyt_quantity, unyt_array, km, s, Mpc, kg, m, Msun, kpc
from unyt.physical_constants import G, kboltz, mp

from velociraptor.observations.objects import ObservationalData
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--modeldir')
parser.add_argument('--simdir')
parser.add_argument('--sim', type=str)
parser.add_argument('--targethalo', type=int)
args = parser.parse_args()

MODELDIR = args.modeldir
SIMDIR = args.simdir
SIM = args.sim
HALO = args.targethalo

SNAPLIST = [151]

# Define common log-r bins globally so all runs use the same binning
common_log_r_bins = np.linspace(-4, np.log10(5), 30)  # log(R/R500)

def gen_observable(SNAPLIST, SNAPDIR, SIM):
    for j in SNAPLIST:
        SNAP = os.path.join(SNAPDIR, f'{SIM}_{j:03d}.hdf5')
        CAESARFILE = os.path.join(SNAPDIR, f'{SIM}_caesar_{j:03d}.hdf5')

        if not os.path.exists(SNAP):
            print(SNAP, "does not exist")
            continue

        snap = yt.load(SNAP)
        data = snap.all_data()
        obj = caesar.load(CAESARFILE)

        output_directory = SNAPDIR

        h0 = obj.simulation.hubble_constant
        omega_matter = snap.cosmology.omega_matter
        cosmology = FlatLambdaCDM(H0=h0*100 * u.km / u.s / u.Mpc, Om0=omega_matter)

        # Define entropy fields
        def gas_entropy(field, data):
            return data["gas", "kT"] * data["gas", "El_number_density"] ** (-2.0 / 3.0)

        def gas_entropy2(field, data):
            return data["gas", "kT"] * data["PartType0", "El_number_density"] ** (-2.0 / 3.0)

        snap.add_field(("gas", "gas_entropy"), units="keV*cm**2", function=gas_entropy, sampling_type="cell")
        snap.add_field(("PartType0", "gas_entropy"), units="keV*cm**2", function=gas_entropy2, sampling_type="particle", force_override=True)

        redshift = obj.simulation.redshift
        halo = obj.halos[HALO]
        halo_center = halo.minpotpos.in_units('Mpc')

        igrm = data.include_above(('PartType0', 'temperature'), 5e5, units='K')
        igrm = igrm.include_below(('PartType0', 'H_nuclei_density'), 0.13, units='cm**-3')
        igrm = igrm.include_equal(('PartType0', 'DelayTime'), 0)

        igrm_pos_x = igrm['gas', 'particle_position_x']
        igrm_pos_y = igrm['gas', 'particle_position_y']
        igrm_pos_z = igrm['gas', 'particle_position_z']

        r_igrm = np.sqrt((igrm_pos_x - halo_center[0])**2 +
                         (igrm_pos_y - halo_center[1])**2 +
                         (igrm_pos_z - halo_center[2])**2)

        igrm_entropy = igrm['gas', 'entropy']
        igrm_mass = igrm['gas', 'mass']

        radial_bins = unyt.unyt_array(np.append([0], np.logspace(-4, np.log10(5), 45)), 'Mpc')
        radial_bin_centres = (radial_bins[:-1] + radial_bins[1:]) / 2.

        entropy_profile_mean = []
        entropy_profile_std = []

        for bin_lo, bin_hi in zip(radial_bins[:-1], radial_bins[1:]):
            mask = (r_igrm >= bin_lo) & (r_igrm < bin_hi)
            entropies = igrm_entropy[mask]
            masses = igrm_mass[mask]

            if len(entropies) > 0 and np.sum(masses) > 0:
                mean = np.average(entropies, weights=masses)
                variance = np.average((entropies - mean)**2, weights=masses)
                std_dev = np.sqrt(variance)
            else:
                mean = np.nan
                std_dev = np.nan

            entropy_profile_mean.append(mean)
            entropy_profile_std.append(std_dev)

        entropy_profile_mean = np.array(entropy_profile_mean)
        entropy_profile_std = np.array(entropy_profile_std)

        # Normalize to K500
        mu = 0.59
        mu_e = 1.14
        omega_b = obj.simulation.omega_baryon
        omega_m = obj.simulation.omega_matter
        f_b = omega_b / omega_m

        H0_si = (100 * h0 * 1e3) / (3.086e22)  # H0 in 1/s
        rho_crit_si = 3 * H0_si**2 / (8 * np.pi * 6.6743e-11)
        rho_crit = unyt_quantity(rho_crit_si, "kg/m**3").to("Msun/kpc**3")

        M500 = halo.virial_quantities['m500c'].to("Msun")
        R500 = halo.virial_quantities['r500c'].to("m")
        ne_500 = (500 * f_b * rho_crit / (mu_e * mp)).to("cm**-3")
        T500 = (G * M500 * mu * mp / R500).to("keV")
        K500 = (T500 * ne_500**(-2/3)).to("keV*cm**2")

        normalized_r = (radial_bin_centres / halo.virial_quantities['r500c'].to("Mpc")).value
        normalized_k = entropy_profile_mean / K500.value
        normalized_k_std = entropy_profile_std / K500.value

        # Take log safely
        normalized_k_safe = np.where(normalized_k > 0, normalized_k, np.nan)
        log_k = np.log10(normalized_k_safe)
        log_k_error = np.log10(1 + normalized_k_std / normalized_k_safe)

        log_r = np.log10(normalized_r)

        # Interpolate both log_k and log_k_error to common bins
        interp_log_k = np.interp(common_log_r_bins, log_r, log_k, left=np.nan, right=np.nan)
        interp_log_k_err = np.interp(common_log_r_bins, log_r, log_k_error, left=np.nan, right=np.nan)

        # Save
        save_data = {
            'log_data': {
                'x': common_log_r_bins,
                'xerr': np.zeros_like(common_log_r_bins),
                'y': interp_log_k,
                'yerr': interp_log_k_err,
                'x_label': r'log($R / R_{500}$)',
                'y_label': r'$\log(K_{\mathrm{mass}} / K_{500})$',
                'name': 'Entropy Profile'
            }
        }

        for key, val in save_data.items():
            output_filename = f"{SIM}_{j:03d}_entropy_profile_common.hdf5"
            x = val['x'] * unyt.dimensionless
            xerr = val['xerr'] * unyt.dimensionless
            y = val['y'] * unyt.dimensionless
            yerr = val['yerr'] * unyt.dimensionless

            processed = ObservationalData()
            processed.associate_x(x, scatter=xerr, comoving=True, description=val['x_label'])
            processed.associate_y(y, scatter=yerr, comoving=True, description=val['y_label'])
            processed.associate_citation("", "")
            processed.associate_name(val['name'])
            processed.associate_comment("Halo3224 - BAL Variations")
            processed.associate_redshift(redshift, redshift, redshift)
            processed.associate_plot_as("line")
            processed.associate_cosmology(cosmology)

            output_path = os.path.join(output_directory, output_filename)
            if os.path.exists(output_path):
                os.remove(output_path)

            processed.write(filename=output_path)
            
    print("\nDONE\n")

gen_observable(SNAPLIST, os.path.join(MODELDIR, SIMDIR), SIM)
