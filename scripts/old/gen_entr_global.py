USE_VERSION = 'v0.2b'
import os
import sys
import yt
import caesar
import unyt
from unyt import physical_constants
import argparse
import numpy as np

import gen_sim_data as gen

from velociraptor.observations.objects import ObservationalData
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

##########################################################################################################################
#### This function generates entropy profile calculations from caesar files for the given simulation in a given directory.
#### Call it using python gen_entr.py --modeldir=path/to/model --simdir=SIMDIR --sim=SIM
##########################################################################################################################


parser = argparse.ArgumentParser()
parser.add_argument('--modeldir')
parser.add_argument('--simdir')
parser.add_argument('--sim', type=str)
args = parser.parse_args()

MODELDIR = args.modeldir
SIMDIR = args.simdir
SIM = args.sim

SNAPLIST = [151]

print(SNAPLIST)

x_units = unyt.keV * unyt.cm**2
length_units = 'Mpccm'
y_units = unyt.Mpc**(-3)

x_label = r'$\log(S \, / \, \mathrm{keV \, cm}^2)$'
y_label = r'$\log(\Phi=dn/d\log S\,[\,\mathrm{dex}^{-1}\,\mathrm{cMpc}^{-3}])$'

min_log_entropy = -2
max_log_entropy = 3

dx = 0.2
minN = 1

def gas_entropy(field, data):
    mu = 0.59  # Mean molecular weight

    # Convert to physical units
    U = data["PartType0", "InternalEnergy"].in_cgs()  # erg/g
    ne_abundance = data["PartType0", "ElectronAbundance"]
    density = data["PartType0", "Density"].in_cgs()  # g/cm³

    # Calculate temperature in Kelvin
    T = (2 * U * physical_constants.mp * mu) / (3 * physical_constants.kboltz)

    # Convert temperature from Kelvin to keV
    T_keV = T * physical_constants.kboltz.to("keV/K")

    # Electron number density in cm^-3
    nH = density / physical_constants.mp
    ne = nH * ne_abundance

    # Entropy in keV * cm² units
    entropy = T_keV * ne ** (-2.0 / 3.0)
    return entropy

def gen_observable(SNAPLIST, SNAPDIR, SIM):
    print(SNAPDIR)
    
    for j in SNAPLIST:
        print(j)
        SNAP = os.path.join(SNAPDIR, '%s_%03d.hdf5' % (SIM, j))
        if not os.path.exists(SNAP):
            print(SNAP, "does not exist")
            continue

        CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%03d.hdf5' % (SIM, j))

        ds = yt.load(SNAP)
        ds.add_field(("gas", "gas_entropy"), function=gas_entropy, units="keV*cm**2", sampling_type="cell")

        Vcom = ds.domain_width.in_units(length_units).prod()

        print('Loading caesar object')
        obj = caesar.load(CAESARFILE)
        redshift = obj.simulation.redshift
        h0 = obj.simulation.hubble_constant
        boxsize = obj.simulation.boxsize.in_units(length_units)

        print("Comoving Volume=%s" % Vcom)
        print('z=%s' % redshift)
        print('h0=%s' % h0)

        omega_matter = ds.cosmology.omega_matter
        cosmology = FlatLambdaCDM(name=r'Flat $\Lambda$CDM', H0=h0*100. * u.km / u.s / u.Mpc, Om0=omega_matter)

        pos_Mpc = np.array([gal.pos.in_units(length_units) for gal in obj.galaxies])

        entropy_values = []

        for gal in obj.galaxies:
            if gal.contamination > 0:
                continue
            gal_pos = gal.pos.in_units("Mpc")
            sp = ds.sphere(gal_pos, (30, 'kpc'))
            entropy = sp["gas", "gas_entropy"]
            mass = sp["gas", "mass"]

            if len(entropy) > 0:
                K_mass = np.average(entropy, weights=mass)
                entropy_values.append(K_mass)
            else:
                entropy_values.append(np.nan)

        entropy_values = np.array(entropy_values)
        entropy_values = entropy_values[np.isfinite(entropy_values)]
        print("Min log10 entropy:", np.min(np.log10(entropy_values)))
        print("Max log10 entropy:", np.max(np.log10(entropy_values)))

        print(f"Number of galaxies with valid entropy: {len(entropy_values)}")

        log_entropy, log_phi, bins, nbins, filter_ = gen.mass_function_with_error(
            entropy_values, pos_Mpc, boxsize Vcom, dlogM=dx,
            min_logM=min_log_entropy, max_logM=max_log_entropy,
            calc_min_logM=False, calc_max_logM=False, minN=minN
        )

        print(f"Number of bins after filtering: {nbins}")

        save_data = {
            'log_data': {
                'x': np.log10(log_entropy),
                'xerr': np.zeros(len(log_entropy)),
                'y': np.log10(log_phi),
                'yerr': np.zeros(len(log_phi)),
                'x_label': r'$\log(S / \mathrm{keV\,cm^2})$',
                'y_label': r'$\log(\Phi)$',
                'name': 'Entropy Function'
            }
        }

        output_directory = SNAPDIR

        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        for key, val in save_data.items():
            output_filename = f"{SIM}_entropy_function_{key}_min{minN}gal_{j:03d}.hdf5"
            comment = f"h-corrected for SWIFT using cosmology: {cosmology}."
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