import numpy as np
import os
import h5py
from velociraptor.observations.objects import ObservationalData
from unyt import unyt_array, dimensionless
from astropy.cosmology import FlatLambdaCDM

# Paths
home_name = "/project/b/babul/aspadawe/data/groups-and-clusters/Oppenheimer+21/"
f_b = 0.157

# Constants
mu_e = 1.14
mu = 0.59
mp = 1.6726e-24
G = 6.672e-08
kboltz = 1.38066e-16
keVperK = 8.616e-08
cmperMpc = 3.086e+24
gperMsol = 1.989e+33

# Mass range filter (log10(M500/Msun))
M500_min = 12.8
M500_max = 13.2

# Define cosmology
cosmology = FlatLambdaCDM(H0=70, Om0=0.3)
redshift = 0.0

# Output path
output_directory = "./"
output_filename = "ACCEPT_entropy_profiles_mass_filtered.hdf5"

def loadACCEPTdata(group_ID):
    """Load ACCEPT cluster entropy data"""
    group_ID_prof = np.loadtxt(home_name + "ACCEPT/ACCEPT.profiles.dat",
                               usecols=(0), dtype=(str), unpack=True)
    r_in, r_out, ne, ne_err, K_itpl, K_flat, K_err = np.loadtxt(
        home_name + "ACCEPT/ACCEPT.profiles.dat",
        usecols=(1,2,3,4,5,6,7), unpack=True
    )
    r_mid = (r_in + r_out) / 2.0
    indexes = np.char.equal(group_ID_prof, group_ID)
    return r_mid[indexes], K_flat[indexes]

def save_ACCEPT_entropy_profiles():
    """Save ACCEPT entropy profiles in velociraptor .hdf5 format (mass-filtered)"""
    accept_file = home_name + "ACCEPT/ACCEPT.main.tab"
    ACCEPT = open(accept_file, "r").readlines()
    nACCEPT = len(ACCEPT)

    all_lR = []
    all_lK = []

    for i in range(2, nACCEPT):  # skip header lines
        cols = ACCEPT[i].split()
        name = cols[0]
        z = float(cols[3])
        Tcl = float(cols[7])  # keV

        # Estimate M500 from Tspec scaling relation (rough)
        M500 = 3e13 * (Tcl/1.0)**1.5   # [Msun]
        logM500 = np.log10(M500)

        # Mass filter
        if logM500 < M500_min or logM500 > M500_max:
            continue

        rho_crit = 8.52e-30 * (1+z)**3
        R500 = (M500 * gperMsol / (rho_crit*500) * (3/(4*np.pi)))**(1/3)

        ne500 = 500*f_b*rho_crit/(mu_e*mp)
        T500 = G*M500*gperMsol*mu*mp / (2*R500*kboltz)
        K500 = T500 * ne500**(-2/3) * keVperK

        # Load data
        R, K = loadACCEPTdata(name)
        lR = np.log10(R / (R500/cmperMpc))
        lK = np.log10(K / K500)

        all_lR.extend(lR[np.isfinite(lK)])
        all_lK.extend(lK[np.isfinite(lK)])

    # Convert to arrays
    all_lR = np.array(all_lR)
    all_lK = np.array(all_lK)

    # Bin the data to get a median profile + scatter
    nbins = 30
    percentiles = np.linspace(0, 100, nbins+1)
    xmed = np.zeros(nbins)
    ymed = np.zeros(nbins)
    ylow = np.zeros(nbins)
    yhi = np.zeros(nbins)

    for i in range(nbins):
        xlo, xhi = np.percentile(all_lR, [percentiles[i], percentiles[i+1]])
        mask = (all_lR >= xlo) & (all_lR <= xhi)
        if np.sum(mask) > 0:
            xmed[i] = np.median(all_lR[mask])
            ymed[i] = np.median(all_lK[mask])
            ylow[i] = np.percentile(all_lK[mask], 16)
            yhi[i] = np.percentile(all_lK[mask], 84)

    # Save in velociraptor observational format
    save_data = {
        'ACCEPT': {
            'x': xmed,
            'xerr': np.zeros_like(xmed),
            'y': ymed,
            'yerr': 0.5*(yhi - ylow),
            'x_label': r'log($R / R_{500}$)',
            'y_label': r'$\log(K / K_{500})$',
            'name': f'ACCEPT Entropy Profiles ({M500_min} ≤ logM500 ≤ {M500_max})'
        },
        'CLoGS':{},
        'sun_2009':{}
    }

    """
    x = save_data['ACCEPT']['x']
    xerr = save_data['ACCEPT']['xerr']
    y = save_data['ACCEPT']['y']
    yerr = save_data['ACCEPT']['yerr']
    x_label = save_data['ACCEPT']['x_label']
    y_label = save_data['ACCEPT']['y_label']
    name = save_data['log_data']['name']

    import matplotlib.pyplot as plt
    # Make the plot
    plt.errorbar(
        x, y, 
        xerr=xerr, yerr=yerr, 
        fmt='o', capsize=3, label=name
    )

    # Labels and title
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(name)
    plt.legend()

    plt.savefig("entropy_profiles.png", dpi=300, bbox_inches="tight")
    plt.close()

    """

    for key, val in save_data.items():
        x = val['x'] * dimensionless
        xerr = val['xerr'] * dimensionless
        y = val['y'] * dimensionless
        yerr = val['yerr'] * dimensionless

        processed = ObservationalData()
        processed.associate_x(x, scatter=xerr, comoving=True, description=val['x_label'])
        processed.associate_y(y, scatter=yerr, comoving=True, description=val['y_label'])
        processed.associate_citation("ACCEPT", "Cavagnolo+2009")
        processed.associate_name(val['name'])
        processed.associate_comment("Filtered ACCEPT clusters by logM500 range")
        processed.associate_redshift(redshift, redshift, redshift)
        processed.associate_plot_as("line")
        processed.associate_cosmology(cosmology)

        output_path = os.path.join(output_directory, output_filename)
        if os.path.exists(output_path):
            os.remove(output_path)

        processed.write(filename=output_path)

    print(f"\nDONE. File written to {output_filename}\n")


if __name__ == "__main__":
    save_ACCEPT_entropy_profiles()
