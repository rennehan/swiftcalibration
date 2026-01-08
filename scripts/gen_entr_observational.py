import numpy as np
import os
import h5py
from velociraptor.observations.objects import ObservationalData
from unyt import unyt_array, dimensionless
from astropy.cosmology import FlatLambdaCDM
import pynbody as pnb
import scaling_relation_functions as scale

catalog = "CLoGS" # CHOOSE FROM CLoGS, ACCEPT, and sun_2009


# Paths
home_name = "/project/b/babul/aspadawe/data/groups-and-clusters/Oppenheimer+21/"
f_b = 0.157

# Constants
mu_e = 1.14
mu = 0.59
mp = 1.6726e-24
G = 6.672e-08
kboltz = 1.38066e-16
kboltz_v2 = 8.617e-8 # keV/K
cmperkpc = 3.086e+21
cmperMpc = 3.086e+24
keVperK = 8.616e-08
gperMsol = 1.989e+33
cmperkpc = 3.086e+21
cmperpc = 3.086e+18

# Mass range filter (log10(M500/Msun))
M500_min = 12.8
M500_max = 13.2

# Define cosmology
cosmology = FlatLambdaCDM(H0=70, Om0=0.3)
redshift = 0.0

# Output path
output_directory = "./"
output_filename = "entropy_profiles_mass_filtered.hdf5"

def flexbin_med_and_onesigmaspread(x,y,nbins):

    xmed = np.zeros(nbins)
    xlow = np.zeros(nbins)
    xhi  = np.zeros(nbins)

    ymed = np.zeros(nbins)
    ylow = np.zeros(nbins)
    yhi  = np.zeros(nbins)
    
    for i in range(nbins):
        xmed[i] = np.percentile(x,(i+0.5)/nbins*100.)
        xlow[i] = np.percentile(x,(i+0.0)/nbins*100.)
        xhi[i] = np.percentile(x,(i+1.0)/nbins*100.)
        indexes_x = np.where((x>=xlow[i]) & (x<=xhi[i]))
        #print("xlo,xhi,index= ",xlow[i],xhi[i],indexes_x)
        ymed[i] = np.percentile(y[indexes_x],50)
        ylow[i] = np.percentile(y[indexes_x],16)
        yhi[i] = np.percentile(y[indexes_x],84)

    return(xmed,xlow,xhi,ymed,ylow,yhi)


#---------------------------------------------------------------------------------------------------------------
def ACCEPT():
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
    
    save_data = {
        'ACCEPT': {
            'x': xmed,
            'xerr': np.zeros_like(xmed),
            'y': ymed,
            'yerr': 0.5*(yhi - ylow),
            'x_label': r'log($R / R_{500}$)',
            'y_label': r'$\log(K / K_{500})$',
            'name': f'ACCEPT Entropy Profiles ({M500_min} ≤ logM500 ≤ {M500_max})'
        }}

    return save_data

#--------------------------------------------------------------------------------
#clogs
def CLOGS():
    def loadCLoGS_entropy(group_ID):
        """Load CLoGS entropy profiles"""
        fname = home_name + "CLoGS_entropy_profiles/LGG" + group_ID + "_entropy.dat"
        if os.path.isfile(fname):
            r_in, r_out, K, Kmin, Kmax = np.loadtxt(fname, usecols=(0, 1, 2, 3, 4), unpack=True)
            r_mid = (r_in + r_out) / 2.
            return r_mid, K
        else:
            return [], []
        

    """Plot averaged CLoGS entropy profiles"""
    CLoGS_lR_coll = []
    CLoGS_ly_coll = []

    CLoGS_list = home_name + "CLoGS_entropy_profiles/CLoGS_Table4.dat"


    CLoGS = open(CLoGS_list, "r").readlines()
    nCLoGS = len(CLoGS)

    for i in range(nCLoGS):
        CLoGS_ID = CLoGS[i].split()[0]
        CLoGS_T = float(CLoGS[i].split()[1])  # keV
        CLoGS_M500 = scale.M500_Tspec_simbac_N1024L100_renier(pnb.array.SimArray(CLoGS_T, units='keV'))
        CLoGS_lM500 = np.log10(CLoGS_M500)
        
        # Only include systems within mass range
        if CLoGS_lM500 >= M500_min and CLoGS_lM500 <= M500_max:
            CLoGS_entropy_exists = int(CLoGS[i].split()[3])
            if CLoGS_entropy_exists:
                R, y = loadCLoGS_entropy(CLoGS_ID)
                
                if len(R) > 0 and len(y) > 0:
                    CLoGS_redshift = 0.0  # Assuming z=0 for simplicity
                    CLoGS_rhocrit = 8.52e-30 * (1 + CLoGS_redshift) ** 3  # g/cm^3
                    CLoGS_rho500 = 500 * CLoGS_rhocrit  # g/cm^3
                    CLoGS_R500 = (CLoGS_M500 * 1.989e+33 / CLoGS_rho500 * (3 / (4 * np.pi))) ** 0.3333  # cm
                    CLoGS_ne500 = 500 * f_b * CLoGS_rhocrit / (mu_e * mp)  # cm^-3
                    CLoGS_T500 = (G * CLoGS_M500 * 1.989e+33 * mu * mp / (2 * CLoGS_R500 * kboltz))
                    CLoGS_K500 = CLoGS_T500 * CLoGS_ne500 ** (-0.6666) * keVperK
                    
                    lR = np.log10(R / (CLoGS_R500 / cmperkpc))
                    ly = np.log10(y / CLoGS_K500)
                    
                    CLoGS_lR_coll.extend(lR)
                    CLoGS_ly_coll.extend(ly)

    # Calculate averaged profile
    CLoGS_lR, CLoGS_lRlo, CLoGS_lRhi, CLoGS_ly, CLoGS_lylo, CLoGS_lyhi = \
        flexbin_med_and_onesigmaspread(np.asarray(CLoGS_lR_coll), np.asarray(CLoGS_ly_coll), 10)
    clogs_xerr = [CLoGS_lR - CLoGS_lRlo, CLoGS_lRhi - CLoGS_lR]  # asymmetric
    clogs_yerr = [CLoGS_ly - CLoGS_lylo, CLoGS_lyhi - CLoGS_ly]  # asymmetric

    save_data={
        'CLoGS':{
        'x': CLoGS_lR,
        'xerr': clogs_xerr,
        'y': CLoGS_ly,
        'yerr': clogs_yerr,
        'x_label': r'log($R / R_{500}$)',
        'y_label': r'$\log(K / K_{500})$',
        'name': f'CLoGS Entropy Profiles ({M500_min} ≤ logM500 ≤ {M500_max})'}
    }

    return save_data

#--------------------------------------------------------------------------------
#sun2009
def SUN2009():
    def load_Sun2009(group_ID, property_base):
        """Load Sun+2009 observational data"""
        fname = home_name + "Sun2009/" + property_base + "/" + group_ID + "_" + property_base
        if os.path.isfile(fname):
            R, y = np.loadtxt(fname, usecols=(0, 1), unpack=True)
            if property_base == 'temperature':
                y /= kboltz_v2
        else:
            R = []
            y = []
        return R, y


    """Plot averaged Sun+2009 entropy profiles"""
    Sun_lR_coll = []
    Sun_ly_coll = []

    Sungroups_list = home_name + "Sun2009/Sun_groups.cat.sort"

    Sungroups = open(Sungroups_list, "r").readlines()
    nSungroups = len(Sungroups)

    for i in range(nSungroups):
        Sun_ID = Sungroups[i].split()[0]
        Sun_T500 = float(Sungroups[i].split()[1])  # keV
        Sun_R500 = float(Sungroups[i].split()[2])
        Sun_M500_13 = float(Sungroups[i].split()[3])
        
        if Sun_M500_13 < 0:
            Sun_M500_13 = 3e+13 * (Sun_T500 / 1.0) ** (1.5) / 1e+13
        
        Sun_lM500 = np.log10(Sun_M500_13) + 13
        
        # Only include systems within mass range
        if Sun_lM500 >= M500_min and Sun_lM500 <= M500_max:
            lR, ly = load_Sun2009(Sun_ID, "K")  # Entropy data
            
            if len(lR) > 0 and len(ly) > 0:
                Sun_lR_coll.extend(lR)
                Sun_ly_coll.extend(ly)

    # Calculate averaged profile

    Sun_lR, Sun_lRlo, Sun_lRhi, Sun_ly, Sun_lylo, Sun_lyhi = \
        flexbin_med_and_onesigmaspread(np.asarray(Sun_lR_coll), np.asarray(Sun_ly_coll), 15)
    sun_xerr = [Sun_lR - Sun_lRlo, Sun_lRhi - Sun_lR]  # asymmetric
    sun_yerr = [Sun_ly - Sun_lylo, Sun_lyhi - Sun_ly]  # asymmetric

    save_data = {
        'sun_2009':{
        'x': Sun_lR,
        'xerr': sun_xerr,
        'y': Sun_ly,
        'yerr': sun_yerr,       
        'x_label': r'log($R / R_{500}$)',
        'y_label': r'$\log(K / K_{500})$',
        'name': f'Sun+2009 Entropy Profiles ({M500_min} ≤ logM500 ≤ {M500_max})'}
    }

    return save_data

#--------------------------------------------------------------------------------

"""
x = save_data['ACCEPT']['x']
xerr = save_data['ACCEPT']['xerr']
y = save_data['ACCEPT']['y']
yerr = save_data['ACCEPT']['yerr']
x_label = save_data['ACCEPT']['x_label']
y_label = save_data['ACCEPT']['y_label']
name = save_data['ACCEPT']['name']

import matplotlib.pyplot as plt
# Make the plot
plt.errorbar(
    x, y, 
    xerr=xerr, yerr=yerr, 
    fmt='o', capsize=3, label=name
)

x = save_data['CLoGS']['x']
y = save_data['CLoGS']['y']
x_label = save_data['CLoGS']['x_label']
y_label = save_data['CLoGS']['y_label']
name = save_data['CLoGS']['name']

# Make the plot
plt.errorbar(
    x, y, 
    fmt='o', capsize=3, label=name
)

x = save_data['sun_2009']['x']
y = save_data['sun_2009']['y']
x_label = save_data['sun_2009']['x_label']
y_label = save_data['sun_2009']['y_label']
name = save_data['sun_2009']['name']

# Make the plot
plt.errorbar(
    x, y, 
    fmt='o', capsize=3, label=name
)

# Labels and title
plt.xlabel(x_label)
plt.ylabel(y_label)
plt.title(f"Entropy Profiles - {M500_min} ≤ logM500 ≤ {M500_max}")
plt.legend()

plt.savefig("entropy_profiles.png", dpi=300, bbox_inches="tight")
plt.close()

"""

if catalog == 'CLoGS':
    save_data = CLOGS()

if catalog == 'ACCEPT':
    save_data = ACCEPT()

if catalog == 'sun_2009':
    save_data = SUN2009()

output_directory = "/home/b/babul/stlock/swift/swiftcalibration/data/"
output_filename = f"{catalog}_entropy_profiles.hdf5"
# Replace the final section with this:
output_path = os.path.join(output_directory, output_filename)

# Remove the file ONCE at the beginning
if os.path.exists(output_path):
    os.remove(output_path)

for key, val in save_data.items():
    x = val['x'] * dimensionless
    xerr = val['xerr'] * dimensionless
    y = val['y'] * dimensionless
    yerr = val['yerr'] * dimensionless

    processed = ObservationalData()
    processed.associate_x(x, scatter=xerr, comoving=True, description=val['x_label'])
    processed.associate_y(y, scatter=yerr, comoving=True, description=val['y_label'])
    
    # Set appropriate citation for each dataset
    if key == 'ACCEPT':
        processed.associate_citation("ACCEPT", "Cavagnolo+2009")
    elif key == 'CLoGS':
        processed.associate_citation("CLoGS", "CLoGS_citation")
    elif key == 'sun_2009':
        processed.associate_citation("Sun+2009", "Sun+2009")
    
    processed.associate_name(val['name'])
    processed.associate_comment(f"Filtered {key} clusters by logM500 range")
    processed.associate_redshift(redshift, redshift, redshift)
    processed.associate_plot_as("line")
    processed.associate_cosmology(cosmology)

    processed.write(filename=output_path)
    print(f"Saved {key} data to {output_path}")

print(f"\nDONE. File written to {output_path}\n")
