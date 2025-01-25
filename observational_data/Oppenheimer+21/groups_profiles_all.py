#import yt
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import mpl_toolkits.axes_grid1 as axgrid
import scipy.stats as scist
import sys
import h5py
from astropy import constants as const
import os.path
import group_data

mu_e = 1.14
mu = 0.59
mp = 1.6726e-24
G = 6.672e-08
kboltz = 1.38066e-16
cmperkpc = 3.086e+21
cmperMpc = 3.086e+24
keVperK = 8.616e-08
gperMsol = 1.989e+33
cmperkpc = 3.086e+21
cmperpc = 3.086e+18

home_name = "/home/ariboppe/groups_review/"

def loadAndersonWang():
    fname = home_name + "data/anderson2015_wang2016.M500_LX.dat"
    if(os.path.isfile(fname)):
        LX,LX_err,M500,M500_err = np.loadtxt(fname,usecols=(6,7,8,9),unpack=True)

    return(M500,M500_err,LX,LX_err)
    
def loadSun2009(group_ID, property_base):
    fname = home_name + "Sun2009/" + property_base + "/" + group_ID + "_" + property_base
    #print("fname= ", fname)
    if(os.path.isfile(fname)):
        R,y = np.loadtxt(fname,usecols=(0,1),unpack=True)
        print("y[-1]= ", y[-1])
    else:
        R = [-99,-99]
        y = [-99,-99]

    return(R,y)

def loadCLoGS_entropy(group_ID):
    r_in,r_out,K, Kmin, Kmax = np.loadtxt(home_name + "CLoGS_entropy_profiles/LGG" + group_ID + "_entropy.dat", usecols=(0,1,2,3,4), unpack=True)
    r_mid = (r_in+r_out)/2.
    return(r_mid,K)                                


def loadACCEPTdata(group_ID,property):

    group_ID_prof = np.loadtxt(home_name + "ACCEPT/ACCEPT.profiles.dat", usecols=(0), dtype=(str), unpack=True)
    r_in,r_out,ne,ne_err,K_itpl,K_flat,K_err  = np.loadtxt(home_name + "ACCEPT/ACCEPT.profiles.dat", usecols=(1,2,3,4,5,6,7),unpack=True)

    #np.genfromtxt(home_name + "ACCEPT/ACCEPT.profiles.dat",dtype=(str,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float,float))#,unpack=True)
    
    r_mid = (r_in+r_out)/2.
    #print("group_ID,group_ID_prof=",group_ID,group_ID_prof[0],group_ID_prof[1],group_ID_prof[2],r_in[0],r_in[1],r_ins[2])
    indexes = np.char.equal(group_ID_prof,group_ID)
    #print(indexes)

    r_mid = r_mid[indexes]
    K_flat = K_flat[indexes]
    #print("rmid, Kitpl= ", r_mid, K_itpl)
    return(r_mid, K_flat)

  #  if(property=="K"):
    #    return(r_mid[indexes],K_flat[indexes])        
    

def loadSun2011pressure():
    fname = home_name + "Sun2011/Sun2011_normalizedpressure.groups.dat"
    R, P_med, P_lo, P_hi = np.loadtxt(fname,usecols=(0,1,2,3),unpack=True)
    return(R,P_med,P_lo,P_hi)

def loadLovisari2015dens():
    fname = home_name + "Lovisari2015/ne_groups.txt"
    R, ne_med, ne_lo, ne_hi = np.loadtxt(fname,usecols=(0,1,2,3),unpack=True)
    return(R,ne_med,ne_lo,ne_hi)

def IntegrateGas(Radius,Mass,OuterRadius):

    TotMass = 0.
    for i in range(len(Radius)):
        if(i>0):
            if(Radius[i]<OuterRadius):
                TotMass += Mass[i]
            if((Radius[i]>OuterRadius) & (Radius[i-1]<OuterRadius)):
                TotMass += Mass[i]*(OuterRadius-Radius[i-1])/(Radius[i]-Radius[i-1])

    return(TotMass)

def IntegrateDensity(Radius,Density,OuterRadius):

    TotMass = 0.
    Router = -99
    for i in range(len(Radius)):
        if((i>0) & (i<len(Radius)-1)):
            Router = (Radius[i] + (Radius[i+1]-Radius[i])/2.)
            Rinner = (Radius[i] - (Radius[i]-Radius[i-1])/2.)
            Volume = 4*np.pi/3*((Router)**3 - (Rinner)**3)*(cmperkpc)**3 # cm^3            
            if(Router < OuterRadius):
                TotMass += 10**Density[i]*(mu_e*mp)*Volume
            if((Router > OuterRadius) & (Rinner < OuterRadius)):
                TotMass += 10**Density[i]*(mu_e*mp)*Volume*(OuterRadius-Rinner)/(Router-Rinner)
            #print("Intdens= ", i,Router,Rinner,OuterRadius,Density[i],Volume,TotMass)
            
    if(Router>OuterRadius):    
        return(TotMass)
    else:
        return(-99)
            
def AbelTransform(Radius,Mass):

    SurfaceDensity = np.zeros(len(Radius))
        
    for i in range(len(Radius)):    
        for j in range(len(Radius)):
            if((Radius[j]>=Radius[i]) & (Radius[j]<Radius[-1])):
                path = (Radius[j+1]-Radius[j])*Radius[j+1]/np.sqrt(Radius[j+1]**2-Radius[i]**2)
                delta_R = Radius[j+1]-Radius[j]
                mean_R = (Radius[j+1]+Radius[j])/2.
                Volume = 4*np.pi/3*((mean_R+delta_R/2.)**3-(mean_R-delta_R/2.)**3)
                #print(i,j,Volume, Mass[j], path)
                SurfaceDensity[i] += 2*path*Mass[j]/Volume
    
    return(SurfaceDensity)

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

def return_most_massive_topbin(x,y,nbins):

    xmin = np.percentile(x,(nbins-0.5)/nbins*100.)
    massive_indexes = np.where(x>=xmin)
    print("xmin,massive_index= ",xmin,massive_indexes)
    return(x[massive_indexes],y[massive_indexes])


def calc_med_and_flexspread(x, y, nbins,bin_lo,bin_hi,percent_lo,percent_hi):
    med = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, 50), bins=nbins, range=[(bin_lo,bin_hi)])[0]
    range_low = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, percent_lo), bins=nbins, range=[(bin_lo,bin_hi)])[0]
    range_high = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=lambda y: np.percentile(y, percent_hi), bins=nbins, range=[(bin_lo,bin_hi)])[0]
    mean = scist.binned_statistic(np.asarray(x), np.asarray(y),statistic=np.mean, bins=nbins, range=[(bin_lo,bin_hi)])[0]
    return(med,range_low,range_high,mean)

file_list = sys.argv[1]
property = sys.argv[2]
coloring = "mass"
legend = int(sys.argv[3])
multi = int(sys.argv[4])

base_name = "/data1/extboppe/"
base_name = "/hpcdata5/users/extboppe/"

title_label = file_list.split(".")[0]

cm = plt.get_cmap('plasma') 
if(multi==1):
    cm = plt.get_cmap('viridis') 


lRlo = -2.0
lRhi = 0.6 # 0.3
nrbins = 26
lRbins = np.linspace(lRlo,lRhi,nrbins)

#if(len(sys.argv)>5):
lhalo_lo = float(sys.argv[5])
lhalo_hi = float(sys.argv[6])
nhalo_bins = int(sys.argv[7])
mhalo_str = '.MH%d_%d_%dbins'%(lhalo_lo*10,lhalo_hi*10,nhalo_bins) 

#else:
#    nhalo_bins = 1
#    mhalo_str = ''

if(multi==1):
    plot_individual_observations = 0
else:
    plot_individual_observations = 0

topborder = 0.93
leftborder = 0.16
rightborder = 0.98
if(multi==2):
    topborder = 0.98

if(coloring=="mass"):
    cbar_lo = lhalo_lo
    cbar_hi = lhalo_hi
    cbarlabel = "log$M_{500}/M_{\odot}$"

if(property=="entropy_norm"):
    lylo = -1.5
    lyhi = 1.0
    if((multi==0) | (multi==2)):
        lylo = -2.0
        lyhi = 0.5
        if(multi==2):
            lylo = -2.2
            
    ylabel = "log $K/K_{500}$" # [keV cm$^{2}$]"
    leftborder = 0.14
if(property=="pressure_norm"):
    lylo = -1.0
    lyhi = 2.0
    ylabel = "log $P/P_{500}$" # [keV cm$^{2}$]"
    leftborder = 0.14
if(property=="pressure"):
    lylo = 1.0
    lyhi = 6.0
    ylabel = "log $P$" # [keV cm$^{2}$]"
    leftborder = 0.14

if(property=="density"):
    lylo = -5.0 #-4.0
    lyhi = -1.0
    ylabel = "log $n_e$ [cm$^{-3}$]"
if(property=="temperature"):
    lylo = -1.0
    lyhi = 1.0
    ylabel = "log $T$ [keV]"
if(property=="metallicity"):
    lylo = -1.0
    lyhi = 0.0
    ylabel = "log $Z/Z_{\odot}$"
    Zsol = 0.0134
if(property=="XraySB_XMMcounts"):
    lylo = -4.0
    lyhi = 0.0
    ylabel = "log SB [cts s$^{-1}$ arcmin$^{-2}$]"
    XMM_erg_to_count_conversion = -39.5
    lRlo = -2.0
    lRhi = 0.6 # 0.3
    nrbins = 9
    lRbins = np.linspace(lRlo,lRhi,nrbins)
if(property=="coolgas"):
    lylo = 3.0
    lyhi = 8.0
    ylabel = "log $\Sigma_{\mathrm{cool}}$ [$M_{\odot}$ kpc$^{-2}$]"
if(property=="warmhotgas"):
    lylo = 4.0
    lyhi = 9.0
    ylabel = "log $\Sigma_{\mathrm{warm-hot}}$ [$M_{\odot}$ kpc$^{-2}$]"
if(property=="dispersionmeasure"):
    ylo = 0
    yhi = 1000
    ylabel = "log DM$_{\mathrm{CGM}}$ [pc cm$^{-3}$]"
    lRlo = 1.0
    lRhi = 3.0 # 0.3
    nrbins = 25
    lRbins = np.linspace(lRlo,lRhi,nrbins)
if(property=="fgas"):
    lylo = 0.00
    lyhi = 0.16
    lMlo = 11.9
    lMhi = 14.7
    cbar_lo = lMlo
    cbar_hi = lMhi
    ylabel = "$f_{\mathrm{gas}}$"
if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    lylo = 10.5
    lyhi = 13.0
    lMlo = 12.9
    lMhi = 14.6
    cbar_lo = lMlo
    cbar_hi = lMhi
    if(property=="m100kpc"):
        ylabel = "log $M_{\mathrm{IGrM,<100 kpc}}$ [$M_{\odot}$]" 
    if(property=="m200kpc"):
        ylabel = "log $M_{\mathrm{IGrM,<200 kpc}}$ [$M_{\odot}$]" 
    if(property=="m400kpc"):
        ylabel = "log $M_{\mathrm{IGrM,<400 kpc}}$ [$M_{\odot}$]" 
if(property=="LX"):
    lylo = 37.0
    lyhi = 44.5
    lMlo = 12.0
    lMhi = 15.0
    cbar_lo = lMlo
    cbar_hi = lMhi
    ylabel = "log $L_{\mathrm{X}}$ [erg s$^{-1}$]" 
if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    fig_all = plt.figure(figsize=(4.5,5.0))
    ax_all = fig_all.add_subplot(111)

    fig_coll = plt.figure(figsize=(4.5,5.0))
    ax_coll = fig_coll.add_subplot(111)
    leftborder = 0.18
else:
    fig_all = plt.figure(figsize=(6.0,5.0))
    ax_all = fig_all.add_subplot(111)

    fig_coll = plt.figure(figsize=(6.0,5.0))
    ax_coll = fig_coll.add_subplot(111)

if(property=="LX"):
    fig_all = plt.figure(figsize=(5.0,5.0))
    ax_all = fig_all.add_subplot(111)

    fig_coll = plt.figure(figsize=(5.0,5.0))
    ax_coll = fig_coll.add_subplot(111)
    leftborder = 0.18

   
f_b = 0.04825/0.307


if(multi):
    sim_list = file_list
    sims = open(sim_list,"r").readlines()
    nsims = len(sims)
else:
    nsims = 1

file_base = file_list

for isim in range(nsims):
    
    lM500_coll = []
    mass_coll = []
    R_coll = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    y_coll = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    y_median = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    y_low = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    y_high = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    y_mean = [[1e-99 for i in range(1)] for i in range(nhalo_bins)]
    count_halos = np.zeros(nhalo_bins)

    if(multi):
        file_list = sims[isim].split()[0]
        sim_color = sims[isim].split()[1]
        sim_label = sims[isim].split()[2]
        sim_linestyle = sims[isim].split()[3]
        sim_zorder = int(sims[isim].split()[4])

        cbar_lo = lhalo_lo
        cbar_hi = lhalo_hi
        if(property=="entropy_norm"):
            title_label = "Entropy"
        if((property=="pressure_norm") | (property=="pressure")):
            title_label = "Pressure"
        if(property=="density"):
            title_label = "Density"
        if(property=="temperature"):
            title_label = "Temperature"            
        if(property=="metallicity"):
            title_label = "Metallicity"            
        if(property=="XraySB_XMMcounts"):
            title_label = "Xray Surface Brightness"            
        if(property=="fgas"):
            title_label = "Gas Fraction inside $R_{500}$"
        if(property=="m100kpc"):
           title_label = "Hot Gas Inside 100 kpc"
        if(property=="m200kpc"):
           title_label = "Hot Gas Inside 200 kpc"
        if(property=="m400kpc"):
           title_label = "Hot Gas Inside 400 kpc"
        if(property=="LX"):
           title_label = "X-ray Luminosity 0.15-1.0$R_{500}$"
        if(property=="dispersionmeasure"):
            title_label = "Dispersion Measure"            


    else:
        sim_color = "blue"
        sim_label = title_label
        sim_linestyle = "-"
        sim_zorder = 5

    halos = open(file_list,"r").readlines()

    nhalos = len(halos)

    for i in range(nhalos):

        file_type = halos[i].split()[0].split(".")[-1]
        if(file_type == "hdf5"):
            halo_hdf5_file = h5py.File(base_name + halos[i].split()[0])

            print(halo_hdf5_file)

            if(property=="LX"):
                #if not "/EventFileAnalysis".attrs["LSoft_01R200"] in halo_hdf5_file:
                if not "/EventFileAnalysis/fR200_LSoftDens" in halo_hdf5_file:
                    continue
            else:
                if not "/RadialProfiles/R_Mass" in halo_hdf5_file:
                    continue

            redshift = 0.0
            lM200 = halo_hdf5_file['Halo'].attrs['lM200']
            lMStar = halo_hdf5_file['Halo'].attrs['lMStar']
            SFR = halo_hdf5_file['Halo'].attrs['SFR']
            if(lM200>20): # SIMBA got is wrong
                lM200= np.log10(lM200)
            if(lMStar>20): # SIMBA is wrong
                lMStar= np.log10(lMStar)
            M500 = 10**(lM200-0.16) # Msol
            lM500 = np.log10(M500)
        else:
            home_name =  "/home/ariboppe/groups_review/"
            ascii_file = home_name + halos[i].split()[0]
            #print(halos[i].split("_")[2][:-5])
            if(property=="entropy_norm"):
                print("halos[i]= ", halos[i])
                if(file_type == "hdf5"):
                    redshift = float(halos[i].split("_")[2][:-5])
                else:
                    redshift = float(halos[i].split("_")[1][:-5])

                M500,entropy,R = np.genfromtxt(ascii_file, usecols=(0,8,9),unpack=True)
                M500 = M500[2]
                entropy = entropy[2:]
                R = R[2:]
                lM500 = np.log10(M500)
            if((property=="density") | (property=="temperature") | (property=="pressure_norm") | (property=="pressure")):
                print(halos[i])
                redshift = float(halos[i].split("_")[1][:-5])
                M500,TkeV,density,R = np.genfromtxt(ascii_file, usecols=(0,8,9,10),unpack=True)
                M500 = M500[2]
                temperature = TkeV[2:]/keVperK
                density = density[2:]
                pressure = density*temperature
                R = R[2:]
                lM500 = np.log10(M500)

                lM200 = lM500+0.16
                lMStar = -99
                #print("M500= ", M500)
                #print("redshift= ", redshift)
        

        #print("hdf5name, SFR,lMStar,LM200=",halo_hdf5_file, SFR,lMStar,lM200)

        if((property=="entropy_norm") | (property=="density") | (property=="temperature") | (property=="pressure_norm") | (property=="pressure") | (property=="metallicity") | (property=="XraySB_XMMcounts") |  (property=="coolgas") |  (property=="warmhotgas") | (property=="dispersionmeasure")):
            if( (lM500>=lhalo_lo) & (lM500<lhalo_hi)):
                hbin =int((lM500-lhalo_lo)/((lhalo_hi-lhalo_lo)/nhalo_bins))
            else:
                continue 
            
        rhocrit = 8.52e-30*(1+redshift)**3 # g/cm^3  
        ne500 = 500*f_b*rhocrit/(mu_e*mp) # cm^-3

        rho500 = 500*rhocrit  # g/cm^3
        R500 = (M500*1.989e+33/rho500*(3/(4*np.pi)))**0.3333   # cm
        
        T500 = (G*M500*1.989e+33*mu*mp/(2*R500*kboltz))
        K500 = T500*ne500**-0.6666*keVperK

        P500 = ne500*T500
    
        print ("ID,M500,R500,T500,ne500,K500= ",halos[i].split()[0],np.log10(M500),np.log10(R500),np.log10(T500),np.log10(ne500),np.log10(K500))

        if(property=="entropy_norm"):
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
                entropy = np.asarray(halo_hdf5_file["RadialProfiles/Entropy_hot"])
                entropy = np.where(np.isnan(entropy), 1e-99, entropy)
            #print(R, entropy)
            R = R*cmperkpc/R500
            y = entropy/K500
            #print(R, y)
        if((property=="pressure_norm") | (property=="pressure")):
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
                density = np.asarray(halo_hdf5_file["RadialProfiles/Density_hot"])
                temperature = np.asarray(halo_hdf5_file["RadialProfiles/Temperature_hot"])
                pressure = density*temperature/(mp*mu_e)
                pressure = np.where(np.isnan(pressure), 1e-99, pressure)
                #pressure = np.asarray(halo_hdf5_file["RadialProfiles/Pressure_hot"])
            R = R*cmperkpc/R500
            if(property=="pressure_norm"):
                y = pressure/P500
            if(property=="pressure"):
                y = pressure
            #print("pressure = " , y)
            #print(R, y)
        if(property=="density"):
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
                density = np.asarray(halo_hdf5_file["RadialProfiles/Density_hot"])
                density = np.where(np.isnan(density), 1e-99, density)
            R = R*cmperkpc/R500
            if(file_type == "hdf5"):
                y = density/(mp*mu_e)
            else:
                y = density
        if(property=="temperature"):
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
                temperature = np.asarray(halo_hdf5_file["RadialProfiles/Temperature_hot"])
                temperature = np.where(np.isnan(temperature), 1e-99, temperature)
            R = R*cmperkpc/R500
            y = temperature*keVperK
        if(property=="metallicity"):
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
                metallicity = np.asarray(halo_hdf5_file["RadialProfiles/Metallicity_hot"])
                metallicity = np.where(np.isnan(metallicity), 1e-99, metallicity)
            R = R*cmperkpc/R500
            y = metallicity/Zsol
        if(property=="XraySB_XMMcounts"):
            XraySB = np.asarray(halo_hdf5_file["EventFileAnalysis/LSoftDens"])
            y = 10**(np.log10(XraySB) + XMM_erg_to_count_conversion)
            R_edges = np.asarray(halo_hdf5_file["EventFileAnalysis/R"])
            R_kpc = np.asarray([10**((np.log10(R_edges[i+1])+np.log10(R_edges[i]))/2.) for i in range(len(R_edges)-1)])
            print("R_kpc,y= ", R_kpc,y)
            R = R_kpc*cmperkpc/R500            
        if(property=="coolgas"):
            if(file_type == "hdf5"):            
                R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
                coolgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass_cool"]) # g
            SurfaceDensity = AbelTransform(R,coolgas)/gperMsol
            R = R*cmperkpc/R500
            y = SurfaceDensity
        if((property=="warmhotgas") | (property=="dispersionmeasure")):
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
                warmhotgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass_warmhot"]) # g
            SurfaceDensity = AbelTransform(R,warmhotgas)
            if(property=="warmhotgas"):
                R = R*cmperkpc/R500
                y = SurfaceDensity/gperMsol # Msol/kpc^2
            if(property=="dispersionmeasure"):
                R = R #kpc
                y = SurfaceDensity/(mp*mu_e) #electrons/kpc^2
                y = y/cmperkpc**2/cmperpc
        if(property=="fgas"):
            R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
                allgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass"]) # g

            Mgas500 = IntegrateGas(R,allgas,R500/cmperkpc)/gperMsol
            fgas500 = Mgas500/M500
            y = fgas500
        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
            R = np.asarray(halo_hdf5_file["RadialProfiles/R"]) # kpc
            if(file_type == "hdf5"):
                R = np.asarray(halo_hdf5_file["RadialProfiles/R_Mass"]) # kpc
                hotgas = np.asarray(halo_hdf5_file["RadialProfiles/Mass_hot"]) # g

            if(property=="m100kpc"):
                y = IntegrateGas(R,hotgas,100)/gperMsol
            if(property=="m200kpc"):
                y = IntegrateGas(R,hotgas,200)/gperMsol
            if(property=="m400kpc"):
                y = IntegrateGas(R,hotgas,400)/gperMsol

            y = np.log10(y)
            
        if(property=="LX"):
            y = np.log10(halo_hdf5_file['EventFileAnalysis'].attrs['LSoft_01R200'])
            ax_all.plot(lM500,y,marker='o',color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))

        if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):
            ax_all.plot(lM500,y,marker='o',color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))
            print("LM500,y= ",lM500,y)
            lM500_coll.extend([lM500])
            mass_coll.extend([y])
        
            #print("M,fgas=", lM500, fgas500, lM_coll,fgas_coll)
                
        else:     # radial profile               
            if((lM500>=cbar_lo) & (lM500<cbar_hi)):
                if(property=="dispersionmeasure"):
                    ax_all.plot(np.log10(R),y,color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))
                else:
                    print("R= ", R, len(R))
                    print("y= ", y, len(y))
                    ax_all.plot(np.log10(R),np.log10(y),color=cm((lM500-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lM500*100))

            if((hbin>=0) & (hbin<nhalo_bins)):
                R_coll[hbin].extend(R)
                y_coll[hbin].extend(y)
                count_halos[hbin] += 1

    if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):
        nflexbins = int(len(mass_coll)/20)
        if(nflexbins>30):
            nflexbins = 30
        M,Mlo,Mhi,y,ylo,yhi=flexbin_med_and_onesigmaspread(np.asarray(lM500_coll),np.asarray(mass_coll),nflexbins)
        M_massive,y_massive = return_most_massive_topbin(np.asarray(lM500_coll),np.asarray(mass_coll),nflexbins)
    
        ax_coll.plot(M, y, color = sim_color ,ls=sim_linestyle,lw=4,label=sim_label,zorder=sim_zorder) #max(y)*100)
        ax_coll.fill_between(M,ylo,yhi, facecolor=sim_color, alpha=0.2,zorder=sim_zorder) #max(y)*100-1)
        ax_coll.plot(M_massive,y_massive,color=sim_color,marker='.',ls='',zorder=sim_zorder) #max(y)*100)
    else:
        for h in range(nhalo_bins):

            lMhalo  = lhalo_lo + (lhalo_hi-lhalo_lo)*(h+0.5)/nhalo_bins
            lMhalo_min = lhalo_lo + (lhalo_hi-lhalo_lo)*(h)/nhalo_bins
            lMhalo_max = lhalo_lo + (lhalo_hi-lhalo_lo)*(h+1)/nhalo_bins

            if(property=="dispersionmeasure"):
                y_median[h], y_low[h], y_high[h], y_mean[h] = calc_med_and_flexspread(np.log10(R_coll[h]),y_coll[h],nrbins,lRlo,lRhi,16,84)
                print("y_median= ",y_median[h])
            else:
                y_median[h], y_low[h], y_high[h], y_mean[h] = calc_med_and_flexspread(np.log10(R_coll[h]),y_coll[h],nrbins,lRlo,lRhi,16,84)
                print("y_median= ",y_median[h])

            #print("calc_med_and_flexspread= ", y_median[h], y_low[h], y_high[h],len(y_coll[h]), sim_label, nrbins, count_halos[h])
            if(multi):
                if(multi==1):
                    if(property=="dispersionmeasure"):                    
                        ax_coll.plot(lRbins, y_median[h], color = sim_color,lw=4,label=sim_label)
                        ax_coll.fill_between(lRbins,y_low[h],y_high[h], facecolor=sim_color, alpha=0.2)
                        ###ax_coll.plot(lRbins, y_mean[h], color = sim_color,lw=4, ls=':')                    
                    else:
                        ax_coll.plot(lRbins, np.log10(y_median[h]), color = sim_color,lw=4,label=sim_label)
                        ax_coll.fill_between(lRbins,np.log10(y_low[h]),np.log10(y_high[h]), facecolor=sim_color, alpha=0.2)
                        ###ax_coll.plot(lRbins, np.log10(y_mean[h]), color = sim_color,lw=4, ls=':')                    
                else: #multi==2
                    if(sim_linestyle=="-"):
                        ax_coll.plot(lRbins, np.log10(y_median[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=sim_linestyle, lw=4) ###,label="log[$M_{500}$]=%4.1f-%4.1f"%(lMhalo_min,lMhalo_max)) #BDO 5/14/21, replacing with colorbar
                        ax_coll.fill_between(lRbins,np.log10(y_low[h]),np.log10(y_high[h]), facecolor=cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)), alpha=0.2)
                        ###ax_coll.plot(lRbins, np.log10(y_mean[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=":", lw=4)
                    else:
                        ax_coll.plot(lRbins, np.log10(y_median[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=sim_linestyle, lw=4)
                        # No fill_between for 2nd simulation.  
                        
            else:
                if(property=="dispersionmeasure"):                    
                    ax_coll.plot(lRbins, y_median[h], color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),lw=4)
                    ax_coll.fill_between(lRbins,y_low[h],y_high[h], facecolor=cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)), alpha=0.2)
                    ###ax_coll.plot(lRbins, y_mean[h], color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=":",lw=4)
                else:
                    ax_coll.plot(lRbins, np.log10(y_median[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),lw=4) ###,label="log[$M_{500}$]=%4.1f-%4.1f"%(lMhalo_min,lMhalo_max)) #BDO 5/14/21, replacing with colorbar
                    ax_coll.fill_between(lRbins,np.log10(y_low[h]),np.log10(y_high[h]), facecolor=cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)), alpha=0.2)
                    ###ax_coll.plot(lRbins, np.log10(y_mean[h]), color = cm((lMhalo-cbar_lo)/(cbar_hi-cbar_lo)),zorder=int(lMhalo*100),ls=":",lw=4)


            ###
            
            data_output_filename_rvir = 'profile_coll_' + sims[isim].split()[2] + '.' + property  + mhalo_str +  '.dat'
            
            data_output_file_rvir = open(data_output_filename_rvir, "w")

            data_output_file_rvir.write("#R/R500 median 16percentile 84percentile nhalos  \n")

            for k in range(len(lRbins)):
                data_output_file_rvir.write("%5.3e %5.2f %5.2f %5.2f %4d\n"%(lRbins[k],np.log10(y_median[h][k]),np.log10(y_low[h][k]),np.log10(y_high[h][k]),count_halos[h]))

            data_output_file_rvir.close()
            ###
            
if(multi==2):
    for isim in range(nsims):
        sim_label = sims[isim].split()[2]
        sim_linestyle = sims[isim].split()[3]
        ax_coll.plot([-99,-99],[-99,-99], color = 'k', ls=sim_linestyle, lw=4, label=sim_label)
                

if(property=="entropy_norm"):
    ax_coll.plot(lRbins, np.log10(1.32/1.09 * (10**(lRbins))**1.1),color="gray",lw=2,ls=":",label="Voit+ (2005)")
    ax_all.plot(lRbins, np.log10(1.32/1.09 * (10**(lRbins))**1.1),color="gray",lw=2,ls=":",label="Voit+ (2005)")

if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    lR = np.linspace(-1.0,3.3,87) # log kpc
    lM500 = np.linspace(12.9,15.1,23) # log Msolar

    Mbase_all = []
    lM500_all = []
    for i in range(len(lM500)):        
        R500 = (10**lM500[i]*1.989e+33/rho500*(3/(4*np.pi)))**0.3333   # cm
        v_c = np.sqrt(G*10**lM500[i]*gperMsol/(R500))
        Density = f_b*v_c**2/(4*np.pi*G*10**lR*cmperkpc) #/(mu_e*mp)
        #Density = (1.32* (((10**lR)/(R500/cmperkpc))**1.1)*K500/T500)**-1.5 # cm**-3
        print("R500,v_c, Density= " , R500,v_c, Density)
        if(property=="m100kpc"):
            Mbase = IntegrateDensity(10**lR,Density,100)/gperMsol
        if(property=="m200kpc"):
            Mbase = IntegrateDensity(10**lR,Density,200)/gperMsol
        if(property=="m400kpc"):
            Mbase = IntegrateDensity(10**lR,Density,400)/gperMsol

        print("LM500,Mbase= ",lM500[i],Mbase)
        Mbase_all.extend([Mbase])
        #lM500_all.extend([lM500])

    print("lM500,Mbase= ",lM500,np.log10(Mbase_all))
    ax_coll.plot(lM500,np.log10(Mbase_all),color="gray",lw=2,ls=":")

    
    #load Sun 2009 data:

plot_only_Lovisari_Sun_groups = True
if((property=="entropy_norm") | (property=="density") | (property=="temperature") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    if(plot_only_Lovisari_Sun_groups):
        group_data.load_and_plot_Lovisari_Sun_groups(property, ax_coll,cbar_lo,cbar_hi)
    else:
        group_data.load_and_plot_Sun2009(property, ax_coll,cbar_lo,cbar_hi)

#Sun_lR_coll = []
#Sun_ly_coll = []
#Sun_lM500_coll = []
#if((property=="entropy_norm") | (property=="density") | (property=="temperature") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
#    Sungroups_list = home_name + "Sun2009/Sun_groups.cat.sort"

#    Sungroups = open(Sungroups_list,"r").readlines()
#    nSungroups = len(Sungroups)

#    Sun_label = 0 

#    for i in range(nSungroups):
#        Sun_ID = Sungroups[i].split()[0]
#        Sun_T500 = float(Sungroups[i].split()[1])
#        Sun_R500 = float(Sungroups[i].split()[2])
#        Sun_M500_13 = float(Sungroups[i].split()[3])
#        if(Sun_M500_13<0):
#            Sun_M500_13 = 3e+13*(Sun_T500/1.0)**(1.5)/1e+13

#        Sun_lM500 = np.log10(Sun_M500_13)+13
#        if(property=="entropy_norm"):
#            lR,ly = loadSun2009(Sun_ID,"K")
#        if(property=="density"):
#            lR,ly = loadSun2009(Sun_ID,"den")
#            lR = lR + np.log10(0.471) # In R2500, got to get to R500
#        if(property=="temperature"):
#            lR,y = loadSun2009(Sun_ID,"T_3d")
#            ly = np.log10(y) + np.log10(Sun_T500)
#        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
#            lR,ly = loadSun2009(Sun_ID,"den")
#            lR = lR + np.log10(0.471) # In R2500, got to get to R500
#            R = 10**lR*Sun_R500
#            if(property=="m100kpc"):
#                m100kpc = IntegrateDensity(R,ly,100)/gperMsol
#                ly = np.log10(m100kpc)
#            if(property=="m200kpc"):
#                m200kpc = IntegrateDensity(R,ly,200)/gperMsol
#                ly = np.log10(m200kpc)
#            if(property=="m400kpc"):
#                m400kpc = IntegrateDensity(R,ly,400)/gperMsol
#                ly = np.log10(m400kpc)
            
#        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
#            ax_all.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",zorder=0)
#            if(i==0):
#                ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,zorder=0,label="Sun+ (2009)")
#            else:
#                ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,zorder=0)

#        else:
#            Sun_lR_coll.extend(lR)
#            Sun_ly_coll.extend(ly)
#            Sun_lM500_coll.extend([Sun_lM500])
#            if(plot_individual_observations):
#                 if(Sun_label == 0):
#                    if((Sun_lM500 > lhalo_lo) & (Sun_lM500 < lhalo_hi)):
#                        #if(multi == 2 ):  # BDO 5/13/21, trying to replace inferno colorbar.   
#                        ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0,label="Sun+ (2009)")
#                        #else:
#                        #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0,label="Sun+ (2009)")
#                        Sun_label = 1
#                else:
#                    if((Sun_lM500 > lhalo_lo) & (Sun_lM500 < lhalo_hi)):
#                        #if(multi == 2):# BDO 5/13/21, trying to replace viridis colorbar.
#                            ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0)
#                        #else:
#                        #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0)
#                            
#            else:
#                ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0)

#    if(plot_individual_observations==0):
        
#        Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi=flexbin_med_and_onesigmaspread(np.asarray(Sun_lR_coll),np.asarray(Sun_ly_coll),15)
#        print("Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,Median_M500=",Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,np.median(np.asarray(Sun_lM500_coll)))
#        Sun_fractional_color = (np.median(np.asarray(Sun_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
#        ax_coll.plot(Sun_lR,Sun_ly,color=cm(Sun_fractional_color),lw=4,ls=":",zorder=0,label="Sun+ (2009)")
#        for si in range(len(Sun_lR)):
#            print("Sun2009:", Sun_lR[si], Sun_ly[si], Sun_lylo[si], Sun_lyhi[si])
#        #ax_coll.plot(Sun_lR,Sun_lylo, color=cm(Sun_fractional_color),ls=":",zorder=0)
#        #ax_coll.plot(Sun_lR,Sun_lyhi, color=cm(Sun_fractional_color),ls=":",zorder=0)
#        #ax_coll.fill_between(Sun_lR,Sun_lylo,Sun_lyhi, facecolor=cm(Sun_fractional_color),alpha=0.2,zorder=0)
        

if(lhalo_hi>14.5):
    
    CLoGS_lR_coll = []
    CLoGS_ly_coll = []
    CLoGS_lM500_coll = []
    if((property=="entropy_norm")):
        CLoGS_list = home_name + "CLoGS_entropy_profiles/CLoGS_Table4.dat"

        CLoGS = open(CLoGS_list,"r").readlines()
        nCLoGS = len(CLoGS)

        CLoGS_label = 0

        for i in range(nCLoGS):
            CLoGS_ID = CLoGS[i].split()[0]
            CLoGS_M500 = float(CLoGS[i].split()[2])*1e+13
            if(property=="entropy_norm"):
                CLoGS_entropy_exists = int(CLoGS[i].split()[3])
                if(CLoGS_entropy_exists):
                    R,y = loadCLoGS_entropy(CLoGS_ID)
                    CLoGS_redshift = 0.0
                    CLoGS_lM500 = np.log10(CLoGS_M500)
                    CLoGS_rhocrit = 8.52e-30*(1+CLoGS_redshift)**3 # g/cm^3  
                    CLoGS_ne500 = 500*f_b*CLoGS_rhocrit/(mu_e*mp) # cm^-3
                    CLoGS_rho500 = 500*CLoGS_rhocrit  # g/cm^3
                    
                    CLoGS_R500 = (CLoGS_M500*1.989e+33/CLoGS_rho500*(3/(4*np.pi)))**0.3333   # cm
                    CLoGS_ne500 = 500*f_b*CLoGS_rhocrit/(mu_e*mp) # cm^-3
                    
                    CLoGS_T500 = (G*CLoGS_M500*1.989e+33*mu*mp/(2*CLoGS_R500*kboltz))
                    CLoGS_K500 = CLoGS_T500*CLoGS_ne500**-0.6666*keVperK
                    lR = np.log10(R/(CLoGS_R500/cmperkpc))
                    ly = np.log10(y/CLoGS_K500)
                    
                    CLoGS_lR_coll.extend(lR)
                    CLoGS_ly_coll.extend(ly)
                    CLoGS_lM500_coll.extend([CLoGS_lM500])
                    if(plot_individual_observations):                
                        if(CLoGS_label == 0):
                            if(CLoGS_lM500 > lhalo_lo):
                                ax_coll.plot(lR,ly,color=cm((CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=5,label="CLoGS Groups")
                                CLoGS_label = 1
                        else:
                            if(CLoGS_lM500 > lhalo_lo):
                                ax_coll.plot(lR,ly,color=cm((CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=5)
                    else:
                        ax_coll.plot(lR,ly,color=cm((CLoGS_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=0)

                else:
                    continue
            
        if(plot_individual_observations==0):
        
            CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi=flexbin_med_and_onesigmaspread(np.asarray(CLoGS_lR_coll),np.asarray(CLoGS_ly_coll),10)
            print("CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi,median_M500=",CLoGS_lR,CLoGS_lRlo,CLoGS_lRhi,CLoGS_ly,CLoGS_lylo,CLoGS_lyhi,np.median(np.asarray(CLoGS_lM500_coll)))
            CLoGS_fractional_color = (np.median(np.asarray(CLoGS_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
            ax_coll.plot(CLoGS_lR,CLoGS_ly,color=cm(CLoGS_fractional_color),lw=4,ls="--",zorder=0,label="CLoGS Groups")
            #ax_coll.plot(CLoGS_lR,CLoGS_lylo, color=cm(CLoGS_fractional_color),ls="--",zorder=0)
            #ax_coll.plot(CLoGS_lR,CLoGS_lyhi, color=cm(CLoGS_fractional_color),ls="--",zorder=0)
            #ax_coll.fill_between(CLoGS_lR,CLoGS_lylo,CLoGS_lyhi, facecolor=cm(CLoGS_fractional_color),alpha=0.2,zorder=0)

if(lhalo_hi>14.5):
    ACCEPT_lR_coll = []
    ACCEPT_ly_coll = []
    ACCEPT_lM500_coll = []
    if((property=="entropy_norm")):
        ACCEPT_list = home_name + "ACCEPT/ACCEPT.main.tab"
        
        ACCEPT = open(ACCEPT_list,"r").readlines()
        nACCEPT = len(ACCEPT)
        
        ACCEPT_label = 0
        ACCEPT_plot_frequency = 5

        for i in range(nACCEPT):
            if(i<2):
                continue # skip header lines.  
            ACCEPT_name = ACCEPT[i].split()[0]
            ACCEPT_redshift = float(ACCEPT[i].split()[3])
            ACCEPT_Tcl = float(ACCEPT[i].split()[7])
            ACCEPT_M500 = 3e+13*(ACCEPT_Tcl/1.0)**(1.5)
            ACCEPT_lM500 = np.log10(ACCEPT_M500)
            if(ACCEPT_lM500<14.0):
                continue
            #print("ACCEPT_name= ",ACCEPT_name)
            if(property=="entropy_norm"):
                R,y = loadACCEPTdata(ACCEPT_name,property)
                ACCEPT_rhocrit = 8.52e-30*(1+ACCEPT_redshift)**3 # g/cm^3  
                ACCEPT_ne500 = 500*f_b*ACCEPT_rhocrit/(mu_e*mp) # cm^-3
                ACCEPT_rho500 = 500*ACCEPT_rhocrit  # g/cm^3
            
                ACCEPT_R500 = (ACCEPT_M500*1.989e+33/ACCEPT_rho500*(3/(4*np.pi)))**0.3333   # cm
                ACCEPT_ne500 = 500*f_b*ACCEPT_rhocrit/(mu_e*mp) # cm^-3
                
                ACCEPT_T500 = (G*ACCEPT_M500*1.989e+33*mu*mp/(2*ACCEPT_R500*kboltz))
                ACCEPT_K500 = ACCEPT_T500*ACCEPT_ne500**-0.6666*keVperK

                #print("ACCEPT name, M500, R500, K500,color= ",ACCEPT_name,ACCEPT_lM500,ACCEPT_R500/cmperMpc,ACCEPT_K500,(ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo))#,R,y)
                lR = np.log10(R/(ACCEPT_R500/cmperMpc))
                ly = np.log10(y/ACCEPT_K500)

            ACCEPT_lR_coll.extend(lR[np.where(np.isnan(ly)==0)])
            ACCEPT_ly_coll.extend(ly[np.where(np.isnan(ly)==0)])
            ACCEPT_lM500_coll.extend([ACCEPT_lM500])
            if(plot_individual_observations):                
                if(ACCEPT_label == 0):
                    if(ACCEPT_lM500 > lhalo_lo):                    
                        if((i-2)%ACCEPT_plot_frequency==0):
                            ax_coll.plot(lR,ly,color=cm((ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="-.",zorder=0,label="ACCEPT Clusters")
                            ACCEPT_label = 1
                else:
                    if(ACCEPT_lM500 > lhalo_lo):
                        if((i-2)%ACCEPT_plot_frequency==0):
                            ax_coll.plot(lR,ly,color=cm((ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="-.",zorder=0)
            else: # Plot a selection of ACCEPT clusters.  
                if((i-2)%ACCEPT_plot_frequency==0):
                    ax_coll.plot(lR,ly,color=cm((ACCEPT_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="-.",zorder=0)
                    
        if(plot_individual_observations==0):
            ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi=flexbin_med_and_onesigmaspread(np.asarray(ACCEPT_lR_coll),np.asarray(ACCEPT_ly_coll),50)
            ACCEPT_fractional_color = (np.median(np.asarray(ACCEPT_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
            print("ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi,median_mass=",ACCEPT_lR,ACCEPT_lRlo,ACCEPT_lRhi,ACCEPT_ly,ACCEPT_lylo,ACCEPT_lyhi,np.median(np.asarray(ACCEPT_lM500_coll)))
            ax_coll.plot(ACCEPT_lR,ACCEPT_ly,color=cm(ACCEPT_fractional_color),lw=4,ls="-.",zorder=10,label="ACCEPT Clusters")
            #ax_coll.fill_between(ACCEPT_lR,ACCEPT_lylo,ACCEPT_lyhi, facecolor=cm(ACCEPT_fractional_color),alpha=0.2,zorder=10)

                    
if(property=="density"):
    if(plot_only_Lovisari_Sun_groups is not True):
        R,y,ylo,yhi = loadLovisari2015dens()
        lR = np.log10(R)
        ly_ne = np.log10(y)
        lylo_ne = np.log10(ylo)
        lyhi_ne = np.log10(yhi)
        ax_coll.plot(lR,ly_ne,color='black',lw=2,ls="-.",zorder=1,label="Lovisari+ (2015)")
        ax_coll.plot(lR,lylo_ne,color='black',lw=1,ls="-.",zorder=1)
        ax_coll.plot(lR,lyhi_ne,color='black',lw=1,ls="-.",zorder=1)

if(property=="pressure_norm"):
    lR,ly_P,lylo_P,lyhi_P = loadSun2011pressure()
#    lR = np.log10(R)
#    ly_P = np.log10(y)
#    lylo_P = np.log10(ylo)
#    lyhi_P = np.log10(yhi)
    ax_coll.plot(lR,ly_P,color='black',lw=2,ls="-.",zorder=1,label="Sun+ (2011)")
    ax_coll.plot(lR,lylo_P,color='black',lw=1,ls="-.",zorder=1)
    ax_coll.plot(lR,lyhi_P,color='black',lw=1,ls="-.",zorder=1)

if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
    L15_M500, L15_Mgas100kpc, L15_Mgas200kpc, L15_Mgas400kpc = np.loadtxt("Mg_LL15.dat",usecols=(0,1,2,3),unpack=True)
    if(property=="m100kpc"):
        indexes = np.where(L15_Mgas100kpc>0)
        L15_lM500 = np.log10(L15_M500[indexes])
        ly = np.log10(L15_Mgas100kpc[indexes])
    if(property=="m200kpc"):
        indexes = np.where(L15_Mgas200kpc>0)
        L15_lM500 = np.log10(L15_M500[indexes])
        ly =  np.log10(L15_Mgas200kpc[indexes])
    if(property=="m400kpc"):
        indexes = np.where(L15_Mgas400kpc>0)
        L15_lM500 = np.log10(L15_M500[indexes])
        ly = np.log10(L15_Mgas400kpc[indexes])
    
    print("Lov15= ", L15_lM500, ly)
    ax_all.plot(L15_lM500,ly,color="black",marker="+",zorder=0)
    ax_coll.plot(L15_lM500,ly,color="black",marker="+",ls="",ms=12,zorder=0,label="Lovisari+ (2015)")
    

if(property=="metallicity"):
    R500_low_Lov, R500_hi_Lov, Arelax_Lov, sigrelax_Lov, Adist_Lov, sigdist_Lov = group_data.load_Lovisari2019()

    lR_500_Lov = np.log10((R500_low_Lov+R500_hi_Lov)/2.)
    lArelax_Lov = np.log10(Arelax_Lov)
    lAdist_Lov = np.log10(Adist_Lov)

    print('Lov= ', lR_500_Lov, lArelax_Lov, lAdist_Lov)
    ax_coll.plot(lR_500_Lov,lArelax_Lov,color='black',marker='o',ms=5,ls='-',label="Lovisari+ (2019) relaxed",zorder=100000)
    ax_coll.plot(lR_500_Lov,lAdist_Lov,color='black',marker='s',ms=5,ls=':',label="Lovisari+ (2019) disturbed",zorder=100000)

    lR_500, Z_3d, Z_3d_upper, Z_3d_lower = group_data.load_Sun2009_1550_Z()

    ax_coll.errorbar(lR_500, np.log10(Z_3d),yerr=[np.log10(Z_3d_lower+Z_3d)-np.log10(Z_3d),-np.log10(Z_3d)+np.log10(Z_3d+Z_3d_upper)],color=cm((np.log10(3.18e+13)-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=10000,label="Sun+ (2009) NGC 1550 ")


if(property=="XraySB_XMMcounts"):
    group_data.load_and_plot_Lovisari_Sun_groups(property, ax_coll,cbar_lo,cbar_hi)
    
if((property=="LX")):
    M500,M500_err,LX,LX_err = loadAndersonWang()
    ax_coll.plot(M500,LX,marker="o",ls='',color='k',label="Wang+ (2016)")
    ax_coll.errorbar(M500,LX,xerr=M500_err,yerr=LX_err,ls='',color='k')
    
if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):

    if(property=="fgas"):
        M500_Xray, fgas500 = np.loadtxt("McCarthy17_fgas.txt", usecols=(0,1), unpack=True)

        ax_all.plot(M500_Xray, fgas500,color="black",marker="x",ls='',label="Observations")
        ax_coll.plot(M500_Xray, fgas500,color="black",marker="x",ls='',label="Observations")

    ax_all.set_title(title_label,fontsize=18)
    ax_all.tick_params(axis='both', which='major', labelsize=12)
    ax_all.set_xlim(lMlo,lMhi)
    ax_all.set_ylim(lylo,lyhi)
    ax_all.set_xlabel("log $M_{500}$ [$M_{\odot}$]",fontsize=16)
    ax_all.set_ylabel(ylabel,fontsize=16)
    if(multi==0):
        fig_all.subplots_adjust(left=leftborder, bottom=0.17,top=0.90,right=0.98)
        fig_all.savefig('profile_all_' + file_base + '.' + property  + mhalo_str + '.png')
else:    
    #ax_all.set_xscale('log')
    #ax_all.set_yscale('log')
    ax_all.set_title(title_label,fontsize=18)
    ax_all.tick_params(axis='both', which='major', labelsize=12)

    ax_all.set_xlim(lRlo,lRhi)
    if(property=="dispersionmeasure"):
        ax_all.set_ylim(ylo,yhi)
        ax_all.set_xlabel("log $R$ [kpc]",fontsize=16)
    else:
        ax_all.set_ylim(lylo,lyhi)
        ax_all.set_xlabel("log $R/R_{500}$",fontsize=16)
    ax_all.set_ylabel(ylabel,fontsize=16)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=cbar_lo, vmax=cbar_hi))
    sm._A = []
    div = axgrid.make_axes_locatable(ax_all)
    cbar = fig_all.colorbar(sm)
    cbar.ax.set_xlabel(cbarlabel,labelpad=8,fontsize=12)

    #if(multi==0):
    fig_all.subplots_adjust(left=leftborder, bottom=0.17,top=0.90,right=0.98)
    fig_all.savefig('profile_all_' + file_base + '.' + property  + mhalo_str + '.png')



#ax_coll.set_xscale('log')
#ax_coll.set_yscale('log')
if(multi < 2):
    ax_coll.set_title(title_label,fontsize=18)
    
if(legend>0):
    if(multi==2):
        ax_coll.legend(loc="best",fontsize=12)
    else:
        ax_coll.legend(loc="best",fontsize=12)
    ax_coll.tick_params(axis='both', which='major', labelsize=12)

if(property=="dispersionmeasure"):
    ax_coll.set_ylim(ylo,yhi)
else:
    ax_coll.set_ylim(lylo,lyhi)

if((property=="fgas") | (property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc") | (property=="LX")):
    plt.locator_params(axis="x", nbins=5)
    ax_coll.set_xlim(lMlo,lMhi)
    ax_coll.set_xlabel("log $M_{500}$ [$M_{\odot}$]",fontsize=16)
else:
    ax_coll.set_xlim(lRlo,lRhi)
    if(property=="dispersionmeasure"):
        ax_coll.set_xlabel("log $R$ [kpc]",fontsize=16)
    else:
        ax_coll.set_xlabel("log $R/R_{500}$",fontsize=16)


    
if((property=="temperature") | ((property=="dispersionmeasure") & (multi!=1)) | ((property=="entropy_norm") & (multi==2))):
    sm = plt.cm.ScalarMappable(cmap=cm, norm=mpl.colors.Normalize(vmin=cbar_lo, vmax=cbar_hi))
    sm._A = []
    div = axgrid.make_axes_locatable(ax_coll)
    cbar = fig_coll.colorbar(sm)
    cbar.ax.set_xlabel(cbarlabel,labelpad=8,fontsize=12)

ax_coll.set_ylabel(ylabel,fontsize=16)
fig_coll.subplots_adjust(left=leftborder, bottom=0.12,top=topborder,right=rightborder)
fig_coll.savefig('profile_coll_' + file_base + '.' + property  + mhalo_str + '.png')

