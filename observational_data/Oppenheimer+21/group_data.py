import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import mpl_toolkits.axes_grid1 as axgrid
import os.path

#home_name = "/home/ariboppe/groups_review/"
home_name = "/home/b/babul/aspadawe/project/data/groups-and-clusters/Oppenheimer+21/"
cm = plt.get_cmap('plasma') 


def load_Sun2009(group_ID, property_base):
    fname = home_name + "Sun2009/" + property_base + "/" + group_ID + "_" + property_base
    #print("fname= ", fname)
    if(os.path.isfile(fname)):
        R,y = np.loadtxt(fname,usecols=(0,1),unpack=True)
        print("y[-1]= ", y[-1])
    else:
        R = [-99,-99]
        y = [-99,-99]

    return(R,y)

def load_Lovisari2015(obs_Lovisari_ID,property):
    
    fname = home_name + "Lovisari2015/" + obs_Lovisari_ID + "_sb.dat"
    
    if(os.path.isfile(fname)):
        R_kpc,y,yerr = np.loadtxt(fname,usecols=(2,4,5),unpack=True)
    else:
        R_kpc = [-99,-99]
        y = [-99,-99]
        yerr = [0,0]

    return(R_kpc,y,yerr)

def load_Sun2009_1550_Z():

    fname = home_name + "Sun2009/N1550_Z_3d"

    R_500_N1550_kpc = 465 

    R_kpc,R_kpc_err,Z_3d,Z_3d_upper,Z_3d_lower = np.loadtxt(fname,usecols=(0,1,2,3,4),unpack=True)

    lR_500 = np.log10(R_kpc/R_500_N1550_kpc)
    
    return(lR_500,Z_3d,Z_3d_upper,Z_3d_lower)
    

def load_and_plot_Sun2009(property,ax_coll,cbar_lo,cbar_hi):

    Sun_lR_coll = []
    Sun_ly_coll = []
    Sun_lM500_coll = []
    
    Sungroups_list = home_name + "Sun2009/Sun_groups.cat.sort"

    Sungroups = open(Sungroups_list,"r").readlines()
    nSungroups = len(Sungroups)

    Sun_label = 0 

    for i in range(nSungroups):
        Sun_ID = Sungroups[i].split()[0]
        Sun_T500 = float(Sungroups[i].split()[1])
        Sun_R500 = float(Sungroups[i].split()[2])
        Sun_M500_13 = float(Sungroups[i].split()[3])
        if(Sun_M500_13<0):
            Sun_M500_13 = 3e+13*(Sun_T500/1.0)**(1.5)/1e+13

        Sun_lM500 = np.log10(Sun_M500_13)+13
        if(property=="entropy_norm"):
            lR,ly = load_Sun2009(Sun_ID,"K")
        if(property=="density"):
            lR,ly = load_Sun2009(Sun_ID,"den")
            lR = lR + np.log10(0.471) # In R2500, got to get to R500
        if(property=="temperature"):
            lR,y = load_Sun2009(Sun_ID,"T_3d")
            ly = np.log10(y) + np.log10(Sun_T500)
        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
            lR,ly = load_Sun2009(Sun_ID,"den")
            lR = lR + np.log10(0.471) # In R2500, got to get to R500
            R = 10**lR*Sun_R500
            if(property=="m100kpc"):
                m100kpc = IntegrateDensity(R,ly,100)/gperMsol
                ly = np.log10(m100kpc)
            if(property=="m200kpc"):
                m200kpc = IntegrateDensity(R,ly,200)/gperMsol
                ly = np.log10(m200kpc)
            if(property=="m400kpc"):
                m400kpc = IntegrateDensity(R,ly,400)/gperMsol
                ly = np.log10(m400kpc)
            
        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
            ax_all.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",zorder=0)
            if(i==0):
                ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,zorder=0,label="Sun+ (2009)")
            else:
                ax_coll.plot(np.log10(Sun_M500_13*1e+13),ly,color="black",marker="x",ls="",ms=12,zorder=0)

        else:

            Sun_lR_coll.extend(lR)
            Sun_ly_coll.extend(ly)
            Sun_lM500_coll.extend([Sun_lM500])
            if(plot_individual_observations):
                if(Sun_label == 0):
                    if((Sun_lM500 > lhalo_lo) & (Sun_lM500 < lhalo_hi)):
                        #if(multi == 2 ):  # BDO 5/13/21, trying to replace inferno colorbar.   
                        ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0,label="Sun+ (2009)")
                        #else:
                        #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0,label="Sun+ (2009)")
                        Sun_label = 1
                else:
                    if((Sun_lM500 > lhalo_lo) & (Sun_lM500 < lhalo_hi)):
                        #if(multi == 2):# BDO 5/13/21, trying to replace viridis colorbar.
                            ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0)
                        #else:
                        #    ax_coll.plot(lR,ly,color='k',lw=1,ls="--",zorder=0)
                            
            else:
                ax_coll.plot(lR,ly,color=cm((Sun_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls=":",zorder=0)

    if(plot_individual_observations==0):
        
        Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi=flexbin_med_and_onesigmaspread(np.asarray(Sun_lR_coll),np.asarray(Sun_ly_coll),15)
        print("Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,Median_M500=",Sun_lR,Sun_lRlo,Sun_lRhi,Sun_ly,Sun_lylo,Sun_lyhi,np.median(np.asarray(Sun_lM500_coll)))
        Sun_fractional_color = (np.median(np.asarray(Sun_lM500_coll))-cbar_lo)/(cbar_hi-cbar_lo)
        ax_coll.plot(Sun_lR,Sun_ly,color=cm(Sun_fractional_color),lw=4,ls=":",zorder=0,label="Sun+ (2009)")
        for si in range(len(Sun_lR)):
            print("Sun2009:", Sun_lR[si], Sun_ly[si], Sun_lylo[si], Sun_lyhi[si])


def load_and_plot_Lovisari_Sun_groups(property,ax_coll,cbar_lo,cbar_hi):

    obs_lR_coll = []
    obs_ly_coll = []
    obs_lM500_coll = []

    obsgroups_list = home_name + "Lovisari2015/Sun_Lovisari_groups.txt"
    
    obsgroups = open(obsgroups_list,"r").readlines()
    nobsgroups = len(obsgroups)

    obs_label = 0 

    for i in range(nobsgroups):
        obs_Sun_ID = obsgroups[i].split()[0]
        obs_Lovisari_ID = obsgroups[i].split()[1]
        obs_T500 = float(obsgroups[i].split()[2])
        obs_R500 = float(obsgroups[i].split()[3])
        obs_M500_13 = float(obsgroups[i].split()[4])
        if(obs_M500_13<0):
            obs_M500_13 = 3e+13*(obs_T500/1.0)**(1.5)/1e+13

        obs_lM500 = np.log10(obs_M500_13)+13
        
        if(property=="XraySB_XMMcounts"):
            print("obs_Lovisari_ID= ", obs_Lovisari_ID)
            R_kpc,y,yerr = load_Lovisari2015(obs_Lovisari_ID,"XraySB_XMMcounts")
            lR = np.log10(R_kpc/obs_R500)
            ly = np.log10(y)
            lyerr = np.log10(yerr)
        if(property=="entropy_norm"):
            lR,ly = load_Sun2009(obs_Sun_ID,"K")
        if(property=="density"):
            lR,ly = load_Sun2009(obs_Sun_ID,"den")
            lR = lR + np.log10(0.471) # In R2500, got to get to R500
        if(property=="temperature"):
            lR,y = load_Sun2009(obs_Sun_ID,"T_3d")
            ly = np.log10(y) + np.log10(obs_T500)
        if((property=="m100kpc") | (property=="m200kpc") | (property=="m400kpc")):
            lR,ly = load_Sun2009(obs_Sun_ID,"den")
            lR = lR + np.log10(0.471) # In R2500, got to get to R500
            R = 10**lR*obs_R500
            if(property=="m100kpc"):
                m100kpc = IntegrateDensity(R,ly,100)/gperMsol
                ly = np.log10(m100kpc)
            if(property=="m200kpc"):
                m200kpc = IntegrateDensity(R,ly,200)/gperMsol
                ly = np.log10(m200kpc)
            if(property=="m400kpc"):
                m400kpc = IntegrateDensity(R,ly,400)/gperMsol
                ly = np.log10(m400kpc)
            

        obs_lR_coll.extend(lR)
        obs_ly_coll.extend(ly)
        obs_lM500_coll.extend([obs_lM500])

        ax_coll.plot(lR,ly,color=cm((obs_lM500-cbar_lo)/(cbar_hi-cbar_lo)),lw=1,ls="--",zorder=0)
