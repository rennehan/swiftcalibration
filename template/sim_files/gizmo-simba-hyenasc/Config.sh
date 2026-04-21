#!/bin/bash            # this line only there to enable syntax highlighting in this file

####################################################################################################
#  Enable/Disable compile-time options as needed: this is where you determine how the code will act
#  From the list below, please activate/deactivate the
#       options that apply to your run. If you modify any of these options,
#       make sure that you recompile the whole code by typing "make clean; make".
#
# This file was originally part of the GADGET3 code developed by
#   Volker Springel (volker.springel@h-its.org). The code has been modified
#   substantially by Phil Hopkins (phopkins@caltech.edu) for GIZMO (to add new modules and clean
#   up the naming conventions and changed many of them to match the new GIZMO conventions)
#
####################################################################################################



####################################################################################################
# --------------------------------------- Boundary Conditions & Dimensions
####################################################################################################
PERIODIC                        # Use this if periodic boundaries are needed (otherwise open boundaries are assumed)
####################################################################################################



####################################################################################################
# --------------------------------------- Hydro solver method
####################################################################################################
HYDRO_MESHLESS_FINITE_MASS      # Lagrangian (constant-mass) finite-volume Godunov method
KERNEL_FUNCTION=3              # Choose the kernel function (2=quadratic peak, 3=cubic spline [default], 4=quartic spline, 5=quintic spline, 6=Wendland C2, 7=Wendland C4)
####################################################################################################



####################################################################################################
# --------------------------------------- Additional Fluid Physics
####################################################################################################
COOLING                        # enables radiative cooling and heating: if GALSF, also external UV background read from file "TREECOOL"
GRACKLE                        # enable GRACKLE: cooling+chemistry package (requires COOLING above; https://grackle.readthedocs.org/en/latest/)
GRACKLE_CHEMISTRY=1            # choose GRACKLE cooling chemistry: (0)=tabular, (1)=Atomic, (2)=(1)+H2+H2I+H2II, (3)=(2)+DI+DII+HD
GRACKLE_SELFSHIELD=3            # Self-shield based on Rahmati+13; see Grackle-3.0 docs for options
####################################################################################################



####################################################################################################
## ------------------------ Gravity & Cosmological Integration Options ---------------------------------
####################################################################################################
# --------------------------------------- TreePM Options (recommended for cosmological sims)
PMGRID=512                     # COSMO enable: resolution of particle-mesh grid		#APB
MULTIPLEDOMAINS=16             # Multi-Domain option for the top-tree level: iso=16,COSMO=64-128
ADAPTIVE_GRAVSOFT_FORALL       # enable adaptive gravitational softening lengths for all particle types
                                # (ADAPTIVE_GRAVSOFT_FORGAS should be disabled). the softening is set to the distance
                                # enclosing a neighbor number set in the parameter file. baryons search for other baryons,
                                # dm for dm, sidm for sidm, etc.

# ------------------------ Zoom Options
PM_PLACEHIGHRESREGION=1+2+16+32 # 2^0+2^1+2^4+2^5 = gas + hi-res dm + stars + black holes
PM_HIRES_REGION_CLIPDM		# split low-res DM particles that enter high-res region (completely surrounded by high-res)
PM_HIRES_REGION_CLIPPING=2500	# for stability: clips particles that escape the hires region in zoom/isolated sims
CHECK_CONTAMINATION
#MAXLEN_OUTPUTLIST=5000
###################################################################################################



####################################################################################################
#------------------ Galaxy formation / Star formation / Supermassive BH Models (with feedback)
####################################################################################################
#---- basic/master switches ---- #
GALSF                          # master switch for galactic star formation model: enables SF, stellar ages, metals, generations, etc.
METALS                         # enable metallicities (with multiple species optional) for gas and stars
GALSF_GENERATIONS=1           # the number of stars a gas particle may spawn (defaults to 1, set otherwise)

##-----------------------------------------------------------------------------------------------------------------------------
#----- old sub-grid models (for large-volume simulations) ---- #
#--------- these are all ultimately variations of the Springel & Hernquist 2005 sub-grid models for the ISM, star formation,
#--------- and stellar winds. their use follows the GADGET-3 use policies. If you are not sure whether you have permission to use them,
#--------- you should contact those authors (Lars Hernquist & Volker Springel)
#------------------------------------------------------------------------------------------------------------------------------
GALSF_JEANS_MIN_T=1           # ISM pressure to resolve M_Jeans with JEANS_MIN_T*N_ngb particles
GALSF_SUBGRID_WINDS            # sub-grid winds ('kicks' as in Oppenheimer+Dave,Springel+Hernquist,Boothe+Schaye,etc)
GALSF_SUBGRID_VARIABLEVELOCITY # winds with velocity scaling based on halo properties (Oppenheimer+Dave); req.GALSF_SUBGRID_WINDS
GALSF_SUBGRID_RECOUPLE=1.0	# recouple when dv<RECOUPLE*csound 
GALSF_SUBGRID_RECOUPLE_STATS	# output recoupling stats in log file
GALSF_SUBGRID_DAA=7.5    		# use eta fit to particle tracking results of DAA+17: -0.35 at M*<5e9, -0.7 above, with eta=9 at 5e9
GALSF_ETA_FUDGE_SLOPE		# Resolution-dependent hi-z suppression of eta: =2 @50/512, =1.5 @25/512, =1 @12.5/512
GALSF_SUBGRID_FIREVEL=1.6     # normalization of velocity scaling; Muratov+15 suggests 0.854
GALSF_SUBGRID_VBOOST=0.25	# boost velocity to this radius (rel. to r200)
GALSF_SUBGRID_HOTWIND=0.3      # fraction of remaining E_SN used to heat wind
GALSF_SUBGRID_HOTWIND_ZSCALE	# add metallicity scaling to E_SN for HOTWIND based on Schaerer 03
GALSF_SUBGRID_HOTWIND_PANDYA	# Pandya+21 hot wind fraction
GALSF_SUBGRID_ISM_RECOUPLE=0.1   # if past this frac of max delay time and in ISM, then immediately recouple
GALSF_SUBGRID_ELIM=1.0		# Limit kinetic wind E to this times E_SN
GALSF_WINDS_RADIAL              # uses POLAR; chooses sign of vXa to be away from gal (req. FOFRAD)

#GALSF_INSTANTANEOUS_METALS     # instantaneously enrich SF particles # APB
#GALSF_SUBGRID_METALLOAD		# winds loaded with extra metals from SNe ejectae. Use only when using INSTANTANEOUS_METALS and AGBFEEDBACK and not when using GALSF_FB_KOBAYASHI # APB
#GALSF_TYPEIA                    # Type Ia enrichment and energy input, for INSTANTANEOUS_METALS # APB
#GALSF_AGBFEEDBACK               # enrichment from AGB stars # APB
#GALSF_AGBWINDHEATING=100        # heating from AGB winds (see Conroy,vanDokkum,Kravtsov 14) # APB

GALSF_FB_KOBAYASHI		# Chemical enrichment of stars using 34 elements (Kobayashi, C. 2007) # APB
GALSF_FB_KOBAYASHI_ALLX  	# ALLX = 34 element if enabled, otherwise its only 3 elements (Z, O, Fe) # APB 
GALSF_FB_KOBAYASHI_AGB  	# APB
#RH_DEBUG_FILE_BAD		# Chem function output.txt file
#RH_DEBUG_LITE
#KroupaIMF                      # Kroupa IMF calculation for the yields in the Chem5 model
ChabrierIMF                     # Chabrier IMF calculation for the yields in the Chem5 model
##GALSF_FB_KOBAYASHI_YIELD_FUDGE=2.0
##RATIO_ENERGY=2
##RH_DEBUG_RATIO_ENERGY


##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#------ PFH physical models for star formation and feedback: these are the FIRE simulation modules (Hopkins et al. 2014) ------ ##
#--------- their use follows the FIRE authorship policy. Any new project using these physics must first be agreed to by all of the
#--------- core development team of the FIRE simulations: P. Hopkins, E. Quataert, D. Keres, C.A. Faucher-Giguere.
#--------- Papers using these modules must offer co-authorship to the members of the FIRE development team.
##-----------------------------------------------------------------------------------------------------
#----------- star formation law ---- #
GALSF_SFR_KMT                   # calc fH2 via the KMT model, form stars if fH2 > 0
GALSF_SFR_KMT_SCALECLUMP=1.0	# scale KMT sub-res clumping factor with this power of epsilon
#----------- dust processes, formation, growth and destruction ---- #
GALSF_DUST                    # enable dust formation and destruction (master switch)
GALSF_DUST_GRAINSIZE=0.1     # grain size [micro m] (to-do: a distribution rather than a fixed value)
GALSF_DUST_PRODUCTION		   # enable dust production via ejecta's condensation
GALSF_DUST_PRODUCTION_BOOST=2.0 # Boost Popping+2017 SN condensation efficiency by this factor
GALSF_DUST_GROWTH			   # enable dust growth
GALSF_DUST_GROWTH_DENSREF=2.3e-24 # [g cm-3] proper
GALSF_DUST_GROWTH_TREF=20.0 # [K]
GALSF_DUST_GROWTH_TAUREF=1.0      # [Gyr] growth-timescale at TREF and DENSREF
GALSF_DUST_DESTRUCTION        # enable dust destruction
GALSF_DUST_DESTRUCTION_EFF=0.3 # efficiency of dust destruction by shock (destroyed mass/shocked mass)
GALSF_DUST_DESTRUCTION_SPUTTERING # dust destruction by (thermal) sputtering


##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#-------------------------------------- SMBH/AGN stuff; also heavily expanded with PFH models
##-----------------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------------
#------ PFH physical models for black hole growth and feedback: these are the FIRE simulation modules, their use follows the same FIRE policy above
#------ The original GADGET-3 BH model (only: BLACK_HOLES,BH_SWALLOWGAS,BH_BONDI,BH_DRAG) follow the GADGET-3 Springel & Hernquist policy above
##-----------------------------------------------------------------------------------------------------
BLACK_HOLES                    # enables Black-Holes (master switch)
BH_HOST_TO_SEED_RATIO=30000     # DAA: The minimum stellar mass for seeding is BH_HOST_TO_SEED_RATIO * All.SeedBlackHoleMass
BH_SEED_FROM_STAR_PARTICLE
BH_FOFRAD			# Use FOFRAD galaxies for seeding and repositioning
BH_SWALLOWGAS                  # enables stochastic accretion of gas particles consistent with growth rate of hole
BH_GRAVACCRETION=0               # Gravitational instability accretion estimator from Hopkins & Quataert 2010
BH_BAL_KICK                    # do BAL winds with stochastic particle kicks at specified velocity (instead of continuous wind solution - requires BH_SWALLOWGAS - )
BH_BAL_KICK_COLLIMATED         # DAA: winds follow the direction of angular momentum within Kernel (only for BH_BAL_KICK winds)
BH_BAL_KICK_MOMENTUM_FLUX=20.0   # DAA: increase the effective mass-loading of BAL winds to reach the desired momentum flux in units of L_bol/c (needs BH_BAL_KICK)
BH_OUTPUT_MOREINFO             # DAA: output additional info to "blackhole_details"
BH_QUENCH_JET			# increases BH kick velocity in quenched halos		#APB: Flag that activates jet
BH_XRAY_FEEDBACK		# Adds X-ray heating following Choi+11, for v_kick>this value	#APB: Flag that activates X-ray feedback
BH_QUENCH_JET_ACTIVATION_ASCALE # Scale jet activation BH masses with expansion factor (APB)
BH_QUENCH_JET_HOTWIND=2000  	# Heat jet particles to Tvir above this wind speed
BH_BONDI_HOT=5                  # log10(Tgas) above which to allow Bondi accretion
BH_SFWIND_SUPPRESS=3.e6		# suppress BH accretion exponentially below this BH mass
BH_JET_MHALO_SCALE=1.e8		# scales max jet velocity with MBH^1/3; velocity is BH_QUENCH_JET at this MBH


##-----------------------------------------------------------------------------------------------------
#--------------------------------------- Smagorinsky Turbulent Eddy Diffusion Model
#---------------------------------------- Users of these modules should cite Hopkins et al. 2017 (arXiv:1702.06148) and Colbrook et al. (arXiv:1610.06590)
#METALS                         # enable metallicities (with multiple species optional) for gas and stars [must be included in ICs or injected via dynamical feedback; needed for some routines]
#TURB_DIFF_METALS               # turbulent diffusion of metals (passive scalars); requires METALS
#TURB_DIFF_ENERGY               # turbulent diffusion of internal energy (conduction with effective turbulent coefficients)
#TURB_DIFF_VELOCITY             # turbulent diffusion of momentum (viscosity with effective turbulent coefficients)
#TURB_DIFF_DYNAMIC              # replace Smagorinsky-style eddy diffusion with the 'dynamic localized Smagorinsky' model from Rennehan et al. (arXiv:1807.11509 and 2104.07673): cite those papers for all methods. more accurate but more complex and expensive.
##-----------------------------------------------------------------------------------------------------


##----------------------------------------------------------------------------------------------------
##--------- Turbulent diffusion/mixing models from Rennehan et al. (2019) ----------------------------
##----------------------------------------------------------------------------------------------------
### Following 6 options activate only metal diffusion
#OUTPUT_VEL_TENSOR
#TURB_DIFF_GRADIENT
#TURB_DIFF_METALS
#TURB_DIFF_GRADIENT_Z
#TURB_DIFF_GRADIENT_HARMONIC
#TURB_DIFF_USE_HSML
#TURB_DIFF_GRADIENT_VISC  # Also activates viscosity (momentum diffusion)


##----------------------------------------------------------------------------------------------------
##--------- Turbulent diffusion/mixing models from Rennehan (2021) -----------------------------------
##----------------------------------------------------------------------------------------------------
#OUTPUT_VEL_TENSOR
#TURB_DIFF_GRADIENT_VISC         # Momentum diffusion (eddy viscosity)
#TURB_DIFF_METALS		# metal diffusion (passive scalar mixing)
#TURB_DIFF_GRADIENT_Z           	# metal diffusion (passive scalar mixing)
#TURB_DIFF_GRADIENT_HARMONIC
#TURB_DIFF_GRADIENT_BALARAC13    # Momentum diffusion (eddy viscosity)
#TURB_DIFF_USE_HSML
#TURB_DIFF_LOW_MACH_FIX
####################################################################################################


####################################################################################################
# --------------------------------------- Output/Input options
####################################################################################################
HAVE_HDF5						# needed when HDF5 I/O support is desired
OUTPUTPOTENTIAL                # forces code to compute+output potentials in snapshots
####################################################################################################



####################################################################################################
# -------------------------------------------- De-Bugging & special (usually test-problem only) behaviors
####################################################################################################
#ID_TRACK_GENERATIONS		# This option is broken, so remove
NO_ISEND_IRECV_IN_DOMAIN
#USE_MPI_IN_PLACE               # MPI debugging: makes AllGatherV compatible with MPI_IN_PLACE definitions in some MPI libraries
DOUBLEPRECISION_FFTW           # FFTW in double precision to match libraries
USE_FFTW3
EVALPOTENTIAL                  # computes gravitational potential
####################################################################################################



####################################################################################################
#---------------------------------------- On the fly FOF groupfinder
#------------------ This is originally developed as part of GADGET-3 (SUBFIND) by V. Springel
#------------------ Use of these modules follows the GADGET-3 policies described above
####################################################################################################
FOF                                # enable FoF output
FOF_PRIMARY_LINK_TYPES=2           # 2^type for the primary dark matter type
FOF_SECONDARY_LINK_TYPES=1+4+8+16+32   # 2^type for the types linked to nearest primaries (gas, low-res dm, stars, bhs)
FOF_GROUP_MIN_LEN=16               # default is 32
FOFRAD=0.0056                      # fast, approximate galaxy (dense gas+stars) finder for winds
FOF_SCALEDEPENDENT		   # sets scale-dep intervals for redoing FOFs
FOF_SAVE_GROUPID		   # save IDs of particles in FOF groups (as opposed to tab files)
LINKLENGTH=0.2                    # Linkinglength for FoF (default=0.2)
###################################################################################################
