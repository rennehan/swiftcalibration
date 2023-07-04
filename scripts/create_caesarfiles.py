USE_VERSION = 'v0.2b'
import os
import sys
import yt
import caesar
import argparse

##########################################################################################################################
#### This function creates ceasar files from the given simulation in a given directory.
#### Call it using python create_caesarfiles.py SNAPDIR MODEL1 MODEL2...
#### Based on Renier Hough's "Create_caesarfiles.py"
##########################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('myArgs', nargs='*')
args = parser.parse_args()
SNAPDIR = args.myArgs[0]
SIM = []
for i in range(1,len(args.myArgs)):
   SIM.append(args.myArgs[i])

#SNAPLIST = [151]
#SNAPLIST = [62,76,105,119]
#SNAPLIST = [36,51,62,105]
#SNAPLIST = [36,51,62,73,78,105,151]
#SNAPLIST = [71,75,83,90,95,100,110,120,130,140]
#SNAPLIST = [36,51,62,78,105,151,20,25,30,40,45,56,66,71,75,83,90,95,100,110,120,130,140]
#SNAPLIST = list(range(0,16))

#SNAPLIST = list(range(0,16))
SNAPLIST = [15]

#SNAPLIST = list(range(0,170))
print (SNAPLIST)
#SNAPLIST = [0,1]
#SIM = 'Kobayashi_NoDust_N128L12_test'

def reduce(SNAPLIST, SNAPDIR, SIM='simba_s12.5n128', FOF6DLOC='Fof6D', NPROC=16):
#  if SIM[0]=="M":
#      SNAPDIR ='/scratch/b/babul/rennehan/cosmo/256/25Mpc/Simba/'+str(SIM)+'/'
#  else:
#      SNAPDIR ='/scratch/b/babul/renierht/output_directory/'+str(SIM)+'/' 
  print (SNAPDIR)
  '''
  SNAPLIST: list of the snapshot numbers
  SIM: Sim name, e.g. m50n512, used in snapshot name (assumes Simba naming convention)
  SNAPDIR = Directory where snapshots are located
  CAESARLOC = Subdirectory to write Caesar catalogs
  FOF6DLOC = Subdirectory to look for / write fof6d files
  '''
  for j in SNAPLIST:
    print (j)
    # path to the file
    SNAP  = os.path.join(SNAPDIR, '%s_%04d.hdf5' % (SIM,j))
    if not os.path.exists(SNAP):
      print(SNAP, "does not exist")
      continue
    
    CAESARFILE = os.path.join(SNAPDIR, '%s_caesar_%04d.hdf5' % (SIM,j))
    #CAESARFILE = '/scratch/b/babul/renierht/output_directory/Kobayashi_NoDust_N256L25_Factor/caesar_%03d.hdf5' % (j)
    FOF6DFILE = os.path.join(SNAPDIR, 'Fof6D/fof_%04d.hdf5'%(j))
    HALOIDFILE = os.path.join(SNAPDIR, FOF6DLOC,  'haloid_%04d.hdf5'%(j))
    HALOID = 'fof' # options atm are 'snap' (Halo ID's from snapshot) or 'fof' (yt's 3DFOF)
    SELFSHIELD  = False # True=compute self-shielded mass; False= use fHI,fH2 from snap
    if 'fh_qr' in SNAPDIR: # Mufasa snaps have no HaloId, need self-shielding correction
      HALOID = 'fof'
      SELFSHIELD = True
    ds = yt.load(SNAP)
    print('Loading caesar object')
    obj = caesar.CAESAR(ds)
    redshift = obj.simulation.redshift
    if USE_VERSION == 'v0.2b':
      #obj.member_search(haloid=HALOID,fof6d_file=FOF6DFILE,nproc=NPROC)
      obj.member_search(haloid=HALOID, blackholes=True, nproc=NPROC)
      #obj.member_search(haloid=HALOID,fsps_bands='uvoir',ssp_model='FSPS',ext_law='composite',fof6d_file=FOF6DFILE,nproc=NPROC)
    elif USE_VERSION == 'v0.1':
      if not os.path.exists(FOF6DFILE): # if no fof6d file, run fof6d and create file
        print('Using caesar v0.1, running fof6d')
        obj.member_search(haloid=HALOID,fof_from_snap=1,fof6d=True,fof6d_outfile=FOF6DFILE,blackholes=True,nproc=NPROC)
      else: # use existing fof6d file
        print('Using caesar_v0.1, inputting fof6d file')
        obj.member_search(fof_from_snap=1,fof6d=True,fof6d_file=FOF6DFILE,blackholes=True,nproc=NPROC,compute_selfshielding=SELFSHIELD,v01_member_search=True)
    else:
      obj.member_search() # just try it and see what happens!
    obj.save(CAESARFILE)

for i in range(0, len(SIM)):
#    reduce(SNAPLIST, SNAPDIR, SIM=SIM[i])    
    reduce(SNAPLIST, os.path.join(SNAPDIR, SIM[i]))

#ssp_table_file='/home/rad/caesar/SSP_Chab_EL.hdf5'
