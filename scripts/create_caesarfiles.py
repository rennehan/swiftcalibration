USE_VERSION = 'v0.2b'
import os
import sys
import yt
import caesar
import argparse

##########################################################################################################################
#### This function creates ceasar files from the given simulation in a given directory.
#### Call it using python create_caesarfiles.py --modeldir=path/to/model --simdir=SIMDIR --sim=SIM
#### Based on Renier Hough's "Create_caesarfiles.py"
##########################################################################################################################

parser = argparse.ArgumentParser()
parser.add_argument('--modeldir')
parser.add_argument('--simdir')
parser.add_argument('--sim', type=str)
args = parser.parse_args()

MODELDIR = args.modeldir
SIMDIR = args.simdir
SIM = args.sim


## Need to set this in here
# SNAPLIST = list(range(0,16))
SNAPLIST = [15]

print (SNAPLIST)


def reduce(SNAPLIST, SNAPDIR, SIM, FOF6DLOC='Fof6D', NPROC=16):
    print (SNAPDIR)
    '''
    SNAPLIST: list of the snapshot numbers
    SIM: Sim name, e.g. m50n512, used in snapshot name (assumes Simba naming convention)
    SNAPDIR = Directory where snapshots are located
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
        FOF6DFILE = os.path.join(SNAPDIR, 'Fof6D/fof_%04d.hdf5'%(j))
        HALOIDFILE = os.path.join(SNAPDIR, FOF6DLOC,  'haloid_%04d.hdf5'%(j))
        HALOID = 'fof' # options atm are 'snap' (Halo ID's from snapshot) or 'fof' (yt's 3DFOF)
        SELFSHIELD  = False # True=compute self-shielded mass; False= use fHI,fH2 from snap
        if 'fh_qr' in SNAPDIR: # Mufasa snaps have no HaloId, need self-shielding correction
            HALOID = 'fof'
            SELFSHIELD = True

        if os.path.exists(CAESARFILE):
            print(CAESARFILE, "exists")
            continue

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
                obj.member_search(fof_from_snap=1,fof6d=True,fof6d_file=FOF6DFILE,blackholes=True,
                                  nproc=NPROC,compute_selfshielding=SELFSHIELD,v01_member_search=True)
        else:
            obj.member_search() # just try it and see what happens!
        obj.save(CAESARFILE)



reduce(SNAPLIST, os.path.join(MODELDIR, SIMDIR), SIM)

#print()
print("DONE")
print()
print()