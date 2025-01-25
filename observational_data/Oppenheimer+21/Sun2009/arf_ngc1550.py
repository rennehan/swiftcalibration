import yt
import numpy as np
import pyxsim
import soxs
import matplotlib.pyplot as plt


arf_fi = soxs.AuxiliaryResponseFile("fi_aim.arf")
ebins_fi = np.linspace(arf_fi.elo[0], arf_fi.ehi[-1], arf_fi.elo.size)#+1)

arf_bi = soxs.AuxiliaryResponseFile("bi_aim.arf")
ebins_bi = np.linspace(arf_bi.elo[0], arf_bi.ehi[-1], arf_bi.elo.size)#+1)

arf_cy0 = soxs.AuxiliaryResponseFile.from_instrument(f"chandra_acisi_cy0")
ebins_cy0 = np.linspace(arf_cy0.elo[0], arf_cy0.ehi[-1], arf_cy0.elo.size)#+1)


print("ebins_bi= ", ebins_bi, len(ebins_bi))
#counts, _ = np.histogram(e, ebins)
#print("arf.eff_area= ", arf.eff_area, len(arf.eff_area))

fig,ax = plt.subplots(ncols=1, nrows=1, figsize=(10,6))
ax.plot(ebins_fi, arf_fi.eff_area, label="arf_fi")
ax.plot(ebins_bi, arf_bi.eff_area, label="arf_bi")
ax.plot(ebins_cy0, arf_cy0.eff_area, label="SOXS Cycle 0")
ax.legend()
ax.set_xlabel("E (kev)")
ax.set_ylabel("Area (cm$^{2}$)")
fig.savefig(f"Chandra_arfs.png")
