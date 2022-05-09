# %%
import os
import sys
import numpy as np
import h5py as h

# %%
if len(sys.argv) != 2:
    exit("You need to only include one command line argument: the HDF5 file name")
filename = sys.argv[1] #"../2022-01-01_Scan/nreo_NR_both_diver_c15_only_Halo_WW/samples/nreo.hdf5"
rawH5 = h.File(filename, "r")
namesplit = filename.split("_")
coupling = int(namesplit[namesplit.index("diver")+1].strip("c"))

# %%
pointID_label = "nreo/pointID"
logL_label = "nreo/LogLike"
dmMass_label = "nreo/#NREO_DiracDM_parameters @NREO_DiracDM::primary_parameters::m"
c0_label = "nreo/#NREO_DiracDM_parameters @NREO_DiracDM::primary_parameters::c0_"+str(coupling)
v0_label = "nreo/#Halo_gNFW_rho0_parameters @Halo_gNFW_rho0::primary_parameters::v0"
vesc_label = "nreo/#Halo_gNFW_rho0_parameters @Halo_gNFW_rho0::primary_parameters::vesc"
rho0_label = "nreo/#Halo_gNFW_rho0_parameters @Halo_gNFW_rho0::primary_parameters::rho0"

pointID = rawH5[pointID_label]
logL = rawH5[logL_label]
dmMass = rawH5[dmMass_label]
c0 = rawH5[c0_label]
v0 = rawH5[v0_label]
vesc = rawH5[vesc_label]
rho0 = rawH5[rho0_label]

# %%
logL_mask = np.array(rawH5[logL_label+"_isvalid"], dtype=int)
# dmMass_mask = np.arramMass_mask, c0_mask))

# %%
pointID_valid = pointID[np.array(logL_mask, bool)]
logL_valid = logL[np.array(logL_mask, bool)]
dmMass_valid = dmMass[np.array(logL_mask, bool)]
c0_valid = c0[np.array(logL_mask, bool)]
v0_valid = v0[np.array(logL_mask, bool)]
vesc_valid = vesc[np.array(logL_mask, bool)]
rho0_valid = rho0[np.array(logL_mask, bool)]

# %%
outputFilename = os.path.splitext(filename)[0]+"-extracted.dat"
with open(outputFilename, mode="w") as file:
    file.write("pointID logL v0 vesc rho0 dmMass c0_i\n")
    for i in range(len(logL)):
        file.write(np.array(pointID_valid).astype(str)[i]+" ")
        file.write(np.array(logL_valid).astype(str)[i]+" ")
        file.write(np.array(v0_valid).astype(str)[i]+" ")
        file.write(np.array(vesc_valid).astype(str)[i]+" ")
        file.write(np.array(rho0_valid).astype(str)[i]+" ")
        file.write(np.array(dmMass_valid).astype(str)[i]+" ")
        file.write(np.array(c0_valid).astype(str)[i]+" ")
        file.write("\n")
