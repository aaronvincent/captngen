# %%
import sys
import numpy as np
import h5py as h

# %%
if len(sys.argv) != 3:
    exit("You need to include exactly two command line arguments: the data file name, and the HDF5 file name")
datafile = sys.argv[1]#"nreo_NR_both_diver_c4_only_Halo_WW.dat"
H5file = sys.argv[2]#"../2022-01-01_REScan/nreo_NR_both_diver_c4_only_Halo_WW/samples/nreo.hdf5"

# %%
totalloglike = "LogLike"
subloglike = ['#CDMSlite_LogLikelihood @DarkBit::CDMSlite_GetLogLikelihood', '#CRESST_II_LogLikelihood @DarkBit::CRESST_II_GetLogLikelihood', '#DarkSide_50_LogLikelihood @DarkBit::DarkSide_50_GetLogLikelihood', '#LUX_2016_LogLikelihood @DarkBit::LUX_2016_GetLogLikelihood', '#PICO_60_2017_LogLikelihood @DarkBit::PICO_60_2017_GetLogLikelihood', '#PandaX_2016_LogLikelihood @DarkBit::PandaX_2016_GetLogLikelihood', '#PandaX_2017_LogLikelihood @DarkBit::PandaX_2017_GetLogLikelihood', '#XENON1T_2018_LogLikelihood @DarkBit::XENON1T_2018_GetLogLikelihood', '#IC79_loglike @DarkBit::IC79_loglike', '#lnL_rho0 @DarkBit::lnL_rho0_lognormal', '#lnL_v0 @DarkBit::lnL_v0_gaussian', '#lnL_vesc @DarkBit::lnL_vesc_gaussian']

# %%
pointID = np.loadtxt(datafile, dtype=int, usecols=0, unpack=True)
logL, captures, newANTARES, newIceCube = np.loadtxt(datafile, dtype=float, usecols=(1,7,9,10), unpack=True)
appendingH5 = h.File(H5file, "a")
points = appendingH5["nreo/pointID"]
length = len(points)

completeANTARES = np.zeros(length, dtype=float)
completeIceCube = np.zeros(length, dtype=float)
modifiedLike = np.zeros(length, dtype=float)
ANTARES_isvalid = np.zeros(length, dtype=int)
IceCube_isvalid = np.zeros(length, dtype=int)
modifiedLike_isvalid = np.zeros(length, dtype=int)
recalcCaps = np.zeros(length, dtype=float)

for i in range(len(points)):
    point = i+1
    if point in pointID:
        indexPoint = np.where(pointID == point)[0][0]
        completeANTARES[i] = newANTARES[indexPoint]
        ANTARES_isvalid[i] = 1
        completeIceCube[i] = newIceCube[indexPoint]
        IceCube_isvalid[i] = 1
        modifiedLike[i] = logL[indexPoint] + newANTARES[indexPoint] + newIceCube[indexPoint]
        modifiedLike_isvalid[i] = 1
        recalcCaps[i] = captures[indexPoint]

# %%
appendingH5.create_dataset(name="nreo/ANTARES_new_LogLike", data=completeANTARES)
appendingH5.create_dataset(name="nreo/ANTARES_new_LogLike_isvalid", data=ANTARES_isvalid)
appendingH5.create_dataset(name="nreo/IceCube_new_LogLike", data=completeIceCube)
appendingH5.create_dataset(name="nreo/IceCube_new_LogLike_isvalid", data=IceCube_isvalid)
appendingH5.create_dataset(name="nreo/Modified_LogLike", data=modifiedLike)
appendingH5.create_dataset(name="nreo/Modified_LogLike_isvalid", data=modifiedLike_isvalid)
appendingH5.create_dataset(name="nreo/Recalculated_Capturerate", data=recalcCaps)
appendingH5.close()
