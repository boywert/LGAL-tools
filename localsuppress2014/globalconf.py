firstfile = 0
lastfile = 127
zlistfile = "/scratch/01937/cs390/data/snap_z2.txt"
config = {}
tau_e = 0.067
delta_tau_e = 0.016
h0 = 0.7
gadgetmass = 1.e10

model_names = ["nore_ori","oka_ori","patchy_ori","nore_infall","oka_infall","patchy_infall"]

model_labels = ["No Suppression, stripping 0","Homogeneous, stripping 0","Patchy suppression, stripping 0", "No Suppression, stripping 1", "Homogeneous, stripping 1", "Patchy Supression, stripping 1"]

model_paths = ["/scratch/01937/cs390/data/outputs/no_reionization/0/","/scratch/01937/cs390/data/outputs/okamoto/0/","/scratch/01937/cs390/Hybrid/sams/720.00/","/scratch/01937/cs390/data/outputs/no_reionization_infall/","/scratch/01937/cs390/data/outputs/okamoto_infall/","/scratch/01937/cs390/Hybrid/sams/1600.00/"]
struct_file = []
for p in model_paths:
    struct_file.append(p+"/inputs/LGalaxyStruct.py")
struct_file[0] = "/scratch/01937/cs390/data/outputs/okamoto/inputs/LGalaxyStruct.py"
struct_file[1] = "/scratch/01937/cs390/data/outputs/okamoto/inputs/LGalaxyStruct.py"
struct_file[2] = "/scratch/01937/cs390/Hybrid/sams/390.00/LGalaxyStruct.py"
struct_file[5] = "/scratch/01937/cs390/Hybrid/sams/1350.00/LGalaxyStruct.py"
struct_file[3] = "/scratch/01937/cs390/data/outputs/no_reionization_infall/inputs/LGalaxyStruct.py"
struct_file[4] = "/scratch/01937/cs390/data/outputs/okamoto_infall/inputs/LGalaxyStruct.py"
#struct_file[6] = "/scratch/01937/cs390/data/outputs/no_reionization/inputs/LGalaxyStruct.py"
model_plot_colors = ['#1b9e77','#d95f02','#7570b3','#1b9e77','#d95f02','#7570b3']
model_plot_patterns = ['-','-','-','--','--','--']
model_xfrac_path = ["","","/scratch/01937/cs390/Hybrid/xfrac/1350.00","","/scratch/01937/cs390/data/SFR2/okamoto_infall/SEMNUM/1350.00/","/scratch/01937/cs390/data/SFR2/no_reionization_infall/SEMNUM/1350.00/",""]
model_fesc = [720,720,720,1600,1600,1600]
use_model = [True,True,True,True,True,True]
xfrac_doubleflag = [0,0,0,0,0,0]
tau_folder = "/scratch/01937/cs390/tauplot/"
