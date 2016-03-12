firstfile = 0
lastfile = 7
zlistfile = "/scratch/01937/cs390/data/snap_z2.txt"
config = {}
tau_e = 0.067
delta_tau_e = 0.016
h0 = 0.7
gadgetmass = 1.e10

model_names = ["nore_ori"]

model_labels = ["test"]

model_paths = ["/scratch/01937/cs390/data/outputs/no_reionization/0/"]
struct_file = []
for p in model_paths:
    struct_file.append(p+"/inputs/LGalaxyStruct.py")
struct_file[0] = "/scratch/01937/cs390/data/outputs/okamoto/inputs/LGalaxyStruct.py"
model_plot_colors = ['#1b9e77','#d95f02','#7570b3','#1b9e77','#d95f02','#7570b3']
model_plot_patterns = ['-','-','-','--','--','--']
model_xfrac_path = ["","","/scratch/01937/cs390/Hybrid/xfrac/1350.00","","/scratch/01937/cs390/data/SFR2/okamoto_infall/SEMNUM/1350.00/","/scratch/01937/cs390/data/SFR2/no_reionization_infall/SEMNUM/1350.00/",""]
model_fesc = [720,720,720,1600,1600,1600]
use_model = [True,True,True,True,True,True]
xfrac_doubleflag = [0,0,0,0,0,0]
tau_folder = "/scratch/01937/cs390/tauplot/"
