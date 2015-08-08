firstfile = 0
lastfile = 127
config = {}

h0 = 0.7
gadgetmass = 1.e10

model_names = ["nosup","nosup_nothres","patchyG","patchyG_nothres"]

struct_file = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization_nothres/inputs/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_5/sams/480.00/LGalaxyStruct.py","/mnt/lustre/scratch/cs390/47Mpc/couple/model_5_nothres/sams/480.00/LGalaxyStruct.py"]

model_labels = ["No Suppression","No Suppression nothreshold","Patchy Suppression (Gradual)","Patchy Suppression (Gradual,nothres)"]

model_paths = ["/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization/","/mnt/lustre/scratch/cs390/47Mpc/outputs/no_reionization_nothres/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_5/sams/480.00/","/mnt/lustre/scratch/cs390/47Mpc/couple/model_5_nothres/sams/480.00/"]

model_plot_patterns = ['r--','g--','b--','y--']


use_model = [False,True,False,True]


