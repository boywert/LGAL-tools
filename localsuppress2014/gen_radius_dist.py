exe = "/scratch/01937/cs390/Size-Analysis-C2-Ray/ZAHN/zahn_bubbles"
label = ["nore_ori",
         "nore_infall",
         "oka_ori",
         "oka_infall",
         "patchy_ori",
         "patchy_infall"]

#30%
folder = "0.3"
filelist = ["/scratch/01937/cs390/data/CSFR/no_reionization/0/SEMNUM/720.00/xfrac3d_9.308.bin",
            "/scratch/01937/cs390/data/CSFR/no_reionization_infall/SEMNUM/1600.00/xfrac3d_9.308.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto/0/SEMNUM/720.00/xfrac3d_9.308.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto_infall/SEMNUM/1600.00/xfrac3d_9.308.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/720.00/xfrac3d_9.026.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/1600.00/xfrac3d_8.892.bin"]
print "mkdir -p "+folder
for i in range(len(filelist)):
    print exe, filelist[i], folder+"/"+label[i], "306 47.0 0.7 0.9 20"
# 50%
folder = "0.5"
filelist = ["/scratch/01937/cs390/data/CSFR/no_reionization/0/SEMNUM/720.00/xfrac3d_8.397.bin",
            "/scratch/01937/cs390/data/CSFR/no_reionization_infall/SEMNUM/1600.00/xfrac3d_8.397.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto/0/SEMNUM/720.00/xfrac3d_8.397.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto_infall/SEMNUM/1600.00/xfrac3d_8.397.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/720.00/xfrac3d_8.172.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/1600.00/xfrac3d_8.064.bin"]
print "mkdir -p "+folder
for i in range(len(filelist)):
    print exe, filelist[i], folder+"/"+label[i], "306 47.0 0.7 0.9 20"

# 70%
folder = "0.7"
filelist = ["/scratch/01937/cs390/data/CSFR/no_reionization/0/SEMNUM/720.00/xfrac3d_7.960.bin",
            "/scratch/01937/cs390/data/CSFR/no_reionization_infall/SEMNUM/1600.00/xfrac3d_7.960.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto/0/SEMNUM/720.00/xfrac3d_7.960.bin",
            "/scratch/01937/cs390/data/CSFR/okamoto_infall/SEMNUM/1600.00/xfrac3d_7.960.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/720.00/xfrac3d_7.760.bin",
            "/scratch/01937/cs390/Hybrid/xfrac/1600.00/xfrac3d_7.391.bin"]
print "mkdir -p "+folder
for i in range(len(filelist)):
    print exe, filelist[i], folder+"/"+label[i], "306 47.0 0.7 0.9 20"


