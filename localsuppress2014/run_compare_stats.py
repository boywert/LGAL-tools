import os
zlist = {}
zlist[0] = "6.00"
zlist[1] = "9.94"
zlist[2] = "9.46"
zlist[3] = "9.03"
zlist[4] = "8.51"
zlist[5] = "7.96"
zlist[6] = "7.48"
zlist[7] = "6.98"
zlist[8] = "6.48"

for i in range(len(zlist)):
    z = zlist[i]
    os.system("python compare_stats.py "+z)
