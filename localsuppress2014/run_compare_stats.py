import os
zlist = {}
zlist[0] = "6.00"

for i in range(len(zlist)):
    z = zlist[i]
    os.system("python compare_stats.py "+z)
