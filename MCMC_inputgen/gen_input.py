import sys
import subprocess
import random
import os
import numpy

def get_template(filename):
    fp = open(filename, "r")
    content = fp.readlines()
    allvars = {}
    vars_order =[]
    for line in content:
        data = line.strip().split()
        if len(data) > 0:
            if (data[0][0] != "%") & (data[0][0] != "-"):
                allvars[data[0]] = data[1]
                vars_order.append(data[0])
    return allvars,vars_order
def make_unique(a):
    ind,indices = numpy.unique(a[:,1],return_index=True)
    return a[indices]
    
def get_mcmc_variables(mcmc_template, output_folder, n_trials):
    fp = open(mcmc_template)
    mcmc_content = fp.readlines()
    fp.close()
    for line in mcmc_content:
        data = line.strip().split()
        if len(data) > 0:
            if (data[0][0] == "%") | (data[0][0] == "#"):
                mcmc_content.remove(line)

    nVars = int(mcmc_content[0].strip())
    mcmc_allvars = {}
    var_order = []
    for i in range(nVars):
        line = mcmc_content[1+i]
        data = line.strip().split()
        if len(data) > 0:
            if data[5] == '1':
                mcmc_allvars[data[0]] = True
                var_order.append(data[0])
    
    p = os.listdir(output_folder)
    sortlist = 1000000*numpy.ones(shape=(2,len(var_order)+2),dtype=numpy.float64)

    for file in p:
        if file.find("senna_gt_10") > -1:
            print file
            listp = numpy.loadtxt(output_folder+"/"+file)
            listp = make_unique(listp)
            listp = numpy.sort(listp,axis=0)[0:n_trials]
            sortlist = numpy.append(sortlist,listp,axis=0)
            sortlist = make_unique(sortlist)
            sortlist = numpy.sort(sortlist,axis=0)[0:n_trials]

    mcmc_set = []
    for i in range(len(sortlist)):
        mcmc_set.append(mcmc_allvars)
        for j in range(len(var_order)):
            key = var_order[j]
            mcmc_set[i][key] = 10.**sortlist[i][j+2]
        print "SfrBurstEfficiency" , mcmc_set[i]['SfrBurstEfficiency']
    return mcmc_set

def gen_input(template,order,mcmc_set,dest_folder,n_trials):
    for i in range(n_trials):
        print "======================================================"
        temp = template.copy()
        for key in order:
            if key in mcmc_set[i]:
                print "haha",i
                print key,mcmc_set[i][key],temp[key]
            else:
                print key,temp[key]    
def main(argv):
    if len(argv) < 5:
        print "Usage python "+argv[0]+" <valid_lgal_input_file> <MCMCParameterPriorsAndSwitches.txt> <MCMC_output_folder> <number_of_trials>"
        exit()
    n_trials =  int(argv[4])
    mcmc_set = get_mcmc_variables(argv[2], argv[3], int(argv[4]))
    template,order = get_template(argv[1])
    #gen_input(template,order,mcmc_set,"",n_trials)
    for i in range(n_trials):
        print mcmc_set[i]['ReheatSlope']
    return 0

if __name__ == "__main__":
    main(sys.argv)
