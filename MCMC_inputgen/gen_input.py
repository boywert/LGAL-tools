import sys
import subprocess
import random
import os
import numpy

def get_template(filename):
    fp = open(filename, "r")
    content = fp.readlines()
    allvars = {}
    for line in content:
        data = line.strip().split()
        if len(data) > 0:
            if (data[0][0] != "%") & (data[0][0] != "-"):
                allvars[data[0]] = data[1]
    return allvars
def make_unique(a):
    ind = numpy.lexsort(a.T)
    a[numpy.concatenate(([True],numpy.any(a[ind[1:]]!= a[ind[:-1]],axis=1)))]
    return a

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
            var_order.append(data[0])
            if data[5] == '1':
                mcmc_allvars[data[0]] = True
            else:
                mcmc_allvars[data[0]] = False
    p = os.listdir(output_folder)
    sortlist = numpy.array([],dtype=numpy.float32)
    for file in p:
        if file.find("senna_gt") > -1:
            print file
            listp = numpy.loadtxt(output_folder+"/"+file)
            listp = make_unique(listp).sort(axis = 1)[0:n_trials]
            print listp
            #numpy.append(sortlist,listp)
            #sortlist = numpy.unique(sortlist).sort(axis=1)[0:n_trials]
            #print sortlist

    # for key in var_order:
    #     if mcmc_allvars[key]:
            
    return mcmc_allvars

def main(argv):
    if len(argv) < 5:
        print "Usage python "+argv[0]+" <valid_lgal_input_file> <MCMCParameterPriorsAndSwitches.txt> <MCMC_output_folder> <number_of_trials>"
        exit()
    
    get_mcmc_variables(argv[2], argv[3], argv[4])
    template = get_template(argv[1])
    return 0

if __name__ == "__main__":
    main(sys.argv)