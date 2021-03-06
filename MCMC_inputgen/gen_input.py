import sys
import subprocess
import random
import os
import numpy
import datetime

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
        if file.find("senna_") > -1:
            print file
            listp = numpy.loadtxt(output_folder+"/"+file)
            listp = make_unique(listp)
            listp = numpy.sort(listp,axis=0)[0:n_trials]
            sortlist = numpy.append(sortlist,listp,axis=0)
            sortlist = make_unique(sortlist)
            sortlist = numpy.sort(sortlist,axis=0)[0:n_trials]

    mcmc_set = []
    lhood = []
    for i in range(len(sortlist)):
        mcmc_set.append({})
        lhood.append(sortlist[i][1])
        for j in range(len(var_order)):
            key = var_order[j]
            mcmc_set[i][key] = 10.**sortlist[i][j+2]
    print lhood
    return mcmc_set,lhood

def gen_input(template,order,mcmc_set,lhood,dest_folder,n_trials):
    for i in range(n_trials):
        fp = open(dest_folder+"/output_"+str(i), "w")
        temp = template.copy()
        temp['OutputDir'] =  temp['OutputDir'].strip()+"/"+str(i)
        os.system("mkdir -p "+temp['OutputDir'])
        print >> fp, "% Sample llhood ="+str(lhood[i])+" generated on "+datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        for key in order:
            if key in mcmc_set[i]:
                print >> fp, key,mcmc_set[i][key]
            else:
                print >> fp, key,temp[key]
        fp.close()
def main(argv):
    if len(argv) < 6:
        print "Usage python "+argv[0]+" <valid_lgal_input_file> <MCMCParameterPriorsAndSwitches.txt> <MCMC_output_folder> <number_of_trials> <dest_folder>"
        exit()
    n_trials =  int(argv[4])
    mcmc_set,lhood = get_mcmc_variables(argv[2], argv[3], int(argv[4]))
    template,order = get_template(argv[1])
    os.system("mkdir -p "+argv[5])
    gen_input(template,order,mcmc_set,lhood,argv[5],n_trials)
    return 0

if __name__ == "__main__":
    main(sys.argv)
