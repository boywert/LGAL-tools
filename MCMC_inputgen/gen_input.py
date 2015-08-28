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
def get_mcmc_variables(mcmc_template, output_folder):
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
    tmpfile = "/tmp/"+str(random.randint(1000000)) 
    exe = "sort -k 2 -n "+output_folder+"/senna_gt_*.txt | head -n 100 >"+tmpfile
    os.system(exe)
    mcmc_trials = open(tmpfile).readlines()
    mcmc_trials = dict(mcmc_trials)
    print mcmc_trials
    # for key in var_order:
    #     if mcmc_allvars[key]:
            
    return mcmc_allvars

def main(argv):
    if len(argv) < 4:
        print "Usage python "+argv[0]+" <valid_lgal_input_file> <MCMCParameterPriorsAndSwitches.txt> <MCMC_output_folder>"
        exit()
    
    get_mcmc_variables(argv[2], argv[3])
    template = get_template(argv[1])
    return 0

if __name__ == "__main__":
    main(sys.argv)
