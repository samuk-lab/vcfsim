import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
import argparse
from IPython.display import SVG, display
from .SimulatorClass import MyVcfSim
import warnings

def multiple_chrom(chromfilename = 'input.txt', seed = 1234, foldername = 'PixyFolder', percentmissing = 0, percentsitemissing = 0, outputfile = 'myvcftest', samp_num = 20, sample_names = None, population_mode = 1, time = 1000, hmm_baseline=None, hmm_multiplier=None, hmm_p_good_to_bad=None, hmm_p_bad_to_good=None):
    
    f = open(chromfilename, "r")

    a = f.readline()
    
    chromlist = []
    ploidylist = []
    lengthlist = []
    nelist = []
    mulist = []

    #check if all lists are equal

    while(a):
        a = a.split()
        if(len(a) != 5):
            print('Error, sample file must only have 5 columns, chromosone, ploidy, sequence length, population size, mutation rate')
            sys.exit(0)
            
        if(not a[1].isdigit()):
            print('Error, ploidy must be a number')
            sys.exit(0)

        try:
            temp = int(eval(a[2]))
        except Exception as e:
            print("Error, sequence length must be an expression or number")
            sys.exit(0)

        try:
            temp = int(eval(a[3]))
        except Exception as e:
            print("Error, effective population must be an expression or number")
            sys.exit(0)
            
        try:
            temp = int(eval(a[4]))
        except Exception as e:
            print("Error, mutation rate must be an expression or number")
            sys.exit(0)
            

        chromlist.append(a[0]), ploidylist.append(int(a[1])), lengthlist.append(int(eval(a[2]))), nelist.append(int(eval(a[3]))),  
        mulist.append(eval(a[4]))
        
        a = f.readline()

    temp = len(chromlist) + len(ploidylist) + len(lengthlist) + len(nelist) + len(mulist)
    if(temp/5 != len(chromlist)):
        print('Error, all columns must be the same length')
        sys.exit(0)

    for i in range(len(chromlist)):
        sim = MyVcfSim(chromlist[i], lengthlist[i], ploidylist[i], nelist[i], mulist[i], percentmissing, percentsitemissing, seed, chromlist[i], samp_num, 'population.txt', 'vcf', sample_names, population_mode, time, hmm_baseline, hmm_multiplier, hmm_p_good_to_bad, hmm_p_bad_to_good)
        
        sim.simulate_vcfs()
    
    with open(outputfile, "w") as f:
        pass

    for i in range(len(chromlist)):
        if i == 0:
            with open(chromlist[i], 'r') as source:
                lines = source.readlines()
            with open(outputfile, 'a') as destination:
                destination.writelines(lines)
        else:
            with open(chromlist[i], 'r') as source:
                lines = source.readlines()
            with open(outputfile, 'a') as destination:
                destination.writelines(lines[6:])
        os.remove(chromlist[i])

def vcf_simulator(chrom = 1, amountofruns = 1, seed = 1234, foldername = 'PixyFolder', sitesize = 10000, ploidy = 2, population = 1700000, mutationrate = 0.0000000055, percentmissing = 0, percentsitemissing = 0, outputfile = 'myvcftest', samp_num = 20, sample_names = None, population_mode = 1, time = 1000, hmm_baseline=None, hmm_multiplier=None, hmm_p_good_to_bad=None, hmm_p_bad_to_good=None):
    
    # when custom sample names are provided we set sample size from the names
    if sample_names is not None:
        try:
            names_len = len(sample_names)
        except Exception as e:
            print('Error, samples must be a list of names')
            sys.exit(0)
        if names_len < 1:
            print('Error, samples list must have at least one name')
            sys.exit(0)
        samp_num = names_len

    for x in range(amountofruns):
        if outputfile is not None:
            outputfilename = outputfile + str(seed) + '.vcf'
        else:
            outputfilename = 'None'
            
        sim = MyVcfSim(chrom, sitesize, ploidy, population, mutationrate, percentmissing, percentsitemissing, seed, outputfilename, samp_num, 'population.txt', 'vcf', sample_names, population_mode, time,  hmm_baseline, hmm_multiplier, hmm_p_good_to_bad, hmm_p_bad_to_good)

        sim.simulate_vcfs()
        seed+=1

def main():
    #argument parser using argparse
    
    parser = argparse.ArgumentParser(description = 'VCFSim is used to accurately simulate VCFs')
    parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')
    optional = parser.add_argument_group('optional arguments')

    optional.add_argument('--chromosome', type=str, nargs = '?', help = 'Chromosome name', required = False)
    optional.add_argument('--replicates', type=int, nargs = '?', help = 'Amount of times for Simulator to run', required = False)
    required.add_argument('--seed', type=int, nargs = '?', help = 'Random seed for VCFSim to use', required = True)
    optional.add_argument('--sequence_length', type=int, nargs = '?', help = 'Size of your site', required = False)
    optional.add_argument('--ploidy', type=int, nargs = '?', help = 'Ploidy for your VCF', required = False)
    optional.add_argument('--Ne', type=int, nargs = '?', help = 'Effective population size of the simulated population', required = False)
    optional.add_argument('--mu',type=float, nargs = '?', help = 'Mutation rate in the simulated population', required = False)  
    required.add_argument('--percent_missing_sites', type=int, nargs = '?', help = 'Percent of rows missing from your VCF', required = False) #percent missing sites   
    required.add_argument('--percent_missing_genotypes', type=int, nargs = '?', help = 'Percent of samples missing from your VCF', required = True) #percent missing genotypes, percent of genotypes that are missing from samples   

    #HMM commands
    optional.add_argument('--hmm_baseline', type=float, nargs='?', help='Baseline site missing probability in good regions', required=False) 
    optional.add_argument('--hmm_multiplier', type=float, nargs='?', help='Multiplier for missingness in bad regions', required=False)
    optional.add_argument('--hmm_p_good_to_bad', type=float, nargs='?', help='Probability of switching from good to bad region', required=False)
    optional.add_argument('--hmm_p_bad_to_good', type=float, nargs='?', help='Probability of switching from bad to good region', required=False)
 
    optional.add_argument('--output_file', nargs = '?', help = 'Filename of outputed vcf, will automatically be followed by seed', required = False)    

    optional.add_argument('--chromosome_file', nargs = '?', help = 'Specified file for multiple chromosome inputs', required = False)

    optional.add_argument('--population_mode', type=int, nargs='?', help='1 = single population (default), 2 = population C splits into A and B at given time', required=False)
    optional.add_argument('--time', type=int, nargs='?', help='Time of split (only used if --population_mode is 2)', required=False)

    #choose one way to specify samples
    #either a numeric sample size or explicit names or a file of names
    sample_group = parser.add_mutually_exclusive_group(required = True)
    sample_group.add_argument('--sample_size', type = int, nargs = '?', help = 'Amount of samples from population in VCF')
    sample_group.add_argument('--samples', nargs = '*', help = 'Custom sample names space separated')
    sample_group.add_argument('--samples_file', nargs = '?', help = 'File with one whitespace separated line of sample names')

    args = parser.parse_args()


    #HMM argument handling
    hmm_args = [args.hmm_baseline, args.hmm_multiplier, args.hmm_p_good_to_bad, args.hmm_p_bad_to_good]

    hmm_provided = any(x is not None for x in hmm_args)

    #if percentmissing is provided and hmm params are also provided, warn and ignore HMM
    if args.percent_missing_sites is not None and hmm_provided:
        warnings.warn("Both --percent_missing_sites and HMM parameters were provided. Using deterministic percent missing. If you want HMM-based site missingness, don't include --percent_missing_sites.")

    if args.percent_missing_sites is None:
        if not hmm_provided:
            raise ValueError("Error: Either --percent_missing_sites must be provided, or all four HMM parameters must be specified.")
        if any(x is None for x in hmm_args):
            raise ValueError("Error: All HMM parameters must be provided together (--hmm_baseline, --hmm_multiplier, --hmm_p_good_to_bad, --hmm_p_bad_to_good).")

    #validate HMM parameter ranges
    if hmm_provided:
        if not (0 <= args.hmm_baseline <= 1):
            raise ValueError("Error: --hmm_baseline must be between 0 and 1")
        if args.hmm_multiplier < 1:
            raise ValueError("Error: --hmm_multiplier must be >= 1")
        if not (0 <= args.hmm_p_good_to_bad <= 1):
            raise ValueError("Error: --hmm_p_good_to_bad must be between 0 and 1")
        if not (0 <= args.hmm_p_bad_to_good <= 1):
            raise ValueError("Error: --hmm_p_bad_to_good must be between 0 and 1")

    #read custom names if provided
    custom_names = None

    if args.samples_file is not None:
        try:
            with open(args.samples_file, 'r') as sf:
                line = sf.readline().strip()

            #Check if it's a comma separated list
            if ',' not in line:
                raise ValueError("Error: Sample names in the file must be comma-separated (e.g., A1,B1,C1)")

            #Convert to space separated string and split into list
            line = line.replace(',', ' ')
            tokens = line.split()

            if len(tokens) < 1:
                print('Error: samples file must contain at least one name')
                sys.exit(0)

            custom_names = tokens

        except Exception as e:
            print(f'Error reading samples file: {e}')
            sys.exit(0)

    elif args.samples is not None and len(args.samples) > 0:
        #Validate that the input is a single comma separated string
        if len(args.samples) == 1 and ',' in args.samples[0]:
            args.samples = args.samples[0].replace(',', ' ').split()
        elif len(args.samples) > 1:
            raise ValueError("Error: --samples must be passed as a single comma-separated string (e.g., A1,B1,C1)")

        custom_names = args.samples

    if args.population_mode is not None:
        
        population_mode = args.population_mode

        if(population_mode == 2):
            if(args.samples is not None and len(args.samples) % 2 != 0):
                raise ValueError("Error: Number of samples must be even when --population_mode is 2")
            if(custom_names is not None and len(custom_names) % 2 != 0):
                raise ValueError("Error: Number of samples must be even when --population_mode is 2")
            if(args.sample_size is not None and args.sample_size % 2 != 0):
                raise ValueError("Error: Number of samples must be even when --population_mode is 2")

    else:
        population_mode = 1

    if args.time is not None:
        time_value = args.time
    else:
        time_value = 1000 #We pass dummy time value if population_mode is 1 to avoid error, however this doesn't effect performance in any way
        if(population_mode == 2):
            raise ValueError("Error: Time must be passed when --population_mode is 2")

    if population_mode != 1 and population_mode != 2:
        print("Error: --population_mode must be either 1 or 2")
        sys.exit(0)


    if args.chromosome_file is not None:
        # when using multiple chromosome file we forward the custom names too
        # sample size must be determined here so we never pass None
        if custom_names is not None:
            effective_samp_num = len(custom_names)
        else:
            effective_samp_num = args.sample_size

        multiple_chrom(chromfilename = args.chromosome_file, seed = args.seed, percentmissing = args.percent_missing_sites, percentsitemissing = args.percent_missing_genotypes, outputfile = args.output_file, samp_num = effective_samp_num, sample_names = custom_names, population_mode = population_mode, time = time_value, hmm_baseline=args.hmm_baseline, hmm_multiplier=args.hmm_multiplier, hmm_p_good_to_bad=args.hmm_p_good_to_bad, hmm_p_bad_to_good=args.hmm_p_bad_to_good)
    
    elif (args.chromosome_file is None and (args.chromosome is None or args.replicates is None or args.sequence_length is None
                                      or args.ploidy is None or args.Ne is None or args.mu is None)):
        print("Error, no parameter file is specified, or missing one of the following arguments: chromosome, replicates, sequence_length, ploidy, population size, or mutation rate")

    elif (args.chromosome_file is not None and (args.chromosome is not None or args.replicates is not None or 
                                           args.sequence_length is not None or args.ploidy is not 
                                           None or args.Ne is not None or args.mu is not None)):
        print("Error,  parameter file is specified, however so are one of the following: chromosome, replicates, sequence_length, ploidy, population size, or mutation rate")
    
    elif args.mu >= 1:
        print("Error: Mutation rate must be less than 1")

    elif type(args.output_file) != str and args.output_file is not None:
        print("Error: Output_file must be a string")
        
    else:
        vcf_simulator(chrom = args.chromosome , amountofruns = args.replicates, seed = args.seed, sitesize = args.sequence_length, ploidy = args.ploidy, population = args.Ne, mutationrate = args.mu, percentmissing = args.percent_missing_sites, percentsitemissing = args.percent_missing_genotypes, outputfile = args.output_file, samp_num = args.sample_size, sample_names = custom_names, population_mode = population_mode, time = time_value,  hmm_baseline=args.hmm_baseline, hmm_multiplier=args.hmm_multiplier, hmm_p_good_to_bad=args.hmm_p_good_to_bad, hmm_p_bad_to_good=args.hmm_p_bad_to_good)
    
    #argument checker after to double check
    
warnings.filterwarnings("ignore")
main()