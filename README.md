# vcfsim <img align="right" width="160" src="https://github.com/user-attachments/assets/228cfba0-bec0-4b74-8010-412d0f184417">
vcfsim is a command-line tool for generating simulated VCF's (variant call format files for encoding genetic data). It combines a coalescent simulation backend (msprime) with clean and efficient postprocessing to produce a wide variety of biologically realistic VCFs, with parameterized levels of missing data. Realistic VCF's can now be easily simulated with just a few command line arguments!

## Authors 
Paimon Goulart (UC Riverside), Kieran Samuk (UC Riverside)

## Installation
First create and activate a conda environment for vcfsim:

```shell
conda create -n vcfsim_env python=3.10
conda activate vcfsim_env
```

vcfsim is currently available on bioconda, and can be installed by using the following command:
```shell
conda install bioconda::vcfsim
```

For more detailed installation instructions, please visit:  
https://bioconda.github.io/recipes/vcfsim/README.html?highlight=vcfsi#package-package%20&#x27;vcfsim&#x27;

## Arguments 
Here is the list of required/optional arguments to run vcfsim

### Required
--seed [SEED] Random seed for vcfsim to use  

--percent_missing_sites [PERCENT_MISSING_SITES] Percent of rows missing from your VCF  

--percent_missing_genotypes [PERCENT_MISSING_GENOTYPES] Percent of samples missing from your VCF  

One of the following three options must also be provided to set the samples:  
- --sample_size [SAMPLE_SIZE] Amount of samples from population in VCF  
- --samples [SAMPLES ...] Custom sample names, space separated (e.g. A1 B1 C1)  
- --samples_file [SAMPLES_FILE] File containing one whitespace separated line of custom sample names  

### Optional
--chromosome [CHROMOSOME] Chromosome name/label  

--replicates [REPLICATES] Number of replicate VCFs to produce (with varying seeds)

--sequence_length [SEQUENCE_LENGTH] Length of the chromosome to be simulated, in basepairs 

--ploidy [PLOIDY] Ploidy for your VCF  

--Ne [NE] Effective population size of the simulated population(s) 

--mu [MU] Mutation rate in the simulated population(s)  

--output_file [OUTPUT_FILE] Filename of outputed vcf, will automatically be followed by seed  

--chromosome_file [CHROMOSOME_FILE] Specified file for multiple chromosome inputs (see below for details)

--population_mode [1|2] Number of populations simulate. 1 = single population (default), 2 = two populations with a shared history (C splits into A & B).

--time [TIME] Split time for population mode 2 (e.g. generations before present). Required if --population_mode 2 is specified.

Instead of specifying --percent_missing_sites [PERCENT_MISSING_SITES], which produces uniform deterministic missingness across sites, we also provide a more advanced option that uses a Hidden Markov Model (HMM). By leaving --percent_missing_sites blank, users may instead provide the following parameters to allow site-level missingness to be spatially clustered (see bwlow for details). In the case that this option is used rather than the --percent_missing_sites parameter, all four HMM parameters must be specified:

--hmm_baseline [HMM_BASELINE] Baseline probability of site missingness when in a low-missing (good) state

--hmm_multiplier [HMM_MULTIPLIER] Multiplier applied to the baseline missingness probability when in a high-missing (bad) state

--hmm_p_good_to_bad [HMM_P_GOOD_TO_BAD] Probability of switching from a good state to a bad state

--hmm_p_bad_to_good [HMM_P_BAD_TO_GOOD] Probability of switching from a bad state to a good state

## Usage
Typical usage for vcfsim is as follows:

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --sample_size 10
```

This will create a VCV with the name "myvcf1234.vcf", i.e. "myvcf" followed by the seed given for the input.  
If input for replicates were requested higher number than 1, 2 for example, then vcfsim will create two output files by the name of myvcf1234 and myvcf1235, adding one to the seed after every run.  

NOTE: An output file doesn't needed to be specified. If no output file is specified, then the vcf will be redirected to STDOUT.

Screenshot of output file:
<img width="1437" height="458" alt="Image" src="https://github.com/user-attachments/assets/11078b68-6a62-44e0-bf8d-87c34544b2a6" />

### Using Hidden Markov Model Parameters
Instead of specifying `--percent_missing_sites`, which produces uniform deterministic missingness across sites, you can instead specify the following Hidden Markov Model (HMM) parameters to introduce more spatially clustered missing data. Note that when missingness is introduced in this manner, the exact proportion of missing sites is stochastic and will vary between runs, but this missingness can still be approximately determined based on the chosen parameter values.

All four HMM parameters must be specified together when using this option.

```shell
vcfsim --chromosome 1 --replicates 1 --seed 4000 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu 0.000001 --percent_missing_genotypes 0 --hmm_baseline 0.05 --hmm_multiplier 6 --hmm_p_good_to_bad 0.002 --hmm_p_bad_to_good 0.005 --output_file myvcf --sample_size 10
```
To do this, sites transition between low-missing (good) and high-missing (bad) states according to a two-state Markov process. Sites in bad regions have a higher probability of being missing.


### Using custom sample names
Instead of `--sample_size`, you can provide explicit sample names:  

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --samples A1 B1 C1 D1
```

This will automatically set the sample size to 4 and label the VCF columns `A1 B1 C1 D1`.  

You can also read the names from a file containing a single whitespace separated line:  

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --samples_file names.txt
```

Where `names.txt` might contain:  
```
A1 B1 C1 D1 E1
```

Otherwise, sample identifiers will default to tsk_0,...,tsk_n


### Simulating a structured population split (population_mode = 2)
To simulate a demographic split between populations A and B from an ancestral population C:

```shell
vcfsim --chromosome 1 --replicates 1 --seed 1234 --sequence_length 10000 --ploidy 2 --Ne 100000 --mu .000001 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --sample_size 10 --population_mode 2 --time 1000
```

### Multiple chromosome inputs
Another way vcfsim can be used is by providing a file for multiple chromosome inputs.  

Your input file should be in the form of a text file, and should be formatted as such:  
<img width="221" alt="Example input file" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/39ca4b31-c58e-4fff-8b52-456849678339">

The columns are in the order of: chromosone, ploidy, sequence length, population size, mutation rate.  
Each row will represent a seperate run of vcfsim, all these runs will be concatenated to the same file in the end.  

The following command should be used when running vcfsim in this way:

```shell
vcfsim --seed 1234 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --sample_size 10 --chromosome_file input.txt
```

You can also combine a param file with custom names:

```shell
vcfsim --seed 1234 --percent_missing_sites 0 --percent_missing_genotypes 0 --output_file myvcf --samples_file names.txt --chromosome_file input.txt
```

When done this way, the output should look like such:  
<img width="1437" height="458" alt="Image" src="https://github.com/user-attachments/assets/11078b68-6a62-44e0-bf8d-87c34544b2a6" />

With the concatenated vcf looking like:  
<img width="1059" alt="ExampleInput" src="https://github.com/Pie115/VCFSimulator-SamukLab/assets/6378028/fb6508eb-34cd-473a-bc19-762858ed4c31">
