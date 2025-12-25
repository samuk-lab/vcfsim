import msprime
import numpy as np
import pandas as pd
import sys
import os
import time
from IPython.display import SVG, display
import io

class MyVcfSim:
    
    def __init__(self, chrom, site_size, ploidy, pop_num, mutationrate, percentmissing, percentsitemissing, randoseed, outputfile, samp_num, samp_file, folder, sample_names = None,
                 population_mode = 2, time = 1000, hmm_baseline = None, hmm_multiplier = None, hmm_p_good_to_bad = None, hmm_p_bad_to_good = None):

        self.chrom = chrom
        self.site_size = site_size
        self.ploidy = ploidy
        self.pop_num = pop_num
        self.mutationrate = mutationrate
        self.percentmissing = percentmissing
        self.percentsitemissing = percentsitemissing
        self.randoseed = randoseed
        self.outputfile = outputfile
        self.samp_num = samp_num + 1
        self.samp_file = samp_file
        self.folder = folder
        self.population_mode = population_mode
        self.time = time
        self.hmm_baseline = hmm_baseline
        self.hmm_multiplier = hmm_multiplier
        self.hmm_p_good_to_bad = hmm_p_good_to_bad
        self.hmm_p_bad_to_good = hmm_p_bad_to_good


        # store custom names if provided
        if sample_names is not None:
            self.sample_names = list(sample_names)
        else:
            self.sample_names = None

        # prebuild small strings used many times
        change_list_zero = []
        for i in range(self.ploidy):
            change_list_zero.append('0')
        joiner_temp = '|'
        self.change_zero = joiner_temp.join(change_list_zero)

        change_list_missing = []
        for i in range(self.ploidy):
            change_list_missing.append('.')
        self.change_missing = joiner_temp.join(change_list_missing)

        self.col_start = None
        self.col_end = None

    def make_site_mask_hmm(self, p_base, multiplier, p_good_to_bad, p_bad_to_good):
        n = self.site_size

        #states are 0 for good and 1 for bad
        states = np.zeros(n, dtype=int)

        #start in good state
        current_state = 0

        for i in range(n):
            states[i] = current_state
            if current_state == 0:
                if np.random.random() < p_good_to_bad:
                    current_state = 1
            else:
                if np.random.random() < p_bad_to_good:
                    current_state = 0

        #now generate the missing data conditioned on the good or bad state
        mask = np.zeros(n, dtype=int)

        for i in range(n):
            if states[i] == 0:
                p = p_base
            else:
                p = p_base * multiplier

            p = min(1.0, p) #if probability is somehow greater than 1 (maybe a large multiplier was put accidentally)
            #we cap this at 1 so no weird stuff

            if np.random.random() < p:
                mask[i] = 1

        return mask

    
    def make_site_mask(self):

        #Generate site mask as we previously were
        if self.percentmissing is not None:
            temp_percent = self.percentmissing / 100
            temp_var = self.site_size * temp_percent

            temp_array = np.zeros(self.site_size, dtype=int)

            k = round(temp_var)
            if k > 0:
                order = np.random.permutation(self.site_size)
                missing_idx = order[:k]
                temp_array[missing_idx] = 1

            return temp_array

        #Generate site mask with our HMM
        else:

            p_base = self.hmm_baseline        #baseline missingness chance when in good state
            multiplier = self.hmm_multiplier #missingness multiplier when in bad state
            p_gb = self.hmm_p_good_to_bad     #chance of going from good state to bad state
            p_bg = self.hmm_p_bad_to_good     #chance of going from bad state to good

            return self.make_site_mask_hmm(p_base, multiplier, p_gb, p_bg)
    
    def row_changes(self, row, vcfdata, tempvcf):
        col_start = self.col_start
        col_end = self.col_end
        altlist = row[col_start:col_end+1].values

        #pick site missing samples
        randomsitemissing = round((self.percentsitemissing / 100) * (self.samp_num - 1))
        randomsites = np.arange(1, self.samp_num)
        np.random.shuffle(randomsites)
        randomsites = randomsites[:randomsitemissing]

        #temporarily set missing samples to all zeros genotype
        change = self.change_zero
        for idx_site in randomsites:
            altlist[idx_site] = change

        #refindex for old reference
        refindex = 0
        if ((1 in randomsites) and ((self.percentsitemissing / 100) != 1)):
            for i in range(1, self.samp_num):
                if i not in randomsites:
                    refindex = i * self.ploidy
                    break

        #Parse ALT as alleles
        alt_field = row['ALT']
        alts = [] if alt_field == "." else alt_field.split(',')

        #Flatten all allele integers from the row depending on ploidy
        allele_codes = []
        if self.ploidy != 1:
            for item in altlist:
                for x in item.split('|'):
                    allele_codes.append(int(x))
        else:
            for x in altlist:
                allele_codes.append(int(x))

        #Determine which allele reference sample is
        referencegenome = allele_codes[refindex]  # 0 = REF, >=1 = ALT index

        #Change REF then remap all GTs
        alleles = [row['REF']] + alts  #0=REF, OR 1,..,n=ALTs

        if referencegenome != 0 and referencegenome < len(alleles):
            new_ref = alleles[referencegenome]
            old_ref = alleles[0]

            #all old ALTs except new ref, plus old REF at end
            new_alts = [alleles[i] for i in range(1, len(alleles)) if i != referencegenome]
            new_alts.append(old_ref)

            #Use disctionary to switch references/alts
            index_map = {0: len(new_alts)}
            index_map[referencegenome] = 0
            for i in range(1, len(alleles)):
                if i in index_map:
                    continue
                index_map[i] = i - 1 if i > referencegenome else i

            mapped_codes = [index_map[c] for c in allele_codes]
            row['REF'] = new_ref
            curr_alts = new_alts
        else:
            mapped_codes = allele_codes[:]
            curr_alts = alts[:]

        #compact ALT indices and shrink ALT list accordingly
        used_alts = sorted({c for c in mapped_codes if c >= 1})
        compact_map = {0: 0}
        for new_i, old_i in enumerate(used_alts, start=1):
            compact_map[old_i] = new_i

        mapped_codes = [compact_map[c] for c in mapped_codes]
        compacted_alts = [curr_alts[i - 1] for i in used_alts]  # i>=1
        row['ALT'] = "." if len(compacted_alts) == 0 else ",".join(compacted_alts)

        #rebuild GT strings
        finaldata = []
        for i in range(0, len(mapped_codes), self.ploidy):
            finaldata.append(mapped_codes[i:i + self.ploidy])
        finaldataoutput = ["|".join(str(i) for i in data) for data in finaldata]

        #apply missing as '.|.,,,,|.' or '.' for ploidy = 1
        for idx_site in randomsites:
            finaldataoutput[idx_site] = self.change_missing

        #edge case of 100% site missing make all red '.'
        if ((self.percentsitemissing / 100) == 1):
            row['REF'] = '.'

        max_alt_index = len(compacted_alts)
        for gt in finaldataoutput:
            if gt == self.change_missing:
                continue
            for a_code in gt.split('|'):
                ai = int(a_code)
                assert 0 <= ai <= max_alt_index, f"Illegal allele {ai} with ALT count {max_alt_index}"

        #write columns back
        row[col_start:col_end+1] = finaldataoutput
        return row


    def make_missing_vcf(self, ts):
        np.random.seed(self.randoseed)
        site_mask = self.make_site_mask()

        buf = io.StringIO()
        ts.write_vcf(buf, site_mask=site_mask)
        buf.seek(0)

        header_lines = []
        data_lines = []

        linenum = 1
        for line in buf:
            if (linenum <= 5):
                header_lines.append(line)
            else:
                if (linenum == 6):
                    line = line[1:]
                data_lines.append(line)
            linenum = linenum + 1

        body_buf = io.StringIO()
        for line in data_lines:
            body_buf.write(line)
        body_buf.seek(0)

        vcfdata = pd.read_csv(body_buf, delimiter = '\t')
        
        a = 'tsk_0'
        self.col_start = vcfdata.columns.get_loc(a)
        self.col_end = vcfdata.columns.get_loc(f"tsk_{self.samp_num - 1}")
        
        tempvcf = vcfdata
        
        vcfdata = vcfdata.apply(self.row_changes, axis = 1, args = (vcfdata, tempvcf))
 
        vcfdata["CHROM"] = self.chrom
        vcfdata["POS"] = vcfdata["POS"] + 1
        vcfdata["ID"] = '.'

        if ("tsk_0" in vcfdata.columns):
            del vcfdata["tsk_0"]

        # rename tsk style columns to custom names when provided
        if (self.sample_names is not None):
            # we rename tsk_1 to first name and so on
            idx_val = 1
            name_index = 0
            while(idx_val < self.samp_num and name_index < len(self.sample_names)):
                old_name = "tsk_" + str(idx_val)
                new_name = self.sample_names[name_index]
                if old_name in vcfdata.columns:
                    vcfdata.rename(columns={old_name: new_name}, inplace=True)
                idx_val = idx_val + 1
                name_index = name_index + 1
        
        if (self.outputfile != 'None'):
            with open(self.outputfile, 'w') as fout:
                
                #Replace the existing ##source line with a single version to display version and command
                for i, line in enumerate(header_lines):
                    if line.startswith("##source="):
                        new_source = '##source=tskit 0.6.4, vcfsim 1.0.24.alpha, ' + ' '.join(sys.argv).replace('\t', ' ').strip() + '\n'
                        header_lines[i] = new_source
                        break

                #Write all header lines including the modified one
                for line in header_lines:
                    fout.write(line)

                #Write the data portion of the VCF
                csv_buf = io.StringIO()
                vcfdata.to_csv(csv_buf, index=False, sep='\t', header=True)
                csv_text = csv_buf.getvalue()
                fout.write('#')
                fout.write(csv_text)
        else:
            display(vcfdata.to_string())
        
    def make_population_file(self):
        np.random.seed(self.randoseed)
        file = open(self.samp_file, "a")

        for x in range(0, self.samp_num):
            file.write("tsk_"+ str(x) + "\t1\n")
        file.close()

    def simulate_vcfs(self):
        np.random.seed(self.randoseed)


        if self.population_mode == 1:
            #same as always
            ts = msprime.sim_ancestry(samples=[msprime.SampleSet(self.samp_num, ploidy=self.ploidy)], population_size=self.pop_num, random_seed=self.randoseed, sequence_length=self.site_size)
        else:
            #two population mode: C (size 2x) splits into A and B at self.time
            demography = msprime.Demography()
            demography.add_population(name="A", initial_size=self.pop_num)
            demography.add_population(name="B", initial_size=self.pop_num)
            demography.add_population(name="C", initial_size=self.pop_num * 2)
            demography.add_population_split(time=self.time, derived=["A", "B"], ancestral="C")
            
            #Split the total number of samples into two groups: A and B
            #We intentionally *do not* simulate samples from C (the ancestral population)
            #This ensures that all samples come from the derived populations A and B

            nA = self.samp_num // 2  #Half of the samples (rounded down) from population A
            nB = self.samp_num - nA  #The remaining samples from population B
            samples = [msprime.SampleSet(nA, population="A", ploidy=self.ploidy), msprime.SampleSet(nB, population="B", ploidy=self.ploidy)]
            
            ts = msprime.sim_ancestry(samples=samples, demography=demography, random_seed=self.randoseed, sequence_length=self.site_size)

        rlg = np.random.default_rng(self.randoseed)
        rchar = rlg.integers(low=0, high=4)

        difference_counter = range(self.site_size)

        tables = ts.dump_tables()

        for x in difference_counter:
            tables.sites.add_row(x, "A")

        ts = tables.tree_sequence()
        ts = msprime.sim_mutations(ts, rate=self.mutationrate, random_seed=self.randoseed)

        self.make_missing_vcf(ts)
