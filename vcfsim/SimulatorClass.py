import msprime
import numpy as np
import sys
import io


# Header lines required by consumers (e.g. pixy) to parse GVCF block records.
GVCF_END_INFO_LINE = (
    '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n'
)
GVCF_NON_REF_ALT_LINE = (
    '##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">\n'
)


class MyVcfSim:

    def __init__(self, chrom, site_size, ploidy, pop_num, mutationrate, percentmissing, percentsitemissing, randoseed, outputfile, samp_num, samp_file, folder, sample_names = None,
                 population_mode = 2, time = 1000, hmm_baseline = None, hmm_multiplier = None, hmm_p_good_to_bad = None, hmm_p_bad_to_good = None, gvcf_output = False, recombination_rate = 0):

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
        self.gvcf_output = bool(gvcf_output)
        self.recombination_rate = recombination_rate

        # store custom names if provided
        if sample_names is not None:
            self.sample_names = list(sample_names)
        else:
            self.sample_names = None

        # prebuild small strings used many times
        joiner_temp = '|'
        self.change_missing = joiner_temp.join(['.' for _ in range(self.ploidy)])

    def make_site_mask_hmm(self, p_base, multiplier, p_good_to_bad, p_bad_to_good, rng):
        n = self.site_size

        #states are 0 for good and 1 for bad
        states = np.zeros(n, dtype=int)

        #start in good state
        current_state = 0

        for i in range(n):
            states[i] = current_state
            if current_state == 0:
                if rng.random() < p_good_to_bad:
                    current_state = 1
            else:
                if rng.random() < p_bad_to_good:
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

            if rng.random() < p:
                mask[i] = 1

        return mask


    def make_site_mask(self, rng):

        #Generate site mask as we previously were
        if self.percentmissing is not None:
            temp_percent = self.percentmissing / 100
            temp_var = self.site_size * temp_percent

            temp_array = np.zeros(self.site_size, dtype=int)

            k = round(temp_var)
            if k > 0:
                order = rng.permutation(self.site_size)
                missing_idx = order[:k]
                temp_array[missing_idx] = 1

            return temp_array

        #Generate site mask with our HMM
        else:

            p_base = self.hmm_baseline        #baseline missingness chance when in good state
            multiplier = self.hmm_multiplier #missingness multiplier when in bad state
            p_gb = self.hmm_p_good_to_bad     #chance of going from good state to bad state
            p_bg = self.hmm_p_bad_to_good     #chance of going from bad state to good

            return self.make_site_mask_hmm(p_base, multiplier, p_gb, p_bg, rng)

    def make_missing_vcf(self, ts):
        rng = np.random.default_rng(self.randoseed)
        site_mask = self.make_site_mask(rng)

        # get VCF header from tskit via a zero-site stub tree sequence.
        # stripping all sites/mutations lets write_vcf() produce only the header
        stub_tables = ts.dump_tables()
        stub_tables.sites.clear()
        stub_tables.mutations.clear()
        stub_ts = stub_tables.tree_sequence()
        hdr_buf = io.StringIO()
        stub_ts.write_vcf(hdr_buf, position_transform="legacy")
        hdr_buf.seek(0)
        raw = list(hdr_buf)
        # raw[-1] is the #CHROM line; raw[:-1] are the ##... meta lines
        col_parts  = raw[-1].rstrip('\n').split('\t')
        fixed_cols = col_parts[:9]   # '#CHROM' … 'FORMAT'
        samp_cols  = col_parts[10:]  # skip tsk_0 (index 9)
        if self.sample_names is not None:
            samp_cols = list(self.sample_names)
        extra_header = ''
        if self.gvcf_output:
            extra_header = GVCF_END_INFO_LINE + GVCF_NON_REF_ALT_LINE
        header_text = ''.join(raw[:-1]) + extra_header + '\t'.join(fixed_cols + samp_cols) + '\n'

        # precompute all per-site missing-sample selections in one NumPy call.
        # using argpartition avoids calling the rando-generator once per site in the loop!
        n_real     = self.samp_num - 1
        k_missing  = round((self.percentsitemissing / 100) * n_real)
        n_included = int(np.sum(site_mask == 0))

        if k_missing > 0 and n_included > 0:
            rand_vals = rng.random((n_included, n_real))
            kth = min(k_missing, n_real - 1)   # kth must be < axis size
            missing_by_site = np.argpartition(rand_vals, kth, axis=1)[:, :k_missing] + 1
            # +1: 0-indexed real sample → 1-indexed (tsk_0 at index 0 is never missing)
        else:
            missing_by_site = np.empty((n_included, 0), dtype=np.intp)

        # stream through variants — one pass gives both genotypes and alleles (its fast!)
        rows = []
        site_idx = 0
        # gVCF collapse state: a maximal run of consecutive invariant sites with
        # an identical per-sample GT and no positional gap is emitted as one
        # block record (ALT=<NON_REF>, INFO=END=<last_pos>). Variant sites and
        # positional gaps (created by site_mask) flush the run.
        gvcf_pending = None  # dict: first_pos, last_pos, ref, gt_fields
        for var in ts.variants():
            if site_mask[var.site.id]:
                continue

            gts       = var.genotypes.astype(np.int16)  # (samp_num * ploidy,)
            alleles   = list(var.alleles)
            n_alleles = len(alleles)
            pos       = int(var.site.position) + 1

            # Allele remapping: tsk_0 haplotype 0 determines the reference allele
            ref_code = int(gts[0])
            if ref_code != 0 and ref_code < n_alleles:
                remap = np.arange(n_alleles, dtype=np.int16)
                remap[0]        = n_alleles - 1
                remap[ref_code] = 0
                for j in range(1, n_alleles):
                    if remap[j] == j:
                        remap[j] = j - 1 if j > ref_code else j
                gts      = remap[gts]
                new_ref  = alleles[ref_code]
                new_alts = [alleles[j] for j in range(1, n_alleles) if j != ref_code] + [alleles[0]]
            else:
                new_ref  = alleles[0]
                new_alts = list(alleles[1:])

            # apply missing: mark selected samples' haplotypes with -1 sentinel
            missing_samps = missing_by_site[site_idx]
            if len(missing_samps) > 0:
                hap_idx = np.concatenate([
                    np.arange(s * self.ploidy, (s + 1) * self.ploidy, dtype=np.intp)
                    for s in missing_samps
                ])
                gts[hap_idx] = -1

            # drop allele codes unused by any non-missing sample
            used_alts = sorted(set(int(g) for g in gts if g >= 1))
            if used_alts:
                max_code = max(used_alts)
                compact  = np.zeros(max_code + 1, dtype=np.int16)
                for new_i, old_i in enumerate(used_alts, 1):
                    compact[old_i] = new_i
                valid = gts >= 0
                gts   = np.where(valid, compact[np.clip(gts, 0, max_code)], np.int16(-1))
                alt_str = ','.join(new_alts[c - 1] for c in used_alts)
            else:
                gts[gts > 0] = 0
                alt_str = '.'

            # build GT strings for real samples (tsk_0 occupies positions 0..ploidy-1)
            real_gts    = gts[self.ploidy:].reshape(n_real, self.ploidy)
            missing_set = set(int(m) for m in missing_samps)
            gt_fields   = [
                self.change_missing if s in missing_set
                else '|'.join(str(int(h)) for h in haps)
                for s, haps in enumerate(real_gts, 1)
            ]

            if self.gvcf_output:
                is_invariant = alt_str == '.'
                if (
                    is_invariant
                    and gvcf_pending is not None
                    and gt_fields == gvcf_pending['gt_fields']
                    and pos == gvcf_pending['last_pos'] + 1
                ):
                    gvcf_pending['last_pos'] = pos
                    site_idx += 1
                    continue
                if gvcf_pending is not None:
                    rows.append(self._gvcf_block_line(gvcf_pending))
                    gvcf_pending = None
                if is_invariant:
                    gvcf_pending = {
                        'first_pos': pos,
                        'last_pos': pos,
                        'ref': new_ref,
                        'gt_fields': gt_fields,
                    }
                    site_idx += 1
                    continue

            rows.append('\t'.join(
                [self.chrom, str(pos), '.', new_ref, alt_str, '.', 'PASS', '.', 'GT']
                + gt_fields
            ) + '\n')
            site_idx += 1

        if self.gvcf_output and gvcf_pending is not None:
            rows.append(self._gvcf_block_line(gvcf_pending))
            gvcf_pending = None

        if self.outputfile != 'None':
            with open(self.outputfile, 'w') as out:
                out.write(header_text)
                out.writelines(rows)
        else:
            sys.stdout.write(header_text + ''.join(rows))

    def _gvcf_block_line(self, pending):
        """
        Render a collapsed run of consecutive invariant sites as a single
        GVCF block record: ALT=<NON_REF>, INFO=END=<last_pos>. A single-site
        "block" is rendered the same way so downstream parsers can rely on
        block records always carrying END.
        """
        return '\t'.join(
            [
                self.chrom,
                str(pending['first_pos']),
                '.',
                pending['ref'],
                '<NON_REF>',
                '.',
                'PASS',
                f"END={pending['last_pos']}",
                'GT',
            ]
            + pending['gt_fields']
        ) + '\n'

    def make_population_file(self):
        np.random.seed(self.randoseed)
        file = open(self.samp_file, "a")

        for x in range(0, self.samp_num):
            file.write("tsk_"+ str(x) + "\t1\n")
        file.close()

    def simulate_vcfs(self):
        if self.population_mode == 1:
            #same as always
            ts = msprime.sim_ancestry(samples=[msprime.SampleSet(self.samp_num, ploidy=self.ploidy)], population_size=self.pop_num, random_seed=self.randoseed, sequence_length=self.site_size, recombination_rate=self.recombination_rate)
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

            nA = (self.samp_num-1) // 2  #Half of the samples (rounded down) from population A
            nB = (self.samp_num-1) - nA  #The remaining samples from population B
            samples = [msprime.SampleSet(nA+1, population="A", ploidy=self.ploidy), msprime.SampleSet(nB, population="B", ploidy=self.ploidy)]

            ts = msprime.sim_ancestry(samples=samples, demography=demography, random_seed=self.randoseed, sequence_length=self.site_size, recombination_rate=self.recombination_rate)

        difference_counter = range(self.site_size)

        tables = ts.dump_tables()

        rlg = np.random.default_rng(self.randoseed)
        ancestral_bases = rlg.choice(np.array(["A", "C", "G", "T"]), size=self.site_size)

        for x in difference_counter:
            tables.sites.add_row(x, ancestral_bases[x])

        ts = tables.tree_sequence()
        ts = msprime.sim_mutations(ts, rate=self.mutationrate, random_seed=self.randoseed)

        self.make_missing_vcf(ts)
