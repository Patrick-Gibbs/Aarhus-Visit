import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from random import randint
from tqdm import tqdm
from bed_reader import open_bed
from scipy.stats import mode
import os
import random
import statistics
import json
import matplotlib.colors as mcolors

# To help with the generation of simulation data #
CHROMOSOME = 0
START = 1
FINISH = 2
NUM_HUMAN_AUTOSOMES = 22

# For the anylisis of simulation data #
CRITICAL_INDEX = -1
LEADING_P_CASUAL = 'leading_p_is_casual'
LEADING_POST_CASUAL = 'highest_post_prob_is_casual'
SIG_P_TP = 'num_true_possitive_p5*10**(-8)'
SIG_P_FP = 'num_false_possitive_p5*10**(-8)'
POST_95_TP = 'num_true_possitive_post95'
POST_95_FP = 'num_false_possitive_post95'
P_CRITICAL = 5*10**(-8)
POST_CRITICAL = 0.95
PROB_COLUMN = 'Probability'
CASUAL_SIZE = 8
SIZE = 3

##############################################
# To help with sampling subsets of the genome #
##############################################

def make_finemap_file(f_name, z_file_base, cors_base, num_samples, samples, model_path):
    """Makes a file for Finemap to run from

    Args:
        `dir` (str): 
            Where the finemap files will be written 
        `z_file_base` (str):
            Where the z files are stored
        `cors_base` (str):
            Where the cors are stored
    """

    f = open(f_name, "w")
    output = []
    output.append(f"z;ld;snp;config;cred;n_samples")
    for i in samples:
        output.append(f"{z_file_base}/sample_{i}.z;{cors_base}/sample_{i}.ld;{model_path}/snp_{i}.snp;{model_path}/config_{i}.config;{model_path}/cred_{i}.cred;{num_samples}")
    f.write('\n'.join(output))

def make_finemap_z_files(dir, gwas_base, binary_base, samples):
    """Makes a file for Finemap to run from

    Args:
        `dir` (str): 
            Where the finemap files will be written 
        `gwas_base` (str):
            Where the gwas files are stored (from LDAK)
        `ukbb_base` (str):
            Where the ukbb binaries are stored
        `samples` (iter):
            The of samples to be generated. 
    """
    if type(samples) == int:
        samples = range(samples)

    for i in samples:
        df = {}
        bim = pd.read_csv(f'{binary_base}/sample_{i}.bim', sep='\t', header=None)
        df['rsid'] = np.array(bim[1])
        df['chromosome'] = np.array(bim[0])
        df['position'] = np.array(bim[3])
        df['allele1'] = np.array(bim[4])
        df['allele2'] = np.array(bim[5])
        gwas_assoc = pd.read_csv(f'{gwas_base}/sample_{i}.assoc', sep=' ')
        df['maf'] = np.array(gwas_assoc['MAF'])
        df['beta'] = np.array(gwas_assoc['Effect'])
        df['se'] = np.array(gwas_assoc['SD'])
        pd.DataFrame(df).to_csv(f'{dir}/sample_{i}.z', sep=' ', index=None)



def random_indivdual_subset(number, fam, output, seed=42):
    """ Creates a plink formatted fam file which specifies a random
    subset of individuals used for a analysis.

    Args:
    `number`: str
        The number of individuals in the subset
    `fam`: str 
        The fam file from which individuals are sampled
    `output`: str
        The location that the new fam file should be written to 
    `seed`: int
        A random seed such that generation of results can occur deterministically.
    Outputs: None
    """
    random.seed(seed)
    accession_numbers = pd.read_csv(fam, header=None, sep='\t')
    index = np.array(sorted(random.sample(
        list(range(len(accession_numbers))), number)))
    # column 0 and column 1
    accession_numbers.iloc[index][[0, 1]].to_csv(
        output, header=None, sep='\t', index=None)


def is_overlap(p_bounds, c_bound):
    """
    Helper function to sample domains of the genome. Tests if a prospective
    bounds to sample the genome overlaps with past samples.
    
    Args:
    `p_bounds`: [(int, int), ..., (int, int)]
        A list of two element tuples which specify the past bounds.
    `c_bound`: (int, int)
        A tuple of the prospective bound
    Returns: bool
        whether `c_bound` overlaps with of in `p_bounds`
    """

    for p_bound in p_bounds:
        if (c_bound[CHROMOSOME] == p_bound[CHROMOSOME]
                and pd.Interval(p_bound[START], p_bound[FINISH]).overlaps(pd.Interval(c_bound[START], c_bound[FINISH]))):
            return True
    return False


def sample_snp(previous_bounds, ukbb_bim, approx_chromosome_size, size=10**6, max_snps=1500):
    """
    Selects a random genome location of `size`. Also ensures that selected region has not before been selected
    and is thus not in previous_bounds.

    Args:
    `previous_bounds`: [(int, int, int), ...]
        A list of tuples which specify the bounds of previously sampled regions
    `ukbb_bim`: pd.DataFrame
        The bim containing the accession ids
    `approx_chromosome_size`: [int, ...]
        A list of the approximate sizes of each chromosome
    `size`: int
        The size of the region to be sampled
    """

    # selected a random snp from the bim
    snp_index = randint(0, len(ukbb_bim))
    snp = ukbb_bim.iloc[snp_index]

    # extract snp infomation
    snp_chromosome = snp['chrom_number']
    snp_pos = snp['chrom_pos']

    

    # extract += 1mb from the randomly seleted SNP
    start = int(snp_pos - size/2)
    end = int(snp_pos + size/2 - 1)

    # get MAF
    

    # if the the bound is outside the chromosome try again
    if (start < 0 or end > approx_chromosome_size[snp_chromosome - 1]
            or is_overlap(previous_bounds, (snp_chromosome, start, end))):
        return sample_snp(previous_bounds, ukbb_bim, approx_chromosome_size, size=size)


    # ensure that the region has at most `max_snps` snps
    upper =  ukbb_bim.iloc[snp_index + max_snps//2]
    lower =  ukbb_bim.iloc[snp_index - max_snps//2]
    if upper['chrom_number'] != snp_chromosome or upper['chrom_pos'] < end:
        return sample_snp(previous_bounds, ukbb_bim, approx_chromosome_size, size=size)
    if lower['chrom_number'] != snp_chromosome or lower['chrom_pos'] > start:
        return sample_snp(previous_bounds, ukbb_bim, approx_chromosome_size, size=size)

    return snp_chromosome, start, end


def sample_genome_domains(sample_size, sample_number, genome_path, output_path, max_snps=1500):
    """
    Samples a random regions of the genome, and write a file speciefyinh the SNPs
    that are in that genomic region. The sample regions are such that they do not
    overlap

    Args
    `sample_size`: int
        The size of the genomic region sampled (in basepairs)
    `sample_number`: int
        The number of genomic samples that should be generated
    `genome_path`: str
        The path to the genome binary (e.g. uk biobank bim)
    `output_path`: str
        The location at which the binary should be output
    Returns: None
    """

    ### sample genome ###
    ukbb_bim = pd.read_csv(f"{genome_path}.bim", sep='\t', header=None)

    # making another column containing the position of each SNP on  a chromosome as an integer
    ukbb_bim['chrom_pos'] = [
        int(item[item.index(':') + 1:item.index('_')]) for item in ukbb_bim[1]]
    ukbb_bim.rename(columns={0: 'chrom_number'}, inplace=True)

    approx_chromosome_size = [max(ukbb_bim[ukbb_bim['chrom_number'] == chrom]['chrom_pos'])
                              for chrom in range(1, NUM_HUMAN_AUTOSOMES + 1)]

    # select 200 1mb regions
    previous_bounds = []
    for i in tqdm(range(sample_number)):
        previous_bounds.append(sample_snp(
            previous_bounds, ukbb_bim, approx_chromosome_size, size=sample_size))
        with open(f'{output_path}_{i}', 'w') as f:
            f.write(
                f'{previous_bounds[i][CHROMOSOME]} {previous_bounds[i][START]} {previous_bounds[i][FINISH]} sample_{i+1}')

def get_bed(bfile, i):
    ukbb_path = f'{bfile}/sample_{i}.bed'
    bed = open_bed(ukbb_path)
    X = bed.read()
    # impute NaN to the mode
    modes = mode(X, nan_policy='omit').mode
    for i in range(X.shape[1]):
        X[np.isnan(X[:, i])] = modes[i]
    return X

def make_ld_matrix(bfile, sample_number, output_path):
    for i in range(sample_number):
        # read plink into python to simulate phenotype
        X = get_bed(bfile, i)
        R = np.corrcoef(X.T)
        # ensure that the matrix is symmetric if not make it such
        R = (R + R.T) / 2
        # round to 5 decimal places
        R = np.round(R, 5)

        np.savetxt(f'{output_path}/sample_{i}.ld', R, delimiter=' ')
    


################################################
# To help with the of simulation of phenotypes #
################################################


def simulate_traits_for_genome_chunks(bfile, phenotypes_output_path, snps_info,
                                      heritability, num_casual, sample):
    """
    Simulates traits for a small sample of the genome.

    Args
    `bfile`: str
        The path to the binary genome file (.bim)
    `phenotypes_output_path`: str
        The path to the path where the phenotypes should be written
    `snps_info` : str
        The path to where information about which SNPs are casual
        is stored
    `heritability` : str
        The Heritability of the trait that is getting simulated
    `num_casual`
        The number of casual variants that drive the trait
    """

    # read in genomic regions
    ukbb_path = f'{bfile}_{sample}'
    bim_df = pd.read_csv(f'{ukbb_path}.bim', delimiter='\t', header=None)
    bim_df.columns = ['Chromosome code', 'Variant identifier',
                      'Position in morgans or centimorgans', 'Base-pair coordinate ', 'Allele 1', 'Allele 2']
    fam_df = pd.read_csv(f'{ukbb_path}.fam', delimiter=' ', header=None)

    # read plink into python to simulate phenotype
    bed = open_bed(f'{ukbb_path}.bed')
    X = bed.read()
    num_markers = X.shape[1]

    # impute NaN to the mode
    modes = mode(X, nan_policy='omit').mode
    for i in range(X.shape[1]):
        X[np.isnan(X[:, i])] = modes[i]

    if num_casual == 1:
        # one snp with effect size one so phenotypes start as being the same as allele state for 1 locus
        casual = np.array([num_markers // 2])
        beta = np.array([1])
        y_g = X[:, casual] @ beta
    elif num_casual > 0:
        casual = np.random.randint(
            num_markers//10, num_markers-num_markers//10, size=num_casual)
        beta = np.random.normal(0, 2, size=num_casual)
        y_g = X[:, casual] @ beta

    if num_casual == 0:
        y = np.random.normal(0, 1, size=(len(X.shape[1])))
    else:
        # add noise and rescale to predefined heritability
        y = (y_g) + np.random.normal(0, statistics.variance(y_g) /
                                    heritability - statistics.variance(y_g), size=(len(y_g)))

    phenotype_df = pd.DataFrame(
        {'acc1': fam_df[0], 'acc2': fam_df[1], 'pheno': y})
    phenotype_df.to_csv(
        f'{phenotypes_output_path}/sample_{sample}', sep=' ', header=None, index=None)

    if num_casual == 0:
        return

    # save the results
    casual_info = bim_df.iloc[casual].copy()
    if len(casual) == 1:
        casual_info = pd.DataFrame(casual_info)

    casual_info['snp_index_in_chunk'] = casual
    casual_info['effect_size'] = beta
    casual_info['total_he'] = heritability

    casual_info.to_csv(
        f'{snps_info}/snp_info_sample_{sample}', sep='\t', header=True)


#######################################
# For the anylisis of simulation data #
#######################################

def get_post_probs(post_path):
    """
    Used to extract the Post. probabilities from a .probs file
    Args:
    `post_path` (str):
        The path to the post. probs file
    
    Returns: np.array
        An array of the post. probs
    """
    df = pd.read_csv(post_path, sep=' ')
    return np.array(df[PROB_COLUMN])


def get_casual_marker(snp_info):
    """
    Used to extract the relevant info from a snp info file

    Args:
    `snp_info` (str):
        The path to the snp info file

    Returns: np.array
        An array of the casual markers
    """
    return np.array((pd.read_csv(snp_info, sep='\t')['snp_index_in_chunk']))


def get_p_value_gwas(gwas_path):
    """"
    Used to extract the P values from a LDAK gwas file
    Args:
    `gwas_path` (str):
        The path to the gwas file
    Returns: np.array
        An array of the P values
    """
    return np.array(pd.read_csv(gwas_path, sep=' ')['P'])


def get_paths_for_chunk_1_snp(p_vals=None, probs=None, casual=None, i=0):
    """
    Generates the paths for the p values, post. probs and casual markers for a single simulation
    
    Args:
        `p_vals` (str): 
            path to p values
        `probs` (str): 
            path to post. probabilities file
        `casual` (str): 
            path to the snp infomation file
        `i` (int): 
            The ith simulation

    Returns:
        tuple: (p_vals, post, casual) containing as numpy arrays the p values, post. probs and casual markers
    """
    p_vals_path = f'{p_vals}/sample_{i}.pvalues'
    post_path = f'{probs}/sample_{i}.probs'
    casual_path = f'{casual}/snp_info_sample_{i}'

    if p_vals is not None:
        p_vals = get_p_value_gwas(p_vals_path)
    if casual is not None:
        post = get_post_probs(post_path)
    if casual is not None:
        casual = get_casual_marker(casual_path)

    return p_vals, post, casual


def compute_stats_chunk_1_snp(p_vals, probs, casual, p_critical=P_CRITICAL, post_critical=POST_CRITICAL):
    """Generates statistics for a single simulation comparing the leading p value and post. prob to the casual marker
    and also counts the number of true and false positives for both the p value and post. prob

    Args:
        `p_vals` (np.array):
            The p values
        `probs` (np.array): 
            The post. probs
        `casual` (np.array): 
            The casual variants
        `p_critical` (float, optional): 
            The p value for a variant to be considered casual. Defaults to P_CRITICAL.
        `post_critical` (float, optional): 
            The post. probabilty for a variant to be consider casual. Defaults to POST_CRITICAL.

    Returns:
        dict: The statistics for the simulation
    """
    return {
        LEADING_P_CASUAL: np.argmin(p_vals) == casual,
        LEADING_POST_CASUAL: np.argmax(probs) == casual,
        SIG_P_TP: len(np.where(p_vals[casual] <= p_critical)),
        SIG_P_FP: len(np.where(np.delete(p_vals, casual) <= p_critical)[0]),
        POST_95_TP: len(np.where(probs[casual] >= post_critical)),
        POST_95_FP: len(np.where(np.delete(probs, casual) >= post_critical)[0])
    }


def mine_all_tp_fp(p_vals_file, probs_file, casual_file, sim_numbers):
    """Finds the number of true and false positives for each simulation

    Args:
        p_vals_file (str): base path to p_values
        probs_file (str): post. probabilities
        casual_file (str): path to the casual variant
        sim_numbers (iterable): The simulations to be mined

    Returns:
        pd.DataFrame: the number of true and false positives for each simulation
    """
    # indivdual TP/FP for each simulation
    sim_numbers = list(sim_numbers)
    stats = {k: [v] for k, v in compute_stats_chunk_1_snp(
        *get_paths_for_chunk_1_snp(p_vals_file, probs_file, casual_file, sim_numbers[0])).items()}

    # gets stats for each chunk
    for i in sim_numbers[1:]:
        one_sim_stats = compute_stats_chunk_1_snp(
            *get_paths_for_chunk_1_snp(p_vals_file, probs_file, casual_file, i))
        for k in stats.keys():
            stats[k].append(one_sim_stats[k])

    return pd.DataFrame(stats)


def count_tp_fp(all_post_vals, all_casual, post_critical):
    """Given a critical value to define casual variants, counts the number of true and false positives

    Args:
        all_post_vals (np.array([int, ...])): 
            The post. probabilities
        all_casual (np.array([int, ...])): 
            The indices of the casual variants
        post_critical (float): 
            The critical value to define casual variants

    Returns:
        (int, int): The number of true and false positives
    """
    num_true_discoveries = len(np.intersect1d(
        np.where(all_post_vals >= post_critical)[0], all_casual))
    num_false_discoveries = len(np.intersect1d(np.where(all_post_vals >= post_critical)[
                                0], np.delete(np.arange(len(all_post_vals)), all_casual)))
    return num_true_discoveries, num_false_discoveries


def restricted_true_discovery_rate(all_post_vals, all_casual, fp_tolerance):
    """ Iteratively finds lowers post. probability threshold that maximises the number of true positives, 
    subject to a false positive tolerance.

    Args:
        `all_post_vals` (np.array([float, ...])): 
            All post. probabilities
        `all_casual` (np.array([int, ...])):
            The indicies of all casual variants
        `fp_tolerance` (float):
            The false positive tolerance

    Returns:
        (TP rate, tolerance): The true possitive rate, and FP tolerance
    """

    arg_sort = np.argsort(all_post_vals)[::-1]
    all_casual = set(all_casual)
    restricted_post_discovery = (None, max(all_post_vals))
    tp, fp = 0, 0
    for index in arg_sort[:int(len(all_casual)/fp_tolerance)+1]:
        tp += index in all_casual
        fp += index not in all_casual

        if tp/(tp + fp) >= fp_tolerance:
            restricted_post_discovery = (
                tp/len(all_casual), all_post_vals[index])
    
    return restricted_post_discovery

def get_power_fdr_series(all_post_vals, all_casual, strict_inc=True):
    arg_sort = np.argsort(all_post_vals)[::-1]
    tolerance = 1e-6 # Set your tolerance value
    disjoint_indices = np.where(~np.isclose(np.diff(all_post_vals[arg_sort]), 0, atol=tolerance))[0]
    if len(disjoint_indices) == 0:
        disjoint_indices = np.array([0])

    probs = []
    power = []
    fdr = []
    tp = 0
    fp = 0
    last_iter = 0

    for disjoint_index in disjoint_indices:
        tp_i = len(np.intersect1d(arg_sort[last_iter:disjoint_index+1], all_casual))

        fp += (disjoint_index - last_iter + 1) - tp_i
        tp += tp_i
        last_iter = disjoint_index + 1
        probs.append(all_post_vals[arg_sort[disjoint_index]])
        power.append(tp/len(all_casual))
        fdr.append(fp/(tp + fp))
    if not strict_inc:
        return probs, power, fdr
    else:
        path = [-1]
        smallest = fdr[-1]
        for i in range(len(fdr)-2, -1, -1):
            if fdr[i] < smallest:
                path.append(i)
                smallest = fdr[i]
        path = np.array(path[::-1])

        return np.array(probs)[path], np.array(power)[path], np.array(fdr)[path]

def get_genome_wide_post_prob(model):
    return get_post_probs(f'{model}/sample_0.probs')


def plot_power_fdr(models_path, casual_info_path, out_path, samples, strict_inc=True, ignore_bad=False, title='', genome_wide_models=False):
    
    lines = []
    for model in os.listdir(models_path):
        all_post_vals, all_casual, _ = combine_genome_sample_r(None, f'{models_path}/{model}', casual_info_path, samples, ignore_bad=ignore_bad)
        _, power, fdr = get_power_fdr_series(all_post_vals, all_casual, strict_inc=strict_inc)
        fdr = [fdr[i] for i in np.repeat(np.arange(len(fdr)), 2)][1:]
        power = [power[i] for i in np.repeat(np.arange(len(power)), 2)][:len(power)*2-1]
        lines.append((model, fdr, power))

    # plot each line in a different color and use model as the legand
    lines = sorted(lines)
    colors = list(mcolors.TABLEAU_COLORS) + list(sns.color_palette('dark'))
    for (model, fdr, power), color in zip(lines, colors):
        print(model)
        plt.plot(fdr, power, label=model, color=color)

    if genome_wide_models:
        print(genome_wide_models)
        for model in os.listdir(genome_wide_models):
            print(model)
            _, fdr, power = get_power_fdr_series(get_genome_wide_post_prob(f'{genome_wide_models}/{model}'), 
                                             casual_info_path)
            plt.plot(fdr, power,label=model, color=color, linestyle='dashed')

    plt.xlabel('FDR')
    plt.ylabel('Power')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')  # place the legend outside
    plt.title(title)
    plt.tight_layout()

    plt.savefig(out_path)


def test_get_power_fdr_series():
    all_post_vals = np.array([0.9, 0.9, 0.9, 0.9, 0.9, 0.8, 0.7, 0.6, 0.6, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.1])
    all_casual = np.array([1, 2, 4, 5, 10])

    



def make_fine_map_probs(base_path, samples):
    for sample in samples:
        df = pd.read_csv(f'{base_path}/snp_{sample}.snp', sep=' ')        
        df = df.sort_values(by='index')
        df = df[['rsid','allele1', 'allele2', 'mean', 'prob']]
        # rename columns to Predictor A1 A2 Mean Probability
        df.columns = ['Predictor', 'A1', 'A2', 'Mean', 'Probability']
        df.to_csv(f'{base_path}/sample_{sample}.probs', sep=' ', index=None)
    




def combine_genome_sample_r(p_vals_file, probs_file, casual_file, sim_numbers, ignore_bad=False):
    """Combine the post. probabilities, casual markers and p values for a set of simulations 
    into a single arrays respectively

    Args:
        `p_vals_file` (str): 
            base path to p values
        `probs_file` (str): 
            base path to post. probabilities
        `casual_file` (str): 
            base path to marker infomation
        `sim_numbers` (iterable): 
            The simulations to be mined
        `ignore_bad` (bool, optional):
            If true disregards simulations with terrible results. Defaults to False.

    Returns:
        (np.array, np.array, np.array): The post. probabilities, casual markers and p values
    """

    all_casual, all_post_vals, all_p_vals = [], [], []
    for i in sim_numbers:
        p_vals, post, casual = get_paths_for_chunk_1_snp(
            p_vals_file, probs_file, casual_file, i)
        if ignore_bad and (sum(p >= 0.99 for p in post)) > 10*len(casual):
            continue
        if casual is not None:
            all_casual.extend(casual + len(all_post_vals))
        if post is not None:
            all_post_vals.extend(post)
        if p_vals is not None:
            all_p_vals.extend(p_vals)

    all_post_vals = np.array(all_post_vals)
    all_casual = np.array(all_casual)
    all_p_vals = np.array(all_p_vals)
    return all_post_vals, all_casual, all_p_vals

def get_stats_chunk_1_snp(pvals_file, probs_file, casual_file, model, sim_numbers=range(20),
                          critical_post_values=[1, 0.99, 0.95, 0.90, 0.80, 0.7], ignore_bad=False):
    """Generates statistics for a model as writes them as a human readeable text file

    Args:
        pvals_file (str): base path to p values
        probs_file (str): base path to post. probabilities
        casual_file (str): base path to the casual variants
        model (str): path to the model used
        sim_numbers (iterable, optional): Simulations to be included. Defaults to range(20).
        critical_post_values (list, optional): critical post. probabilities. 
        ignore_bad (bool, optional): Whether to include simulations with terrible results. Defaults to False.
    """

    # indivdual TP/FP for each simulation
    df = mine_all_tp_fp(pvals_file, probs_file,
                        casual_file, sim_numbers=sim_numbers)
    df.to_csv(f'{model}/_results.csv')

    # get 95% discovery rate
    all_post_vals, all_casual, all_p_vals = combine_genome_sample_r(pvals_file, probs_file,
                                                                    casual_file, sim_numbers, ignore_bad=ignore_bad)

    post_discoverys = [restricted_true_discovery_rate(
        all_post_vals, all_casual, post_start) for post_start in critical_post_values]
    tp_fp_list = [count_tp_fp(all_post_vals, all_casual, critical)
                  for critical in critical_post_values]

    ignore_bad_suffix = '_ignore_bad' if ignore_bad else ''
    with open(f'{model}/_summary_results{ignore_bad_suffix}.csv', 'w') as f:
        f.write("Number of times leading p_val is casual " +
                str(sum(df[LEADING_P_CASUAL])) + '\n')
        f.write("Number of times leading post. is casual " +
                str(sum(df[LEADING_POST_CASUAL])) + '\n')
        for critical, tp_fp, post_discovery in zip(critical_post_values, tp_fp_list, post_discoverys):
            percent_var_detected, cutoff = post_discovery
            if percent_var_detected == None:
                percent_var_detected = 0
            f.write(f"{(critical*100):>6.2f}% true positive discovery rate: {percent_var_detected*100:>3.2f}% of variants detected with post. cut of off {cutoff:>.2f}\n")
        for critical, tp_fp, post_discovery in zip(critical_post_values, tp_fp_list, post_discoverys):
            symbol = '=' if critical == 1 else '>'
            f.write(
                f"treating post {symbol} {critical:<3.2f} as casual: TP {tp_fp[0]:3d} FP {tp_fp[1]:3d}\n")



def get_json_summary(gwas_path, models_dir, output_path, real_snp_base,
                     critical_post_values=[1, 0.99, 0.95, 0.90, 0.80, 0.7, 0.5, 0.2], sim_numbers=range(200), ignore_bad=False):
    """
    Creates a json file containing statistics for each model

    Args:
    `gwas_path` (str):
        The path to the gwas files
    `models_dir` (str):
        The path to the models
    `output_path` (str):
        The path to the output json file
    `real_snp_base` (str):
        The base path to the real snp locations
    `critical_post_values` (iterable, optional):
        The critical post. values to be used.
    `sim_numbers` (iterable, optional):
        The simulations to be mined
    `ignore_bad` (bool, optional):
        If true disregards simulations with terrible results. Defaults to False. 

    Returns: None   
    """

    models = os.listdir(models_dir)

    json_dict = {}
    for model in models:
        post_base = f'{models_dir}/{model}'
        all_post_vals, all_casual, all_pvals \
            = combine_genome_sample_r(gwas_path, post_base, real_snp_base, sim_numbers, ignore_bad=ignore_bad)

        post_discoverys \
            = [restricted_true_discovery_rate(all_post_vals, all_casual, post_start)
               for post_start in critical_post_values]

        tp_fp_list \
            = [count_tp_fp(all_post_vals, all_casual, critical)
               for critical in critical_post_values]

        json_dict[model] = ({f'{restriction}_d_r': rate[0] for rate, restriction, in zip(
            post_discoverys, critical_post_values)})
        json_dict[model].update({f'{restriction}_tp': tp_fp[0]
                                for tp_fp, restriction, in zip(tp_fp_list, critical_post_values)})
        json_dict[model].update({f'{restriction}_fp': tp_fp[1]
                                for tp_fp, restriction, in zip(tp_fp_list, critical_post_values)})

    gwas_discoverys \
        = [restricted_true_discovery_rate(-np.log(all_pvals), all_casual, post_start)
            for post_start in critical_post_values]
    json_dict['gwas'] = ({f'{restriction}_d_r': rate[0] for rate, restriction, in zip(
        gwas_discoverys, critical_post_values)})

    with open(output_path, 'w') as f:
        f.write(json.dumps(json_dict))


def model_compare_figures(json_path, output, statistics=None):
    """
    From a json file containing the results of multiple models, creates a bar chart for each statistic

    Args:
    `json_path` (str):
        The path to the json file
    `output` (str):
        The output path
    `statistics` (iterable, optional):
        The statistics to be plotted. Defaults to None, in which case all statistics are plotted.
    
    Returns: None
    """

    json_dict = json.load(open(json_path, 'r'))
    models = list(json_dict.keys())

    if statistics == None:
        # assuming all models have the same keys
        statistics = set.intersection(
            *[set(json_dict[model].keys()) for model in models])
    for statistic in statistics:
        model_values = []
        for model in models:
            model_values.append(json_dict[model][statistic])

        model_values = [0 if x is None else x for x in model_values]
        plt.bar(x=models, height=model_values)
        plt.xlabel('models')
        plt.ylabel(statistic)
        plt.xticks(rotation=300)
        plt.tight_layout()
        plt.savefig(f'{output}_{statistic}.png')
        plt.clf()


def plot_post_probs(ax, probs, casual_i, title=None):
    """Generates a plot of the post. probabilities with the casual marker highlighted

    Args:
        ax (plt.ax): matplotlib axis
        probs (np.array): Array of post. probabilities
        casual_i (np.array): The indicies of the casual marker
        title (str, optional): Figure title. Defaults to None.
    """

    ax.scatter(np.arange(len(probs))[
               casual_i], probs[casual_i], s=CASUAL_SIZE, marker='*', c='red')
    ax.scatter(np.delete(np.arange(len(probs)), casual_i),
               np.delete(probs, casual_i), s=SIZE)
    ax.set_ylabel("Post. Prob")

    if title is not None:
        ax.set_title(title)


def plot_gwas(ax, p_vals, casual_i, p_critical):
    """
    Plot the GWAS results.

    Parameters:
        ax (matplotlib.axes.Axes): The axes object to plot on.
        p_vals (numpy.ndarray): The p-values for each marker.
        casual_i (int): The index of the causal marker.
        p_critical (float): The critical p-value threshold.

    Returns:
    None
    """

    inf_positions = np.isinf(p_vals)
    # Find the maximum of non-infinite elements
    max_value = np.max(p_vals[~inf_positions])
    # Replace inf with max_value + 100
    p_vals[inf_positions] = max_value + 100

    ax.scatter([casual_i], p_vals[casual_i],
               s=CASUAL_SIZE, marker='*', c='red')
    ax.scatter(np.arange(len(p_vals))[
               inf_positions], p_vals[inf_positions], s=CASUAL_SIZE, marker='s', c='green')


    left_over = np.concatenate([casual_i, np.arange(
        len(p_vals))[inf_positions]])
    ax.scatter(np.delete(np.arange(len(p_vals)), left_over), np.delete(p_vals, left_over), s=SIZE)
    ax.axhline(y=-np.log(p_critical), color='black', linestyle='-')
    ax.set_xlabel("Marker")
    ax.set_ylabel("-Log p wald")


def make_plots_for_genome_region(gwas_path, probs_path, output_path, real_snp_location_df, p_critical=5*10**(-8)):
    """
    Used for fine mapping where one genomic region is plotted at a time. Makes two subplots, one for the posterior probabilities 
    and one for the GWAS.

    Parameters:
        gwas_path (str): 
            The file path to the GWAS data.
        probs_path (str): 
            The file path to the posterior probabilities data.
        output_path (str):
            The file path to save the generated plots.
        real_snp_location_df (pandas.DataFrame): 
            The DataFrame containing the real SNP locations.
        p_critical (float, optional):
            The critical p-value threshold for plotting GWAS results. Defaults to 5*10**(-8).

    Returns: None
    """
    casual_i = get_casual_marker(real_snp_location_df)
    df = pd.read_csv(probs_path, sep=' ')
    probs = np.array(df[PROB_COLUMN])

    # Post. Prob Plot
    fig, axs = plt.subplots(2, 1)
    plot_post_probs(axs[0], probs, casual_i)

    # Basic GWAS
    p_vals = -np.log(get_p_value_gwas(gwas_path))
    assert len(p_vals) == len(probs)
    plot_gwas(axs[1], p_vals, casual_i, p_critical)

    plt.savefig(f'{output_path}.png')


def make_stacked_plots(gwas_path, models_dir, output_path, real_snp_location_base, i, subplot_height=4, subplot_width=10):
    """
    Used for fine mapping where the whole genome is plotted. Makes a stacked plot for each model, where the top plots are the
    posterior probabilities and the bottom plot is the GWAS.

    Parameters:
        gwas_path (str): 
            The file path to the GWAS data.
        models_dir (str): 
            The directory containing the models.
        output_path (str):
            The file path to save the generated plots.
        real_snp_location_base (str): 
            The base path to the real SNP locations.
        i (int): 
            The index of the simulation to plot.
        subplot_height (int, optional):
            The height of each subplot. Defaults to 4.
        subplot_width (int, optional):
            The width of each subplot. Defaults to 10.

    Returns: None
    """
    
    models = os.listdir(models_dir)

    # Create figure with dynamic height based on the number of models
    fig, axs = plt.subplots(
        len(models) + 1, 1, figsize=(subplot_width, subplot_height * (len(models) + 1)))

    casual_i = get_casual_marker(
        f'{real_snp_location_base}/snp_info_sample_{i}')

    for j, model in enumerate(sorted(models)):
        post_path = f'{models_dir}/{model}/sample_{i}.probs'
        post = get_post_probs(post_path)
        plot_post_probs(axs[j], post, casual_i, title=model)

    p_vals = p_vals = - \
        np.log(get_p_value_gwas(f'{gwas_path}/sample_{i}.pvalues'))
    plot_gwas(axs[-1], p_vals, casual_i, p_critical=P_CRITICAL)
    plt.savefig(f'{output_path}/sample_{i}')


def plot_region(center_marker, gwas_p_path, bayes_post_path, output_location, markers_examined=1000, point_size=4):
    """Used for genome wide fine mapping where plot the whole genome as well as one specific genomeic region"""

    fig, axs = plt.subplots(3, 1, figsize=(17, 8))
    casual_size = point_size + 5  # Larger size for casual point

    # Find Peak from GWAS
    df = pd.read_csv(gwas_p_path, sep=' ')
    p_vals = -np.log(np.array(df['P']))
    num_p_vals = len(p_vals)

    domain_start = int(max(0, center_marker - markers_examined/2))
    domain_end = int(min(len(p_vals), center_marker + markers_examined/2))

    # Basic GWAS
    p_critical = 5*10**(-8)
    marker_chrom = [int(loci[:loci.index(':')]) for loci in df['Predictor']]
    colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

    axs[0].scatter(np.arange(len(p_vals)), p_vals, c=colors, s=0.5)
    axs[0].axhline(y=-np.log(p_critical), color='black', linestyle='-')
    axs[0].set_ylabel("-Log p wald")
    axs[0].axvspan(domain_start, domain_end, color='green', alpha=0.5)

    p_vals = p_vals[domain_start:domain_end]
    center_marker_index = int((domain_end - domain_start)/2)
    # GWAS Plot
    axs[1].scatter(np.delete(domain_start+np.arange(len(p_vals)), center_marker_index),
                   np.delete(p_vals, center_marker_index), s=point_size)
    axs[1].scatter([center_marker], [p_vals[center_marker_index]],
                   s=casual_size, marker='*', c='red')
    axs[1].axhline(y=-np.log(p_critical), color='black', linestyle='-')
    axs[2].set_xlabel("Marker")
    axs[1].set_ylabel("-Log p wald")

    #### Post. Prob Plot####
    df = pd.read_csv(bayes_post_path, sep=' ')
    num_post_prob = len(df)
    df = df[domain_start:domain_end]
    probs = np.array(df[PROB_COLUMN])
    axs[2].scatter(np.delete(domain_start+np.arange(len(probs)), center_marker_index),
                   np.delete(probs, center_marker_index), s=point_size)
    axs[2].scatter([center_marker], [probs[center_marker_index]],
                   s=casual_size, marker='*', c='red')
    axs[2].set_ylabel("Post. Prob")

    plt.savefig(output_location)







#### combining genomes ###

def combine_phenotypes(phenotypes_path, output_path, samples):
    """
    Combines the phenotypes of multiple samples into a single file

    Args:
    `phenotypes_path` (str):
        The path to the phenotypes
    `output_path` (str):
        The path to the output file
    `samples` (iterable):
        The samples to be combined
    `num_phenotypes` (int, optional):
        The number of phenotypes to be combined. Defaults to 1.

    Returns: None
    """

    samples = list(samples)

    # accession ids
    output_df = pd.read_csv(f'{phenotypes_path}/sample_{samples[0]}', sep = ' ', header=None)[[0, 1]]




    # column 2 is the phenotype
    pheno_combine = np.zeros(len(output_df))
    for sample in samples:
        pheno = np.array(pd.read_csv(f'{phenotypes_path}/sample_{sample}', sep=' ', header=None)[2])
        # normalise pheno
        pheno_combine += (pheno - np.mean(pheno))/np.std(pheno)

    pheno_combine /= len(samples)

    output_df[3] = pheno_combine

    output_df.to_csv(output_path, sep=' ', header=None, index=None)