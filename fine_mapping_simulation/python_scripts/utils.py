import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pandas as pd
from random import randint
from tqdm import tqdm
from bed_reader import open_bed, sample_file
from scipy.stats import mode
import os, sys, getopt, argparse, random, statistics


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
POST_95_TP =  'num_true_possitive_post95'
POST_95_FP =  'num_false_possitive_post95'
P_CRITICAL = 5*10**(-8)
POST_CRITICAL = 0.95



##############################################
# To help with sampling subsets of the genom #
##############################################

def random_indivdual_subset(number,fam,output,seed=42):
    """
    Creates a plink formated fam file which specifies a random
    subset of individuals used for a anylisis.
    
    Inputs:
    `number`: str
        The number of individuals in the subset
    `fam`: str 
        The fam file from which indivduals are sampled
    `output`: str
        The location that the new fam file should be written to 
    `seed`: int
        A random seed such that generation of results can occur deterministically.
    Outputs: None
    """
    random.seed(seed)
    accession_numbers = pd.read_csv(fam,header=None,sep='\t')
    indx = np.array(sorted(random.sample(list(range(len(accession_numbers))), number)))
    # column 0 and column 1
    accession_numbers.iloc[indx][[0,1]].to_csv(output, header=None,sep='\t', index=None)

def is_overlap(p_bounds, c_bound):
    """
    Helper function to sample domains of the genome. Tests if a prospective
    bounds to sample the genome overlaps with past samples.
    Input
    `p_bounds`: [(int, int), ..., (int, int)]
        A list of two element tuples which specify the past bounds.
    `c_bound`: (int, int)
        A tuple of the prospective bound
    Output: bool
        whether `c_bound` overlaps with of in `p_bounds`
    """
    
    for p_bound in p_bounds:
        if (c_bound[CHROMOSOME] == p_bound[CHROMOSOME]
            and pd.Interval(p_bound[START], p_bound[FINISH]).overlaps(pd.Interval(c_bound[START], c_bound[FINISH]))):
            return True
    return False 

def sample_snp(previous_bounds, ukbb_bim, approx_chromosome_size, size=10**6):
    """
    Selects a random genome location of `size`. Also ensures that selected region has not before been selected
    and is thus not in previous_bounds
    """

    # selected a random snp from the bim
    snp = ukbb_bim.iloc[randint(0,len(ukbb_bim))]
    
    # extract snp infomation
    snp_chromosome = snp['chrom_number']
    snp_pos = snp['chrom_pos']

    # extract += 1mb from the randomly seleted SNP
    start = int(snp_pos - size/2)
    end = int(snp_pos + size/2 - 1)


    # if the the bound is outside the chromosome try again
    if (start < 0 or end > approx_chromosome_size[snp_chromosome - 1] 
            or is_overlap(previous_bounds, (snp_chromosome, start, end))):
        return sample_snp(previous_bounds, size=size)
    

    return snp_chromosome, start, end

def sample_genome_domains(sample_size, sample_number, genome_path, output_path):
    """
    Samples a random regions of the genome, and write a file speciefyinh the SNPs
    that are in that genomic region. The sample regions are such that they do not
    overlap

    Inputs:
    `sample_size`: int
        The size of the genomic region sampled (in basepairs)
    `sample_number`: int
        The number of genomic samples that should be generated
    `genome_path`: str
        The path to the genome binary (e.g. uk biobank bim)
    `output_path`: str
        The location at which the binary should be output
    Outputs: None
    """

    ### sample genome ###
    ukbb_bim = pd.read_csv(f"{genome_path}.bim", sep='\t', header=None)

    # making another column containing the position of each SNP on  a chromosome as an integer
    ukbb_bim['chrom_pos'] = [int(item[item.index(':') + 1:item.index('_')]) for item in ukbb_bim[1]]
    ukbb_bim.rename(columns = {0:'chrom_number'}, inplace = True)

    approx_chromosome_size = [max(ukbb_bim[ukbb_bim['chrom_number'] == chrom]['chrom_pos']) for chrom in range(1,NUM_HUMAN_AUTOSOMES + 1)]

    # select 200 1mb regions
    previous_bounds = []
    for i in tqdm(range(sample_number)):
        previous_bounds.append(sample_snp(previous_bounds, ukbb_bim, approx_chromosome_size, size=sample_size))
        with open(f'{output_path}_{i}', 'w') as f:
            f.write(f'{previous_bounds[i][CHROMOSOME]} {previous_bounds[i][START]} {previous_bounds[i][FINISH]} sample_{i+1}')
    
################################################
# To help with the of simulation of phenotypes #
################################################
def simulate_traits_for_genome_chunks(bfile, phenotypes_output_path, snps_info, 
                                      heritability, num_casual, sample):
    """
    Simulates traits for a small sample of the genome.

    Inputs:
    `bfile`: str
        The path to the binary genome file (.bim)
    `phenotypes_output_path`: str
        The path to the path where the phenotypes should be written
    `snps_info` : str
        The path to where infomation about which SNPs are casual
        is stored
    `heritability` : str
        The Heritability of the trait that is getting simulated
    `num_casual`
        The number of casual variants that drive the trait
    """

    # read in genomic regions
    ukbb_path = f'{bfile}_{sample}'
    bim_df = pd.read_csv(f'{ukbb_path}.bim', delimiter='\t', header = None)
    bim_df.columns = ['Chromosome code', 'Variant identifier', 'Position in morgans or centimorgans', 'Base-pair coordinate ', 'Allele 1', 'Allele 2']
    fam_df = pd.read_csv(f'{ukbb_path}.fam', delimiter=' ', header = None)

    # read plink into python to simulate phenotype

    bed = open_bed(f'{ukbb_path}.bed')
    X = bed.read()
    num_markers = X.shape[1]

    # impute NaN to the mode
    modes = mode(X, nan_policy='omit').mode
    for i in range(X.shape[1]):
        X[np.isnan(X[:,i])] = modes[i]


    if num_casual == 1:
        # one snp with effect size one so phenotypes start as being the same as allele state for 1 locus
        casual = np.array([num_markers // 2])
        beta = np.array([1])
        y_g =  X[:,casual] @ beta
    else:
        casual = np.random.randint(num_markers//10,num_markers-num_markers//10, size=num_casual)
        beta = np.random.normal(0,2,size=num_casual)
        y_g = X[:,casual] @ beta
        
    # add noise and rescale to predermined heritiability
    y = (y_g) + np.random.normal(0,statistics.variance(y_g)/heritability - statistics.variance(y_g), size= (len(y_g)))

    # save the results
    casual_info = bim_df.iloc[casual].copy()
    if len(casual) == 1:
        casual_info = pd.DataFrame(casual_info)

    casual_info.to_csv(f'{snps_info}/snp_info_sample_{sample}', sep='\t', header = True)
    phenotype_df = pd.DataFrame({'acc1': fam_df[0], 'acc2': fam_df[1], 'pheno': y})
    phenotype_df.to_csv(f'{phenotypes_output_path}/pheno_{sample}', sep=' ', header=None, index=None)


#######################################
# For the anylisis of simulation data #
#######################################

def get_post_probs(p_val_path):
    """Used to extract the Post. probabilties from a .probs file"""
    df = pd.read_csv(p_val_path, sep=' ')
    return np.array(df['Probability'])

def get_casual_marker(snp_info):
    """Used to extract the releveant info from a snp info file"""
    return np.array((pd.read_csv(snp_info, sep='\t')['snp_index_in_chunk']))

def get_p_value_gwas(gwas_path):
    """"Used to extract the P values from a LDAK gwas file"""
    return np.array(pd.read_csv(gwas_path, sep=' ')['P'])

def get_paths_for_chunk_1_snp(simulation_dir, i):
    path_to_pvals_dir = f'{simulation_dir}/fine_map_gwas'
    path_to_probs_dir = f'{simulation_dir}/fine_map_bayesr'
    path_to_casual_dir = f'{simulation_dir}/casual_snp_info'

    pvals_path = f'{path_to_pvals_dir}/gwas_sample_{i}.pvalues'
    post_path = f'{path_to_probs_dir}/bayesr_sample_{i}.probs'
    casual_path = f'{path_to_casual_dir}/snp_info_sample_{i}'


    pvals = get_p_value_gwas(pvals_path)
    post = get_post_probs(post_path)
    casual = get_casual_marker(casual_path)


    return pvals, post, casual

def compute_stats_chunk_1_snp(pvals,probs,casual, p_critical=P_CRITICAL, post_critical=POST_CRITICAL):
    return {
        LEADING_P_CASUAL: np.argmin(pvals) == casual,
        LEADING_POST_CASUAL: np.argmax(probs) == casual,
        SIG_P_TP: len(np.where(pvals[casual] <= p_critical)),
        SIG_P_FP: len(np.where(np.delete(pvals, casual) <= p_critical)[0]),
        POST_95_TP: len(np.where(probs[casual] >= post_critical)),
        POST_95_FP: len(np.where(np.delete(probs, casual) >= post_critical)[0])
    }

def count_tp_fp(all_post_vals, all_casual, post_critical):
    num_true_discoveries = len(np.intersect1d(np.where(all_post_vals >= post_critical)[0], all_casual))
    num_false_discoveries = len(np.intersect1d(np.where(all_post_vals >= post_critical)[0], np.delete(np.arange(len(all_post_vals)), all_casual))) 
    return num_true_discoveries, num_false_discoveries

def restricted_true_discovery_rate(all_post_vals, all_casual, post_critical):
        post_start = max(all_post_vals)
        num_true_discoveries, num_false_discoveries = count_tp_fp(all_post_vals, all_casual, post_start)
        restricted_post_discovery=(None, post_start)
        print(post_start)
        while post_start > 0.01:
            num_true_discoveries, num_false_discoveries = count_tp_fp(all_post_vals, all_casual, post_start)
            if num_true_discoveries/(num_true_discoveries + num_false_discoveries) >= post_critical:
                restricted_post_discovery=(num_true_discoveries/len(all_casual), post_start)
            post_start *= 0.98
        return restricted_post_discovery

def get_stats_chunk_1_snp(simulation_dir, output_dir, sim_numbers=range(20), critical_post_values=[1, 0.99, 0.95, 0.90, 0.80, 0.7], ignore_bad=False):
    sim_numbers = list(sim_numbers)

    stats = {k:[v] for k,v in compute_stats_chunk_1_snp(*get_paths_for_chunk_1_snp(simulation_dir, sim_numbers[0])).items()}

    # gets stats for each chunk
    for i in sim_numbers[1:]:
        one_sim_stats = compute_stats_chunk_1_snp(*get_paths_for_chunk_1_snp(simulation_dir, i))
        for k in stats.keys(): stats[k].append(one_sim_stats[k])


    # get 95% discovery rate
    all_casual = []
    all_post_vals = []
    all_pvals = []
    for i in sim_numbers:
        pvals, post, casual = get_paths_for_chunk_1_snp(simulation_dir, i)
        if ignore_bad and (sum (p >= 0.99 for p in post)) > 10*len(casual):
            continue
        all_casual.extend(casual + len(all_pvals))
        all_post_vals.extend(post)
        all_pvals.extend(pvals)

    # using 95% discovery rate
    all_post_vals = np.array(all_post_vals)
    all_casual = np.array(all_casual)
    

    post_discoverys = [restricted_true_discovery_rate(all_post_vals, all_casual, post_start) for post_start in critical_post_values]
    tp_fp_list = [count_tp_fp(all_post_vals, all_casual, critical) for critical in critical_post_values]


    df = pd.DataFrame(stats)
    df.to_csv(f'{simulation_dir}/results.csv')

    ignore_bad_suffix = '_ignore_bad' if ignore_bad else ''
    with open(f'{simulation_dir}/summary_results{ignore_bad_suffix}.csv', 'w') as f:
        f.write("Number of times leading p_val is casual " + str(sum(df[LEADING_P_CASUAL]))+ '\n')
        f.write("Number of times leading post. is casual " + str(sum(df[LEADING_POST_CASUAL])) + '\n')
        for critical, tp_fp, post_discovery in zip(critical_post_values, tp_fp_list, post_discoverys):
            percent_var_detected, cutoff = post_discovery
            if percent_var_detected == None:
                percent_var_detected = 0
            f.write(f"{(critical*100):>6.2f}% true possitive discovery rate: {percent_var_detected*100:>3.2f}% of variants detected with post. cut of off {cutoff:>.2f}\n")
        for critical, tp_fp, post_discovery in zip(critical_post_values, tp_fp_list, post_discoverys):
            symbol = '=' if critical == 1 else '>'
            f.write(f"treating post {symbol} {critical:<3.2f} as casual: TP {tp_fp[0]:3d} FP {tp_fp[1]:3d}\n")

def make_plots_for_genome_region(gwas_path,bayesr_path,output_path,real_snp_location_df):
    """
    Used for fine mapping where one geneomic region is plotted at a time
    """
    casual_i = get_casual_marker(real_snp_location_df)
    size = 3
    casual_size = 8

    df = pd.read_csv(bayesr_path, sep=' ')
    probs = np.array(df['Probability'])

    marker_chrom = [int(loci[:loci.index(':')]) for loci in df['Predictor']]
    colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

    #### make a figure with GWAS and Post.Prob ####

    # Post. Prob Plot 
    fig, axs = plt.subplots(2, 1)
    axs[0].scatter(np.arange(len(probs))[casual_i], probs[casual_i], s=casual_size, marker = '*', c='red')
    axs[0].scatter(np.delete(np.arange(len(probs)), casual_i), np.delete(probs, casual_i), s=size)

    axs[0].set_ylabel("Post. Prob")


    # Basic GWAS
    p_critical = 5*10**(-8)
    p_vals = -np.log(get_p_value_gwas(gwas_path))

    inf_positions = np.isinf(p_vals)

    # Step 2: Find the maximum of non-infinite elements
    max_value = np.max(p_vals[~inf_positions])

    # Step 3: Replace inf with max_value + 100
    p_vals[inf_positions] = max_value + 100

    axs[1].scatter([casual_i],p_vals[casual_i], s=casual_size, marker = '*', c='red')
    axs[1].scatter(np.delete(np.arange(len(probs)), casual_i),np.delete(p_vals, casual_i), s=size)
    axs[1].axhline(y=-np.log(p_critical), color='black', linestyle='-')

    axs[1].set_xlabel("Marker")
    axs[1].set_ylabel("-Log p wald")
    plt.savefig(f'{output_path}.png')

def plot_region(center_marker, gwas_p_path, bayes_post_path, output_location, markers_examined = 1000, point_size = 4):


    """Used for genome wide fine mapping where plot the whole genome as well as one specific genomeic region"""


    fig, axs = plt.subplots(3, 1, figsize=(17,8))
    casual_size = point_size + 5 # Larger size for casual point

    # Find Peak from GWAS
    df = pd.read_csv(gwas_p_path, sep=' ')
    p_vals = -np.log(np.array(df['P']))
    num_p_vals = len(p_vals)


    domain_start = int(max(0, center_marker - markers_examined/2))
    domain_end = int(min(len(p_vals), center_marker + markers_examined/2))

    # Basic GWAS
    p_critical = 5*10**(-8)
    marker_chrom = [int(loci[:loci.index(':')]) for loci in  df['Predictor']]
    colors = ['blue' if chrom % 2 == 0 else 'red' for chrom in marker_chrom]

    axs[0].scatter(np.arange(len(p_vals)),p_vals, c=colors, s=0.5)
    axs[0].axhline(y=-np.log(p_critical), color='black', linestyle='-')
    axs[0].set_ylabel("-Log p wald")
    axs[0].axvspan(domain_start, domain_end, color='green', alpha=0.5)

    p_vals = p_vals[domain_start:domain_end]
    center_marker_index = int((domain_end - domain_start)/2)
    # GWAS Plot
    axs[1].scatter(np.delete(domain_start+np.arange(len(p_vals)), center_marker_index), np.delete(p_vals,center_marker_index), s=point_size)
    axs[1].scatter([center_marker], [p_vals[center_marker_index]], s=casual_size, marker = '*', c='red')
    axs[1].axhline(y=-np.log(p_critical), color='black', linestyle='-')
    axs[2].set_xlabel("Marker")
    axs[1].set_ylabel("-Log p wald")

    #### Post. Prob Plot####
    df = pd.read_csv(bayes_post_path, sep=' ')
    num_post_prob = len(df)
    df = df[domain_start:domain_end]
    probs = np.array(df['Probability'])
    axs[2].scatter(np.delete(domain_start+np.arange(len(probs)), center_marker_index),  np.delete(probs,center_marker_index), s=point_size)
    axs[2].scatter([center_marker], [probs[center_marker_index]], s=casual_size, marker = '*', c='red')
    axs[2].set_ylabel("Post. Prob")

    print(num_post_prob, num_p_vals)
    plt.savefig(output_location)