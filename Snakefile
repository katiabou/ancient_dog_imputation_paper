import pandas as pd

##### load config file #####
configfile: 'config.yaml'

##### Define total chromosome numbers ##### 
CHROM = [f"chr{i}" for i in range(1,int(config['chromosome_number'])+1)]

##### set wildcard constraints ##### 
wildcard_constraints:
    chrom="chr\d+",
    chrom_con="chr\d+",
    info_cutoff="\d+.\d+",
    maf="\d+.\d+",
    maf_cutoff="\d+.\d+",
    info="\d+.\d+",
    canid_subset="\w+",
    cov_cutoff="\d+.\d+",
    hom_win_het="\d+"

#Bams used for the validation
samples_df = pd.read_table(config['bam_targets'], dtype=str, delimiter="\t").set_index("Original_ID", drop=False)
SAMPLE = list(samples_df['Original_ID'])

#Bams which will be imputed (might need to be updated from the online spreadsheet)
bams_df = pd.read_table(config['bam_imputation'], dtype=str, delimiter="\t").set_index("Sample", drop=False)
BAM = list(bams_df['Sample'])

#Define chromosome for accuracy check
CHROM_CON = ['chr1']

#Define coverage values to downsample target bams to
COVERAGE_VAL = ['2', '1', '0.5', '0.2', '0.1', '0.05']

#Define INFO score cutoffs for filtering imputed sites
INFO_CUTOFF = ['0.0','0.8', '0.9', '0.95']

#datasets to use 
CANID_SUBSET = ['dogs', 'wolves', 'dogwolf']

#paths to programmes
glimpse_chunk = config['glimpse_chunk']
glimpse_impute = config['glimpse_impute']
glimpse_ligate = config['glimpse_ligate']
glimpse_sample = config['glimpse_sample']
glimpse_concordance = config['glimpse_concordance']

##### Rules to include #####

#Rule for benchmarking
# include: "rules/ref_panel.smk"
# include: "rules/ref_panel_concordance_only_dogs.smk" 
include: "rules/GLIMPSE_concordance.smk" 
# include: "rules/GLIMPSE_concordance_only_dogs.smk" 
# include: "rules/GLIMPSE_concordance_HC.smk" 
# include: "rules/GLIMPSE_concordance_transversions.smk" 
# include: "rules/GLIMPSE_concordance_estimate_sites_per_MAF_bin.smk" 
# include: "rules/filter_MAF_INFO_concordance.smk" 
# include: "rules/validation_filtering_pseudohaploid.smk"
# include: "rules/smartpca_concordance_pseudohaploid_HC_genotyped.smk" 
# include: "rules/smartpca_concordance_pseudohaploid_no_filters.smk"  
# include: "rules/plink_ROH_validation.smk"  
# include: "rules/plink_ROH_concordance.smk" 


#Full imputation rules
# include: "rules/GLIMPSE_impute.smk"
# include: "rules/filter_MAF_INFO_imputation.smk"
# include: "rules/smartpca_imputed_v2.smk"
# include: "rules/plink_ROH_phased.smk"


rule all:
    input:
        #### output files for imputation benchmarking
        ### ROHs with ROHan
        #expand('{path}/output/GLIMPSE_concordance/plots/ROH_concordance_ROHan/{sample}_allchrom_INFO_{info}_MAF_{maf}_all_sites_hom_win_het_{hom_win_het}_accuracy_imputed_ROHan.png', sample=SAMPLE, coverage_val=COVERAGE_VAL, info=config['info_cutoff'], maf=config['maf_cutoff'], hom_win_het=config['roh']['homozyg_window_het'], path=config['path']),
        #new concordance plots:
        expand('{path}/output/GLIMPSE_concordance/plots/glimpse_concordance/rsquare_accuracy_allchrom_filtered-{sample}.png', sample=SAMPLE, path=config['path']),
        #expand('{path}/output/GLIMPSE_concordance/concordance_INFO_filtered_transversions/concordance_{sample}_allchrom_{coverage_val}x-INFO_{info_cutoff}_filtered_transversions.rsquare-mod.grp.txt.gz', sample=SAMPLE, coverage_val=COVERAGE_VAL, info_cutoff=INFO_CUTOFF,  path=config['path']),
        #expand('{path}/output/GLIMPSE_concordance/plots/glimpse_concordance_tranversions/concordance_allchrom_0.5x_1x_filtered-all_sites_transversions.png', path=config['path']),
        #expand('{path}/output/GLIMPSE_concordance_only_dogs/plots/glimpse_concordance/concordance_allchrom_0.5x_1x_filtered-all_sites_only_dogs.png', path=config['path'])
