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
    #sample="[A-Za-z_0-9]+",

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

MAF = config['maf']

#datasets to use 
CANID_SUBSET = ['dogs', 'wolves', 'dogwolf']

#admixture Ks to run
#ADMIX_K = range(2,10)

#paths to programmes
glimpse_chunk = config['glimpse_chunk']
glimpse_impute = config['glimpse_impute']
glimpse_ligate = config['glimpse_ligate']
glimpse_sample = config['glimpse_sample']
glimpse_concordance = config['glimpse_concordance']

##### Rules to include #####

#Concordance rules
include: "rules/ref_panel_concordance.smk"
include: "rules/ref_panel_concordance_only_dogs.smk"
include: "rules/GLIMPSE_concordance.smk"
include: "rules/GLIMPSE_concordance_only_dogs.smk"
include: "rules/GLIMPSE_concordance_HC.smk"
include: "rules/GLIMPSE_concordance_transversions.smk"
include: "rules/GLIMPSE_concordance_estimate_sites_per_MAF_bin.smk"
include: "rules/filter_MAF_INFO_concordance.smk"
include: "rules/validation_filtering_pseudohaploid.smk"
include: "rules/smartpca_concordance_pseudohaploid_HC_genotyped.smk"  
include: "rules/plink_ROH_validation.smk"
include: "rules/plink_ROH_concordance.smk"


#Full imputation rules
include: "rules/GLIMPSE_impute.smk"
include: "rules/filter_MAF_INFO_imputation.smk"
include: "rules/smartpca_imputed.smk"
include: "rules/plink_ROH_phased.smk"


#no filter concordance PCA:
#include: "rules/smartpca_concordance_pseudohaploid_no_filters.smk"

#maf specific concordance:
#include: "rules/GLIMPSE_concordance_per_pop_maf.smk"



#full imputation steps
#include: "rules/admixture_imputed.smk"



rule all:
    input:
        expand('{path}/output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_corr_{canid_subset}_accuracy_HC.png', sample=SAMPLE, chrom_con=CHROM_CON, info=config['info'], maf=MAF, canid_subset=CANID_SUBSET, path=config['path']),
        expand('{path}/output/GLIMPSE_concordance/plots/smartpca_concordance_PH_HC_genotyped/merged_ph_called_modern_imputed.{sample}_{chrom_con}_INFO_{info}_MAF_{maf}_corr_{canid_subset}.png', sample=SAMPLE, chrom_con=CHROM_CON, info=config['info'], maf=MAF, canid_subset=CANID_SUBSET, path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/plots/PCA/merged_phased_ref_panel.all_chr_MAF_{maf_cutoff}_recalibrated_INFO_{info}_dogs.png', maf_cutoff=config['maf_cutoff'], info=config['info'], path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_1_dogs_only.png', maf_cutoff=config['maf_cutoff'], info=config['info'], path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_1_dogs_only.png', maf_cutoff=config['maf_cutoff'], info=config['info'], path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_1_dogs_coeff_tests.tsv', info=config['info'], maf_cutoff=config['maf_cutoff'], path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_1_dogs_coeff_tests.tsv', info=config['info'], maf_cutoff=config['maf_cutoff'], path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_1_wolves_count_length_all_ROHs.png', info=config['info'], maf_cutoff=config['maf_cutoff'], path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/plots/ROHs/merged_phased_ref_panel.allchrom_MAF_{maf_cutoff}_INFO_{info}_transversions_hom_win_het_1_wolves_count_length_all_ROHs.png', info=config['info'], maf_cutoff=config['maf_cutoff'], path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_1_{canid_subset}_islands.bed', maf_cutoff=config['maf_cutoff'], info=config['info'], hom_win_het=config['roh']['homozyg_window_het'], canid_subset=CANID_SUBSET, path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_1_{canid_subset}-genes-unique_islands.bed', maf_cutoff=config['maf_cutoff'], info=config['info'], hom_win_het=config['roh']['homozyg_window_het'], canid_subset=CANID_SUBSET, path=config['path']),
        expand('{path}/output/GLIMPSE_imputation/ROH_islands_deserts/imputed_modern_window_bed_{maf_cutoff}_INFO_{info}_all_sites_hom_win_het_1_{canid_subset}_deserts_GO_terms.txt', maf_cutoff=config['maf_cutoff'], info=config['info'], hom_win_het=config['roh']['homozyg_window_het'], canid_subset=CANID_SUBSET, path=config['path']),

### SOS:  only going with het=1 for this figure, otherwise I get an error (since it does not find long ROHs for het=0 in the ancient samples!)
               