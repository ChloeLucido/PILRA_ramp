#!/usr/local/bin/bash

#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=128G
#SBATCH --mail-user=gageblack@byu.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

PLINK_BIN=/fslhome/gagesb/fsl_groups/fslg_KauweLab/compute/src/plink-1.9/plink
SOURCE_DIR=/fslhome/gagesb/fsl_groups/fslg_KauweLab/compute/Heritability_2019/data
SOURCE_FILE=$SOURCE_DIR/adgc_hrc_merged_unrelated
COVAR=covar_noAPOE.txt
OUT=gwas_noAPOE

$PLINK_BIN --noweb --bfile $SOURCE_FILE  \
	--logistic --r2 dprime --adjust \
	--pheno pheno.txt \
	--pheno-name AD_status \
	--covar $COVAR \
	--covar-name sex,age,pc1_merged,pc2_merged,pc3_merged,pc4_merged,pc5_merged,pc6_merged,pc7_merged,pc8_merged,pc9_merged,pc10_merged \
	--out $OUT
