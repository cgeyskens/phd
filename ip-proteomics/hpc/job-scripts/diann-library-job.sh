#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --job-name="diann-spec_lib_zenotof_1miscleavage"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=00:10:00
#SBATCH --account=lp_big_wice_cpu
#SBATCH --partition=dedicated_big_bigmem
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cydric.geyskens@kuleuven.be

echo JOB STARTED

cd /scratch/leuven/357/vsc35768/ip-proteomics/

apptainer exec --bind /scratch/leuven/357/vsc35768/ip-proteomics:/ip-proteomics \
    diann-2.2.0.sif /diann-2.2.0/diann-linux \
    --threads 72 \
    --verbose 1 \
    --qvalue 0.01 \
    --out-lib "/ip-proteomics/mouse_ref_proteome_zenotof_1miscleavage_OxM_N-terAc_20251006_diann_lib.tsv" \
    --gen-spec-lib \
    --predictor \
    --fasta "/ip-proteomics/crap_fasta_20250304.fasta" \
    --cont-quant-exclude cRAP- \
    --fasta "/ip-proteomics/UP000000589_10090_20250930.fasta" \
    --fasta-search \
    --met-excision \
    --min-pep-len 7 \
    --max-pep-len 30 \
    --min-pr-mz 350 \
    --max-pr-mz 1250 \
    --min-pr-charge 1 \
    --max-pr-charge 4 \
    --cut K*,R* \
    --missed-cleavages 1 \
    --unimod4 \
    --var-mods 1 \
    --var-mod UniMod:35,15.994915,M \
    --var-mod UniMod:1,42.010565,*n \
    --peptidoforms \
    --relaxed-prot-inf 