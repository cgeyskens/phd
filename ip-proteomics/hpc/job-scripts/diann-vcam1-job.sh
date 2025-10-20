#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --job-name="diann-v2.2.0_exp17_zenotof_1miscleavage-wo-gpr37l1"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --time=3:00:00
#SBATCH --account=lp_big_wice_cpu
#SBATCH --partition=dedicated_big_bigmem
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cydric.geyskens@kuleuven.be

echo JOB STARTED

cd /scratch/leuven/357/vsc35768/ip-proteomics/

apptainer exec --bind /scratch/leuven/357/vsc35768/ip-proteomics:/ip-proteomics \
    --bind /scratch/leuven/357/vsc35768/ip-proteomics/exp17-zenotof-diann-output-2miscleavage:/out \
    diann-2.2.0.sif /diann-2.2.0/diann-linux \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_3.mzML" \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_4.mzML" \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_5.mzML" \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_6.mzML" \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_7.mzML" \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_8.mzML" \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_9.mzML" \
    --f "/ip-proteomics/exp17-ms-zenotof-convert-output/2ul_CG_10.mzML" \
    --lib "/ip-proteomics/mouse_ref_proteome_zenotof_1miscleavage_OxM_N-terAc_20251006_diann_lib.predicted.speclib" \
    --threads 72 \
    --verbose 1 \
    --out "/out/exp17-zenotof-diann-output-1miscleavage.tsv" \
    --qvalue 0.01 \
    --matrices \
    --out-lib "/out/exp17-zenotof-diann-output-1miscleavage.parquet" \
    --temp "/out/" \
    --xic \
    --fasta "/ip-proteomics/crap_fasta_20250304.fasta" \
    --cont-quant-exclude cRAP- \
    --fasta "/ip-proteomics/UP000000589_10090_20250930.fasta" \
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
    --mass-acc 20.0 \
    --mass-acc-ms1 12 \
    --peptidoforms \
    --reanalyse \
    --relaxed-prot-inf \
    --rt-profiling
