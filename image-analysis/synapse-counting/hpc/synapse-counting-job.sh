#!/bin/bash

#SBATCH --cluster=wice
#SBATCH --job-name="synapse_counting_pipeline"
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --time=04:00:00
#SBATCH --account=lp_big_wice_cpu
#SBATCH --partition=dedicated_big_bigmem
#SBATCH --cluster=wice
#SBATCH --array=1-4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=cydric.geyskens@kuleuven.be

echo JOB STARTED

# loading nextflow and checking version
module load Nextflow
nextflow -version

cd $VSC_SCRATCH/synapse-counting

# define parameters for each run based on the SLURM_ARRAY_TASK_ID
case $SLURM_ARRAY_TASK_ID in
    1)
        PROTEIN_AND_SYNAPTIC_MARKER="VCAM1_LacZ_VGLUT1_PSD95"
        INPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95"
        OUTPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp9_VCAM1_VGLUT1-PSD95_output_data"
        ;;
    2)
        PROTEIN_AND_SYNAPTIC_MARKER="VCAM1_LacZ_VGAT_GEPH"
        INPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp12_VCAM1_VGAT-GEPH"
        OUTPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp12_VCAM1_VGAT-GEPH_output_data"
        ;;
    3)
        PROTEIN_AND_SYNAPTIC_MARKER="GPR37L1_LacZ_VGLUT1_PSD95"
        INPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp10_GPR37L1_VGLUT1-PSD95"
        OUTPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp10_GPR37L1_VGLUT1-PSD95_output_data"
        ;;
    4)
        PROTEIN_AND_SYNAPTIC_MARKER="GPR37L1_LacZ_VGAT_GEPH"
        INPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp11_GPR37L1_VGAT-GEPH"
        OUTPUT_DIR="/scratch/leuven/357/vsc35768/synapse-counting/IHC_Exp11_GPR37L1_VGAT-GEPH_output_data"
        ;;
esac

# run the pipeline
nextflow run /scratch/leuven/357/vsc35768/synapse-counting/pipeline/run-pipeline.nf \
    --protein_and_synaptic_marker "$PROTEIN_AND_SYNAPTIC_MARKER" \
    --input_dir "$INPUT_DIR" \
    --output_dir "$OUTPUT_DIR" \
    --preprocessing_params "/scratch/leuven/357/vsc35768/synapse-counting/pipeline/preprocess_params_202511.json" \
    --optimization_params "/scratch/leuven/357/vsc35768/synapse-counting/pipeline/optimization_param_ranges_local_peaks_202511.json" \
    --name_segments "0 1 2 3 6 8 9" \
    --presynapse_channel "0" \
    --postsynapse_channel "1" \
    --nr_of_optimization_trials "100" \
    -profile hpc

echo JOB ENDED
