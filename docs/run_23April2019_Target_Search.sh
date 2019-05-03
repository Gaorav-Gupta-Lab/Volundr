#!/bin/bash
#Parameter file to run Völundr Synthetic Lethal module Target_Search
#File generated 2019-04-24 15:11:45

# SLURM Commands
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=07-00:00:00
#SBATCH --mem=80g

python3 /nas/longleaf/home/dennis/scripts/Volundr/Volundr.py --options_file /pine/scr/d/e/dennis/23April2019/run_23April2019_Target_Search.sh
exit


--Target_Search	True
--Statistics	False

--FASTQ1	/pine/scr/d/e/dennis/23April2019/23April2019.fastq.gz
--Index_File	/pine/scr/d/e/dennis/23April2019/23April2019.bed
--Target_File	/nas/longleaf/home/dennis/Reference_Files/Targets/CRISPR/mDDR_sgRNA.bed
--Master_Index_File	/nas/longleaf/home/dennis/Reference_Files/Indicies/Ion_Indices.bed
--Working_Folder	/pine/scr/d/e/dennis/23April2019/

--Verbose	DEBUG
--Job_Name	23April2019
--Spawn	19
--RevComp	False
--Species	Mouse
--Analyze_Unknowns	False
--Delete_Demultiplexed_FASTQ	True
--Compress	False
--Target_Mismatch	1
--Min_Length	120
--Target_Length	20
--Target_Start	20
--Index_Mismatch	1
--Target_Padding	2
--Expected_Position	112
--AnchorSeq	AAACACCG
--AnchorMismatch	1
--AnchorStart	100
--AnchorStop	125
