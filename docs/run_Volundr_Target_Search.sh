#!/bin/bash
#Parameter file to run Völundr Synthetic Lethal module Target_Search
#File generated 2019-12-12

# SLURM Commands
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --time=07-00:00:00
#SBATCH --mem=60g

python3 /full/path/to/Volundr.py --options_file /full/path/to/run_Volundr_Target_Search.sh
exit

--Target_Search	True
--Statistics	False

--Working_Folder	/full/path/to/working/folder/
--FASTQ1	/full/path/to/fastq_file.fq.gz
--Target_File	/full/path/to/target_file.bed
--Master_Index_File	/full/path/to/Master_Index_File.bed
--SampleManifest	/full/path/to/SampleManifest.txt


--Verbose	INFO
--Job_Name	ASW1
--Spawn	
--Species	
--Platform	Illumina
--Analyze_Unknowns	False
--Delete_Demultiplexed_FASTQ	True
--Compress	False
--RevComp	False
--Target_Mismatch	1
--Min_Length	120
--Target_Length	20
--Target_Start	20
--Index_Mismatch	2
--Target_Padding	2
--Expected_Position	61
--AnchorSeq	AAACACCG
--AnchorMismatch	1
--AnchorStart	35
--AnchorStop	65
