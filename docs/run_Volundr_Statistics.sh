#!/bin/bash
#Parameter file to run Völundr Synthetic Lethal module Statistics
#File generated 2019-12-22

python3 /full/path/to/Volundr.py --options_file /full/path/to/run_Volundr_Statistics.sh
exit

--TargetSearch	False
--Statistics	True

--Target_File	/full/path/to/Target_File.bed
--Working_Folder	/full/path/to/Working_Folder
--DataFiles	/full/path/to/Counts_Files/
--Master_Index_File	/full/path/to/Master_Indices_File.bed
--SampleManifest	/full/path/to/SampleManifest.txt
--Verbose	INFO
--Job_Name	ASW1
--Species	Mouse

--Control_Sample	Name of control sample from SampleManifest file
--Library_Control	Name of sample used to control for sgRNA diversity
--Bad_sgRNA_Percentile	2.5

# Null distribution upper and lower cutoffs and permutation count
--UpperPercentile	97.5
--LowerPercentile	2.5
--PermutationCount	100000

# False discovery for multiple sample correction.
--Alpha	0.1

# How many mismatches to use from Counts Files
--Target_Mismatch	1

# Choose which additional output files you would like
--Write_TDnorm_Log2_sgRNA_Control_File	False
--Write_TDnorm_Log2_sgRNA_Sample_File	False
--Write_Log2_sgRNA_File	False
--Write_Permuted_Log2_Data_File	False
