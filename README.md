# UPDhmm


Description: ------------------
Steps for running UPDhmm package.
1. Read vcf file. IMPORTANT:
   a) Is QC have been already done?
   b) Specify samples names (father="x,mother="y,proband="z)
2. Use the function calcualte_events inside the pacakge.

This function will return a list with dataframes, one per chromosome with the blocks f predicted states. 

