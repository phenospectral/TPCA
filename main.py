"""
Supplementary Methods (Python code) associated with generation of Supplementary Data for Plant Cell paper
Zhang et al. (2024) on multi-omic analysis of the long-term (2-week) heat stress response of the plant species,
Setaria viridis (Poaceae): a popular model species for the study of C4 photosynthesis that is closely related to
the common cereal crop, Setaria italica known colloquially as "foxtail millet".

The following script was written by Dr Adam James Carroll of the Australian National University
for the purpose of analysing correlations between responses of mRNA transcripts and the proteins they
encode across different functional categories of genes/proteins.

The script could be readily re-used for other datasets by changing input files
"""

from tpca import TPCA

"""
read in the Excel file containing protein group fold change data (proteinGroups.txt output from
MaxQuant software with an added column of log base2 fold changes) the index_col here needs to 
contain a unique value for each row that matches exactly to the index_col value at the 
corresponding row in the transcript data table below
"""

tpca = TPCA()

tpca.get_protein_data('data/Peng_Setaria_Protein_Data_20240521.xlsx', 0)
tpca.get_transcript_data('data/Peng_Setaria_Transcript_Data_20240521.xlsx', 1)
tpca.prepare()

options = {'write_html': True}

results = tpca.analyse_TP_fbins(["all"], options)
app = tpca.plot(options)












