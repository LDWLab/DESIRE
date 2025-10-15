#How to process the mitoribosomal protein data

To process the mitoribosomal protein data, make sure you have the following files:

```
-Files from 1_ALL_HITS (most recent version is v1.5)
-BINF_Project.py
-compound_organism_master_table.csv
-Assembly_Exeption_List_240523.csv
-Exeption_List_070823.csv
-Location_Checker_1.csv
```

Make sure you have the necessary files in the same directory

Remove headers for all csv files

Missing mS27,mS33,mS34 aln files


To batch run the files use the proteo_script.sh file which will process all the files at once, it will take around 2-3 hours

To add the files into the DESIRE SQL database, make sure you have the following files:

```
-upload_accession_hmm.py
-upload_aln_asp1.py
```

Log in to the DESIRE database and navigate to the alignment directory

Run the following lines of code (making sure to replace uL02m with the relevant mito-protein):

```
python3 upload_accession_hmm.py -c ./CSV/uL02m.csv
python3 upload_aln_asp1.py uL02m_aligned.fas e
```
