Howden et al. TOMATO NUCLEAR PROTEOMICS
The two Galaxy pipelines can be accessed via 
http://ppserver/galaxy/u/remco/h/nuclear-proteomics ##### NEED UPDATING 
http://ppserver/galaxy/u/remco/h/nuclear-proteomics-v2 ##### NEED UPDATING

## Pt 1. Data preprocessing

1. List with protein IDs manually extracted from Proteomics output file "proteinGroups.txt"
	- Shortened ID created by taking only gene number, not description
	e.g. from long IDs in column replace "\s" for "\t" and select first part, remove rest  
	- All Phyca11, "CON" and "REV" sequences were removed
	- Saved as Step1_proteomics_List.text

2. Edit supplied fasta file
	- Galaxy pipeline pt1 steps 1 to 4 show how the fasta file was: 
	transformed to table, descriptions removed, turned back into fasta
	- File saved as Step2_Tom_Pc.fasta
	- first entry was wrong, manually altered after downloading
	
3. Create fasta file with all sequences from proteomics data
	- Run create_subset.py
	specify List and Fasta file
	- Step3_Proteomics.fasta created, 3653 sequences
	
	
## Pt 2. Identification of nuclear proteins	
	
4. Upload Step3_Proteomics.fasta and look for nuclear proteins, galaxy pipline pt2
	- Run on Galaxy: 
	Wolf PSort (plant), NLStradamus, Predict NLS and NoD ran to detect nuclear localisations
	- Concatenate lists, remove duplicate lines using count
	- Saved as 140716_Nucl-IDs.txt 2534 genes
	
5. Create Venn diagram, using Venny
	- Sequence IDs are copied from lists in Galaxy pipeline and entered in Venny software
	  Oliveros (2007)
	

## Pt 3. Data Filtering

6. Clear up proteinGroups file. 
	- Use R script: Filter_proteinGroups.R
	requires Step1_proteomics_List.txt (file with additional shortened IDs), 
	         Step4_Nuclear_IDs.txt
	- Run script to: 
		Select only Nuclear IDs: output written to Step6_proteinGroups_ID_nucl.txt
		>> This file is also supplied as "Supplementary file 5.txt"
	    Simplify the table, keep LFQ columns only: written to Step6_proteinLFQs.txt
	    additional files written out: 
	    Step6_proteinGroups_ID_nonucl.txt  
	    >> This file is also supplied as "Supplementary file 6.txt",
	    Step6_proteinLFQ_nonucl.txt

	NOTE: These data files were identified in Perseus (for proteins without missing values)
	or as described below (for proteins with missing values)	    

7. Detect Proteins that are Present/Absent between treatments
	- Use R script: Detect_PresenceAbsence.R
	requires Step6_proteinLFQs.txt 
	- Creates four files proteins UP/DOWN by presence absence at 8 and 24 hour
	Step7_eightDOWN.txt, Step7_eightUP.txt, Step7_twfrDOWN.txt, Step7_twfrUP.txt
	- Create Venn Diagram using protein names for each list.

8. Detect Proteins that are 2 fold upregulated in at least one repeat and show presence/
   absence in the other
   - Use R script: CalculateRatios.R
   requires Step6_proteinLFQs.txt 
	- Creates four files proteins UP/DOWN by at 8 and 24 hour
	Step8_eightUP.txt, Step8_eightDOWN.txt, Step8_twfrUP.txt, Step8_twfrDOWN.txt
	- Create Venn Diagram using protein names for each list (added to list of Step7)
	
	
	>>> files for step 7 and step 8: 	8 hours are combined in "Supplementary file 9.txt"
										24 hours are combined in "Supplementary file 10.txt"
										
## Pt 4. Data filtering, non-nuclear comparisons

9.  Comparisons for non nuclear proteins, as done in 7 
	- Use R script: Detect_PresenceAbsence.R
	requires Step6_proteinLFQs_nonucl.txt 
	- Creates four files proteins UP/DOWN by presence absence at 8 and 24 hour
	Step9_eightDOWN.txt, Step9_eightUP.txt, Step9_twfrDOWN.txt, Step9_twfrUP.txt

10. Detect Proteins that are 2 fold upregulated in at least one repeat and show presence/
   absence in the other, as done in 8
   - Use R script: CalculateRatios.R
   requires Step6_proteinLFQs_nonucl.txt 
	- Creates four files proteins UP/DOWN by at 8 and 24 hour
	Step10_eightUP.txt, Step10_eightDOWN.txt, Step10_twfrUP.txt, Step10_twfrDOWN.txt
	
	>> files for step 9 and step 10: 	8 hours combined in "Supplementary file 17.txt"
										24 hours combined in "Supplementary file 18.txt"
										
										
## Pt 5. Data filtering NI8 vs NI24 and I8 vs I24

11.  Comparisons as done in 7 
	- Use R script: Detect_PresenceAbsence.R
	requires Step6_proteinLFQs.txt, Step6_proteinLFQs_nonucl.txt  
	- Creates eight files proteins UP/DOWN by presence absence at within I and within NI for
	predicted nuclear proteins and four for not predicted nuclear
	Step11_NI_DOWN.txt, Step11_NI_UP.txt, Step11_I_DOWN.txt, Step11_I_UP.txt	
	Step11_NI_DOWN_nonucl.txt, Step11_NI_UP_nonucl.txt, Step11_I_DOWN_nonucl.txt, Step11_I_UP_nonucl.txt

12. Detect Proteins that are 2 fold upregulated in at least one repeat and show presence/
   absence in the other, as done in 8
   - Use R script: CalculateRatios.R
   requires Step6_proteinLFQs.txt, Step6_proteinLFQs_nonucl.txt  
	- Creates eight files proteins UP/DOWN by Infected / Noninfected and nuclear/not nuclear.
	Step12_NI_DOWN.txt, Step12_NI_UP.txt, Step12_I_DOWN.txt, Step12_I_UP.txt
	Step12_NI_DOWN_nonucl.txt, Step12_NI_UP_nonucl.txt, Step12_I_DOWN_nonucl.txt, Step12_I_UP_nonucl.txt
	
	>> Files for Step 11 and 12:	Nuclear:
									Infected samples: combined in "Supplementary File 11.txt"
									Non-Infected: combined in "Supplementary File 12.txt" 
									Non Nuclear: 
									Infected samples combined in "Supplementary File 19.txt"
									Non-Infected: combined in "Supplementary File 20.txt"