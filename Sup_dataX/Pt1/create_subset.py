##This script extracts sequences from a supplied fasta file.
##Sequence IDs are supplied in a *.text file
##Selected secquences will be writen to a new fasta file


def ExtractFastaSeq (List, Fasta, OutFile):
    """List must be a list of names matching the names in Fasta.
    Fasta is a Fasta file, sequences from Fasta will be extracted if they are in List
    OutFile is a newly created .fasta file"""
    from Bio import SeqIO 
    FastaFh=open(Fasta)
    Fh = open(List)
    SeqList = []
    for names in [l.strip() for l in Fh.readlines()]:
		SeqList.append(names)
    print SeqList #gives a dictionary with all names as keys and all clusters as values
    
    # now loop write out the sequence for each fasta
    Towrite=[seq for seq in SeqIO.parse(FastaFh,'fasta')
    			if seq.name in SeqList]
    with open(OutFile, 'w') as f:
        SeqIO.write(Towrite, f, 'fasta')
        print '%s sequences added to fasta'%(len(Towrite))
    	print '%s sequences missing'%(len(SeqList)-len(Towrite))
    
    	

List='Step1_Proteomics_List.text'
Fasta='Step2_Tom_Pc.fasta'
OutFile='Step3_Proteomics.fasta'

ExtractFastaSeq(List, Fasta, OutFile)