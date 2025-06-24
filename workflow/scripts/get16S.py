from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord
import sys
import os

gb_file = open("prokaryotes.gbff", "r")
orig_stdout = sys.stdout
outfile=open('18S.fa','w')
sys.stdout = outfile
for gb_record in SeqIO.parse(gb_file, "genbank"):
    # now do something with the record
    for feat in gb_record.features:
        if feat.type == "rRNA":
            #print(feat.gene)
            product = feat.qualifiers['product'][0]
            if '16S ribosomal RNA' in product:
                sequence = feat.extract(gb_record.seq)
                print('>'+gb_record.name+"_length_"+str(len(gb_record)))
                print(sequence)
sys.stdout = orig_stdout
outfile.close()

os.system("blastn -db /p/work1/gschuler/blastdb/SILVA_138.2_SSURef_NR99_tax_silva.fasta -query ./18S.fa -outfmt '6 qseqid stitle evalue' -max_target_seqs 1 >prok_18S_blast.out")
