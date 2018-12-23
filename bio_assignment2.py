#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Sat Nov 11 14:38:42 2017

@author: kotsabo
"""

# import required Biophython functions 
from Bio import Entrez
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import MuscleCommandline

# 1. answer

Entrez.email = 'A.N.Other@example.com'

# example is based on the blue cone opsin, human
my_protein = 'NP_000890.1' 

handle = Entrez.efetch(db = "protein", id = my_protein, rettype = "gb", retmode = "text")
record = SeqIO.read(handle, "genbank")
handle.close()

# show the sequence record
# here we have choosen the human KITLG gene
#print(record)

result_handle = NCBIWWW.qblast('blastp', 'swissprot', record.seq)
# This may take some time to run

# parse the results
result_handle.seek(0)
blast_record = NCBIXML.read(result_handle)

print('Gene name\te-value\tscore\taccession number')
for a in blast_record.alignments:
    if a.hsps[0].expect > 1e-1:
        break
    print(a.title.split('|')[4].split(' ')[0] + '\t' + str(a.hsps[0].expect) + '\t' + str( a.hsps[0].score) + '\t' + a.accession) 
    
# read the full ncbi entries
a = blast_record.alignments[0]
sp_ids = []
for a in blast_record.alignments:
    sp_ids.append(a.title.split('|')[3])
handle = Entrez.efetch(db = "protein", id = ",".join(sp_ids), retmode = "xml")
data = Entrez.read(handle)
species = []
print('Alignment score\tSpecies')
for i,d in enumerate(data):
    species.append(d['GBSeq_source'])
    print(str(blast_record.alignments[i].hsps[0].score) + '\t' + d['GBSeq_source'])
  

# now work with all results with e-value below this value:
E_VALUE_THRESH = 1e-6

# the following will write all results into a FASTA file for the MSA 
def get_seqrecs(alignments, threshold):
    # a little helper function to get the sequence records
    for i,aln in enumerate(alignments):
        for hsp in aln.hsps:
            if hsp.expect < threshold:
                id = aln.title.split('|')[4].split(' ')[0].split('_')[0]+'_'+species[i].replace(' ','_')[:9]
                print(id)
                yield SeqRecord(Seq(hsp.sbjct), id = id)
                break
 
best_seqs = get_seqrecs(blast_record.alignments, E_VALUE_THRESH)
# write out to a fasta file
SeqIO.write(best_seqs, 'family_alignment.fasta', 'fasta')


# run Muscle MSA
#cmdline = MuscleCommandline('./muscle3.8.31_i86linux64', input = 'family_alignment.fasta', 
#                            out = 'family_alignment.aln', clw = True)
#cmdline()

# 2

from Bio import AlignIO

alignment = AlignIO.read('family_alignment.aln', 'clustal')
print(alignment)

from Bio.Phylo.TreeConstruction import DistanceCalculator
calculator = DistanceCalculator('identity')
dm = calculator.get_distance(alignment)
print(dm)


from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
constructor = DistanceTreeConstructor(calculator, 'upgma')
tree = constructor.build_tree(alignment)
print(tree)


from Bio import Phylo
import matplotlib.pyplot as plt
fig = plt.figure(figsize = (12, 12))
ax = plt.subplot(111)
Phylo.draw(tree, axes = ax)
fig.savefig('tree.pdf')




