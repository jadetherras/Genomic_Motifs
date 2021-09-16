#Genomic_Motifs
a programm to read and reconstruct some motifs in a genome sequence. 

Two fonctionnalities :
- find all the sequences with a score upper than a given level, when the user give the position-weight matrix.
- recreate the position-weight matrix when the user give a list of positions, the genome and the length of the motif.

Two extensions :
- compare sequences find with fonctionnality 1 with real data, using a bedgraph file given real affinity of the transcription factor with a special nucleotide.
- optimize the matrix find with fonctionnality 2 by iterate in the best score's sequence of the given sequences of interest, until the matrix stagnate. 

example fonctionnality 1 :

./Genomic_Motifs -F -f ../test/promoters.fasta -m ../test/DBP.mat -T 0 (-B ../test/test_bedgraph.bedgraph)
of
./Genomic_Motifs -F -f ../chr7.fa -m ../test/DBP.mat -T 1 (-B ../BMAL1_ZT06_selection.bedgraph)

(need our relative way to chr7.fa and BMAL1_ZT06_selection.bedgraph)
optionnal : -o to set the output file

exemple fonctionnality 2 :
./Genomic_Motifs -M -l 6 -b ../BMAL1_chr7.bed -f ../chr7.fa (-E)

optionnal : -W for the output of the sequence of interest and -s to set the name of the output file 


