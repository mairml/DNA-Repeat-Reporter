# DNA-Repeat-Reporter
Find exact direct and inverted repeats in a given FASTA sequence

Note: This script was written and run using Python version 2.7.3. The program was additionally tested with and supported by Python 3.5.

For full explaination and example use of this script in finding large, inverted repeats, please see attached pdf.

==============Program Summary==============

This script finds exact direct or inverted nucleotide repeats in a input FASTA file. The program then outputs two files: a repeat report, and a FASTA file containing repeat regions.

The repeat report contains the following information:
1.	Size of the motif sequence
2.	Index or indices at which the motif occurs
3.	Sequence of the motif
4.	Index or indices at which the reverse complement sequence of the motif occurs, if it occurs
5.	Sequence of reverse complement motif, if applicable.

The second output is a FASTA file containing the sequences of all detected repeats, within the context of filler ‘N’ nucleotides to denote sequences not repeated within the genome. This was done to allow visualization of sequence motifs against the original sequence using simple text viewers, without the need for sequence alignment. 

NOTE: This program does not only return maximal repeats, but maximal repeats and their repeated substring elements. Please see attached pdf for further explanation. 

==============Program Use==============

From your command prompt:

path/to/script/repeat-reporter.py input_sequence.fasta min_window report_outfile.txt seq_outfile.txt

Where min_window is an integer value designating the minimum repeat size you would like to detect (example: 5)
