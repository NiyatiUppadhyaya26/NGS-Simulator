# NGS-Simulator
The program takes an input fasta file; first it will divide input file randomly according to the number provided by user and again these subsequences are divided into 50bp length. Then file is aligned with reference genome  using NGS Tools.

Run script providing inputs as: 
python NGS_simulator.py hg.genome 1000 1000_file.fasta 50bp_file.fasta 1000_file.fastq Pan_troglodytes.genome
