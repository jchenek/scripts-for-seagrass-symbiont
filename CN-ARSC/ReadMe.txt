#usage:
#path_to_dir: A directory composed of CDS (faa) files output by Prodigal.
perl C_N_amino_acid_calculation.pl <IN C_N_amino_acid_info.txt> <IN path_to_dir> > <OU C_N_ARSC>

#example
perl C_N_amino_acid_calculation.pl C_N_amino_acid_info.txt ./example_dir/ > example_CN_ARSC.tsv