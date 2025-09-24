#!usr/bin/perl
use warnings;
#usage perl .pl <IN amino acid info> <IN dir path> > <OU >
($amino_info, $dir_path) = @ARGV;

open IN, "<$amino_info";;
while (<IN>) {
#	print "$_";
	chomp;
	($full,$one,$three,$carbon_total,$nitrogen_total,$carbon,$nitrogen) = split "\t";
	$carbon_total ++; #useless, just avoid warning
	$nitrogen_total ++; #useless, just avoid warning
	$three ++; #useless, just avoid warning
	$full ++; #useless, just avoid warning
	$carbon{$one} = $carbon;
	$nitrogen{$one} = $nitrogen;
}

print "genome\tC-ARSC\tN-ARSC\n";
$dir_path =~ s/\/$//;
$dir = "$dir_path/*";
my @filelist = glob( $dir );
foreach $file (@filelist) {
if($file =~ m/faa/){
#	print "$file\n";
	chomp $file;
	$carbon = 0;
	$nitrogen = 0;
	$total = 0;
	open IN, "$file";
	while (<IN>) {
		chomp;
		s/\*//g;
		if (m/>/) {next}
		@array = split "";
		foreach $char (@array) {
			$total ++;
			if(exists $carbon{$char}){
			$carbon = $carbon + $carbon{$char};
			}
			if(exists $nitrogen{$char}){
			$nitrogen = $nitrogen + $nitrogen{$char};
			}
		}
	}
close IN;
	$carbon /= $total;
	$nitrogen /= $total;
	$file =~ s{.+\/(.*)\.faa}{$1};
	print  "$file\t$carbon\t$nitrogen\n";
}
}

