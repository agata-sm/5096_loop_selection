#!c:/perl/bin/perl.exe


# script to obtain observed read counts for selected regions invloved in chromatin loops

# in
# bedgraph from hicDetectLoops
# chr1	18360000	18370000	chr1	18980000	18990000	5.362075554241188e-46



# out
# not interacting loops; i.e. pairs of loops which both regions do not interact with any of the regions which are included in the file (i.e. have another interactor)
# single interactors (within each chr)
# loop coordinates chr-tab-start-tab-end
# loop p values


# 2 scenarios

# (i)
# region 1.1 - region 2.1
# region 1.1 - region 2.2


# (ii)
# region 1.1 - region 2.1
# region 2.1 - region 2.2


# proj 5096
# 20 x 2020

# v 1.0. output format change to tab delimited
# 18xii2020

# author: Agata Smialowska


# v 0.4
# handles the double p value pval1,pval2 in files output by loops_compare.v.0.1.pl
#	#FILE ES_ICE and ES_KR#contig	start	end	contig	start	end	pval_file1,pval_file2
#	chr8	60500000	60510000	chr8	60590000	60600000	3.40849550800013e-20,2.071745244240804e-23



use warnings;
use strict;
use diagnostics;
use Getopt::Long;
use List::MoreUtils qw(any uniq);


#use Data::Dumper;

my $script_name="loop_selection.v.1.0.pl";


if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--loops: /path/to/loops.bedgraph\n";
 	print "--outfile: /path/to/outfile\n";
 }


else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'loops=s'		=>	\(my $infile_loops),
		'outfile=s'		=>	\(my $outfile)
	) or die "Error in command line arguments";


	################
	################
	### loops
	################

	my %loops; # HoH
	my %loop_regions;

	#this used if the loops to be compared are from the file whcih is a result of comparison of loops from different matix balancing methods
	my $file1_pref;
	my $file2_pref;

	#get first line
	open (INFILE_LOOPS, "<","$infile_loops") or die "Cannot open input file $infile_loops: $!"; 
	my $firstLine = <INFILE_LOOPS>; 
	close INFILE_LOOPS;


	#file from comparison of loops from different matix balancing methods by loops_compare.v.0.1.pl
	if($firstLine=~m/FILE\s(\w+)\sand\s(\w+)\#contig./){ #FILE ES_ICE and ES_KR#contig	start	end	contig	start	end	pval_file1,pval_file2 
		($file1_pref,$file2_pref)=($1,$2);
	}
	#loop only in one one of the files by loops_compare.v.0.1.pl
	elsif($firstLine =~m/FILE\s(\w+)\#contig./){
		$file1_pref=$1;
		$file2_pref="na";
	}
	elsif($firstLine =~m/\d+|chr/){
		print "$firstLine\n";
		$file1_pref=$infile_loops;
		$file2_pref="na";
	}



	#parse the loops


	open (INFILE_LOOPS, "<","$infile_loops") or die "Cannot open input file $infile_loops: $!"; 

	while(<INFILE_LOOPS>){

		chomp $_; 
		unless($_=~m/^\#/){
			my @line=split /\t/; #chr1	18360000	18370000	chr1	18980000	18990000	5.362075554241188e-46

			my $reg1="$line[0]:$line[1]::$line[2]";
			my $reg2="$line[3]:$line[4]::$line[5]";
			my $pval="$line[6]";
			my $reg_1_2="$reg1\t$reg2";

			$loop_regions{$reg1}=$reg2;
			$loop_regions{$reg2}=$reg1;

			$loops{$reg1}{$reg2}=$pval;		#$HoH{$who}{$key} = $value; 

		}
	}

	close(INFILE_LOOPS);



	#print Dumper \%loops;




	################
	################
	### get non-interacting loops
	################

	open (OUTFILE, ">", "$outfile") or die "Cannot open output file $outfile: $!";

	my $outfile_multi="$outfile\_nonunique_interactors.txt";
	open (OUTFILE_MULTI, ">", "$outfile_multi") or die "Cannot open output file $outfile_multi: $!";


	my $header;
	if ($file2_pref eq qw /na/){
		$header="\#chr-reg1\tstart-reg1\tend-reg1\tchr-reg2\tstart-reg2\tend-reg2\tpval";

	}elsif($file2_pref){
		my $pval_file1_header="pval_$file1_pref";
		my $pval_file2_header="pval_$file2_pref";
		$header="\#chr-reg1\tstart-reg1\tend-reg1\tchr-reg2\tstart-reg2\tend-reg2\t$pval_file1_header,$pval_file2_header";
	}
	print OUTFILE "$header\n";
	

	foreach my $region1 (sort keys %loops){

		my $interacting_regions2=scalar keys %{ $loops{$region1} };
		#print "elements: $interacting_regions2\n";


		if ($interacting_regions2 ==1){
			foreach my $region2 (keys %{ $loops{$region1} } ) {

				unless (exists $loops{$region2}) {
					#print "UNIQUE: $region1\t$region2: $loops{$region1}{$region2}\n";

					my $region1p = $region1 =~ s/:{1,2}/\t/gr;
					my $region2p = $region2 =~ s/:{1,2}/\t/gr;

					
					my $line="$region1p\t$region2p\t$loops{$region1}{$region2}";

					print OUTFILE "$line\n";

				}else{
					#print "interacting elsewhere: $region1\t$region2: $loops{$region1}{$region2}\n";

					my $region1p = $region1 =~ s/:{1,2}/\t/gr;
					my $region2p = $region2 =~ s/:{1,2}/\t/gr;

					print OUTFILE_MULTI "$region1\t$region2\t$loops{$region1}{$region2}\n";

				}
				
			}	
		}

		
	}



	close(OUTFILE);
	close(OUTFILE_MULTI);

exit;
}

