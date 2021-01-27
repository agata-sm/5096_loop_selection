#!c:/perl/bin/perl.exe


# 5096

# script to filter loops based on genomic distance and TAD / compartment / and CTCF


# input: output by annotate_loops.v.2.0.pl

# #region1:chr:start-end	region2:chr:start-end	pval_ES_ICE,pval_ES_KR	TAD_frag1	compartment_frag1	TAD_frag2	compartment_frag2	CTCF_frag1	-log10qval_frag1	-log10qval_frag2
#chr15:21110000..21120000	chr15:21400000..21410000	1.9145855385978147e-18,1.8089309717493912e-10	chr15:15246000:21400000	B
#	chr15:15246000:21400000	B
#	na	na	chr15:21407310:21407728	5.23457

# v 2.0
# tab delimited input
# includes filtering by distance to the nearest loop

##chr-reg1	start-reg1	end-reg1	chr-reg2	start-reg2	end-reg2	pval	loop_frag_distance	dist_to_chr_5end	dist_to_chr_3end	TAD_frag1	compartment_frag1	TAD_frag2	compartment_frag2	CTCF_frag1	CTCF-log10qval_frag1	CTCF_frag2	CTCF-log10qval_frag2	gene_frag1_coordinates	tgene_frag1_dist	gene_frag1_gene_id	gene_frag1_gene_name	gene_frag1_gene_biotype	gene_frag2_coordinates	gene_frag2_dist	gene_frag2_gene_id	gene_frag2_gene_name	gene_frag2_gene_biotype	distance_to_closestLoop_frag1	distance_to_closestLoop_frag2
#chr10	100870000	100880000	chr10	102150000	102160000	2.0244898783815442e-07	1270000	100870000	28534993	chr10:100597000:102370000	A	chr10:100597000:102370000	A	na	na	chr10:102154908:102155204	3.33570	chr10:100707534:100741635	-128366	ENSMUSG00000112104	Gm35722	lincRNA	chr10:101681486:102391469	0	ENSMUSG00000019888	Mgat4c	protein_coding	-120001 250001




use warnings;
use strict;
use diagnostics;
use Getopt::Long;


#use Data::Dumper;

my $script_name="filter_loops.v.2.0.pl";




if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--infile: /path/to/infile output by annotate_loops.pl\n";
	print "--outfile: /path/to/outfile.txt\n";
	print "--distance minimal genomic distance to separate interacting regions\n";
	print "--loop_dist minimal genomic distance to another detected loop\n";
	print "--compartment A/B compartment of interacting regions; possible format A:B (region 1 in A and region 2 in B) or A (both regions in A)\n";
	print "Please note that the filtering options are hierarchical, i.e. distance needs to be defined to filter by loop_dist and distance and loop_dist need to be defined to filter by compartment.\n";
}else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'infile=s'		=>	\(my $infile),
		'outfile=s'		=>	\(my $outfile),
		'compartment=s'	=>	\(my $compartment),
		'distance=s'		=>	\(my $distance),
		'loop_dist=s'		=>	\(my $loop_dist)		
	) or die "Error in command line arguments";


	open(OUTFILE, ">",$outfile) or die "Cannot open outfile $outfile: $!";


	open(INFILE, "<",$infile) or die "Cannot open infile $infile: $!";

	while(<INFILE>){
		chomp $_;


		if($_=~m/\#/){
			print OUTFILE "$_\n";
		}else{
#chr10	100870000	100880000	chr10	102150000	102160000	2.0244898783815442e-07	1270000	100870000	28534993	chr10:100597000:102370000	A	chr10:100597000:102370000	A	na	na	chr10:102154908:102155204	3.33570	chr10:100707534:100741635	-128366	ENSMUSG00000112104	Gm35722	lincRNA	chr10:101681486:102391469	0	ENSMUSG00000019888	Mgat4c	protein_coding	-120001250001
			my $line_=$_;
			my @line=split/\t/;

			my $chr_r1=$line[0];
			my $start_r1=$line[1];
			my $end_r1=$line[2];
			my $chr_r2=$line[3];
			my $start_r2=$line[4];
			my $end_r2=$line[5];

			#print "$line_\n";

			if(defined $distance){

				if( ($start_r2-$start_r1>= abs($distance) ) | ($start_r2-$end_r1>= abs($distance) ) ){

					my $foo1=$start_r2-$start_r1;
					my $foo2=$start_r2-$end_r1;

					#print "$line_\n";
					#print "$foo1\t\t$foo2\n";

					if(defined $loop_dist){ #distance_to_closestLoop_frag1	distance_to_closestLoop_frag2 (fields 29 and 30)

						if( ( abs($line[28]) >= $loop_dist ) && ( abs($line[29]) >= $loop_dist ) ){

							#print "$line[28]\t\t$line[29]\t\t$loop_dist\n\n";



							if(defined $compartment){#TAD_frag1	compartment_frag1	TAD_frag2	compartment_frag2 (fields 11-14)

								if ( ($line[11] eq $compartment ) && ($line[13] eq $compartment) ){

									#print "$line[11]\t$line[13]\n";
										
									print OUTFILE "$line_\n";	
								}

							} else{
								print OUTFILE "$line_\n";
							}

						}

					}else{
						print OUTFILE "$line_\n";
					}

					


				}	
			}else{
				print OUTFILE "$line_\n";
			}

		}
	}



}
exit;


# #region1:chr:start-end	region2:chr:start-end	pval_ES_ICE,pval_ES_KR	TAD_frag1	compartment_frag1	TAD_frag2	compartment_frag2	CTCF_frag1	-log10qval_frag1	-log10qval_frag2
#chr15:21110000..21120000	chr15:21400000..21410000	1.9145855385978147e-18,1.8089309717493912e-10	chr15:15246000:21400000	B
#	chr15:15246000:21400000	B
#	na	na	chr15:21407310:21407728	5.23457

