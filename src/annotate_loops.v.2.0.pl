#!c:/perl/bin/perl.exe


# script to annotated chromatin loops (bedgraph) with A/B compartment (tab delimited chrom	start	end	Compartment) and CTCF peaks (narrowPeak)

# proj 5096

# v 2.0 added info on distance to the nearest region involved in loop formation

use warnings;
use strict;
use diagnostics;
use Getopt::Long;


#use Data::Dumper;

my $script_name="annotate_loops.v.2.0.pl";


sub check_overlap{
	my $region_loop=shift;
	my $region_to_test=shift;

	my ($chr_r,$coords)=split/:/,$region_loop;
	my ($start_r,$end_r)=split/\.{2,}/,$coords;
	my($tst_chr,$tst_start,$tst_end)=split/:/,$region_to_test;


	my $overlap="na";

	if ($chr_r eq $tst_chr){

		if ( ($start_r >= $tst_start) && ($end_r <= $tst_end) ){
			$overlap="yes";
		}elsif(  ($start_r <= $tst_start) && ($end_r >= $tst_start) ){
			$overlap="yes";
		}elsif( ($start_r <= $tst_end) && ($end_r >= $tst_end)  ){
			$overlap="yes";
		}

	}

	return($overlap);	
}




if ($ARGV[0]eq qw '-h'){
	print "please provide arguments for: \n perl $script_name [arguments]\n";
	print "arguments:\n";
	print "--loops: /path/to/loops.bedgraph\n";
 	print "--outfile: /path/to/outfile\n";
 	print "--compartments: /path/to/A_B.txt: tab delimited: chrom start end Compartment\n";
 	print "--ctcf: /path/to/ctcf.narrowPeak\n";
 	print "--genes_bedops: /path/to/output of BEDOPS closest-features --dist --closest against annotated genes\n";
 	print "--loops_bedops: /path/to/output of BEDOPS closest-features --no-overlaps --dist --closest against all regions forming loops\n";
 	print "--chrom_sizes: /path/to/chrom.sizes UCSC format\n";
 }


else{
	my $parameters=join(' ', @ARGV);
	print "parameters: $parameters\n";
	
	#commandline parsing for parameters
	GetOptions(
		'loops=s'		=>	\(my $infile_loops),
		'outfile=s'		=>	\(my $outfile),
		'compartments=s'	=>	\(my $infile_compartments),
		'ctcf=s'		=>	\(my $infile_ctcf),
		'genes_bedops=s'		=>	\(my $infile_genes_bedops),	
		'loops_bedops=s'		=>	\(my $infile_loops_bedops),	
		'chrom_sizes=s'		=>	\(my $infile_chrom_sizes)		
	) or die "Error in command line arguments";



	#####################
	#####################
	## chrom sizes
	#####################

	my %chrom_sizes;


	open (INFILE_CHS, "<","$infile_chrom_sizes") or die "Cannot open input file $infile_chrom_sizes: $!"; 

	while(<INFILE_CHS>){

		chomp $_; 
		my ($chr,$len)=split/\t/;
		$chrom_sizes{$chr}=$len;
	}

	close(INFILE_CHS);

	#####################
	#####################
	## A/B
	#####################

	my %compartments;

	open (INFILE_AB, "<","$infile_compartments") or die "Cannot open input file $infile_compartments: $!"; 

	while(<INFILE_AB>){

		#chomp $_; 
		$_ =~ s/[\r\n]+$//; #the file comes from exporting from MS Excel thus different EOL
		
		unless($_=~m/^chrom/){ #chrom	start	end	Compartment
			
			my @line=split /\t/; #chr1	0	4328000	B

			my $region="$line[0]:$line[1]:$line[2]";
			$compartments{$region}=$line[3];
		}


	}
	close(INFILE_AB);


	#####################
	#####################
	## ctcf
	#####################

	my %ctcf_peaks;

	open (INFILE_CTCF, "<","$infile_ctcf") or die "Cannot open input file $infile_ctcf: $!"; 

	while(<INFILE_CTCF>){

		chomp $_; 
		
		#chr start end name score strand signalValue pValue qValue peak
		my @line=split /\t/; #chr9	96136358	96136758	.	1000	.	12.54216	-1.00000	0.17286	200

		my $region="$line[0]:$line[1]:$line[2]";	
		$ctcf_peaks{$region}=$line[8]; #-log10qvalue pos 9	

	}
	close(INFILE_CTCF);

	
	#####################
	#####################
	## closest genes
	## chr8	18560000	18570000|chr8	18595130	18803189	gene_id "ENSMUSG00000039842"; gene_version "15"; gene_name "Mcph1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000020492"; havana_gene_version "3";	.	+|25131
	#####################

	my %closest_genes;

	open (INFILE_GENES, "<","$infile_genes_bedops") or die "Cannot open input file $infile_genes_bedops: $!"; 

	while(<INFILE_GENES>){

		chomp $_; 
		
		my @line=split /\|/;

		my @line_loop=split /\t/,$line[0];
		my @line_gene=split /\t/,$line[1];
		my $dist=$line[2];

		#print "@line_loop\n";

		my $loop_region="$line_loop[0]:$line_loop[1]..$line_loop[2]";
		my $gene_coords="$line_gene[0]:$line_gene[1]:$line_gene[2]";
		#my $feature_location=$line_gene[4];

		my @gene_info=split/;\s/,$line_gene[3];

		my $geneid=$gene_info[0];
		$geneid=~s/gene_id\s+//;
		$geneid=~s/"//g;

		my $gene_biotype=$gene_info[4];
		$gene_biotype=~s/gene_biotype\s+//;
		$gene_biotype=~s/"//g;
		$gene_biotype=~s/;//;

		my $gene_name=$gene_info[2];
		$gene_name=~s/gene_name\s+//;
		$gene_name=~s/"//g;

		my $closest_gene="$gene_coords\t$dist\t$geneid\t$gene_name\t$gene_biotype";

		$closest_genes{$loop_region}=$closest_gene;

	}
	close(INFILE_GENES);

	####################
	####################
	## closest loops
	## ES_mapq30_10kb.ICE.single_interactor.loop_regions.sorted.loopregion.bed
	## chr1	3630000	3640000|chr1	3660000	3670000|20001

	my %closest_loops;

	open (INFILE_LOOPS, "<","$infile_loops_bedops") or die "Cannot open input file $infile_loops_bedops: $!"; 

	while(<INFILE_LOOPS>){

		chomp $_; 
		
		my @line=split /\|/;

		my @line_loop=split /\t/,$line[0];
		my @line_loop_closest=split /\t/,$line[1];
		my $dist=$line[2];

		my $loop_region="$line_loop[0]:$line_loop[1]..$line_loop[2]";
		
		my $loop_closest_coords="$line_loop_closest[0]:$line_loop_closest[1]:$line_loop_closest[2]";

		my $closest_loop=$dist;


		$closest_loops{$loop_region}=$closest_loop;

	}
	close(INFILE_LOOPS);


	#####################
	#####################
	## loops
	#####################

	open (OUTFILE, ">","$outfile") or die "Cannot open output file $outfile: $!";


	open (INFILE_LOOPS, "<","$infile_loops") or die "Cannot open input file $infile_loops: $!"; 

	while(<INFILE_LOOPS>){

		chomp $_; 
			
		#header region1:chr:start-end	region2:chr:start-end	pval_ES_ICE,pval_ES_KR
		if ($_=~m/\#/){
			my $header_loops=$_;
			my $header1="#chr_frag1\tstart_frag1\tend_frag1\tchr_frag2\tstart_frag2\tend_frag2\tpval\tloop_frag_distance";
			my $header2="dist_to_chr_5end\tdist_to_chr_3end";
			my $header3="TAD_frag1\tcompartment_frag1\tTAD_frag2\tcompartment_frag2\tCTCF_frag1\tCTCF-log10qval_frag1\tCTCF_frag2\tCTCF-log10qval_frag2";
			my $header4="gene_frag1_coordinates\tgene_frag1_dist\tgene_frag1_gene_id\tgene_frag1_gene_name\tgene_frag1_gene_biotype";
			my $header5="gene_frag2_coordinates\tgene_frag2_dist\tgene_frag2_gene_id\tgene_frag2_gene_name\tgene_frag2_gene_biotype";
			my $header6="distance_to_closestLoop_frag1\tdistance_to_closestLoop_frag2";
			print OUTFILE "$header1\t$header2\t$header3\t$header4\t$header5\t$header6\n";

		}else{

			##chr-reg1	start-reg1	end-reg1	chr-reg2	start-reg2	end-reg2	pval
			#chr10	102470000	102480000	chr10	103130000	103140000	1.1588345613282038e-14
			my @line=split /\t/;  

			my $chr_r1=$line[0];
			my $start_r1=$line[1];
			my $end_r1=$line[2];

			my $chr_r2=$line[3];
			my $start_r2=$line[4];
			my $end_r2=$line[5];

			my $loop_pval=$line[6];

			my $line_loop1= join ("\t", @line);

			my $loop_frag_distance=$start_r2-$end_r1;

			my $chrlen=$chrom_sizes{$chr_r1};

			my $dist_toL=$start_r1;
			my $dist_toR=$chrlen-$end_r2;



			my $line_loop="$line_loop1\t$loop_frag_distance\t$dist_toL\t$dist_toR";

			my $region1="$chr_r1:$start_r1..$end_r1";
			my $region2="$chr_r2:$start_r2..$end_r2";


			#for each region find its TAD / compartment and ctcf binding status
			########################
			# TAD

			my $TAD_frag1;
			my $TAD_frag2;
			while ( (my ($TAD_region,$compartment)) = each (%compartments) ){

				my($tad_chr,$tad_start,$tad_end)=split/:/,$TAD_region;


				if ($tad_chr eq $chr_r1){
					my $overlap_1=check_overlap($region1,$TAD_region); #loop first because different delimiter

					if ($overlap_1 eq qw /yes/){
 						$TAD_frag1="$TAD_region\t$compartment";

					}
					
				}
				if ($tad_chr eq $chr_r2){
					my $overlap_2=check_overlap($region2,$TAD_region); #loop first because different delimiter

					if ($overlap_2 eq qw /yes/){
						#print "$overlap_2\t$TAD_region\t$region2\t";
 						$TAD_frag2="$TAD_region\t$compartment";

					}	
				}
			}

			unless(defined $TAD_frag1){
				$TAD_frag1="na\tna";
			}
			unless(defined $TAD_frag2){
				$TAD_frag2="na\tna";
			}

			my $line_TAD="$TAD_frag1\t$TAD_frag2";


			########################
			# CTCF

			my $CTCF_frag1;
			my $CTCF_frag2;
			while ( (my ($CTCF_region,$minuslog10qvalue)) = each (%ctcf_peaks) ){
				my($peak_chr,$peak_start,$peak_end)=split/:/,$CTCF_region;


				if ($peak_chr eq $chr_r1){
					my $overlap_1=check_overlap($region1,$CTCF_region); #loop first because different delimiter

					if ($overlap_1 eq qw /yes/){
 						$CTCF_frag1="$CTCF_region\t$minuslog10qvalue";

					}
					
				}
				if ($peak_chr eq $chr_r2){
					my $overlap_2=check_overlap($region2,$CTCF_region); #loop first because different delimiter

					if ($overlap_2 eq qw /yes/){
 						$CTCF_frag2="$CTCF_region\t$minuslog10qvalue";

					}	
				}
			}

			unless(defined $CTCF_frag1){
				$CTCF_frag1="na\tna";
			}
			unless(defined $CTCF_frag2){
				$CTCF_frag2="na\tna";
			}


			my $line_CTCF="$CTCF_frag1\t$CTCF_frag2";


			########################
			# annotated gene
			
			my $gene_f1=$closest_genes{$region1};
			my $gene_f2=$closest_genes{$region2};

			my $gene_frag1=$gene_f1;
			my $gene_frag2=$gene_f2;

			my $line_gene="$gene_frag1\t$gene_frag2";


			########################
			# closest loop

			my $closest_loop_f1=$closest_loops{$region1};
			my $closest_loop_f2=$closest_loops{$region2};

			my $line_closest_loop="$closest_loop_f1\t$closest_loop_f2";


			########################
			# final line

			my $line="$line_loop\t$line_TAD\t$line_CTCF\t$line_gene\t$line_closest_loop";
			
			#$line=~ y/\n//d;

			print OUTFILE "$line\n";
		}

	}
	close(INFILE_LOOPS);
	close(OUTFILE);



}
exit;


