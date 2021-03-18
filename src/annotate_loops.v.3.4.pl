#!c:/perl/bin/perl.exe


# script to annotated chromatin loops (bedgraph) with A/B compartment (tab delimited chrom	start	end	Compartment) and CTCF peaks (narrowPeak)

# proj 5096

# v 2.0 added info on distance to the nearest region involved in loop formation

# v 3.0 add both closest genes 5' and 3' to each fragment
# v 3.1 add both closest loops 5' and 3' to each fragment
# v 3.2 correct treatment of cases where annotated gene is NA
# v 3.3 compare the closest gene distance to the loop length to correctly output the outside / inside the loop status
# v 3.4 add distance status (inside / outside the loop) annoation for closest loop

use warnings;
use strict;
use diagnostics;
use Getopt::Long;


#use Data::Dumper;

my $script_name="annotate_loops.v.3.4.pl";


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

	print "Chromosome sizes from file $infile_chrom_sizes read in.\n";

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

	print "TADs from file $infile_compartments read in.\n";


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

	print "CTCF peaks from file $infile_ctcf read in.\n";

	
	#####################
	#####################
	## closest genes
	## chr8	18560000	18570000|chr8	18595130	18803189	gene_id "ENSMUSG00000039842"; gene_version "15"; gene_name "Mcph1"; gene_source "ensembl_havana"; gene_biotype "protein_coding"; havana_gene "OTTMUSG00000020492"; havana_gene_version "3";	.	+|25131
	#####################

	my %closest_genes1;
	my %closest_genes2;

	open (INFILE_GENES, "<","$infile_genes_bedops") or die "Cannot open input file $infile_genes_bedops: $!"; 

	while(<INFILE_GENES>){

		chomp $_; 
		
		my @line=split /\|/;

		my @line_loop=split /\t/,$line[0];

		#in case there is no annotated gene, i.e. NA
		my @line_gene1;
		my $dist1=$line[2];
		if($line[1] eq "NA"){
			@line_gene1=qw(NA NA NA NA);
		}else{
			@line_gene1=split /\t/,$line[1];
		}


		my @line_gene2;
		my $dist2=$line[4];
		if ($line[3] eq "NA"){
			@line_gene2=qw(NA NA NA NA);
		}else{
			@line_gene2=split /\t/,$line[3];
		}



		#my @line_gene1=split /\t/,$line[1];
		#my $dist1=$line[2];
		#my @line_gene2=split /\t/,$line[3];
		#my $dist2=$line[4];

		#print "@line_loop\n";

		my $loop_region="$line_loop[0]:$line_loop[1]..$line_loop[2]";
		my $gene1_coords="$line_gene1[0]:$line_gene1[1]:$line_gene1[2]";
		my $gene2_coords="$line_gene2[0]:$line_gene2[1]:$line_gene2[2]";

		#my $feature_location=$line_gene[4];

		my $geneid1;
		my $gene1_biotype;
		my $gene1_name;

		if($line_gene1[3] eq "NA"){
			$geneid1="NA";
			$gene1_biotype="NA";
			$gene1_name="NA";
		}else{
			my @gene1_info=split/;\s/,$line_gene1[3];
			
			$geneid1=$gene1_info[0];
			$geneid1=~s/gene_id\s+//;
			$geneid1=~s/"//g;

			$gene1_biotype=$gene1_info[4];
			$gene1_biotype=~s/gene_biotype\s+//;
			$gene1_biotype=~s/"//g;
			$gene1_biotype=~s/;//;

			$gene1_name=$gene1_info[2];
			$gene1_name=~s/gene_name\s+//;
			$gene1_name=~s/"//g;
		}


		my $closest_gene1="$gene1_coords\t$dist1\t$geneid1\t$gene1_name\t$gene1_biotype";

		$closest_genes1{$loop_region}=$closest_gene1;


		my $geneid2;
		my $gene2_biotype;
		my $gene2_name;

		if($line_gene2[3] eq "NA"){
			$geneid2="NA";
			$gene2_biotype="NA";
			$gene2_name="NA";
		}else{
			my @gene2_info=split/;\s/,$line_gene2[3];
			
			$geneid2=$gene2_info[0];
			$geneid2=~s/gene_id\s+//;
			$geneid2=~s/"//g;

			$gene2_biotype=$gene2_info[4];
			$gene2_biotype=~s/gene_biotype\s+//;
			$gene2_biotype=~s/"//g;
			$gene2_biotype=~s/;//;

			$gene2_name=$gene2_info[2];
			$gene2_name=~s/gene_name\s+//;
			$gene2_name=~s/"//g;
		}

		my $closest_gene2="$gene2_coords\t$dist2\t$geneid2\t$gene2_name\t$gene2_biotype";

		$closest_genes2{$loop_region}=$closest_gene2;

	}
	close(INFILE_GENES);

	print "Closest genes from file $infile_genes_bedops read in.\n";


	####################
	####################
	## closest loops
	## ES_mapq30_10kb.ICE.single_interactor.loop_regions.sorted.loopregion.bed
	## chr1	3630000	3640000|chr1	3660000	3670000|20001

	my %closest_loops1;
	my %closest_loops2;

	open (INFILE_LOOPS, "<","$infile_loops_bedops") or die "Cannot open input file $infile_loops_bedops: $!"; 

	while(<INFILE_LOOPS>){

		chomp $_; 
		
		my @line=split /\|/;
		#print "@line\n";
		#print "$line[1]\n$line[2]\n\n";

		my @line_loop=split /\t/,$line[0];


		#print "line[1]='$line[1]'\n";

		my @line_loop_closest1;
		my $dist1=$line[2];
		if($line[1] eq "NA"){
			@line_loop_closest1=qw(NA NA NA);
		}else{
			@line_loop_closest1=split /\t/,$line[1];
		}


		my @line_loop_closest2;
		my $dist2=$line[4];
		if ($line[3] eq "NA"){
			@line_loop_closest2=qw(NA NA NA);
		}else{
			@line_loop_closest2=split /\t/,$line[3];
		}


		# my @line_loop_closest1;
		# my $dist1=$line[2];
		# if($line[1] ne /NA/){
		# 	@line_loop_closest1=split /\t/,$line[1];
		# }elsif($line[1] eq "NA"){
		# 	#print "loop1 line[1]='$line[1]'\n";
		# 	@line_loop_closest1=qw(NA NA NA);
		# 	#print "loop1: @line_loop_closest1\n";
		# }
	
		# my @line_loop_closest2;
		# my $dist2=$line[4];
		# if ($line[3] ne "NA"){
		# 	@line_loop_closest2=split /\t/,$line[3];
		# }elsif($line[3] eq "NA"){
		# 	#print "loop2 line[3]='$line[3]'\n";
		# 	@line_loop_closest2=qw(NA NA NA);
		# 	#print "loop2: @line_loop_closest2\n";

		# }


		my $loop_region="$line_loop[0]:$line_loop[1]..$line_loop[2]";
		
		my $loop_closest1_coords="$line_loop_closest1[0]:$line_loop_closest1[1]:$line_loop_closest1[2]";
		my $closest_loop1="$loop_closest1_coords\t$dist1";

		my $loop_closest2_coords="$line_loop_closest2[0]:$line_loop_closest2[1]:$line_loop_closest2[2]";
		my $closest_loop2="$loop_closest2_coords\t$dist2";


		$closest_loops1{$loop_region}=$closest_loop1;
		$closest_loops2{$loop_region}=$closest_loop2;

	}
	close(INFILE_LOOPS);

	print "Closest loops from file $infile_loops_bedops read in.\n";


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
			my $header4_1="gene5_frag1_coordinates\tgene5_frag1_dist\tgene5_frag1_gene_id\tgene5_frag1_gene_name\tgene5_frag1_gene_biotype\tgene5_frag1_position_status";
			my $header4_2="gene3_frag1_coordinates\tgene3_frag1_dist\tgene3_frag1_gene_id\tgene3_frag1_gene_name\tgene3_frag1_gene_biotype\tgene3_frag1_position_status";

			my $header5_1="gene5_frag2_coordinates\tgene5_frag2_dist\tgene5_frag2_gene_id\tgene5_frag2_gene_name\tgene5_frag2_gene_biotype\tgene5_frag2_position_status";
			my $header5_2="gene3_frag2_coordinates\tgene3_frag2_dist\tgene3_frag2_gene_id\tgene3_frag2_gene_name\tgene3_frag2_gene_biotype\tgene3_frag2_position_status";

			my $header6_1="coordinates_closestLoop5_frag1\tdistance_to_closestLoop5_frag1\tposition_status_closestLoop5_frag1\tcoordinates_closestLoop3_frag1\tdistance_to_closestLoop3_frag1\tposition_status_closestLoop3_frag1";
			my $header6_2="coordinates_closestLoop5_frag2\tdistance_to_closestLoop5_frag2\tposition_status_closestLoop5_frag2\tcoordinates_closestLoop3_frag2\tdistance_to_closestLoop3_frag2\tposition_status_closestLoop3_frag2";
			
			print OUTFILE "$header1\t$header2\t$header3\t$header4_1\t$header4_2\t$header5_1\t$header5_2\t$header6_1\t$header6_2\n";


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
			
			my $gene_f1_g1=$closest_genes1{$region1};
			my $gene_f1_g2=$closest_genes2{$region1};

			my $gene_f2_g1=$closest_genes1{$region2};
			my $gene_f2_g2=$closest_genes2{$region2};




			my @closest_gene1f1=split/\t/,$gene_f1_g1;
			my $dist_g1f1=$closest_gene1f1[1];

			my $position_status_g1f1="outside_the_loop";
			if($dist_g1f1 eq "NA"){
				$position_status_g1f1="NA";
			}else{
				if($dist_g1f1==0){
					$position_status_g1f1="overlapping_the_loop_contact";
				}
			}


			my @closest_gene2f1=split/\t/,$gene_f1_g2;
			my $dist_g2f1=$closest_gene2f1[1];

			my $position_status_g2f1="outside_the_loop";
			if($dist_g2f1 eq "NA"){
				$position_status_g2f1="NA";
			}else{
				if($dist_g2f1==0){
					$position_status_g2f1="overlapping_the_loop_contact";
				}elsif($dist_g2f1>0){
					if($dist_g2f1>$loop_frag_distance){
						$position_status_g2f1="outside_the_loop";
					}else{
						$position_status_g2f1="inside_the_loop";
					}
				}
			}

			my $gene_frag1="$gene_f1_g1\t$position_status_g1f1\t$gene_f1_g2\t$position_status_g2f1";
			

			my @closest_gene1f2=split/\t/,$gene_f2_g1;
			my $dist_g1f2=$closest_gene1f2[1];

			my $position_status_g1f2="outside_the_loop";
			if($dist_g1f2 eq "NA"){
				$position_status_g1f2="NA";
			}else{
				if($dist_g1f2==0){
					$position_status_g1f2="overlapping_the_loop_contact";
				}elsif($dist_g1f2<0){
					if(abs($dist_g1f2)>$loop_frag_distance){
						$position_status_g2f1="outside_the_loop";
					}else{
						$position_status_g1f2="inside_the_loop";
					}
				}
			}


			my @closest_gene2f2=split/\t/,$gene_f2_g2;
			my $dist_g2f2=$closest_gene2f2[1];

			my $position_status_g2f2="outside_the_loop";
			if($dist_g2f2 eq "NA"){
				$position_status_g2f2="NA";
			}else{
				if($dist_g2f2==0){
					$position_status_g2f2="overlapping_the_loop_contact";
				}
			}

			my $gene_frag2="$gene_f2_g1\t$position_status_g1f2\t$gene_f2_g2\t$position_status_g2f2";


			my $line_gene="$gene_frag1\t$gene_frag2";


			########################
			# closest loops

			
			my $closest_loop_f1_1=$closest_loops1{$region1};
			my $closest_loop_f1_2=$closest_loops2{$region1};

			my $closest_loop_f2_1=$closest_loops1{$region2};
			my $closest_loop_f2_2=$closest_loops2{$region2};







			my @closest_l1f1=split/\t/,$closest_loop_f1_1;
			my $dist_l1f1=$closest_l1f1[1];

			my $position_status_l1f1="outside_the_loop";
			if($dist_l1f1 eq "NA"){
				$position_status_l1f1="NA";
			}else{
				if(abs($dist_l1f1)==1){
					$position_status_l1f1="adjacent_to_the_loop_contact";
				}elsif($dist_l1f1==0){
					$position_status_l1f1="overlapping_the_loop_contact";
				}
			}


			my @closest_l2f1=split/\t/,$closest_loop_f1_2;
			my $dist_l2f1=$closest_l2f1[1];

			my $position_status_l2f1="outside_the_loop";
			if($dist_l2f1 eq "NA"){
				$position_status_l2f1="NA";
			}else{
				if(abs($dist_l2f1)==1){
					$position_status_l2f1="adjacent_to_the_loop_contact";
				}elsif($dist_l2f1==0){
					$position_status_l2f1="overlapping_the_loop_contact";
				}elsif($dist_l2f1>0){
					if($dist_l2f1>$loop_frag_distance){
						$position_status_l2f1="outside_the_loop";
					}else{
						$position_status_l2f1="inside_the_loop";
					}
				}
			}

			my $closest_loop_f1="$closest_loop_f1_1\t$position_status_l1f1\t$closest_loop_f1_2\t$position_status_l2f1";
			

			my @closest_l1f2=split/\t/,$closest_loop_f2_1;
			my $dist_l1f2=$closest_l1f2[1];

			my $position_status_l1f2="outside_the_loop";
			if($dist_l1f2 eq "NA"){
				$position_status_l1f2="NA";
			}else{
				if(abs($dist_l1f2)==1){
					$position_status_l1f2="adjacent_to_the_loop_contact";
				}elsif($dist_l1f2<0){
					if(abs($dist_l1f2)>$loop_frag_distance){
						$position_status_l1f2="outside_the_loop";
					}elsif($dist_l1f2==0){
						$position_status_l1f2="overlapping_the_loop_contact";
					}else{
						$position_status_l1f2="inside_the_loop";
					}
				}
			}


			my @closest_l2f2=split/\t/,$closest_loop_f2_2;
			my $dist_l2f2=$closest_l2f2[1];

			my $position_status_l2f2="outside_the_loop";
			if($dist_l2f2 eq "NA"){
				$position_status_l2f2="NA";
			}else{
				if(abs($dist_l2f2)==1){
					$position_status_l2f2="adjacent_to_the_loop_contact";
				}elsif($dist_l2f2==0){
					$position_status_l2f2="overlapping_the_loop_contact";
				}
			}

			my $closest_loop_f2="$closest_loop_f2_1\t$position_status_l1f2\t$closest_loop_f2_2\t$position_status_l2f2";


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


#chr_frag1	start_frag1	end_frag1	chr_frag2	start_frag2	end_frag2	pval	loop_frag_distance	dist_to_chr_5end	dist_to_chr_3end	TAD_frag1	compartment_frag1	TAD_frag2	compartment_frag2	CTCF_frag1	CTCF-log10qval_frag1	CTCF_frag2	CTCF-log10qval_frag2	gene5_frag1_coordinates	gene5_frag1_dist	gene5_frag1_gene_id	gene5_frag1_gene_name	gene5_frag1_gene_biotype	gene5_frag1_position_status	gene3_frag1_coordinates	gene3_frag1_dist	gene3_frag1_gene_id	gene3_frag1_gene_name	gene3_frag1_gene_biotype	gene3_frag1_position_status	gene5_frag2_coordinates	gene5_frag2_dist	gene5_frag2_gene_id	gene5_frag2_gene_name	gene5_frag2_gene_biotype	gene5_frag2_position_status	gene3_frag2_coordinates	gene3_frag2_dist	gene3_frag2_gene_id	gene3_frag2_gene_name	gene3_frag2_gene_biotype	gene3_frag2_position_status	coordinates_closestLoop5_frag1	distance_to_closestLoop5_frag1	position_status_closestLoop5_frag1	coordinates_closestLoop3_frag1	distance_to_closestLoop3_frag1	position_status_closestLoop3_frag1	coordinates_closestLoop5_frag2	distance_to_closestLoop5_frag2	position_status_closestLoop5_frag2	coordinates_closestLoop3_frag2	distance_to_closestLoop3_frag2	position_status_closestLoop3_frag2
#chr10	102490000	102500000	chr10	102820000	102830000	1.3715597984223548e-36	320000	102490000	27864993	chr10:102494000:102829000	B	chr10:102494000:102829000	B	chr10:102491741:102492208	5.23457	chr10:102823463:102823832	5.23457	chr10:102481755:102490486	0	ENSMUSG00000019890	Nts	protein_coding	overlapping_the_loop_contact	chr10:102512221:102549736	12222	ENSMUSG00000044921	Rassf9	protein_coding	inside_the_loop	chr10:102819319:102838659	0	ENSMUSG00000112046	Gm5175	transcribed_processed_pseudogene	overlapping_the_loop_contact	chr10:102998706:103030215	168707	ENSMUSG00000036602	Alx1	protein_coding	outside_the_loop	chr10:102150000:102160000	-330001	outside_the_loop	chr10:102820000:102830000	320001	outside_the_loop	chr10:102490000:102500000	-320001	outside_the_loop	chr10:102870000:102880000	40001	outside_the_loop

