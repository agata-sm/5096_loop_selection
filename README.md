Scripts for selection of chromatin loops from bedgraph files. Selected loops are formed by regions which are the sole interactors of each other.

These scripts can be used to select chromatin loops which comprise two regions which are the sole interactors of each other (i.e. do not interact with any other region).
The input is `bedgraph` file as output by `hicDetectLoops` from [hicexplorer](https://hicexplorer.readthedocs.io/en/latest/content/tools/hicDetectLoops.html).
The script to annotate the selected loops with functional features requires output of [BEDOPS](https://bedops.readthedocs.io/en/latest/).

The annotation script has been developed to include annotations available for data from GSE96107 (PMC5651218).


## Workflow description

Starting from the "raw" result of loop detection by  `hicDetectLoops`:

1. Loops in which both interacting regions are the sole interactors of each other are selected;

2. Selected loops from (1) are annotated using the following information: CTCF peaks, TAD, A/B compartment, closest genomic feature, distance to the closest region involved in chromatn loop formation.

3. The annotated list can be further filtered by genomic distance and compartment.


## Installation 


Firstly, you need `git` on your computer. A tutorial on how to install it can be found [here](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) or [here](https://github.com/git-guides/install-git).

To install the scripts simply type

```
git clone https://github.com/agata-sm/5096_loop_selection
```

Each script displays help message when option `-h` is used instead of other command line paramaters.

The scripts were tested under `perl (v5.26.2)`. There is no guarantee they will work under any other `perl` version, especially the older ones.


## Usage

### Loop selection

To select loops with single interactors:

```
Scripts="/path/to/src"

Loops="/path/to/sample.hicDetectLoops.bedgraph"
SingleInteractorLoops="/path/to/sample.single_interactor.txt"

perl $Scripts/loop_selection.v.1.0.pl --loops $Loops --outfile $SingleInteractorLoops
```

### Loop annotation

First, prepare the features for annotation:

* CTCF peaks: in `narrowPeak` format

* TAD and compartment: tab-delimited chromosome, start, end, compartmemt

* closest features (i.e. genes) **sample.single_interactor.loop_regions.sorted.gene.bed**

```
perl -ne 'chomp $_; unless ($_ =~ m/#/) {@cols = split /\t/; print "$cols[0]\t$cols[1]\t$cols[2]\n$cols[3]\t$cols[4]\t$cols[5]\n";}' $SingleInteractorLoops > sample.single_interactor.loop_regions.bed

sort-bed sample.single_interactor.loop_regions.bed >sample.single_interactor.loop_regions.sorted.bed

Genes_mm10="/path/to/Mus_musculus.GRCm38.102.chr.gene.sorted.bed"

closest-features --dist --closest sample.single_interactor.loop_regions.sorted.bed $Genes_mm10> sample.single_interactor.loop_regions.sorted.gene.bed
```

`Mus_musculus.GRCm38.102.chr.gene.sorted.bed` was created from the gtf downloaded from Ensembl:

```
perl -ne 'chomp $_; unless ($_ =~ m/#/) {@cols = split /\t/; $start=$cols[3]-1; my $end=$cols[4]; if( $cols[2] eq qw /gene/ ){print "chr$cols[0]\t$start\t$end\t$cols[8]\t$cols[5]\t$cols[6]\n";} }' Mus_musculus.GRCm38.102.chr.gtf  > Mus_musculus.GRCm38.102.chr.gene.bed

sort-bed Mus_musculus.GRCm38.102.chr.gene.bed >Mus_musculus.GRCm38.102.chr.gene.sorted.bed
```

* closest loop **sample.single_interactor.loop_regions.sorted.loopregion.bed**

```
perl -ne 'chomp $_; unless ($_ =~ m/#/) {@cols = split /\t/; print "$cols[0]\t$cols[1]\t$cols[2]\n$cols[3]\t$cols[4]\t$cols[5]\n";}' $Loops >sample.hicDetectLoops.loopregion.bed

sort-bed sample.hicDetectLoops.loopregion.bed >sample.hicDetectLoops.loopregion.sorted.bed

closest-features --no-overlaps --dist --closest sample.single_interactor.loop_regions.bed sample.hicDetectLoops.loopregion.sorted.bed> sample.single_interactor.loop_regions.sorted.loopregion.bed
```

* `mm10.chrom.sizes` downloaded from [UCSC](https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes)


Finally, to annotate loops with functional features:

```
Ctcf="/path/to/CTCF.IDR0.05.filt.narrowPeak"
Compart="/path/to/TAD.txt"

Chroms="/path/to/UCSC/mm10.chrom.sizes"


Loops="sample.ICE.single_interactor.txt"
LoopsBedops="sample.single_interactor.loop_regions.sorted.loopregion.bed"
Genes="sample.single_interactor.loop_regions.sorted.gene.bed"

Outfile="sample.single_interactor.annotated.txt"

perl $Scripts/annotate_loops.v.2.0.pl --loops $Loops --outfile $Outfile --compartments $Compart  --ctcf $Ctcf --genes_bedops $Genes --loops_bedops $LoopsBedops --chrom_sizes $Chroms
```

## Filtering


The annotated loops in files suffixed `single_interactor.annotated.txt` can be filtered using perl script filter_loops.v.2.0.pl.

The usage is:

```
Infile="/path/to/sample.single_interactor.annotated.txt"

Outfile="/path/to/sample.single_interactor.annotated.filtered.txt"

perl $Scripts/filter_loops.v.2.0.pl --infile $Infile --outfile $Outfile <filtering options>
```

The filtering options available are:


* `--distance` internal loop distance, i.e. the minimal distance between the two interacting fragments in bp; max 2000000 as this was the cutoff used in loop identification;

* `--loop_dist` minimal distance to other fragments involved in chromatin loop formation in bp;

* `--compartment` the identity of the chromatin A/B compartment of both interacting regions;


Please note that the filtering options are hierarchical, i.e.`distance` needs to be defined to filter by `loop_dist` and `distance` and `loop_dist` need to be defined to filter by `compartment`.

For example, to obtain a list of loops with intra-loop distance 1 Mb, distance to other regions involved in loop formation of minimum 30 kb and both loop regions within the closed chromatin compartment B:

```
perl $Scripts/filter_loops.v.2.0.pl --infile $Infile --outfile $Outfile --distance 1000000 --loop_dist 30000 --compartment B
```

