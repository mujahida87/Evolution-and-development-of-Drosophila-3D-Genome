# Evolution-and-development-of-Drosophila-3D-Genome
To study the role of transposable elements (TEs) in 3D genome organization in Drosophila, we measure TE expression using RNA-seq data across development in Drosophila Pseudoobscura.
To measure the TE-expression we used standard tools/software along with custom code to produce RPKM values.

Step1:
01.Star_index.sh
module load star/2.7.3a

STAR --runMode genomeGenerate --runThreadN 4 --genomeDir ./ --genomeFastaFiles Dpse.fa --sjdbGTFfile Dpse.gtf

Step2:
mapping using star
STAR --runMode alignReads --runThreadN 8 --genomeDir ./ --readFilesIn $fq1 $fq2 --sjdbGTFfile Dpse.gtf --readFilesCommand gunzip -c --outFileNamePrefix $out. --outSAMtype BAM Unsorted --winAnchorMultimapNmax 100 --outFilterMultimapNmax 100 --outFilterScoreMinOverLread 0.3 --outFilterMatchNminOverLread 0.3

Step3:
module load stringtie/2.1.1
stringtie -e -B -p 4 -A tissue.gene_abund -o tissue.gtf -G Dpse.gtf tissueAligned.sort.bam

Step4:
feature count
featureCounts -a Dpse.TEs.gtf -o $out $bam -t exon -f --largestOverlap -M -p --countReadPairs -T 4

Step5:

#!/usr/bin/perl
use strict;
use List::Util qw/max min/;
use Statistics::Descriptive;

open(IN1,"$ARGV[0]") || die "Can't open IN1!\n";
open(IN2,"$ARGV[1]") || die "Can't open IN2!\n";
open(IN3,"$ARGV[2]") || die "Can't open IN3!\n";
open(OUT,">$ARGV[3]") || die "Can't open OUT!\n";

my $total_reads;
while(<IN1>){
        chomp;
        next if(/^Status/);
        my @tmp=split /\s+/,$_;
        $total_reads+=$tmp[-1];
}
print OUT "#total_mapped_reads: $total_reads\n";

my %hash;
while(<IN2>){
        chomp;
        next if(/^#/);
        s/\"|\;//g;
        my @tmp=split /\s+/,$_;
        $hash{$tmp[0]}{$tmp[3]}="$tmp[9]\t$tmp[13]\t$tmp[15]";
}

while(<IN3>){
        chomp;
        next if(/^#/);
        if(/^Geneid/){
                print OUT "TE_id\tfamily_id\tclass_id\tchr\tstart\tend\tstrand\tlength\treadcount\tRPKM\n";
        }
        else{
                my @tmp=split /\s+/,$_;
                next unless(exists $hash{$tmp[1]}{$tmp[2]});
                print OUT "$hash{$tmp[1]}{$tmp[2]}\t$tmp[1]\t$tmp[2]\t$tmp[3]\t$tmp[4]\t$tmp[5]\t$tmp[6]\t";
                my $rpkm=($tmp[-1]*1000*1000000)/($total_reads*$tmp[5]) if($total_reads>0 && $tmp[5]>0);
                print OUT "$rpkm\n";
        }
}


then run

perl TE_RPKM.pl tissue.TE_FC.summary Dpse.TEs.gtf Hm.TE_FC tissue.TE_FC.RPKM
