# PiNSiR

This protocol has been developed for calculating Pi_n and Pi_s values from SNPs called using GATK pipeline. 

The R code requires three files:
1. A VCF file with SNP calls also from monomorphic sites.

2. A file with coordinates of deleterious SNPs.

3. A (gzipped) file of theta values computed with Angsd, contains site-wise pi.

Getting each of the files, step by step:

1. A VCF file with SNP calls also from monomorphic sites.

-Run GATK according to best practices protocols, including the options

   GATK4.x: --output-mode EMIT_ALL_CONFIDENT_SITES
   
   GATK3.x: --includeNonVariantSites
   
-After getting this VCF, carry out filtering. There will be two steps,
 the second one is in protocol step 2. In step one, you will prepare a
 file containing all good quality positions, also monomorphic ones -
 you will need this to assess the total number of called sites
 (coverage). To do this apply only filters according to the coverage
 and SNP call quality. Don't include filters involving heterozygosity
 (won't work in monomorphic positions).

2. A VCF file with SNP annotations from SnpEff, use this to get a text
file with coordinates for deleterious poisitions.

-Filter the VCF file to include only variant positions. You may still
 keep the indels, since they contribute to the total nucleotide
 diversity. Notice however that they are still handled as single locus
 polymorphisms in subsequent calculations.

-Run SnpEff to get the SNP impact predictions.

-Get the positions with high or moderate impact:

zcat myvcf.gz | grep "HIGH\|MODERATE" |  awk '{print $1"\t"$2}' > HIGH_MODERATE.coord 

-Needs to be done only once for the whole VCF.

3. Select high quality gene models for Pi_N calculations, output will be coordinates of exons.

Use the filter.hq.genes function. For this you will need a gff file
with gene predictions and a peptide fasta file of the predicted genes
(will explore whether it's necessary to have the fasta - if not this
then we still need the genome fasta and translation to proteins in
R). For now, you can use for example gffread to produce a peptide
fasta file from gff3 and genome sequence.

4. Prepare a (gzipped) file of theta values computed with ANGSD,
contains site-wise pi. NOTE: It is highly recommended to use ANGSD for
this, the vcftools has this option but it doesn't use genotype
likelihoods as ANGSD.

An alternative is PIXY
(https://www.biorxiv.org/content/10.1101/2020.06.27.175091v1.full) but
I haven't yet tried this out. One more option is to calculate these
values inside R (to be implemented later).

For best results you need to have an ancestral genome to get an
unfolded sfs.

How to make ancestral state calls (easy way):

Collect the .bam files of the outgroup individuals into a text file, say ancs.txt, and run ANGSD to get state calls:

   angsd -b ancs.txt -doFasta 2 -doCounts 1 -P 40 -out ANCESTOR.fasta

How to get SFS:

angsd -vcf-pl my.vcf.gz -doSaf 1 -anc ANCESTOR.fasta.fa.gz -out angsd.out
/opt/angsd/misc/realSFS angsd.out.saf.idx -P 12 > angsd.out.saf.sfs

Then calculate site-wise pi:

/opt/angsd/misc/realSFS saf2theta angsd.out.saf.idx -sfs angsd.out.saf.sfs -outname angsd.out.theta

/opt/angsd/misc/thetaStat print angsd.out.theta.thetas.idx > angsd.out.theta.thetas.txt

To save space you can zip it:

pigz angsd.out.theta.thetas.txt

5. When these files are ready you can run the pi summarizations with R.
There are two functions: Pi_N.par and Pi_S.par

Example usage:

library(PiNSiR)

#1. Get HQ gene models

gff.file="CA_v0.6_v1_BTI_annotation.fix.gff3"
peptide.fasta="coffea_arabica_v0.6.faa"
gene.coord=filter.hq.genes(gff.file,peptide.fasta)

#Get SNP positions with high and moderate effect
zerofold.pos=read.delim("Arabica_sgC.HIGH_MODERATE.coord",header=F,as.is=T)

#Get site-wise theta - ANGSD needs to be run for each population separately
TH=read.delim("Arabica_sgC.NoRep.wild.theta.thetas.txt.gz",as.is=T)

#Calculate Pi_N
#Pi_N
#nSNP=read.table("Arabica_sgC.filtered.all.exons.wild.snp.count")[,1]*2/3
#Pi_0=Pi_N.par(TH,gene.coord$HQ,zerofold.pos,nSNP)
Pi_0=Pi_N.par(TH,gene.coord$HQ,zerofold.pos)

#Calculate Pi_S
#Pi_S
#nSNP=read.table("Arabica_sgC.filtered.all.intergenic.wild.snp.count")[,1]
#Pi_int=Pi_S.par(TH,gene.coord$ALL,nSNP)
Pi_int=Pi_S.par(TH,gene.coord$ALL)
