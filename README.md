# PiNSiR

This protocol has been developed for calculating Pi_n and Pi_s values from SNPs called using GATK pipeline.

***NOTE:The code is still very much experimental, but I'm determined to improve it with time.***

Required software:
 1. GATK for SNP calling
 2. Installation of ANGSD, versions supporting vcf-pl input are recommended, here I've used v0.933 available through conda.
 3. snpEff v4.3t
 
The R code for Pi_n and Pi_s calculations requires three files:
   1. A VCF file with SNP calls also from monomorphic sites.
   2. A file with coordinates of deleterious SNPs.
   3. A (gzipped) file of theta values computed with Angsd, contains site-wise pi.

Getting each of the files, step by step:

1. A VCF file with SNP calls also from monomorphic sites.

-Run GATK according to best practices protocols, including the options:

   GATK4.x: --output-mode EMIT_ALL_CONFIDENT_SITES
   
   GATK3.x: --includeNonVariantSites
   
-After getting the VCF, carry out filtering. You will need to do it in two stages.
 In stage one, you will prepare a vcf file containing all good quality positions, also monomorphic ones. 
 This file is used to assess the total number of called sites (=genome coverage). To do this, apply only 
 filters according to the coverage and SNP call quality. Don't include filters involving heterozygosity
 (you'll just drop out the monomorphic positions). 
 
 Instead of using the VCF containing all positions, it is also possible to work with "normal" 
 VCF file containing only the heterozygous loci. However, this causes an ascertainment bias
 which needs to be corrected, and I haven't yet implemented the correction. However,
 you may still proceed at your own risk . Then the estimated diversities are biased towards 
 zero, because the code assumes zero diversity in all non-reported positions.

2. A VCF file with SNP annotations from SnpEff, use this to get a text
file with coordinates for deleterious positions.

-Filter the VCF file to include only variant positions and remove the indels, 
 because ANGSD doesn't support them. For the same reason you also need to remove deletion characters "\*"
 (it's not pretty but I eventually ended up removing them with sed 's\\,\*\\\\g'). 
 
-Run SnpEff to get the SNP impact predictions.

-Get the positions with high or moderate impact, you only need the coordinates:

zcat myvcf.gz | grep "HIGH\\|MODERATE" |  awk '{print $1"\t"$2}' > HIGH_MODERATE.coord 

-This step needs to be done only once for the whole VCF. 

3. Prepare a (gzipped) file of theta values computed with ANGSD,
contains site-wise pi. 

NOTE: It is highly recommended to use ANGSD for pi calculations, the vcftools has this option but it doesn't use genotype
likelihoods as ANGSD.

An alternative is PIXY
(https://www.biorxiv.org/content/10.1101/2020.06.27.175091v1.full) but
I haven't yet tried this out. One more option is to calculate these
values inside R (to be implemented later).

For best results you need to have an ancestral genome to get an
unfolded sfs.

How to make ancestral state calls (easy way):

Collect the .bam files of the outgroup individuals into a text file, say ancs.txt, 
and run ANGSD to get state calls:

   angsd -b ancs.txt -doFasta 2 -doCounts 1 -P 40 -out ANCESTOR.fasta

How to get SFS:

angsd -vcf-pl my.vcf.gz -doSaf 1 -anc ANCESTOR.fasta.fa.gz -out angsd.out

realSFS angsd.out.saf.idx -P 12 > angsd.out.saf.sfs

Then calculate site-wise pi:

realSFS saf2theta angsd.out.saf.idx -sfs angsd.out.saf.sfs -outname angsd.out.theta

thetaStat print angsd.out.theta.thetas.idx > angsd.out.theta.thetas.txt

To save space you can zip it.

Now we can finally do the calculations!

4. To get an accurate calculation of Pi_N values, we first need to find out which 
gene models are of good quality. I wrote a filter.hq.genes function to do this, 
it tests that the combined exon length is divisible with three and whether the 
protein starts with a methionine. To find these high quality gene models you 
will need a gff file with gene predictions and a peptide fasta file of the 
predicted genes. 
(will explore whether it's necessary to have the fasta - if not this,
then we still need the genome fasta and translation to proteins in
R which will take more time. For now, you can use for example gffread to produce a peptide
fasta file from gff3 and genome sequence.)

5. For the actual calculations there are two functions: Pi_N.par and Pi_S.par (par for parallel)

Example usage:

library(PiNSiR)

#1. Get HQ gene models

gff.file="my.genes.gff3"

peptide.fasta="my.protein.faa"

gene.coord=filter.hq.genes(gff.file,peptide.fasta)

#Get SNP positions with high and moderate effect

zerofold.pos=read.delim("snpEff.HIGH_MODERATE.coord",header=F,as.is=T)

#Get site-wise theta - ANGSD needs to be run for each population separately

TH=read.delim("myvcf.mygroup.theta.thetas.txt.gz",as.is=T)

#Calculate Pi_N
#Pi_N
#nSNP=read.table("myvcf.filtered.all.exons.mygroup.snp.count")[,1]*2/3
#Pi_0=Pi_N.par(TH,gene.coord$HQ,zerofold.pos,nSNP)

Pi_0=Pi_N.par(TH,gene.coord$HQ,zerofold.pos)

#Calculate Pi_S
#Pi_S
#nSNP=read.table("myvcf.filtered.all.intergenic.mygroup.snp.count")[,1]
#Pi_int=Pi_S.par(TH,gene.coord$ALL,nSNP)

Pi_int=Pi_S.par(TH,gene.coord$ALL)
