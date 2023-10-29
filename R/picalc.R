filter.ninth=function(x,ID){
  idvec=vector("character",length(x))
  for (i in 1:length(x)){
    y=strsplit(x[i],";")[[1]]
    idvec[i]=strsplit(y[grep(ID,y)],"=")[[1]][2]
    idvec[i]=gsub(" ","",idvec[i])
  }
  return(idvec)
}

#' Identify zero-fold degenerate positions (ie non-synonymous sites) from the gene models.
#' @param gff.file The input gff3 file with gene models. Recommended input: high quality gene models from filter.hq.genes
#' @return A data frame containing coordinates of non-synonymous positions.
#' 
#' @export
folds.0=function(gff){
  genes=unique(gff[,4])
  for (i in 1:length(genes)){
     gff1=gff[gff[,4]==genes[i],]
     exo=gff1[,3]-gff1[,2]+1
     if (gff1[1,5]=="-"){
        posvec=rep(c(3,2,1),length.out=sum(exo))
     }else{
        posvec=rep(c(1,2,3),length.out=sum(exo))
     }
     coordvec=c()
     for (j in 1:nrow(gff1))
        coordvec=c(coordvec,gff1[j,2]:gff1[j,3])
     zero.fold=coordvec[posvec==1 | posvec==2]
     df.0fold=data.frame(Chr=gff1[1,1],pos=zero.fold,ID=gff1[1,4])
     if (i==1){
       all.0fold=df.0fold
     }else{
       all.0fold=rbind(all.0fold,df.0fold)
     }
  }
  return(data.frame(zero.fold=all.0fold))
}     
 

get.1st2nd.codon=function(gff1){
  exo=gff1[,3]-gff1[,2]+1
  if (gff1[1,5]=="-"){
     posvec=rep(c(3,2,1),length.out=sum(exo))
  }else{
     posvec=rep(c(1,2,3),length.out=sum(exo))
  }
  coordvec=c()
  for (j in 1:nrow(gff1))
     coordvec=c(coordvec,gff1[j,2]:gff1[j,3])
  zero.fold=coordvec[posvec==1 | posvec==2]
  df.0fold=data.frame(Chr=gff1[1,1],pos=zero.fold,ID=gff1[1,4],stringsAsFactors=F)
  return(df.0fold)
}

folds.0.par=function(gff,cores=2){
  genes=unique(gff[,4])
  genes.split=split(gff,gff[,4])
  codon.list=parallel::mcmapply(get.1st2nd.codon,genes.split,mc.cores=min(cores,length(genes.split)),SIMPLIFY=F)
  codon.df=plyr::rbind.fill(codon.list)
  return(codon.df)
}     

#' Identify four-fold degenerate positions (4dtv) from the gene models.
#' @param gff.file The input gff3 file with gene models. Recommended input: high quality gene models from filter.hq.genes
#' @param peptide.fasta Path to the input fasta file with protein sequences of the gene models.
#' @return A data frame containing coordinates of 4dtv sites.
#' 
#' @export
folds.4=function(gff,peptide.fasta){
  fourfold=c("V","S","P","T","A","G")
  #read peptide fasta
  pep=seqinr::read.fasta(peptide.fasta,seqtype="AA")
  pep=pep[unique(gff[,4])]
  ngenes=length(pep)
  for (i in 1:ngenes){
     gff1=gff[gff[,4]==names(pep)[i],]
     exo=gff1[,3]-gff1[,2]+1
     if (gff1[1,5]=="-"){
        posvec=rep(c(3,2,1),length.out=sum(exo))
        protvec=cumsum(posvec==3)
     }else{
        posvec=rep(c(1,2,3),length.out=sum(exo))
        protvec=cumsum(posvec==1)
     }
     coordvec=c()
     for (j in 1:nrow(gff1))
        coordvec=c(coordvec,gff1[j,2]:gff1[j,3])
     zero.fold=coordvec[posvec==1 | posvec==2]
     #Val (V),Ser (S),Pro (P),Thr (T),Ala (A),Gly (G)
     pepvec=which(pep[[i]] %in% fourfold)
     four.fold=coordvec[protvec %in% pepvec & posvec==3]
     df.4dtv=data.frame(Chr=gff1[1,1],pos=four.fold,ID=gff1[1,4])     
     if (i==1){
       all.4dtv=df.4dtv
     }else{
       all.4dtv=rbind(all.4dtv,df.4dtv)
     }
     message(i)
  }
  return(four.fold=all.4dtv)
}     

#' Filter high quality genes
#'
#' This function selects high quality genes from the set of annotated genes in a given gff3 file.
#' It checks that the gene length is divisible by three (so that codons are intact) and that the
#' translated protein starts with a methionine. 
#'
#'
#' @param gff.file Path to the input gff3 file with gene models.
#' @param peptide.fasta Path to the input fasta file with protein sequences of the gene models.
#' @param bed Logical, whether to return coordinates in a bed format (True; 0-based) or gff3 coordinates (False, default; starts from 1).
#' @param gene.identifier The unique identifier for matching exons with gene models in the gff3 file. Default is "Parent"
#' @return A list with two items: HQ - high quality gene models and ALL - all gene models. Each list item is a reduced gff3 file.
#' @export
filter.hq.genes=function(gff.file,peptide.fasta,bed=F,region.identifier="exon",gene.identifier="Parent"){

  #read in gff
  gff=read.delim(gff.file,comment.char="#",header=F,as.is=T)
  gff=gff[gff[,3]==region.identifier,c(1,4,5,9,7)]
  #Prepare a gff file with exons, 4th column is the parent gene ID.
  gff[,4]=filter.ninth(gff[,4],gene.identifier)
  # split according to gene ID
  gff.split=split(gff,gff[,4])

  #read peptide fasta
  pep=seqinr::read.fasta(peptide.fasta,seqtype="AA")

  #Check that the protein begins with M
  pep.OK=which(sapply(pep,function(x) x[1]=="M"))

  #Check that the cds is a multiple of 3 (intact codons)
  gene.OK=which(sapply(gff.split,function(x) sum(x[,3]-x[,2]+1)%%3)==0)
  #Whole protein is OK - both checks passed
  prot.OK=intersect(names(pep.OK),names(gene.OK))

  #High quality gene models.
  hq=gff[gff[,4] %in% prot.OK,]
  #Bed is zero based
  if (bed){
    hq[,2]=hq[,2]-1
    hq[,3]=hq[,3]-1
  }
  #All exons 
  all=gff
  if (bed){
    all[,2]=all[,2]-1
    all[,3]=all[,3]-1
  }
  return(list(HQ=hq,ALL=all))
}


contig.Pin=function(TH,exs,ex.mut,datatype){
  th.sum=0;Nsnp=nrow(ex.mut)
  accept.vec=vector("logical",nrow(ex.mut))
  for (i in 1:nrow(ex.mut)){
    #find if mutation is within HQ exons
    pos=as.numeric(ex.mut[i,2])
    nd=which(exs[,2] <= pos & exs[,3] >= pos)    
    if (length(nd)>0)
      accept.vec[i]=T
   }
   ii=which(TH[,2] %in% ex.mut[accept.vec,2])
   th.sum=sum(exp(TH[ii,4]))
   #For now, assume no missing data
   if (datatype=="SNP")
      N.nonsy=2*sum(exs[,3]-(exs[,2]+3)+1)/3
   #option2: assume all positions are called and if missing, it was filtered out
   if (datatype=="full")
      N.nonsy=length(ii)
   return(list(sum.theta=th.sum,N=N.nonsy))
}

#' Calculate Pi_N for loci with a high impact or moderate impact mutation using parallel computation
#' @param TH Site-wise theta calculated using ANGSD
#' @param exon.sel Selected gene models
#' @param exon.mut Coordinates of the non-synonymous positions.
#' @param nSNP If given, divides the total Pi with the given number of SNPs. 
#' @param datatype "full" assumes that all positions are called Pi is a mean of observed values, "SNP" assumes only variant positions and the Pi is divided by the total number of non-synonymous positions in the selected gene models.
#' @param cores Number of cores for parallel computation
#' @return a list with Pi, chromosome-wise cumulative sums and numbers of SNPs.
#' @export
Pi_N.par=function(TH, exon.sel, exon.mut,nSNP=NA,datatype="full",cores=15){
  message("Calculating Pi_N - parallel")
  TH.split=split(TH,TH[,1])
  exon.mut.list=split(exon.mut,exon.mut[,1])
  exon.sel.list=split(exon.sel,exon.sel[,1])
  sel.contigs=intersect(names(TH.split),names(exon.sel.list))
  TH.split=TH.split[sel.contigs]
  exon.mut.list=exon.mut.list[sel.contigs]
  exon.sel.list=exon.sel.list[sel.contigs]
  contig.pi=parallel::mcmapply(contig.Pin,TH.split,exon.sel.list,exon.mut.list,datatype,mc.cores=min(cores,length(TH.split)))
  if (is.na(nSNP)){
     res=list(Pi=sum(unlist(contig.pi["sum.theta",]))/sum(unlist(contig.pi["N",])),chrom.sum=contig.pi["sum.theta",],chrom.N=contig.pi["N",])
     message("The number of called nucleotides is not known - now using total exon length")
     message("This assumes that all uncalled positions are monomorphic.")
  }else{
     res=list(Pi=sum(unlist(contig.pi["sum.theta",]))/nSNP,chrom.sum=contig.pi["sum.theta",],chrom.N=contig.pi["N",])
  }
  return(res)
}

Pi_N=function(TH, exon.sel, exon.mut){
  # Calculate cumulative pi_N
  message("Calculating Pi_N")
  TH.split=split(TH,TH[,1])

  ind.list=vector("numeric",0)
  th.sum=0;Nsnp=nrow(exon.mut)
  for (i in 1:nrow(exon.mut)){
    x=as.character(exon.mut[i,1])
    pos=as.numeric(exon.mut[i,2])
    nd=which(exon.sel[,1]==x & exon.sel[,2]+3 <= pos & exon.sel[,3] >= pos)
    if (length(nd)>0){
      ii=which(TH.split[[x]][,2]==pos)
      if (length(ii)>0){
        th.sum=th.sum+exp(TH.split[[x]][ii,4])
        ind.list=c(ind.list,i)
#        message(i,"/",Nsnp)
     }
   }
  }
  N.nonsy=2*sum(exon.sel[,3]-(exon.sel[,2]+3)+1)/3
  return(th.sum=th.sum,N=N.nonsy)
}

pi.sum=function(x,b){
    ex=which(b[,2]<x[nrow(x),2])
    b=b[ex,]
    # Drop out exons
    exclude=vector("integer",0)
    for (i in 1:nrow(b)){
       TH.ex=which(x[,2]>=b[i,2] & x[,2]<=b[i,3])
       exclude=c(exclude,TH.ex)
    }
    sum.theta=sum(exp(x[setdiff(1:nrow(x),exclude),4]),na.rm=T)
    N.neutral=(x[nrow(x),2]-x[1,2]+1)-sum((b[,3]-b[,2]+1))
    return(list(sum.theta=sum.theta,N.neutral=N.neutral))
}

#' Calculate Pi_S based on intergenic loci.
#' This is obtained by excluding all exonic loci.
#' @param TH Site-wise theta calculated using ANGSD
#' @param exon.all Coordinates of all gene models
#' @param nSNP If given, divides the total Pi with the given number of SNPs. 
#' @param cores Number of cores for parallel computation
#' @return a list with Pi from intergenic regions, chromosome-wise cumulative sums and numbers of SNPs.

#' @export
Pi_S.par=function(TH,exon.all,nSNP=NA,cores=15){
  message("Calculating Pi_s - parallel")
  TH.split=split(TH,TH[,1])
  exon.split=split(exon.all,exon.all[,1])
  sel.contigs=intersect(names(TH.split),names(exon.split))
  TH.split=TH.split[sel.contigs]
  exon.split=exon.split[sel.contigs]
  # calculate neutral pi. If SNP is missing, assume diversity to be zero.
  contig.pi=parallel::mcmapply(pi.sum,TH.split,exon.split,mc.cores=min(cores,length(TH.split)))
  if (is.na(nSNP)){
    res=list(Pi=sum(unlist(contig.pi["sum.theta",]))/sum(unlist(contig.pi["N.neutral",])),chrom.sum=contig.pi["sum.theta",],chrom.N=contig.pi["N.neutral",])
     message("The number of called nucleotides is not known - now using the range of called SNPs minus exon length.")
     message("This assumes that all uncalled positions are monomorphic (pi=0).")
  }else{
    res=list(Pi=sum(unlist(contig.pi["sum.theta",]))/nSNP,chrom.sum=contig.pi["sum.theta",],chrom.N=contig.pi["N.neutral",])
  }
  return(res)
}

#' Non-parallel Pi_S based on intergenic loci.
Pi_S=function(TH.split,coord){
  coord.split=split(coord,coord[,1])
  message("Calculating Pi_s")
  # calculate neutral pi. If SNP is missing, assume diversity to be zero.
  sum.theta=vector("numeric",length(TH.split))
  N.neutral=vector("numeric",length(TH.split))
  # Go through the contigs
  for (k in 1:length(TH.split)){
    x=TH.split[[k]]
    b=coord.split[[k]]
    ex=which(b[,2]<x[nrow(x),2])
    b=b[ex,]
    # Drop out exons
    exclude=vector("integer",0)
    for (i in 2:nrow(b)){
       TH.ex=which(x[,2]>=b[i,2] & x[,2]<=b[i,3])
       exclude=c(exclude,TH.ex)
#     if (length(TH.ex)>0)
#       x=x[TH.ex[length(TH.ex)]:nrow(x),]
#       message(i,"/",nrow(b))
    }
    sum.theta[k]=sum(exp(x[setdiff(1:nrow(x),exclude),4]))
    N.neutral[k]=(x[nrow(x),2]-x[1,2]+1)-sum((b[,3]-b[,2]+1))
    message(names(TH.split)[k])
  }
  res=list(Pi=sum(sum.theta)/sum(N.neutral),chrom.sum=sum.theta,chrom.N=N.neutral)
  return(res)
}

