filter.ninth=function(x,ID="Parent"){
  idvec=vector("character",length(x))
  for (i in 1:length(x)){
    y=strsplit(x[i],";")[[1]]
    idvec[i]=strsplit(y[grep(ID,y)],"=")[[1]][2]
    idvec[i]=gsub(" ","",idvec[i])
  }
  return(idvec)
}

#' @export
filter.hq.genes=function(gff.file,peptide.fasta,bed=F){

  gff=read.delim(gff.file,comment.char="#",header=F,as.is=T)
  gff=gff[gff[,3]=="exon",c(1,4,5,9)]
  #Prepare a gff file with exons, 4th column is the parent gene ID.
  gff[,4]=filter.ninth(gff[,4])
  # split according to gene ID
  gff.split=split(gff,gff[,4])

  #Check that the cds is a multiple of 3 (intact codons)
  gene.OK=which(sapply(gff.split,function(x) sum(x[,3]-x[,2]+1)%%3)==0)

  pep=seqinr::read.fasta(peptide.fasta,seqtype="AA")
  #Check that the protein begins with M
  pep.OK=which(sapply(pep,function(x) x[1]=="M"))
  prot.OK=intersect(names(pep.OK),names(gene.OK))

  # In case we want to get a more accurate estimate, look at codons to define N nonsyn positions
  # Not yet implemented
  #nuc=read.fasta("/data2/jtsaloja/Lychee/Assembly/Latest/HiC_Correct_Final.40Xhic_gmap.cds",seqtype="DNA")
  #nuc.OK=nuc[prot.OK]
  #nuc.codons=sapply(nuc.OK,splitseq)
  #codon.table=table(unlist(nuc.codons))

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


contig.Pin=function(TH,exs,ex.mut){
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
   N.nonsy=2*sum(exs[,3]-(exs[,2]+3)+1)/3
   return(list(sum.theta=th.sum,N=N.nonsy))
}

#' Calculate Pi_N for loci with a high impact or moderate impact mutation
#' @export
Pi_N.par=function(TH, exon.sel, exon.mut,nSNP=NA,cores=15){
  message("Calculating Pi_N - parallel")
  TH.split=split(TH,TH[,1])
  exon.mut.list=split(exon.mut,exon.mut[,1])
  exon.sel.list=split(exon.sel,exon.sel[,1])
  sel.contigs=intersect(names(TH.split),names(exon.sel.list))
  TH.split=TH.split[sel.contigs]
  exon.mut.list=exon.mut.list[sel.contigs]
  exon.sel.list=exon.sel.list[sel.contigs]
  contig.pi=parallel::mcmapply(contig.Pin,TH.split,exon.sel.list,exon.mut.list,mc.cores=min(cores,length(TH.split)))
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
  N.nonsy=sum(exon.sel[,3]-(exon.sel[,2]+3)+1)/3
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

