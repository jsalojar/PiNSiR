
# Exonic SNPs
# bed.pres=list()
# TH.int=TH[1,]
# for (i in 1:nrow(TH)){
#   kk=(bed.split[[TH[i,1]]][,1]==TH[i,1] & bed.split[[TH[i,1]]][,2] <= TH[i,2] & bed.split[[TH[i,1]]][,3] >= TH[i,2])
#   if (!any(kk))
#      TH.int=rbind(TH.int,TH[i,])
#   else
#      bed.pres[[TH[i,1]]]=unique(rbind(bed.pres[[TH[i,1]]],bed.split[[TH[i,1]]][which(kk),]))
#   message(i)
#}
