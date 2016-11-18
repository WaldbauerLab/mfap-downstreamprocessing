
## REQUIRED USER INPUTS 
# Path to working directory where MFAP .PSMs.tsv file is located and output will be written
workdir<-'D:/Orbitrap data/Jul2016_OneStep/Controls/MFAP_Results/Syn 9'
# Name of MFAP-generated .PSMs.tsv file WITHOUT extensions
filename<- 'Nov16-vSyn9-B16T3'
# Maximum q-value (%) to consider for cutoff point
top<-1.0
# Minimum q-value (%) to consider for cutoff point
bottom<-0.1

## LOADING LIBRARIES
library(BioPhysConnectoR)
library(plyr)

## LOADING INITIAL INPUT
setwd(workdir)
file1<-paste0(filename,'_cutoff_info.tsv')
file2<-paste0(filename,'_cutoff_plot.pdf')
file3<-paste0(filename,'_protein_output.tsv')
file4<-paste0(filename,'_summary.tsv')
file5<-paste0(filename,'_spectra_dist.pdf')
input<-paste0(filename,'.PSMs.tsv')

Input_Table<-read.delim(input, header=TRUE)

## CLEANING INITIAL INPUT
Input_Table<-lapply(Input_Table,as.vector)
nametag<-names(Input_Table)
Input_Table<-mapply(c,Input_Table)
colnames(Input_Table)<-nametag

## SORTING BY Q VALUE AND REMOVING ALL POINTS WITH Q VALUE > TOP
Numbermat<-matrix(NA,nrow=nrow(Input_Table),ncol=3)
Numbermat[,1]<-as.numeric(Input_Table[,'Morpheus.Score'])
Numbermat[,2]<-as.numeric(Input_Table[,'Q.Value....'])
Numbermat<-mat.sort(Numbermat,2)
Cutoff<-min(which(Numbermat[,2]>=top))
Numbermat<-Numbermat[1:Cutoff,]
Numbermat<-mat.sort(Numbermat,1)
positions<-which(Input_Table[,'Q.Value....'] %in% Numbermat[,2])
Cleaned_Table<-Input_Table[positions,]


## CREATING FIRST QPLOT
Qplotvals<-matrix(NA,nrow=nrow(Cleaned_Table),ncol=2)
Qplotvals[,1]<-1:nrow(Qplotvals)
Numbermat[,3]<-1:nrow(Numbermat)
Numbermat<-mat.sort(Numbermat,1,decreasing = TRUE)
Qplotvals[,2]<-Cleaned_Table[,'Q.Value....']
# plot(x=Qplotvals[,1],y=Qplotvals[,2],xlab='Spectrum Number',ylab='Q value'
#     ,main='Q Values Ordered By Morpheus Score')

## SELECT CUTOFF Q VALUE TO USE ((REMEMBER Q VALUES ARE ALREADY PERCENTS))
## FINDING ALL SPECTRA WITH Q VALUE BETWEEN BOTTOM AND TOP
LastZero<-as.numeric(max(which(Cleaned_Table[,'Q.Value....']<=bottom)))
zerozero<-as.numeric(max(which(Cleaned_Table[,'Q.Value....']==0)))
plotters<-Qplotvals[(zerozero):nrow(Qplotvals),]
targets<-Qplotvals[(LastZero+1):nrow(Qplotvals),]
plotters<-apply(plotters,2,as.numeric)
targets<-apply(targets,2,as.numeric)
plotQ<-plotters[,2]
Q.Value<-targets[,2]
targets[,2]<-targets[,2]*(targets[nrow(targets),1]/targets[nrow(targets),2])
# plot(targets,main="Cleaned Q Values Sorted by Ascending M Score",
#     xlab="Spectrum Number",ylab="Scaled Q Value")
specs<-nrow(targets)
z<-1
lengthvec<-c(NA)
lengthcount<-0
for(i in 2:specs){
  dif<-targets[i,2]-targets[i-1,2]
  if(dif<=0){
    lengthcount<-lengthcount+1
  }
  if(dif>0){
    lengthvec[z]<-lengthcount
    z<-z+1
    lengthcount<-0
  }
}
lengthvec<-lengthvec+1
cumlength<-cumsum(lengthvec)
rankvec<-rank(lengthvec,ties.method ='last')
rec<-which(rankvec==(max(rankvec)-1))
rec<-cumlength[rec]
pos<-which(rankvec>=(max(rankvec)-9))
spikes<-cumlength[pos]
indices<-Q.Value[(spikes)]

spikeinfo<-matrix(NA,nrow=length(indices),ncol=9)
colnames(spikeinfo)<-c('Legend_Key','FDR','Num_Points_Included','Increase_in_Points','Increase_in_Points_cum',
                       'Increase_in_FDR','Increase_in_FDR_cum','FDR/point(x10^4)','FDR/point_cum(x10^4)')
spikeinfo[,2]<-indices
spikeinfo[,3]<-spikes
spikeinfo[2:length(spikeinfo[,4]),4]<-diff(spikeinfo[,3])
spikeinfo[2:length(spikeinfo[,5]),5]<-cumsum(spikeinfo[2:nrow(spikeinfo),4])
spikeinfo[2:length(spikeinfo[,6]),6]<-diff(spikeinfo[,2])
spikeinfo[2:nrow(spikeinfo),7]<-cumsum(spikeinfo[2:nrow(spikeinfo),6])
spikeinfo[,8]<-(spikeinfo[,6]/spikeinfo[,4])*10^4
spikeinfo[,9]<-(spikeinfo[,7]/spikeinfo[,5])*10^4
spikeinfo[,1]<-1:nrow(spikeinfo)
write.table(spikeinfo,file=file1)

##CREATING PLOT SHOWING CHOICES FOR Q CUTOFF VALUES AND THEIR LEGEND KEY FOR CUTOFF OUTPUT FILE
colors<-c('red2','maroon1','darkorange','gold','lawngreen','green','lightseagreen','royalblue','royalblue4','mediumorchid4')
plot(x=plotters[,1],y=plotQ,xlab='Spectra Number',ylab='Q Value (%)',
     main='Spectra in Order of Descending Morpheus Score w/ Cutoff Points Highlighted')
abline(v=targets[spikes,1],col=colors)
points(targets[spikes,1],indices,col=colors,pch=19)
#abline(v=targets[1,1]+rec-1,lty=3,lwd=2)
legend(LastZero,max(Q.Value),c(spikeinfo[,1]),lty=c(rep(1,nrow(spikeinfo))),
       lwd=c(rep(1,nrow(spikeinfo))),col=c(colors,'black'),pt.cex=1,cex=0.75,bty='n')
pdf(file=file2)
plot(x=plotters[,1],y=plotQ,xlab='Spectra Number',ylab='Q Value (%)',
     main='Spectra in Order of Descending Morpheus Score w/ Cutoff Points Highlighted')
abline(v=targets[spikes,1],col=colors)
points(targets[spikes,1],indices,col=colors,pch=19)
#abline(v=targets[1,1]+rec-1,lty=3,lwd=2)
legend(LastZero,max(Q.Value),c(spikeinfo[,1]),lty=c(rep(1,nrow(spikeinfo))),
       lwd=c(rep(1,nrow(spikeinfo))),col=c(colors,'black'),pt.cex=1,cex=0.75,bty='n')
dev.off()


#####################################################################################################################
######     RUN ABOVE SECTIONS FIRST, THEN CHOOSE CUTOFF POINT BASED ON PLOT
######     INPUT CHOSEN CUTOFF POINT NUMBER BELOW, THEN RUN NEXT SECTIONS
#####################################################################################################################


####SELECT QCutoff
CutoffPoint<-4
QCutoff<-spikeinfo[CutoffPoint,'FDR']
point<-which(Numbermat[,2]==QCutoff)
point<-Numbermat[point,1]
Passed<-Input_Table[which(as.numeric(Input_Table[,'Morpheus.Score'])>=point),]
Passed<-mat.sort(Passed,'Protein.Description')
Proteins<-unique(Passed[,'Protein.Description'])


## CREATE PROTEIN DATA OUTPUT FILE
Output_list<-list(NA)
for(i in 2:length(Proteins)+1){
  lower<-match(Proteins[(i-1)],Passed[,'Protein.Description'])
  upper<-match(Proteins[i],Passed[,'Protein.Description'])
  if(is.na(upper)==TRUE){
    upper=lower
  }
  Indices<-c(lower:(upper-1))
  Num_Spectra<-as.numeric(length(Indices))
  Num_Peptides<-as.numeric(length(unique(Passed[Indices,'Base.Peptide.Sequence'])))
  Median_Log2<-median(as.numeric(Passed[Indices,'log2.18O.16O.Ratio']))
  SD_Log2<-sd(as.numeric(Passed[Indices,'log2.18O.16O.Ratio']))
  SE_Log2<-SD_Log2/((Num_Peptides)^(1/2))
  Output_list[[i-1]]<-list(name=Proteins[i-1],Num_Spectra=Num_Spectra,Num_Peptides=Num_Peptides,Median_Log2=Median_Log2,SD_Log2=SD_Log2,SE_Log2=SE_Log2)
}
Output<-do.call(rbind,Output_list)
Output<-Output[2:nrow(Output),]
Output2 = apply(Output, 2, unlist)
Output_Table = as.data.frame(Output2,stringsAsFactors= F)
Output_Table[,2:6]<- lapply(Output_Table[,2:6], as.numeric)

Output_Table<- cbind(Output_Table, decoy = grepl("DECOY",Output_Table$name))

Good_Output <- Output_Table[Output_Table$decoy == F,]
norm_factor <- median(Good_Output$Median_Log2)
Normalized_Log2 = Good_Output$Median_Log2 - norm_factor

Good_Output <- cbind(filename,Good_Output,norm_factor, Normalized_Log2)
write.table(Good_Output[,-c(8)],file3,sep="\t",col.names=NA)

## WRITE OUT SUMMARY STATS FILE
num_Proteins = nrow(Output_Table[Output_Table$decoy == F,])
Protein_FDR = nrow(Output_Table[Output_Table$decoy == T,])/nrow(Good_Output)
Multispectra_FDR = nrow(Output_Table[Output_Table$decoy == T & Output_Table$Num_Spectra > 1,])/nrow(Good_Output[Good_Output$Num_Spectra > 1,])
Num_spectra = nrow(Passed)

Sum_Stats = c(num_Proteins, Num_spectra, norm_factor,Protein_FDR, Multispectra_FDR)
names(Sum_Stats)<- c("num_Proteins","num_Spectra", "norm_factor","Protein_FDR", "Multispectra_FDR")
Spectra_counts = table(Output_Table$Num_Spectra)

#library(lattice)
#histogram(~Num_Spectra, data = Output_Table[Output_Table$Num_Spectra<= 20,], nint = 20, type = "count",xlab = "spectra per peptide",
#	main = "histogram of spectra/peptide, n <=20")
#pdf(file = file5)
#histogram(~Num_Spectra, data = Output_Table[Output_Table$Num_Spectra<= 20,], nint = 20, type = "count",xlab = "spectra per peptide",
#	main = "histogram of spectra/peptide, n <=20")
#dev.off()

options(scipen=999)
Summ_Output <- list(Sum_Stats, Spectra_counts)
names(Summ_Output)<- c("", "Spectral counts")
Summ_Output
sink(file4)
Summ_Output
sink()

