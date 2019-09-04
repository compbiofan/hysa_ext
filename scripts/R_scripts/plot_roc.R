png("roc.png",width=2,height=2,units="in",res=600,pointsize=3)
par(
  mar      = c(6, 6, 3, 3),
  xaxs     = "i",
  yaxs     = "i",
  cex.axis = 1.5,
  cex.lab  = 1.5,
  cex.main = 1.5
)
delly<-read.table('delly.roc.csv')
plot(delly$V1,delly$V2,xlim=c(0,1),ylim=c(0,1.1),xlab="Sensitivity",ylab="Specificity",main="Comparison of ROC on Large Deletions (NA12878)", col="green3",pch=18)
ours<-read.table('hysa.roc.csv')
points(ours$V1, ours$V2,col="purple",pch=20)
pbhoney<-read.table('pbhoney.roc.csv')
points(pbhoney$V1, pbhoney$V2, col="red", pch=21)
svclassify<-read.table('svclassify.roc.csv')
points(svclassify$V1, svclassify$V2, col="blue", pch=22)
custom<-read.table('custom.roc.csv')
points(custom$V1, custom$V2, col="magenta", pch=23)
legend("topright",c("HySA (hybrid)","Delly (Illumina)", "PBHoney (Pacbio)", "custom (Pacbio)", "svclassify (hybrid)"), col=c("purple","green3", "red", "magenta", "blue"), pch=c(20, 18, 21, 23, 22),cex=1.2,pt.cex=1.2)
dev.off()
