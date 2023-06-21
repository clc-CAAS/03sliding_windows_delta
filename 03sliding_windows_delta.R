#!/usr/bin/env Rscript
argv<-commandArgs(TRUE) 
pdf(argv[2],10,3)
par(mar=c(3,4.5,1.5,4.5),oma=c(0,1,0,1),mgp=c(3.2,0.65,0),cex.axis=1.15,las=1,cex.lab=1.3)
a <- read.table(argv[1],header=F)
###Chr len for Pos infor
len <- c(0)
a[,1] <- as.factor(a[,1])
org = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)
seq1 = c("Chr1" , "Chr2",  "Chr3" , "Chr4" , "Chr5" , "Chr6" , "Chr7" , "Chr8" , "Chr9", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15" ,"Chr16", "Chr17", "Chr18")
for (i in org) {
  snp <- subset(a,a$V1==seq1[i])
  len <- c(len,snp[,2][length(snp[,2])])
}
plot(1,1,type='n',xaxs='i',yaxs='i',cex=0.1,col='steelblue1',axes=F,main='',xlab='',ylab=expression(Delta(SNP-index)),xlim=c(0,sum(len)),ylim=c(-1,1))
par(new=T)
len_sum <- 0
len_sum_d <- c(0)
for (i in 1:18){
  d <- subset(a,V1 == paste('Chr',i,sep=''))
  len_sum =len_sum+ len[i]
  len_sum_d <-c(len_sum_d,len_sum)
  pos <- (d$V2)
  plot(len_sum+pos,d$V8,type='h',xaxs='i',yaxs='i',col='gray81',axes=F,main='',xlab='',ylab='',xlim=c(0,sum(len)),ylim=c(0,10000))  # snp number distribution
  par(new=T)
  for (n in 1:length(d[,2])){
    plot(1,1,type='n',xaxs='i',yaxs='i',cex=0.1,col='blue',axes=F,main='',xlab='',ylab='',xlim=c(0,sum(len)),ylim=c(-1,1)) # add layer out
    points(len_sum+d[n,2],d[n,6],type='p',xaxs='i',yaxs='i',cex=0.2,col='blue',pch=20) # plot point
    points(len_sum+d[n,2],d[n,7],type='p',xaxs='i',yaxs='i',cex=0.2,col='red',pch=20) # plot point
    par(new=T)
  }
  par(new=T)
}
len_sum_d <-c(len_sum_d,sum(len))
len_sum_d1<-len_sum_d[2:20]

axis(4,tcl=-0.3,lwd=2,lwd.ticks=2,at=seq(-1,1,0.4),labels=seq(0,10000,2000))
mtext('Number of variations', side = 4, line = 3.2,cex=1.3,las=0)
axis(2,tcl=-0.3,lwd=2,lwd.ticks=2,at=seq(-1,1,0.5),labels=seq(-1,1,0.5))
axis(1,at=c(len_sum_d1[1],len_sum_d1[2],len_sum_d1[3],len_sum_d1[4],len_sum_d1[5],len_sum_d1[6],len_sum_d1[7],len_sum_d1[8],len_sum_d1[9],len_sum_d1[10],len_sum_d1[11],len_sum_d1[12],len_sum_d1[13],len_sum_d1[14],len_sum_d1[15],len_sum_d1[16],len_sum_d1[17],len_sum_d1[18],len_sum_d1[19]),labels=F,lwd=2,lwd.ticks=2,tcl=-0.3)
axis(1,tick=F,at=c(len_sum_d1[2]/2,len_sum_d1[2]+(len_sum_d1[3]-len_sum_d1[2])/2,len_sum_d1[3]+(len_sum_d1[4]-len_sum_d1[3])/2,len_sum_d1[4]+(len_sum_d1[5]-len_sum_d1[4])/2,len_sum_d1[5]+(len_sum_d1[6]-len_sum_d1[5])/2,len_sum_d1[6]+(len_sum_d1[7]-len_sum_d1[6])/2,len_sum_d1[7]+(len_sum_d1[8]-len_sum_d1[7])/2,len_sum_d1[8]+(len_sum_d1[9]-len_sum_d1[8])/2,len_sum_d1[9]+(len_sum_d1[10]-len_sum_d1[9])/2,len_sum_d1[10]+(len_sum_d1[11]-len_sum_d1[10])/2,len_sum_d1[11]+(len_sum_d1[12]-len_sum_d1[11])/2,len_sum_d1[12]+(len_sum_d1[13]-len_sum_d1[12])/2,len_sum_d1[13]+(len_sum_d1[14]-len_sum_d1[13])/2,len_sum_d1[14]+(len_sum_d1[15]-len_sum_d1[14])/2,len_sum_d1[15]+(len_sum_d1[16]-len_sum_d1[15])/2,len_sum_d1[16]+(len_sum_d1[17]-len_sum_d1[16])/2,len_sum_d1[17]+(len_sum_d1[18]-len_sum_d1[17])/2,len_sum_d1[18]+(len_sum_d1[19]-len_sum_d1[18])/2),labels=c('A01','A02','A03','A04','A05','A06','A07', 'A08','A09','A10','B01','B02','B03','B04','B05','B06','B07','B08'))

for (i in 2:18){
  segments(len_sum_d1[i],-1,len_sum_d1[i],1,lty=2,lwd=1.35,col=rgb(0,0,0,0.5))
}
abline(h=0,lty=2,lwd=1.6,col='black')
legend('topright',legend='95% confidence interval',col='red',bty='n',lty=1,lwd=2,cex=1)


