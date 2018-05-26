
ds<-read.table('dS_value.txt')
ds<-subset(ds,(V1<=2) & V1>0.000)
hist(ds$V1,breaks = 150)
ds$V1<-log10(ds$V1)
hist(ds$V1,breaks = 150)
quit()
