install.packages("ggplot2")
library(ggplot2)

install.packages("tidyr")
library(tidyr)

setwd("/homes/fzhao01/Bureau/quantlib/IMT2018-master/project4")
my_data <- read.table("data.txt",col.names=c("timesteps","NPV1","NPV2"), sep = "", dec = "\t")

ggplot(my_data)+geom_line(aes(x=timesteps, y=NPV1, group=1), color="chartreuse3")+
geom_line(aes(x=timesteps, y=NPV2, group=1),color="red")+labs(title="Différence")

my_dataNew <- gather(my_data, key=type, value= dataNPV, NPV1, NPV2)
my_dataNew[,2] <- as.factor(my_dataNew[,2])
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()
my_dataNew[,2] <- factor(my_dataNew[,2],levels=c("NPV2","NPV1"))
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()+labs(title="Comparaison between BinomialVanillaEngine and BinomialVanillaEngine_2", subtitle = "Cox-Ross-Rubinstein")

my_dataJr <- read.table("dataJr.txt",col.names=c("timesteps","NPV1","NPV2"), sep = "", dec = "\t")

ggplot(my_dataJr)+geom_line(aes(x=timesteps, y=NPV1, group=1), color="chartreuse3")+
  geom_line(aes(x=timesteps, y=NPV2, group=1),color="red")+labs(title="Différence")

my_dataNew <- gather(my_dataJr, key=type, value= dataNPV, NPV1, NPV2)
my_dataNew[,2] <- as.factor(my_dataNew[,2])
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()
my_dataNew[,2] <- factor(my_dataNew[,2],levels=c("NPV2","NPV1"))
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()+labs(title="Comparaison between BinomialVanillaEngine and BinomialVanillaEngine_2", subtitle = "Jarrow-Rudd")

my_dataAe <- read.table("dataAe.txt",col.names=c("timesteps","NPV1","NPV2"), sep = "", dec = "\t")

ggplot(my_dataAe)+geom_line(aes(x=timesteps, y=NPV1, group=1), color="chartreuse3")+
  geom_line(aes(x=timesteps, y=NPV2, group=1),color="red")+labs(title="Différence")

my_dataNew <- gather(my_dataAe, key=type, value= dataNPV, NPV1, NPV2)
my_dataNew[,2] <- as.factor(my_dataNew[,2])
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()
my_dataNew[,2] <- factor(my_dataNew[,2],levels=c("NPV2","NPV1"))
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()+labs(title="Comparaison between BinomialVanillaEngine and BinomialVanillaEngine_2", subtitle = "Additive equiprobabilities")


my_dataBt <- read.table("dataBt.txt",col.names=c("timesteps","NPV1","NPV2"), sep = "", dec = "\t")

ggplot(my_dataBt)+geom_line(aes(x=timesteps, y=NPV1, group=1), color="chartreuse3")+
  geom_line(aes(x=timesteps, y=NPV2, group=1),color="red")+labs(title="Différence")

my_dataNew <- gather(my_dataBt, key=type, value= dataNPV, NPV1, NPV2)
my_dataNew[,2] <- as.factor(my_dataNew[,2])
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()
my_dataNew[,2] <- factor(my_dataNew[,2],levels=c("NPV2","NPV1"))
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()+labs(title="Comparaison between BinomialVanillaEngine and BinomialVanillaEngine_2", subtitle = "Binomial Trigeorgis")

my_dataTian <- read.table("dataTian.txt",col.names=c("timesteps","NPV1","NPV2"), sep = "", dec = "\t")

ggplot(my_dataTian)+geom_line(aes(x=timesteps, y=NPV1, group=1), color="chartreuse3")+
  geom_line(aes(x=timesteps, y=NPV2, group=1),color="red")+labs(title="Différence")

my_dataNew <- gather(my_dataTian, key=type, value= dataNPV, NPV1, NPV2)
my_dataNew[,2] <- as.factor(my_dataNew[,2])
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()
my_dataNew[,2] <- factor(my_dataNew[,2],levels=c("NPV2","NPV1"))
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()+labs(title="Comparaison between BinomialVanillaEngine and BinomialVanillaEngine_2", subtitle = "Binomial Tian")
