install.packages("ggplot2")
library(ggplot2)
setwd("/homes/fzhao01/Bureau/quantlib/IMT2018-master/project4")
my_data <- read.table("data.txt",col.names=c("timesteps","NPV1","NPV2"), sep = "", dec = "\t")

ggplot(my_data)+geom_line(aes(x=timesteps, y=NPV1, group=1), color="chartreuse3")+
geom_line(aes(x=timesteps, y=NPV2, group=1),color="red")+labs(title="Différence")

install.packages("tidyr")
library(tidyr)
my_dataNew <- gather(my_data, key=type, value= dataNPV, NPV1, NPV2)
my_dataNew[,2] <- as.factor(my_dataNew[,2])
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()
my_dataNew[,2] <- factor(my_dataNew[,2],levels=c("NPV2","NPV1"))
ggplot(my_dataNew, aes(x=timesteps, y=dataNPV, color=type, group=type))+geom_line()+labs(title="Comparaison between BinomialVanillaEngine and BinomialVanillaEngine_2")
