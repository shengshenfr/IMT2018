setwd("/homes/fzhao01/Bureau/quantlib/IMT2018-master/project4")
my_data <- read.table("data.txt",col.names=c("timesteps","NPV"), sep = "", dec = "\t")
plot(my_data$timesteps, my_data$NPV, type = "l", xlab="timesteps",ylab="NPV")
my_dataOp <- read.table("dataBs.txt",col.names=c("timesteps","NPV"), sep = "", dec = "\t")
plot(my_dataOp$timesteps, my_dataOp$NPV, type = "l", xlab="timesteps",ylab="NPV")
