# This code uses a Dirchlet distribution to draw successive samples
#for creating control random data for the neutral model 


install.packages("MCMCpack")
library(MCMCpack)
#file names
metadata_file<- "colonization.metadata.final2.txt"
shared_file <- "2colonization.subsample.shared"

#read in files
metadata<- read.table(file=metadata_file, header=T, stringsAsFactors = FALSE)
shared<- read.table(file=shared_file, header=T, row.names=2)
shared <- shared[, c(-1,-2)] #remove unneeded columns from shared file

#subset shared file based on rownames for day 1 from metadata 
day1_ids <- metadata$id[metadata$day=='D01']
day1_shared <- subset(shared, rownames(shared) %in% day1_ids)

#convert to abundance
day1_abund <- (day1_shared/2400)*100

#init empty df for storing data
npops <- length(day1_abund) #number of OTUs
full_gen <- matrix(ncol=npops)

#loop to create distribution 

#loop through each mouse 
for(m in day1_ids){
init_distribution <- as.numeric(day1_abund[m,]) #pull out abund for one mouse at a time 
generations <- 21
strength <- 2400 #number of sequences, analogous to Ntm
gen <- matrix(rep(0, generations*npops), ncol=npops)
gen[1,] <- round(strength * init_distribution)	#set the initial frequencies 

#loop through each generation to create random distribution for each day
for(i in 2:generations){
  gen[i,] <- round(rdirichlet(1, unlist(gen[i-1,]) )* strength)
  gen <- as.data.frame(gen)
  rownames(gen)[i-1] <- paste(m,"day",i-1, sep="_")
}
rownames(gen)[generations] <- paste(m, "day", generations, sep="_")
full_gen <- rbind(full_gen, gen)
}
#clean up df
full_gen <- full_gen[-1,]
colnames(full_gen) <- colnames(day1_abund)


#need to reformat rownames so they match day ids in metadata to use in nm script 

full_gen_renamed <- rownames_to_column(full_gen, var="id")
full_gen_renamed <- separate(full_gen_renamed, id, into=c("name", "day"), sep="_d", remove=TRUE)
full_gen_renamed <- separate(full_gen_renamed, day, into=c("dayB","dayC"), sep="_", remove=TRUE)
full_gen_renamed <- full_gen_renamed[-2]
full_gen_renamed <- separate(full_gen_renamed, name, into="id", sep=".D01", remove=TRUE)
full_gen_renamed$dayC <- as.numeric(full_gen_renamed$dayC)

full_gen_renamed$dayC <- ifelse(full_gen_renamed$dayC < 10, paste("D0", full_gen_renamed$dayC, sep=""), paste("D", full_gen_renamed$dayC, sep=""))

full_gen_renamed <- unite(full_gen_renamed, "sample_id", c(id, dayC), sep=".")
rownames(full_gen_renamed) <- NULL 
full_gen_renamed <- column_to_rownames(full_gen_renamed, var="sample_id")


