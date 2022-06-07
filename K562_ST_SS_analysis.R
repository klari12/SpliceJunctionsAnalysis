## K562 sample
K562_AJ <- read.table("AnnoJunctionsK562.txt", header = TRUE)

#Libraries
library(dplyr)
library(ggplot2)
library(magrittr)
library(scales)
library(hrbrthemes)
# For creating a csv file that can be read by python's pandas:

#write.csv(K562_AJ,"AnnoJunctions_K562.csv", row.names= FALSE)

K562_AJ <- as.data.frame(K562_AJ)
dim(K562_AJ)
View(K562_AJ)
summary(K562_AJ$score)


## We will look out for some individual genes to observe how many splice forms are found:
wash7p_k <- K562_AJ[K562_AJ$gene_ids=="ENSG00000227232.5",]
wash7p_k <- na.omit(wash7p_k)
nrow(wash7p_k)
View(wash7p_k)


apol1k <- K562_AJ[K562_AJ$gene_names=="APOL1",]
apol1k <- na.omit(apol1)
nrow(apol1k)

ift27k <- K562_AJ[K562_AJ$gene_names=="IFT27",]
ift27k <- na.omit(ift27k)
nrow(ift27k)


tmbim1 <- K562_AJ[K562_AJ$gene_names=="TMBIM1",]
tmbim1 <- na.omit(tmbim1)
nrow(tmbim1)

# Number of spliced forms when the score count is filtered >0, >5, and >10
K562_0 <- read.csv("number_splice_forms_K562_0.csv")
head(K562_0)
View(K562_0)


K562_5 <- read.csv("number_splice_forms_K562_5.csv")
head(K562_5)
View(K562_5)

K562_10 <- read.csv("number_splice_forms_K562_10.csv")
head(K562_10)
View(K562_10)

help(tabulate)
freqs <- tabulate(K562_10$number_spliced_forms)
probs <- freqs / sum(freqs)

#Visualization of the number of the frequency of the number of splice forms per splice junction
ggplot(K562_10, aes(x = number_spliced_forms)) + geom_point(stat = "count") +
  scale_y_log10() + 
  scale_x_log10()

# Frequency plot with log transformed y axis
ggplot() + 
  geom_point(data=K562_0, aes(x= number_spliced_forms), color="red", stat="count")+
  geom_point(data=K562_5, aes(x= number_spliced_forms),color="blue", stat="count") + 
  geom_point(data=K562_10, aes(x= number_spliced_forms),color="green", stat="count")+
  scale_y_log10("Frequency", breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  scale_x_continuous("Number of spliced forms", 0:25) 

## In order to compare samples not just visually, the exponential fit and decay rate should be calculated.


#### Shannon entropy analysis
# We will analyse AS events using Shannon entropy to observe how evenly are different AS states used in the sample.

## First type of filtering: after the sum of the common scores of a splice junction
normalized_entropy_K562_n <- read.csv("new_version_normalized_entropy_groups_K562_10.csv")
head(normalized_entropy_K562_n)
View(normalized_entropy_K562_n)


normalized_entropy_K562_n$X <- NULL
head(normalized_entropy_K562_n)

normalized_entropy_K562_n <- normalized_entropy_K562_n %>% 
  rename(
    Norm_entropy = X0
  )
head(normalized_entropy_K562_n)
mean(normalized_entropy_K562_n$Norm_entropy)
summary(normalized_entropy_K562_n)
dim(normalized_entropy_K562_n)
library(ggplot2)

ggplot(normalized_entropy_K562_10, aes(Norm_entropy)) +
  geom_histogram(aes(y = stat(count)),color="black", fill="white", bins=30) + 
  scale_y_continuous("Frequency")  +
  scale_x_continuous("Normalized entropy value")

ggplot(normalized_entropy_K562_10, aes(x = Norm_entropy)) + 
  geom_freqpoly(bins=50)

boxplot(normalized_entropy_K562_10)


## Second type of filtering: before joining the junctions
# >= 0 threshold
normalized_entropy_K562_0 <- read.csv("normalized_entropy_groups_K562_0.csv")
head(normalized_entropy_K562_0)
View(normalized_entropy_K562_0)


normalized_entropy_K562_0$X <- NULL
head(normalized_entropy_K562_0)

normalized_entropy_K562_0 <- normalized_entropy_K562_0 %>% 
  rename(
    Norm_entropy = X0
  )
head(normalized_entropy_K562_0)
mean(normalized_entropy_K562_0$Norm_entropy)
summary(normalized_entropy_K562_0)
dim(normalized_entropy_K562_0)

# Visualization
ggplot(normalized_entropy_K562_0, aes(Norm_entropy)) +
  geom_histogram(aes(y = stat(count)),color="black", fill="white", bins=30) + 
  scale_y_continuous("Frequency")  +
  scale_x_continuous("Normalized entropy value")

ggplot(normalized_entropy_K562_0, aes(x = Norm_entropy)) + 
  geom_freqpoly(bins=50)

boxplot(normalized_entropy_K562_0)

# >= 5 threshold
normalized_entropy_K562_5 <- read.csv("normalized_entropy_groups_K562_5.csv")
head(normalized_entropy_K562_5)
View(normalized_entropy_K562_5)


normalized_entropy_K562_5$X <- NULL
head(normalized_entropy_K562_5)

normalized_entropy_K562_5 <- normalized_entropy_K562_5 %>% 
  rename(
    Norm_entropy = X0
  )
head(normalized_entropy_K562_5)
mean(normalized_entropy_K562_5$Norm_entropy)
summary(normalized_entropy_K562_5)
dim(normalized_entropy_K562_5)

ggplot(normalized_entropy_K562_5, aes(Norm_entropy)) +
  geom_histogram(aes(y = stat(count)),color="black", fill="white", bins=30) + 
  scale_y_continuous("Frequency")  +
  scale_x_continuous("Normalized entropy value")

ggplot(normalized_entropy_K562_5, aes(x = Norm_entropy)) + 
  geom_freqpoly(bins=50)

boxplot(normalized_entropy_K562_5)

# >= 10 threshold

#### Not normalized Shannon entropy analysis
# This data allows us to compare our results from those produced by Whippet.
not_normalized_entropy_K562 <- read.csv("not_normalized_entropy_K562_5.csv")
head(not_normalized_entropy_K562)
View(not_normalized_entropy_K562)
dim(not_normalized_entropy_K562)

not_normalized_entropy_K562$X <- NULL
head(not_normalized_entropy_K562)

not_normalized_entropy_K562 <-not_normalized_entropy_K562 %>% 
  rename(
    Norm_entropy = X0
  )
head(not_normalized_entropy_K562)
mean(not_normalized_entropy_K562$Norm_entropy)
summary(not_normalized_entropy_K562)


ggplot(not_normalized_entropy_K562, aes(Norm_entropy)) +
  geom_histogram(aes(y = stat(count)), color="black", fill="white") +   scale_y_continuous("Frequency")  +
  scale_x_continuous("Entropy value") 

# To know the percentage of high-entropy splice junctions and save them in a dataset for further analysis
summary(not_normalized_entropy_K562)
high_entropy <- length(not_normalized_entropy_K562[((not_normalized_entropy_K562$Norm_entropy) >= 1),])
high_entropy/length(not_normalized_entropy_K562$Norm_entropy)
dataset_he_k562 <- not_normalized_entropy_K562[((not_normalized_entropy_K562$Norm_entropy) >= 1),]

### Alternative splicing probability 
SDprob <- read.csv("alt_splicing_prob_K562_5.csv")
head(SDprob)
mean(SDprob)
View(SDprob)
dim(SDprob)

SDprob$X <- NULL
head(SDprob)

SDprob <-SDprob %>% 
  rename(
    Alternative_splicing_probability = X0
  )
head(SDprob)
mean(SDprob$Alternative_splicing_probability)
summary(SDprob)

## This plot has two major peaks: at very low and at very high numbers. It shows 
# that there is a tendency to have very low, close to 0 probability values (minor forms)
# while having another predominant, major form.
p <- ggplot(SDprob, aes(x=Alternative_splicing_probability)) + 
  geom_density() + scale_x_continuous("Alternative splicing probability")
p


ggplot(SDprob, aes(Alternative_splicing_probability)) +
  geom_histogram(aes(y = stat(count)), bins=30, binwidth=0.05, fill="#404080", color="#404080", alpha=0.6) +
  scale_y_log10("Density", 0:1)  +
  scale_x_log10("Alternative splicing probability")
  