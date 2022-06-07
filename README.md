# SpliceJunctionsAnalysis
Scripts to analyse splice junctions. Based on the work by Osterreich et al. (2016)
This paper explores the question of how much coding potential is gained by alternative splicing (AS), and until what extend does AS affect the expansion for the genome.
It uses a splice site centric approach and Shannon entropy to analyse how different transcriptomes use different splice forms for its splice junctions.
As in the following diagram, a splice junction would be <em>pi</em> , a splice form would be <em>p1</em> , and the AS probability, <em>pj</em>:

<img src="splicejunction1.jpg">
![alt text](https://github.com/klari12/SpliceJunctionAnalysis/blob/main/splicejunction1.jpg?raw=true)

## Samples
To run Regtools, the program that will identify the splice junctions, we need BAM files as an input.
In this case BAM files from samples belonging to the ENCONDE project (https://www.encodeproject.org/) were used.
It is recommended to use paired-end reads, with sequences that are 100bp or longer to be able to detect as many splice junctions (new and annotated) as possible (Chhangawala, Rudy, Mason and Rosenfeld, 2015).

## Workflow
### Regtools
To carry out this analysis, we will use the output produced by Regtools (https://regtools.readthedocs.io/en/latest/). We will use its functions <em>junctions extract</em>  and <em>junctions annotate</em> to obtain the splice junctions from the BAM files of our samples of interest.
The output of the last function is a table that includes all splice forms including:
- Chromosome position
- Number of reads supporting the splice form
- Type of AS event
- Transcripts and genes that overlap the splice junction according to the input annotation (GTF) file
### Python script
With the python script included in this project we will be able to analyse the splice junctions:
- Filter splice junctions to keep only the ones that are supported by a minimum number of reads
- Obtain the number of splice forms present in each splice junction
- Analyse the Shannon entropy of all splice junctions
- Analyse the alternative splicing probability 
### R script
The R script attached will allow us to 
1) create the csv that is necessary to analyse Regtools output file with Python (pandas)
2) further analyse and create plots with the data generated by the python script
3) hypothesis testing in case we want to compare different samples between each other

