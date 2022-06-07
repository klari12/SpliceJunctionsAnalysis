# SpliceJunctionsAnalysis
Scripts to analyze splice junctions. Based on the work by <a href="https://www.biorxiv.org/content/10.1101/048124v1" target="_blank" >Osterreich el al. (2016)</a> .
This paper explores the question of how much coding potential is gained by alternative splicing (AS), and until what extent does AS affect the expansion for the genome.
It uses a splice site centric approach and Shannon entropy to analyze how different transcriptomes use different splice forms for its splice junctions.
As in the following diagram, a splice junction would be <em>pi</em> , a splice form would be <em>p1</em> , and the AS probability, <em>pj</em>:

<img src="https://github.com/klari12/SpliceJunctionsAnalysis/blob/main/splicejunction1.png">

## Samples
To run Regtools, the program that will identify the splice junctions, we need BAM files as an input.
In this case BAM files from samples belonging to the ENCONDE project (https://www.encodeproject.org/) were used.
It is recommended to use paired-end reads, with sequences that are at least 100bp or longer to be able to detect as many splice junctions (new and annotated) as possible (Chhangawala, Rudy, Mason and Rosenfeld, 2015).
The ENCODE project used as an example here is a K562 cell line sample (https://www.encodeproject.org/experiments/ENCSR071ZMO/) and a particular BAM file of one of the replicates (https://www.encodeproject.org/files/ENCFF201FDK/)

## Workflow
### Regtools
To carry out this analysis, we will use the output produced by Regtools (https://regtools.readthedocs.io/en/latest/). We will use its functions <em>junctions extract</em>  and <em>junctions annotate</em> to obtain all splice junctions, new and annotated, from the BAM files of our samples of interest.
The last function's output is a table that includes all splice forms and the following information for each one of them:
- Chromosome position
- Number of reads supporting the splice form
- Type of AS event
- Transcripts and genes that overlap the splice junction according to the input annotation (GTF) file
### Python script
The python script included in this project can be used to analyze the file that contains splice junctions created by Regtools. It can:
- Filter splice junctions to keep only the ones that are supported by a minimum number of reads in two different ways: 1) filter by the total number of reads that support a splice junctions (<em>new_version</em>) or 2) filter by the number of reads that support each splice form, disregarding the splice junction where they belong (<em>normalized_...</em>)
- Obtain the number of splice forms present in each splice junction
- Obtain Shannon entropy values of all splice junctions
- Obtain the alternative splicing probability for every splice junction.
It will create four output files including:
- Number of splice forms for each splice junctions
- Normalized Shannon entropy values for each splice junction 
- Not-normalized Shannon entropy values
- Alternative splicing probabilities for each splice junction
### R script
The R script attached will allow us to:
1) Create the csv that is necessary to analyze the Regtools output file with Python (pandas)
2) Further analyse and create plots using the files generated by the python script
3) Carry out hypothesis testing in case we want to compare the percentage of high entropy events between samples.

