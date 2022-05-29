library(Rsubread)

## Data Import
# We can search for all .fastq.gz files in the data directory using the list.files command. The pattern argument takes a regular expression. In this case we are using the $ to mean the end of the string, so we will only get files that end in “.fastq.gz”
fastq.files <- list.files(path = "./data", pattern = ".fastq.gz$", full.names = TRUE)
fastq.files




## Alignment
#Build Index
buildindex(basename="hg38",reference="./genome/hg38.analysisSet.fa", 
           indexSplit = TRUE,
           gappedIndex = TRUE,
           memory=4000)



# Aligning reads to  reference genome: Now we can align our 12 fastq.gz files using the align command.
align(index="hg38",readfile1=fastq.files)

## List BAM Files
bam.files <- list.files(path = "./data", pattern = ".BAM$", full.names = TRUE)
bam.files


props <- propmapped(files=bam.files)
props


## Quality Controlm
qs <- qualityScores(filename="./data/046_S1_L001_R1_001.fastq.gz",nreads=100)
# Check dimension of qs
dim(qs)
# Check first few elements of qs with head
head(qs)


boxplot(qs)

fc <- featureCounts(bam.files, annot.inbuilt="hg38",
                    annot.ext="./annotation/hsa.gff3",
                    isGTFAnnotationFile=TRUE,
                    GTF.featureType="miRNA", GTF.attrType="ID")

fc$stat

## Take a look at the dimensions to see the number of genes
dim(fc$counts)

## Take a look at the first 6 lines
head(fc$counts)
head(fc$annotation)
fc_counts <- as.data.frame(fc$counts)
names(fc_counts) <- gsub("_L001_R1_001.fastq.gz.subread.BAM","", names(fc_counts))
names(fc_counts) <- gsub("_L002_R1_001.fastq.gz.subread.BAM","", names(fc_counts))
names(fc_counts) <- gsub(".*_","",names(fc_counts))
save(fc_counts,fc, props,qs, file="miRNA-seq.Rdata")




