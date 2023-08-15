# Differential-expression-and-functional-enrichment-analysis-for-RNA-sequence
In this tutorial we will do analysis for RNA seq data by using different packages in R and get our volcano plot and differentially expressed gene and then get functional enrichment gene n the pathway  
# Loading relevant libraries 
library(readr)
library(rhdf5)
library(tidyverse)
library(tximport)
library(ensembldb) 
library(EnsDb.Hsapiens.v86) 
library(beepr) 
library(matrixStats)
library(edgeR)
library(cowplot)

#set a directory
setwd("C:/Users/ahmed.DESKTOP-2322KN8/Downloads/AI-Biology")

#create a study design 

sample_name=list.files("C:/Users/ahmed.DESKTOP-2322KN8/Downloads/AI-Biology")
cell_type=c(rep("neuroplastoma",39), rep("non_neuroplastoma" ,2))
study.design=data.frame(sample_name,cell_type)

#creat a path file 
abundance.path=file.path(study.design$sample_name,"abundance.tsv")

#chek the file path
all(file.exists(abundance.path))

# Getting annotation files ----
Tx <- transcripts(EnsDb.Hsapiens.v86, columns=c("tx_id", "gene_name"))

# creat an object named with Tx to get annotation file, specify two columns to get with you (transcript_id, gene_name)
Tx <- as_tibble(Tx) #to transfer data to tidytable
#note: 
#The tbl_df class is a subclass of data.frame, created in order to have different default behaviour. The colloquial term "tibble" refers to a data frame that has the tbl_df class. Tibble is the central data structure for the set of packages known as the tidyverse, including dplyr, ggplot2, tidyr, and readr.
#need to change first column name to 'target_id' for further matching 
Tx <- dplyr::rename(Tx, target_id = tx_id)
#transcrip ID needs to be the first column in the dataframe (select transcript ID and gene ID only)
Tx <- dplyr::select(Tx, "target_id", "gene_name")

#importing data of kalisto using teximport----
Txi_gene <- tximport(abundance.path, #files
                     type = "kallisto", 
                     tx2gene = Tx, #transcript and gene
                     txOut = FALSE, #How does the result change if this =FALSE vs =TRUE? #False means output will be gene while True means transcript
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE) # to ignore version of transcript

# Examine your data up to this point -----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts

#restore the sample names 
colnames(myCounts)=sample_name

myDGEList <- DGEList(myCounts) 

#Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts. Next we’ll use box/violn plots to check the distribution of the read counts on the log2 scale. We can use the cpm function to get log2 counts per million, which are corrected for the different library sizes
cpm <- cpm(myDGEList) #to extract count per milion

#convert into log 
log2.cpm <- cpm(myDGEList, log=TRUE)

# to work on tidyversy you need dataframe so transfer into dataframe
log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID") # to take care about row names (gene id)

# use the tidy package to 'pivot' (organize and make it more handy) your dataframe (from wide to long)
log2.cpm.df.pivot <- pivot_longer(log2.cpm.df,
                                  cols = "neuroblastoma1":"non-neuroblastoma2",
                                  names_to = "samples", 
                                  values_to = "expression") 
#DO violin plot
ggplot(log2.cpm.df.pivot) +
    aes(x=samples, y=expression, fill=samples) + # fill the color based on sample
    geom_violin(trim = FALSE, show.legend = FALSE) + # shape of the plot
    stat_summary(fun = "median", # to show summary
                 geom = "point", 
                 shape = 95, 
                 size = 10, 
                 color = "black", 
                 show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
         title="Log2 Counts per Million (CPM)",
         subtitle="unfiltered, non-normalized",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw() 

# is there is any genes that has 0 count ?
table(rowSums(myDGEList$counts==0)==41)

#FALSE  #TRUE 
#33013  #2357 

#The line below is important! This is where the filtering starts
# Be sure to adjust this cutoff for the number of samples in the smallest group of comparison.
# how many gene that has greater than 1 
cpm #datamatrix
keepers <- rowSums(cpm>1)>=5 # we write it 5 not 10 to not cause biase
keepers
# now use base R's simple subsetting method to filter your DGEList based on the logical produced above
myDGEList.filtered <- myDGEList[keepers,] #to subsetting data (give me the rows that has True value)
myDGEList.filtered
myDGEList
dim(cpm)  #35370    #41
dim(myDGEList.filtered) # to check dimension  #17715    #41

log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
log2.cpm.filtered.df

# pivot this FILTERED data, just as you did earlier
log2.cpm.filtered.df.pivot <- pivot_longer(log2.cpm.filtered.df, # dataframe to be pivoted
                                           cols = "neuroblastoma1":"non-neuroblastoma2", # column names to be stored as a SINGLE variable
                                           names_to = "samples", # name of that new variable (column)
                                           values_to = "expression") # name of new variable (column) storing all the values (data)

#ploting filtered data
ggplot(log2.cpm.filtered.df.pivot) +
    aes(x=samples, y=expression, fill=samples) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 shape = 95, 
                 size = 10, 
                 color = "black", 
                 show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
         title="Log2 Counts per Million (CPM)",
         subtitle="filtered, non-normalized",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()

# Normalize your data ----
#from edger use calcNOrmfactor to normalize ourdata #it us TMM normalization that normalize between sample
myDGEList.filtered
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM") #Trimmed mean of M 
myDGEList.filtered.norm

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE) # to extract count
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID") #to transform intodatafram

## pivot this NORMALIZED data, just as you did earlier
log2.cpm.filtered.norm.df.pivot <- pivot_longer(log2.cpm.filtered.norm.df, # dataframe to be pivoted
                                                cols = "neuroblastoma1":"non-neuroblastoma2", # column names to be stored as a SINGLE variable
                                                names_to = "samples", # name of that new variable (column)
                                                values_to = "expression") # name of new variable (column) storing all the values (data)

#ploting for normalized and filered data
ggplot(log2.cpm.filtered.norm.df.pivot) +
    aes(x=samples, y=expression, fill=samples) +
    geom_violin(trim = FALSE, show.legend = FALSE) +
    stat_summary(fun = "median", 
                 geom = "point", 
                 shape = 95, 
                 size = 10, 
                 color = "black", 
                 show.legend = FALSE) +
    labs(y="log2 expression", x = "sample",
         title="Log2 Counts per Million (CPM)",
         subtitle="filtered, TMM normalized",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()

# Set up your design matrix ----
group <- factor(targets$group) #catgory for data as we did before
design <- model.matrix(~0 + group) #use model matrix function to create model matrix that will fit linnear model, use 0 for no intercept and group for the data 
colnames(design) <- levels(group)

# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design) # creat linear model by function called lmfit

# Contrast matrix ----
contrast.matrix <- makeContrasts(difference = non_neuroplastoma -  neuroplastoma, # creat matrix that have disease and healty based on model matrix
                                 levels=design)

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix) # to extact linerar model based on matrix

#get bayesian stats for your linear model fit
ebFit <- eBayes(fits) # to get basic stastistcs
write.fit(ebFit, file="lmfit_results.txt") # to save file of stastistics 

# TopTable to view DEGs -----
myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, sort.by="logFC") #represent data in a table

# convert to a tibble
myTopHits.df <- myTopHits %>%
    as_tibble(rownames = "geneID") # to transform into datafram

# now plot
ggplot(myTopHits.df) +
    aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
    geom_point(size=2) +
    geom_hline(yintercept = -log10(0.05), linetype="longdash", colour="grey", size=1) +
    geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
    geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
    annotate("rect", xmin = 1, xmax = 12, ymin = -log10(0.05), ymax = 7.5, alpha=.2, fill="#BE684D") +
    annotate("rect", xmin = -1, xmax = -12, ymin = -log10(0.01), ymax = 7.5, alpha=.2, fill="#2C467A") +
    labs(title="Volcano plot",
         subtitle = "Cutaneous leishmaniasis",
         caption=paste0("produced on ", Sys.time())) +
    theme_bw()

#extract genes -----
# decideTests to pull out the DEGs and make Venn Diagram ----
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=1)

# take a look at what the results of decideTests looks like
head(results)
summary(results)
vennDiagram(results, include="up") 

# retrieve expression data for your DEGs ---
diffGenes <- v.DEGList.filtered.norm$E[results[,1] !=0,] #filtering 
dim(diffGenes)
#308  #41
head(diffGenes)
dim(diffGenes)

#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

#choosing the top 10 p value
head(arrange(myTopHits.df,adj.P.Val),10)
thresh=head(arrange(myTopHits.df,adj.P.Val),10)$adj.P.Val[10]
myTopHits.df$delabel[myTopHits.df$adj.P.Val <=thresh]= (myTopHits.df$geneID[myTopHits.df$adj.P.Val <=thresh])

#Choosing the up and down regulated genes
myTopHits.df$diffexpressed=NA #add coulmn called "diffexpressed"
myTopHits.df$diffexpressed[myTopHits.df$logFC> 0.6 &myTopHits.df$adj.P.Val< 0.05]="UP"
myTopHits.df$diffexpressed[myTopHits.df$logFC< -0.6 &myTopHits.df$adj.P.Val< 0.05]="DOWN"

#enhanced plot
ggplot(data = myTopHits.df, aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed, label= delabel)) +
    geom_vline(xintercept = c(-0.6, 0.6), col = "red", linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), col = "red", linetype = 'dashed') +
    geom_point(size = 2) +
    
    scale_color_manual(values = c("#00AFBB", "grey", "red"), # to set the colours of our variable
                       labels = c("Downregulated", "Not significant", "Upregulated")) + # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
    ggtitle("neuroplastoma vs non-neuroplastoma") + #to add my title
    geom_text_repel(max.overlaps = Inf)

#functional enrichment
library(tidyverse)
library(gplots) #for heatmaps
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources

myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=17715 , sort.by="logFC")

#determint the threshold 
myTopHits.df.mod=myTopHits.df[abs(myTopHits.df$logFC)>=(2) & myTopHits.df$adj.P.Val< 0.05,]

#use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost.res <- gost(myTopHits.df.mod$geneID, organism = "hsapiens", correction_method = "fdr")
data=gost.res$result

#determine the terms
data_GO_BP=data[data$source=="GO:BP"&data$p_value<0.05,]

#table of the determine pathway
publish_gosttable(
    data_GO_BP,
    highlight_terms = NULL,
    use_colors = TRUE,
    show_columns = c("source", "term_name", "term_size", "intersection_size"),
    filename = NULL,
    ggplot=TRUE)

#gene name by using evcodes
gost.res.mod <- gost(myTopHits.df.mod$geneID, organism = "hsapiens", correction_method = "fdr", evcodes = TRUE)

#choosing GO
data_GO_BP.mod=data.mod[data.mod$source=="GO:BP"&data$p_value<0.05,]

#Reusut in table 
my_table=data_GO_BP.mod[,c(4,5,6,9,10,11,16,3)]
write.csv(my_table, file = "my_table.csv")































