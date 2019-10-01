# Cell Type Modularity in Cell Autonomous Immunity

This is a progress report on my Biology thesis as of September 2019.


### Background
==========

Traditional views on immunity to infection promote the idea of a defense
repertoire composed solely of specialized immune cells. Yet most
nucleated cells are endowed with a capacity to defend themselves against
infection. This extended form of host defense -- termed cell-autonomous
immunity -- encompasses both immune and non-immune cell lineages that
protect the mammalian host from infection and is a greatly understudied
area of research. Nearly all of these lineages have immunoreceptors for
signaling host defense programs including the interferon (IFN), tumor
necrosis factor (TNF), interleukin-1 (IL-1), and Toll-like receptor
(TLR) families, which mobilize the transcription of hundreds of gene
products involved in resistance to important human pathogens.
These gene products exhibit a wide variety of activities, inhibiting
pathogens at almost every step of their replicative life cycle. For
example, some immune proteins block pathogen entry into the host cell,
while others may block microbial proliferation after pathogen uptake or
invasion. Remarkably, only a few of these host defense proteins
and their accompanying activities have been analyzed in detail, and most
of the genes elicited by these immunoreceptors still remain
uncharacterized. One particular aspect of this lack of characterization
involves exactly how this type of host immunity varies across cell
types: do all cell types respond similarly to infection, or are certain
responses cell type specific? This question has repercussions both in
understanding exactly which genes elicited in these responses are
responsible for these processes and in understanding how exactly cells
across different types evolved these defenses.

This exploration of the variation of cell autonomous immune processes
among cell types is the focus of this research. I hypothesized that
cells activated by the immunoreceptors typically secreted during
infection tend to suppress their specialized functions to concentrate
towards host defense, which will occur through broadly the same
processes (and therefore suppress genes encoding cell type specific gene
products and activate genes that encode products with similar functions,
or even activate the same proteins). The preliminary results gathered
seems to support this idea, but much more analysis into the individual
functions of the stimulated genes is required. This progress report will
show the initial partial results obtained, as well as describe and
comment on the adopted methodology to see what can be improved upon
going forward. Finally I will discuss the future work that will have to
be completed in order for this research to be satisfactory.

### Preliminary Results
===================

Initial steps of this research involved gathering enough raw genomic
data from different cell types for it to be possible to be looking
across a satisfactory amount of cell type differences. I looked at 19
different experiments analyzing the genomic response of different cell
types to IFN$\gamma$ in different dosages and after different time
periods. While I wanted to use only RNAseq data, there were not enough
experiments done on IFN$\gamma$ stimulation to gather enough cell types
for our purposes (this will be discussed in section 3). Therefore, I
decided to focus solely on microarray data to ensure that the datasets
could be reliably compared. After gathering the data, I found 8
experiments done on 8 different cell types where the time between
stimulation and microarray analysis averaged at 28.5 hours with a mode
of 24 hours (only 3 experiments diverged from 24 hours). Unfortunately,
IFN$\gamma$ dosage was much harder to control and much more variable. In
these eight experiments however, the average dosage was 40 ng/ml with a
mode of 10 ng/ml. These 8 experiments were therefore picked for all
subsequent analysis. The eight cell types they involved were the
following: A549 cells, bronchial epithelial cells
(referred to as AALEB), retinal epithelial cells (referred
to as ARPE), fibroblasts, hepatocytes, keratinocytes, microglial cells,
and thyroid epithelial cells.

After selection, I proceeded to gather the data from the GEO Accession
Browser . While most of it tended to have been previously
normalized for processing, as well as being processed to extract the
fold changes of the individual genes before and after stimulation, some
of it was only submitted using raw data. Therefore, I processed this
data in R using Bioconductor , specifically with the *affy*
package to process Affymetrix microarray data and the *lumi*
package to process Illumina microarray data. I used the BioMart
database and its associated Bioconductor package to fetch gene names and
standardize their usage across all 8 datasets, thus obtaining comparable
fold-change data indexed by Ensembl gene ID. I then filtered for
genes that either increased in expression more twice fold (which I
labeled as activated), or decreased by more than 0.5 fold (i.e. their
expression was diminished by a factor of less than 0.5, I labeled these
genes as suppressed). Once this labeling was completed, I generated two
27533 row and 8 column binary matrices indexed by cell types and gene
names. In one matrix, values of '1' imply the gene is activated in that
cell type upon stimulation. In the other, values of '1' imply the gene
is suppressed for that cell type upon stimulation. This data is the crux
of the analysis that follows.

Using the two binary matrices, I could calculate the number of times a
particular gene was activated or suppressed among the different cell
types, and genes are more 'shared' in
activation than they are in suppression. This seems to
point towards the idea that cells suppress different groups of genes but
activated similar ones when stimulated. I ran a Gene Ontology
enrichment analysis over both the activated and suppressed genes for all
8 cell types, concentrating on the Biological function annotation
tag. This was done using the R package *TopGO*. I gathered all
the statistically significant biological function terms that there were
for each cell type (for both activated and suppressed) and obtained a
count of how many times these terms would be shared among the different
cell types. The distributions between activated and suppressed are very different, the
terms associated with the activated genes being much more shared than
the terms associated with the suppressed genes. This seems to imply the
activated genes are also much similar in function. Finally, I just selected all of the genes whose annotation tag was
'transcription factor' and again checked their distributions depending
on whether they appeared in the activated or the suppressed groups.
Again, the activated transcription factors tend to be more shared among
cell types than the suppressed ones. This qualitative analysis also has
statistical power, I ran an un-paired, unequal sample size ANOVA on all
three activated/suppressed distributions and received a p-value of less
than 0.01 for all three of them, showing the difference in these
distributions is very significant.

The final kind of analysis I did after was to start integrating this
data, which as of now involved only isolated genes, into an interaction
network. This still needs to be explored but I will outline the process
that was followed doing so. This was done following what I learned by
reading Dmitry Zinoviev's _Complex Network Analysis in
Python_. Using the *Newtorkx* package
in combination with the "interactome" data that can be obtained through
the use of the BioPlex 2.0 database,[^1] I generated a Python
network object with 10961 nodes and 56553 edges and proceeded to
annotate it with the activation and suppression data I obtained
previously. I then performed two forms of very preliminary analysis. In
the first, I ranked the various nodes depending on various measures of
centrality (i.e. measures of how 'close' to each other these nodes are).
I then compared the values of the suppressed genes to the activated
genes, with the hypothesis that there would be differences between the
two groups of nodes. This approach did not support my hypothesis, 
because there was no difference, but I believe that might be due to looking at the data the wrong way.

### Issues in Methodology
=====================

The most salient issues in the methodology largely can be divided into
two main points.

The first is the use of these strong outliers in the data, the A549
suppressed genes and the AALEB genes. The A549 problem can be traced
back to the quality of the raw data as it contains only one trial making
the establishment of p-values for the fold change of the genes
impossible, possibly inflating the number of genes that seem to be
suppressed. Unfortunately, there are not a lot of way to correct this.
The processed data provided by the original experimenters, which might
have included more trials, provides absolute values of fold change and
therefore does not not distinguish between suppressed and activated
genes. Unfortunately, the way to fix this is
limited to observing whether the same conclusions that are obtained from
the data that includes the A549 genes are also supported when that data
is removed. This will be discussed in the next section. On the other
hand, once I had noticed the extremely low number of stimulated genes
for the AALEB cell type, I went back to observe whether the
normalization of the data had been done correctly. It resulted that
while I cannot correct the A549, the outlier nature of the AALEB data is
actually due to an error in the code that I had written to normalize
this data. This error resulted in the normalization not occurring, thus
making the data used to calculate fold changes skewed. In contrast with
the A549 data, this can easily and quickly be corrected.

The other methodological issue is more general. It is related to using
microarray data in comparison to RNAseq data. Microarrays delivered
comparatively much more background signal, from which it is harder to
extract expression variation, and microarrays only provide useful data
for known transcripts. RNAseq can in fact provide absolute
compared to relative values of expression, which can help get expression
values even for transcripts that are expressed at very small levels, and
thus determine how their expression changes upon stimulation more
reliably. The ideal data to perform this analysis on would be data that
results from independently running RNAseq expreiments on the stimulation
of a chosen number of cell types. Unfortunately, because of both time
constraint as well as financial constraints, this is not feasible.
Perhaps as RNAseq becomes more affordable this research can be validated
using that data.


### Future Work
===========

There is a lot of future work that has to be done. This can be divided
in to three general steps: first the correction of data, second the
analysis of function, and third the integration (and subsequent
verification) with the interactomic network. For the first point, the
data analysis workflow needs to be run again with corrected data, to see
if the results above stay more or less the same. This corrected data
consists of the removal of the A549 data and the renormalization of the
AALEB data. Therefore, there will be only 7 cell types this analysis
will be performed on. As well as seeing whether one obtains the same
results as above for the steps already completed, I do not believe one
should use the A549 data going forward unless the previous results are
confirmed (and even then, any step including that data should be
followed by the same step confirming the previous results by not
removing that data). The AALEB data on the other hand just simply needs
to be renormalized, and if the stimulated gene counts after
renormalization are not as skewed, the data can be retained for
analysis.

The second step is a deeper analysis of function of the genes, both
those activated and those suppressed. The data up until now has simply
focused on whether the genes (or their functional annotation) were more
or less shared among different cell types. These results have shown that
activated genes tend to be more shared both individually and
functionally than suppressed genes. However, there still has not been an
in depth about what these genes *do*! To fully verify the hypothesis
that cells types suppress their cell type specific activity to activate
what is mostly the same collection of genes that activate cell
autonomous immunity processes one needs to explore the Gene Ontology
enrichment analysis cell type by cell type, comparing which annotated
terms tend to be most activated or most suppressed not simply by raw
count but by keeping in mind the functions of the individual cell types.
Insight into what exactly the functions of the suppressed genes are can
tell us what activities the cell starts de-prioritizing to focus on
combating possible pathogen invasion.

Finally, the third step consists in integrating and verifying the
insight obtained by performing the second step into the interactome
network. The clusters we obtained lend themselves very easily to a kind
of functional analysis. If we identify a strongly activated or
suppressed gene among different cell types (or even a family of genes
with connected functional annotations), we can then look at the clusters
which involve this gene and see the function of the proteins it itself
interacts with! This could produce results such as perhaps whole
clusters of genes are stimulated in certain cell types and not others.
Possible rationalizations of why this happen can lead to further
downstream bench experiments. I have not been satisfied with the closeness 
measures analysis of the nodes in the network. This is probably because I
have not yet found a way to quantify satisfactorily how stimulated a
particular gene is across cell types (for example, if one gene is
activated in two cell types and suppressed in two others, the resulting
stimulation value, according to the way it has been measured up to now
will be zero, and yet this gene is obviously more significant to analyze
than a gene that is not stimulated in any direction in any cell type,
yet that will have a stimulation value of zero too). This could require
some more research in methods to measure combined state values, but will
only happen after the first part that was outlined in the third step.


[^1]: This database is a network of a large quantity (about 11,000) of
    the interacting proteins in the HEK293T cell type, ideal for my
    purposes.