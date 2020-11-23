# camponotus_genomes1
Co-evolution of Camponotus and their endosymbionts

### 01 - Process to assemble the Camponotus de novo genome
    00 - Convert the pacbio bam files
    01 - Run Canu for initial assembly
    02 - Several steps to run the 3d-dna pipeline for scaffolding with Hi-C data
    02b - Run blobtools to identify and filter potential contamination
    03 - Process to annotate repetitive and transposable elements
    04 - Process to annotate genes
    05 - Steps to calculate the Camponotus molecular clock
    06 - summarize and plot genome summaries in 100 kbp windows

### 02 - Assembly and analysis of Blochmannia genomes
    01 - rename files
    02 - run the minys pipeline
    03 - notes about processing the minys assemlies in bandage and geneious
    04 - run pgap annotation pipeline
    05 - extract coding regions for each genome
    06 - blast all cds to Blochmannia floridanus reference
    07 - pull out all seqs with matching blast hits into a single file
    08 - align sequences in each file
    09 - estimate a phylogenetic tree for each alignment
    10 - combine all the trees into a single file
    11 - run species tree analyses
    12 - tests for selection in each gene using hyphy
    13 - analyses looking at the pangenome (e.g., gene loss)
    14 - exploratory look at diversity in Blochmannia genomes (not used)
    15 - compare host and endosymbionts phylogenies

### 03 - Camponotus resequencing and bioinformatics
    01 - rename files for batch scripts
    02 - index the reference genome
    03 - trim and align all individuals to reference
    04 - R script to create batch cluster genotyping scripts
    04b,c - calculate and plot sequencing coverage for each individual
    05 - prep vcfs for filtering
    06 - R script to create batch cluster filtering scripts
    07 - subset all vcf files into 100 kbp windows
    08 - setup step(s)
    09 - calculate stats for each window
    10 - prepare files for phylogenetic analyses
    11 - phylogenetic analyses
    12 - combine all phylogenies
    13 - combine all window stats calculations
    14 - exploratory analysis of runs of homozygosity (not used)
  
### 04 - Camponotus demographic analyses
    01 - Prep vcf files
    02 - Create bootstrap files
    03 - Run MSMC
    04 - plot output and make pop. size calculations
