library(SwabSeqModular)
#install.packages('readxl')
library(readxl)
#if Biostrings isn't installed
#install.packages("BiocManager")
#BiocManager::install("Biostrings")
library(Biostrings)
#conda-forge ggplot2 is broken, need to install scales package manually 
#install.packages('scales')
library(ggplot2)

#local-----------------------
#setwd('~/.basespace/')
#system(paste0('bs -c swabseqPersonal.cfg list runs'))
#runID='284756474'
#bcl.save.dir='/data3/swabseq/240827/'
#dir.create(bcl.save.dir)
#status=system(paste0("bs -c swabseqPersonal.cfg  download run --id ", runID, " -o ", bcl.save.dir))
#setwd(bcl.save.dir)
#sampleSheet=paste0(bcl.save.dir, 'SampleSheet.csv')
#status=as.character(system(paste("bcl-convert --bcl-input-directory . --output-directory out/ --strict-mode false --no-lane-splitting true --force --sample-sheet", sampleSheet),intern=F))
#path.to.index='/home/jbloom/Dropbox/code/SwabSeqModular/projects/240827_swabseqv1test_TJ/'

# hoffman2 -----------------------------------------------------
# so far the gotchas I've seen are that you need to put bs and bcl-convert in your PATH
# the basespace cli config file needs to be in the ~/.basespace directory or it can't be read ...
# also you need to edit your sample sheet to have it listed as index instead of index1

# always define this variable, recommend downloading and bcl-convert in scratch/
bcl.save.dir='/u/home/j/jsbloom/scratch/swabseq/240827/'

# this chunk only needs to run once per sequencing run------------------------------------------
setwd('/u/home/j/jsbloom/.basespace/')
system(paste0('bs -c swabseqPersonal.cfg list runs'))
runID='284756474'
dir.create(bcl.save.dir)
status=system(paste0("bs -c swabseqPersonal.cfg  download run --id ", runID, " -o ", bcl.save.dir))
setwd(bcl.save.dir)
sampleSheet=paste0(bcl.save.dir, 'SampleSheet.csv')
status=as.character(system(paste("bcl-convert --bcl-input-directory . --output-directory out/ --strict-mode false --no-lane-splitting true --force --sample-sheet", sampleSheet),intern=F))
# -----------------------------------------------------------------------------------------------

#load in index file and rejigger --------------------------------------
path.to.index='/u/home/j/jsbloom/project-kruglyak/swabseq/SwabSeqModular/projects/240827_swabseqv1test_TJ/'
index.key=read_excel(paste0(path.to.index, 'SwabSeqV1_Index.xlsx'))
names(index.key)[4:5]=c('index', 'index2')

index.key=index.key[,c(1,2,3,5,4)]
names(index.key)[4:5]=c('index', 'index2')

index.key$index=as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(index.key$index)))
index.key$index2=as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(index.key$index2)))
#------------------------------------------------------------------------
    
amplicons=list(
  S2 =          'TATCTTCAACCTAGGACTTTTCTATT',
  S2_spike000=  'ATAGAACAACCTAGGACTTTTCTATT',
  S2_spike001=  'GTGTATCTCACGAAGCGACCCTTTGG',
  S2_spike002=  'CCTCGCTAGGACGTCGCTATGACGCC',
  S2_spike003=  'AGCACGACTTGATCTAACTGACACTA',
  S2_spike004=  'TAAGTAGGACTTCGATTGGATGGAAT',
  RPP30=        'TCTGACCTGAAGGCTCTGCGCGGACT',
  N=            'CCCCGCATTACGTTTGGTGGACCCTC',
  N_spike000=   'GGGGCGTAATGCAAACCACCTGGCTC',
  fluA=         'TCGCTCACTGGGCACGGTGAGCGTGA',
  fluA_spike=   'TCGCTCACTTGGCACGGTGAGCGTGA',
  fluB=         'ACTCTTCGAGCGTTTTGATGAAGGAC',
  fluB_spike000='TGAGAAGCTCGCAAAACTACTTCCTG')


in.con=gzfile(paste0(bcl.save.dir, 'out/Undetermined_S0_R1_001.fastq.gz'), 'rt')
#line.buffer is 4x the number of reads you want to read at a time 
#use line.buffer to control memory usage 
#load everything
#note, running this function closes the connection in.con
#in.con=file(paste0(bcl.save.dir, 'out/Undetermined_S0_R1_001.fastq'), 'rt')
# read the whole thing
countTables=countAmplicons(in.con, index.key, amplicons, line.buffer=5e7, nthreads=4)
saveRDS(countTables, file = '/home/jbloom/Dropbox/code/SwabSeqModular/projects/240827_swabseqv1test_TJ/countTables.RDS')

#for debugging, load first 2.5e6 reads only for example
countTables2=countAmplicons(in.con, index.key, amplicons, line.buffer=40e6/4, max.lines=10e6/4, nthreads=4)

#note, replace countTables2 with countTables to visualize the full data set 
x=data.table::rbindlist(countTables2$count.tables)
x$Col=as.factor(gsub('^.', '', x$Sample_Well))
x$Row=factor(gsub('..$', '', x$Sample_Well), levels=rev(toupper(letters[1:16])))
x %>%  ggp2::ggplot(ggplot2::aes(x=Col, y=Row, fill=log10(Count))) + 
  ggplot2::geom_raster() +
  ggplot2::coord_equal() +
  ggplot2::facet_grid(amplicon~Plate_ID) +
  ggplot2::scale_fill_viridis_c(option = 'plasma')+ggplot2::ggtitle('') 
