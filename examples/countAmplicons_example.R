dir.path='/home/jbloom/Dropbox/code/SwabSeqModular/'
#source(paste0(dir.path,'/R/coreFunctions.R'))

amplicons=list(
            S2=         'TATCTTCAACCTAGGACTTTTCTATT',
            S2_spike000='ATAGAACAACCTAGGACTTTTCTATT',
            S2_spike001='GTGTATCTCACGAAGCGACCCTTTGG',
            S2_spike002='CCTCGCTAGGACGTCGCTATGACGCC',
            S2_spike003='AGCACGACTTGATCTAACTGACACTA',
            S2_spike004='TAAGTAGGACTTCGATTGGATGGAAT',
            RPP30      ='CGCAGAGCCTTCAGGTCAGAACCCGC'
        )

cfg=list()
cfg$i7_plate_key_file =  paste0(dir.path,"/data/s2_r.csv") #, package="swabseqr")
cfg$i5_plate_key_file =  paste0(dir.path,"/data/s2_f.csv") 
#system.file("keys", "s2_f.csv", package="swabseqr")

index.key=generateExpectedIndices(cfg)
#head(index.key)
#   Plate_ID Sample_Well  Sample_ID      index     index2
#     <char>      <char>     <char>     <char>     <char>
#1:   Plate1         A01 Plate1-A01 GATACAGATG CAAGTCTCAT
#2:   Plate1         A02 Plate1-A02 TCTAAGCCTC AGATAGCCTC
#3:   Plate1         A03 Plate1-A03 GAACGGACAG TCTTCAACGT
#4:   Plate1         A04 Plate1-A04 ACGATGTTAC CAGAGAATCT
#5:   Plate1         A05 Plate1-A05 ACCAAGGTGG TCTTGCCGTT
#6:   Plate1         A06 Plate1-A06 ACTCTTGTAG TTAGATCGCT

in.con=gzfile('/home/jbloom/Dropbox/FileTransfer/swabseq/MINISEQ_Undetermined_S0_R1_001.fastq.gz', 'rt')
#line.buffer is 4x the number of reads you want to read at a time 
#use line.buffer to control memory usage 
countTables=countAmplicons(in.con, index.key, amplicons, line.buffer=1e7)





#example visualization 
x=data.table::rbindlist(countTables$count.tables)
x$Col=as.factor(gsub('^.', '', x$Sample_Well))
x$Row=factor(gsub('..$', '', x$Sample_Well), levels=rev(toupper(letters[1:16])))
x %>%  ggplot2::ggplot(ggplot2::aes(x=Col, y=Row, fill=log10(Count))) + 
  ggplot2::geom_raster() +
  ggplot2::coord_equal() +
  ggplot2::facet_grid(amplicon~Plate_ID) +
  ggplot2::scale_fill_viridis_c(option = 'plasma')+ggplot2::ggtitle('S2') 




#IndexSet=as.character(params$dlong$Plate_ID)
#IndexSet[IndexSet %in% c('Plate1', 'Plate2', 'Plate3', 'Plate4')]='Winter'
#IndexSet[IndexSet %in% c('Plate5', 'Plate6', 'Plate7', 'Plate8')]='Spring'
#IndexSet[IndexSet %in% c('Plate9', 'Plate10', 'Plate11', 'Plate12') ]='Summer'
#IndexSet[IndexSet %in% c('Plate13', 'Plate14', 'Plate15', 'Plate16')]='Fall'
#
#
#IndexSetW=as.character(params$dwide$Plate_ID)
#IndexSetW[IndexSetW %in% c('Plate1', 'Plate2', 'Plate3', 'Plate4')]='Winter'
#IndexSetW[IndexSetW %in% c('Plate5', 'Plate6', 'Plate7', 'Plate8')]='Spring'
#IndexSetW[IndexSetW %in% c('Plate9', 'Plate10', 'Plate11', 'Plate12') ]='Summer'
#IndexSetW[IndexSetW %in% c('Plate13', 'Plate14', 'Plate15', 'Plate16')]='Fall'



