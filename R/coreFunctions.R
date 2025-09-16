#seqinr
#XML
#plater
#dplyr
#tidyr
#data.table
#Biobase
#S4Vectors
#ggplot 
#optional
#seqtrie


#' @title Reverse complement a sequence
#' @description `revcomp` returns the reverse complement of a sequence
#' @param x A character vector
#' @return reverse complement of x
revcomp2=function (x) {    toupper(seqinr::c2s(rev(seqinr::comp(seqinr::s2c(x)))))}



#generate expected indices, hard code semi-combinatorial indexing here ---------------------------------------------------

#' @title Generate expected indices
#' @description Generate expected i7 and i5 index sequences given indices in plater format
#' @param cfg A list with i5_plate_key_file and i7_plate_key_file locations
#' @param diri (optional) directory containing Illumina RunParameters.xml file
#' @return a data.table with Plate_ID, Sample_Well, Sample_ID, index, and index2  ... where index is i7 and index2 is i5
#' @export
generateExpectedIndices=function(cfg, diri=NULL) {
    #how to check run chemistry from xml file  ------------------------------------------
    chemistry=NULL
    if(!is.null(diri)){
    xmlinfo=XML::xmlToList(XML::xmlParse(paste0(diri, '/RunParameters.xml')))
    chemistry=xmlinfo$Chemistry
    }
        #MiniSeq High / MiniSeq Rapid High / NextSeq Mid / NextSeq High
        #don't reverse comp i5 for miniseq rapid or miseq 
    i5RC.toggle=TRUE
        #addition to handle lack of chemistry tag for nextseq2000 11/19/21
    if(!is.null(chemistry)) {     if(chemistry=="MiniSeq Rapid High" | chemistry=="MiSeq") {i5RC.toggle=F}  }
#----------------------------------------------------------------------------------------

    i7s=plater::read_plates(cfg$i7_plate_key_file, well_ids_column="Sample_Well")
    i7s=tidyr::gather(i7s, Plate_ID, index, 3:ncol(i7s))
       
    i5s=plater::read_plates(cfg$i5_plate_key_file, well_ids_column="Sample_Well")
    i5s=tidyr::gather(i5s, Plate_ID, index2, 3:ncol(i5s))
      
        #r=i7
        #f=i5
    semi.key.i7=c(c(1,2,3,4),c(2,3,4,1),c(3,4,1,2), c(4,1,2,3))
    semi.key.i5=rep(seq(1,4),4)

    i7ss=split(i7s, i7s$Plate_ID)
    i5ss=split(i5s, i5s$Plate_ID)


    i7ss=data.table::rbindlist(i7ss[semi.key.i7])
    i5ss=data.table::rbindlist(i5ss[semi.key.i5])

    i7ss$Plate_ID=rep(seq(1,16),each=384)
    i5ss$Plate_ID=rep(seq(1,16),each=384)

    i7ss$Plate_ID=paste0('Plate', i7ss$Plate_ID)
    i7ss$Sample_ID=paste0(i7ss$Plate_ID,'-', i7ss$Sample_Well)
    i7ss$index=as.vector(sapply(i7ss$index, revcomp2))


    i5ss$Plate_ID=paste0('Plate', i5ss$Plate_ID)
    i5ss$Sample_ID=paste0(i5ss$Plate_ID,'-', i5ss$Sample_Well)
    if(i5RC.toggle) { i5ss$index2=as.vector(sapply(i5ss$index2, revcomp2)) }


    i5ss= i5ss %>% dplyr::select(Sample_ID, index2)
    #i7s$bc_set='N1_S2_RPP30'
    index.key=dplyr::right_join(i7ss,i5ss,by='Sample_ID') %>% dplyr::select(-Plate) %>% dplyr::select(Plate_ID,Sample_Well,Sample_ID,index,index2)
    #----------------------------------------------------------------------------------------------------------------
}


#' @title initialize amplicon count table 
#' @param index.key index.key from generateExpectedIndices
#' @param amplicons, a list containing names and expected sequences of amplicons
#' @return table per amplicon to hold amplicon counts per expected index
#' @export
initAmpliconCountTables=function(index.key, amplicons) {
    #Munginging sample sheet-------------------------------------------------------------------
    ss=index.key
    ss$mergedIndex=paste0(ss$index, ss$index2)

    # this code would be obviated if indices designate wells, for most analyses here there are different indices for s2/s2spike and rpp30
    # subset of indices for S2/S2 spike
#    if(sum(grepl('-1$', ss$Sample_ID))==0){
#        ssS=ss
#        ssR=ss
#    } else {
#        ssS=ss[grep('-1$', ss$Sample_ID),]
#        #subset of indices for RPP30
#        ssR=ss[grep('-2$', ss$Sample_ID),]
#    }

    #initalize output count tables ------------------------------------------------------------
    count.tables=list()
    for(a in names(amplicons)){
  #      if(grepl('^S',a)){    count.tables[[a]]=ssS   } 
  #      if(grepl('^R',a)){    count.tables[[a]]=ssR   } 
            count.tables[[a]]=ss
            count.tables[[a]]$Count=0
            count.tables[[a]]$amplicon=a
    }
    return(count.tables)
}

#' @title make all hamming 1 distant sequences for a given sequence 
#' @param x a sequences
#' @return the set of all unique hamming distance 1 sequences from a given sequence
#' @export
make_hamming1_sequences=function(x) {
    #eseq=seqinr::s2c(x)
    eseq= base::strsplit(x, '')[[1]]
    #eseqs=c(seqinr::c2s(eseq))
    eseqs=x
    for(i in 1:length(eseq)){
        eseq2=eseq
        eseq2[i]='A'
        eseqs=c(eseqs,paste0(eseq2, collapse=''))
        #eseqs=c(eseqs,seqinr::c2s(eseq2))
        eseq2[i]='C'
        eseqs=c(eseqs,paste0(eseq2, collapse=''))
        #eseqs=c(eseqs,seqinr::c2s(eseq2))
        eseq2[i]='T'
        eseqs=c(eseqs,paste0(eseq2, collapse=''))
        #eseqs=c(eseqs,seqinr::c2s(eseq2))
        eseq2[i]='G'
        eseqs=c(eseqs,paste0(eseq2, collapse=''))
        #eseqs=c(eseqs,seqinr::c2s(eseq2))
        eseq2[i]='N'
        eseqs=c(eseqs,paste0(eseq2, collapse=''))
       # eseqs=c(eseqs,seqinr::c2s(eseq2))

    }
    eseqs=unique(eseqs)
    return(eseqs)
}

#' @title error correct index reads and assign amplicon counts 
#' @param rid indices of sequences matching an expected amplicon 
#' @param count.table table per amplicon to hold amplicon counts per expected index
#' @param ind1 index1 sequences
#' @param ind2 index2 sequences
#' @return count.table the set of all unique hamming distance 1 sequences from a given sequence
#' @export
errorCorrectIdxAndCountAmplicons=function(rid, count.table, ind1,ind2){
    # get set of unique expected index1 and index2 sequences
    index1=unique(count.table$index)
    index2=unique(count.table$index2)

    i1h=lapply(index1, make_hamming1_sequences)
    names(i1h)=index1
    ih1=Biobase::reverseSplit(i1h)
    ih1.elements=names(ih1)
    ih1.indices=as.vector(unlist(ih1))
    i1m=ih1.indices[S4Vectors::match(ind1[rid],ih1.elements)]

    i2h=lapply(index2, make_hamming1_sequences)
    names(i2h)=index2
    ih2=Biobase::reverseSplit(i2h)
    ih2.elements=names(ih2)
    ih2.indices=as.vector(unlist(ih2))
    i2m=ih2.indices[S4Vectors::match(ind2[rid],ih2.elements)]
    
    idm=paste0(i1m,i2m)
    tS2=table(match(idm, count.table$mergedIndex))
    tbix=match(as.numeric(names(tS2)), 1:nrow(count.table))
    count.table$Count[tbix]=as.vector(tS2)+count.table$Count[tbix]
    return(count.table)
}

#unused code that could enable Levenshtein distance matching for amplicon sequences 
seqtrie_match=function(#character vector of observed sequences
                       seq, 
                       #character vector of expected sequences
                       uexp, 
                       #max levenshtein distance (lv) to search for a match 
                       maxDist=0, 
                       #max parallel threads for match operation
                       lv.nthreads=1,
                       #name to prepend to match table 
                       name='') {
                
                useq=unique(seq)
                toggle=F

                if(length(useq)==1 ) {
                   if( useq=="") {
                       toggle=T
                   }
                }
                if(toggle) {
                    match.df=data.frame(obs=seq,exp=as.character(NA),expInd=as.integer(NA),dist=as.integer(NA))
                   
                } else{
                    seqtrie.res   = seqtrie::dist_search(useq,uexp, max_distance=maxDist, nthreads=lv.nthreads) 
                    seqtrie.res.s = split(seqtrie.res, seqtrie.res$query)

                    seqtrie.res.closestMatch=sapply(seqtrie.res.s, function(x) x$target[which.min(x$distance)])
                    seqtrie.res.closestDist=sapply(seqtrie.res.s, function(x) min(x$distance))
                    seqtrie.res.df=data.frame(obs=names(seqtrie.res.s),
                                              exp=seqtrie.res.closestMatch,
                                              expInd=match(seqtrie.res.closestMatch, uexp),
                                              dist=seqtrie.res.closestDist)
                    seqtrie.res.df=seqtrie.res.df[match(useq, seqtrie.res.df$obs),]
      
                    match.df = seqtrie.res.df[match(seq,seqtrie.res.df$obs),]
                    nm=which(is.na(match.df$obs))
                    match.df$obs[nm]=seq[nm]
                }
                match.df = match.df %>% dplyr::rename_with(~paste0(name, .))
                rownames(match.df)=NULL
                return(match.df)
            }

#' @title Count expected amplicons
#' @description Count expected amplicons (per sample) and also return total reads per sample.
#' @param in.con gzfile/connection to R1 FASTQ (trimmed if needed)
#' @param index.key data.table with Plate_ID, Sample_Well, Sample_ID, index (i7), index2 (i5)
#' @param amplicons named list of expected R1 amplicon sequences
#' @param line.buffer number of lines to read per chunk (not reads; 4 lines/read)
#' @param max.lines optional hard stop on total lines read (for debugging)
#' @param nthreads threads for stringfish (if available)
#' @return list(count.tables=..., amp.match.summary.table=..., total.count.table=...)
#' @export
countAmplicons <- function(in.con, index.key, amplicons,
                           line.buffer = 5e6, max.lines = NULL, nthreads = 1) {

  lines_read <- 0L

  # your original per-amplicon tables
  count.tables <- initAmpliconCountTables(index.key, amplicons)

  # running amplicon summary (add a "no_align" bin)
  amp.match.summary.table <- setNames(integer(length(amplicons) + 1L),
                                      c(names(amplicons), "no_align"))

  # build the Hamming-1 neighborhood for amplicons (your original approach)
  amph1=lapply(amplicons, make_hamming1_sequences)
  amph1=Biobase::reverseSplit(amph1)
  amph1.elements=names(amph1)
  amph1.indices=as.vector(unlist(amph1))

  # totals table that uses the SAME sample order as your per-amplicon tables 
  # we copy the first amplicon's table structure to guarantee row order alignment
  first_amp <- names(count.tables)[1]
  if (length(first_amp) == 0L) stop("initAmpliconCountTables returned no tables.")
  total_tbl_work <- data.table::copy(data.table::as.data.table(count.tables[[first_amp]]))
  # zero Count without using :=
  if (!"Count" %in% names(total_tbl_work)) {
    total_tbl_work$Count <- 0L
  } else {
    data.table::set(total_tbl_work, j = "Count", value = as.integer(0))
  }
    
  if ("amplicon" %in% names(total_tbl_work)) data.table::set(total_tbl_work, j = "amplicon", value = NULL)

  repeat {
    chunk <- readLines(in.con, n = line.buffer)
    n_all <- length(chunk)
    if (n_all == 0L) break
    if (n_all %% 4L != 0L)
      stop("FASTQ chunk not multiple of 4 lines (truncated input?).")

    lchunk <- n_all %/% 4L
    lines_read <- lines_read + n_all

    # FASTQ parsing
    nlines <- seq_len(n_all)
    mod4   <- nlines %% 4L
    header <- chunk[mod4 == 1L]
    rd1    <- chunk[mod4 == 2L]

    # indices from header (everything after last colon)
    if (nthreads > 1 && requireNamespace("stringfish", quietly = TRUE)) {
      tmp <- stringfish::sf_gsub(header, ".*:", "", nthreads = nthreads)
    } else {
      tmp <- gsub(".*:", "", header, perl = TRUE)
    }
    ind1 <- substring(tmp, 1L, 10L)
    ind2 <- substring(tmp, 12L, 21L)

    # ----- amplicon classification (R1) -----
    amp.match=amph1.indices[S4Vectors::match(rd1, amph1.elements)]
    no_align=sum(is.na(amp.match) )

    #summarize amplicon matches
    amp.match.summary=table(amp.match)
    amp.match.summary=amp.match.summary[match(names(amplicons),names(amp.match.summary))]
    amp.match.summary=c(amp.match.summary, no_align)
    names(amp.match.summary) <- c(names(amp.match.summary[-length(amp.match.summary)]),"no_align")
    amp.match.summary.table=amp.match.summary.table+amp.match.summary

    # ----- your original per-amplicon counting path -----
    per.amplicon.row.index=lapply(names(amplicons), function(x) which(amp.match==x))
    names(per.amplicon.row.index)=names(amplicons)

    for (a in names(count.tables)) {
      if (!length(per.amplicon.row.index[[a]])) next
      count.tables[[a]] <- errorCorrectIdxAndCountAmplicons(
        per.amplicon.row.index[[a]],
        count.tables[[a]],
        ind1, ind2
      )
    }

    # totals bump = run the SAME index EC/counting on *all reads* 
    # this assigns every read to a sample (if indices EC) regardless of amplicon match
    total_tbl_work <- errorCorrectIdxAndCountAmplicons(
      seq_along(ind1),  # all reads in this chunk
      total_tbl_work,
      ind1, ind2
    )

    # optional early stop (max.lines is in *lines*, not reads)
    if (!is.null(max.lines) && lines_read >= max.lines) {
      break
    }
  }

  close(in.con)

  # collapse totals to one row per Sample_ID (in case Sample_IDs repeat)
  tt <- as.data.frame(total_tbl_work)
  total.count.table <- stats::aggregate(Count ~ Plate_ID + Sample_Well, data = tt, FUN = sum)
  names(total.count.table)[3] <- "Total_Count"
  
  list(
    count.tables            = count.tables,
    amp.match.summary.table = amp.match.summary.table,
    total.count.table       = total.count.table
  )
}

