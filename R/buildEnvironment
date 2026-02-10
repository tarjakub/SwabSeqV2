usethis::use_pipe(export =TRUE)
####@importFrom magrittr "%>%"
#'
#' Build SwabSeq RVP environment variable
#'
#' This function creates an env() containing information about directory structures and key run parameters
#'
#' @export
buildEnvironment=function( remote.dir, localmirror.dir, bcl.dir,
                          threads=8,
                          lbuffer=30e6,
                          readerBlockSize=1e8,
                          fastqcr=F,
                          i7_plate_key_file=NULL, i5_plate_key_file=NULL) {
    # for 16GB of ram consider lbuffer=5e6 and readerBlockSize=1e6
    cfg = new.env(parent=emptyenv())
    path_to_bs=tryCatch(system('command -v bs', intern=T), error=function(e) {return(NULL)})
    path_to_bcl2fastq=tryCatch(system('command -v bcl2fastq', intern=T), error=function(e) {return(NULL)})
    bs_config_present=file.exists(paste0(Sys.getenv("HOME"), '/.basespace/default.cfg'))
    if(is.null(path_to_bs)) {print('bs CLI tool not found in PATH')}
    if(is.null(path_to_bcl2fastq)) {print('bcl2fastq not found in PATH')}
    if(!(bs_config_present)) {print('basespace cfg file not found at ~/.basespace/default.cfg')}
    
    #i7_plate_key_file= paste0(basedir.dir, 'reference/PCR2_r.csv')
    #i5_plate_key_file= paste0(basedir.dir, 'reference/PCR2_f.csv')
    
    if(is.null(i7_plate_key_file)){
    i7_plate_key_file = system.file("keys", "PCR2_r.csv", package="SwabSeqV2")
    }
    if(is.null(i5_plate_key_file)){
    i5_plate_key_file = system.file("keys", "PCR2_f.csv", package="SwabSeqV2")
    }
    
    cfg$i7_plate_key_file=i7_plate_key_file
    cfg$i5_plate_key_file=i5_plate_key_file

  # Directory structure on samba share ------------------------------------------------------
    # seq.dir is the path to location of yaml, summary stats for each run, and reports
    seq.dir=paste0(remote.dir, 'seq/')
    seq.run.dir=paste0(remote.dir, 'seq/runs/')

    #path to swabseq remote file share sample tracking folder
    sampleTracking.dir=paste0(remote.dir, 'swabseqsampletracking/')

    #location of tracking water tubes
    waterTubesKeyDir=paste0(seq.dir, 'water_tubes/')
    yaml.cfg.file=paste0(seq.dir, 'config.yaml')

  # Directory structures for local mirroring bits of samba share -------------------------------
    localMirrorSeq.dir=paste0(localmirror.dir, 'seq/')
    localMirrorSeq.run.dir=paste0(localmirror.dir, 'seq/runs/')
    localTracking.dir=paste0(localmirror.dir, '/swabseqsampletracking/')

    cfg$remote.dir=remote.dir
    cfg$localmirror.dir=localmirror.dir
    cfg$seq.dir=seq.dir
    cfg$seq.run.dir=seq.run.dir
    cfg$bcl.dir=bcl.dir
    cfg$tracking.dir=sampleTracking.dir
    cfg$localtracking.dir=localTracking.dir
    cfg$localMirrorSeq.dir=localMirrorSeq.dir
    cfg$localMirrorSeq.run.dir=localMirrorSeq.run.dir
    cfg$yaml.cfg.file=yaml.cfg.file
    cfg$waterTubesKeyDir=waterTubesKeyDir

    cfg$path_to_bs=path_to_bs
    cfg$path_to_bcl2fastq=path_to_bcl2fastq
    cfg$bs_config_present=bs_config_present

    cfg$coreVars=setAnalysisVariables(threads=threads,lbuffer=lbuffer,readerBlockSize=readerBlockSize,fastqcr=fastqcr)

    #-------------------------------------------------------------------------------------------
    return(cfg)
  }
  #' Set Analysis Variables
  #' 
  #' This function is called by buildEnvironment() and sets key analysis parameters, exposed here for user overwriting
  #'
  #' @export
  setAnalysisVariables=function(versi=3, lbuffer=30000000, readerBlockSize=1e8, threads=8, fastqcr=FALSE, Ratio=NULL){
      # versi 1 for extracted RNA samples
      # versi 2 for extraction free samples

      flagLowPositive=FALSE
      AmpliconTotal=500
      Rpp=10

      # define respiratory viral panel
      targets=c("CoV2", "FluA_H1", "FluA_H3", "H1N1", "H3N2", "H5N1", "FluB", "RSV")

      if (versi == 1) {
        default_ratio=c(CoV2=0.05, FluA_H1=0.05, FluA_H3=0.05, H1N1=0.05, H3N2=0.05, FluB=0.05, H5N1=0.05, FluB=0.05, RSV=0.05)
      }
      if (veresi == 2) {
        default_ratio=c(CoV2=0.1, FluA_H1=0.1, FluA_H3=0.1, H1N1=0.1, H3N2=0.1, FluB=0.1, H5N1=0.1, FluB=0.1, RSV=0.1)
      }
      
      if(is.null(Ratio)){
        Ratio=default_ratio
      } else if (length(Ratio) == 1 && is.null(names(Ratio))) {
          Ratio=stats::setNames(rep(as.numeric(Ratio), length(targets)), targets)
        } else {
          nm=names(Ratio)
          Ratio_num=as.numeric(Ratio)
          out=stats::setNames(rep(default_ratio, length(targets)), targets)
          if (!is.null(nm)) out[nm] <- Ratio_num
          Ratio <- out
        }
      }

      # Expected amplicon sequences
      amplicons <- list(
        CoV2 = "ATTGGCTACTACCGA",
        CoV2_spike = "TAACCGATGATGGCT",
        FluA_H1 = "TCTATCGTTCCATCA",
        FluA_H3 = "TTCTATCATCCCGTC",
        FluA_spike = "AGATAGTAGGGCAGT",
        H1N1 = "ATTTGAAAGGTTTGA",
        H1N1_spike = "TAAACTTTCCAAACT",
        H3N2 = "ACAACGCGGAGCTTC",
        H3N2_spike = "TGTTGCGCCTCGAAG",
        H5N1 = "GCAAGTTCCCTAGCA",
        H5N1_spike = "CGTTCAAGGGATCGT",
        FluB = "GAGGTGGGTCCGGGA",
        FluB_spike = "CTCCACCCAGGCCCT",
        RSV = "ACTACTTCCCACCAA",
        RSV_spike = "TGATGAAGAGTAGTT",
        RPP30 = "GAAGGCTCTGCGCGG"
      )
      # Spike Control Map
      spike_map <- c(
        CoV2 = "CoV2_spike",
        FluA_H1 = "FluA_spike",
        FluA_H3 = "FluA_spike",
        H1N1 = "H1N1_spike",
        H3N2 = "H3N2_spike",
        H5N1 = "H5N1_spike",
        FluB  = "FluB_spike",
        RSV = "RSV_spike"
      )
      
      if(fastqcr) {
        if(!file.exists(path.expand("~/bin/FastQC"))) { fastqcr::fastqc_install() }
      }

      return(list(
        versi=versi,
        lbuffer=lbuffer,
        readerBlockSize=readerBlockSize,
        threads=threads,
        fastqcr=fastqcr,
        flagLowPositive=flagLowPositive,
        AmpliconTotal=AmpliconTotal,
        Rpp=Rpp,
        Ratio=Ratio,
        amplicons=amplicons,
        spike_map=spike_map,
        empty_well_set=c('', 'TBET', 'No Tube', 'NO TUBE', 'Empty', 'EMPTY', ' ', 'NA')
        ))
      }   
