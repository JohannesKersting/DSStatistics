#import librarys
library(data.table)
library(ggplot2)
library(ggbeeswarm)
library(viridis)
library(corrplot)
library(tidyr)
library(GenomicRanges)
library(GenomicFeatures)

############ groundtruth ##############
#######################################

read_info <- function(path){
  data <- fread(path)
  data <- separate(data,col=transcriptid,into=c("gene_id","type"),sep="_")
  short_data <- data[,.(diff_spliced = (.SD[,any(DEstatus.V1)]|.SD[,any(DEstatus.c2)]),type=.SD$type[1]),by=gene_id]
  return(short_data)
}

read_annotation <- function(path){
  data <- fread(path)
  es <- get_es(data[event_annotation=="es"])
  mee <- get_mee(data[event_annotation=="mee"])
  mes <- get_mes(data[event_annotation=="mes"])
  afe <- get_afe(data[event_annotation=="afe"])
  ale <- get_ale(data[event_annotation=="ale"])
  ir <- get_ir(data[event_annotation=="ir"])
  a3 <- get_a3(data[event_annotation=="a3"])
  a5 <- get_a5(data[event_annotation=="a5"])
  return(list(mee=mee,mes=mes,es=es,afe=afe,ale=ale,ir=ir,a3=a3,a5=a5))
}

get_a5 <- function(data){
  gene_ids <- data[,gsub("_a5","",variant)]
  
  genomic_ranges <- IRanges(start=as.integer(data$genomic_start),
                            end=as.integer(data$genomic_end)
  )
  
  tx_ranges <- IRanges(start=as.integer(data$transcriptomic_start),
                       end=as.integer(data$transcriptomic_end)
  )
  names(genomic_ranges) <- gene_ids
  names(tx_ranges) <- gene_ids
  
  return(list(genomic_ranges=genomic_ranges,tx_ranges=tx_ranges))
}

get_a3 <- function(data){
  gene_ids <- data[,gsub("_a3","",variant)]
  
  genomic_ranges <- IRanges(start=as.integer(data$genomic_start),
                            end=as.integer(data$genomic_end)
  )
  
  tx_ranges <- IRanges(start=as.integer(data$transcriptomic_start),
                       end=as.integer(data$transcriptomic_end)
  )
  names(genomic_ranges) <- gene_ids
  names(tx_ranges) <- gene_ids
  
  return(list(genomic_ranges=genomic_ranges,tx_ranges=tx_ranges))
}


get_ir <- function(data){
  gene_ids <- data[,gsub("_ir","",variant)]
  
  genomic_ranges <- IRanges(start=as.integer(data$genomic_start),
                            end=as.integer(data$genomic_end)
  )
  
  tx_ranges <- IRanges(start=as.integer(data$transcriptomic_start),
                       end=as.integer(data$transcriptomic_end)
  )
  names(genomic_ranges) <- gene_ids
  names(tx_ranges) <- gene_ids
  
  return(list(genomic_ranges=genomic_ranges,tx_ranges=tx_ranges))
}

get_ale <- function(data){
  
  genomic_ranges <- IRanges(start=as.integer(data$genomic_start),
                            end=as.integer(data$genomic_end),
                            type = data$type
  )
  
  tx_ranges <- IRanges(start=as.integer(data$transcriptomic_start),
                       end=as.integer(data$transcriptomic_end),
                       type = data$type
  )
  
  names(genomic_ranges) <- data$variant
  names(tx_ranges) <- data$variant
  return(list(genomic_ranges=genomic_ranges,tx_ranges=tx_ranges))
}

get_afe <- function(data){
  
  genomic_ranges <- IRanges(start=as.integer(data$genomic_start),
                            end=as.integer(data$genomic_end),
                            type = data$type
  )
  
  tx_ranges <- IRanges(start=as.integer(data$transcriptomic_start),
                       end=as.integer(data$transcriptomic_end),
                       type = data$type
  )
  
  names(genomic_ranges) <- data$variant
  names(tx_ranges) <- data$variant
  return(list(genomic_ranges=genomic_ranges,tx_ranges=tx_ranges))
}

get_es <- function(data){
  gene_ids <- data[,gsub("_es","",variant)]
  
  genomic_ranges <- IRanges(start=as.integer(data$genomic_start),
                            end=as.integer(data$genomic_end)
  )
  
  tx_ranges <- IRanges(start=as.integer(data$transcriptomic_start),
                       end=as.integer(data$transcriptomic_end)
  )
  names(genomic_ranges) <- gene_ids
  names(tx_ranges) <- gene_ids
  
  return(list(genomic_ranges=genomic_ranges,tx_ranges=tx_ranges))
}

get_mee <- function(data){
  
  genomic_ranges <- IRanges(start=as.integer(data$genomic_start),
                            end=as.integer(data$genomic_end),
                            type = data$type
  )
  
  tx_ranges <- IRanges(start=as.integer(data$transcriptomic_start),
                       end=as.integer(data$transcriptomic_end),
                       type = data$type
  )
  
  names(genomic_ranges) <- data$variant
  names(tx_ranges) <- data$variant
  return(list(genomic_ranges=genomic_ranges,tx_ranges=tx_ranges))
}

get_mes <- function(data){
  
  gene_ids <- data[,gsub("_mes","",variant)]
  
  genomic_ranges <- IRangesList(lapply(1:nrow(data), function(i){
    genomic_start <- data[i]$genomic_start
    genomic_end <- data[i]$genomic_end
    genomic_starts <- strsplit(genomic_start,",")[[1]]
    genomic_ends <- strsplit(genomic_end,",")[[1]]
    return(IRanges(start=as.integer(genomic_starts),end=as.integer(genomic_ends)))
  }))
  
  tx_ranges <- IRangesList(lapply(1:nrow(data), function(i){
    tx_start <- data[i]$transcriptomic_start
    tx_end <- data[i]$transcriptomic_end
    tx_starts <- strsplit(tx_start,",")[[1]]
    tx_ends <- strsplit(tx_end,",")[[1]]
    return(IRanges(start=as.integer(tx_starts),end=as.integer(tx_ends)))
  }))

  names(genomic_ranges) <- gene_ids
  names(tx_ranges) <- gene_ids
  
  mes = list(genomic_ranges = genomic_ranges, tx_ranges = tx_ranges)
  return(mes)
}

load_introns <- function(annotation, txdb){
  transcript_introns <- intronsByTranscript(txdb,use.names=T)
  
  #es
  es_ranges <- annotation$es$genomic_ranges
  es_introns <- transcript_introns[paste(names(es_ranges),"es",sep="_")]
  names(es_introns) <- gsub("_es","",names(es_introns))
  
  annotation$es$introns <- get_introns(es_ranges,es_introns)

  
  #mes
  mes_ranges <- annotation$mes$genomic_ranges
  mes_merged_ranges <- IRanges(start= min(start(mes_ranges)),
                               end = max(end(mes_ranges)))
  names(mes_merged_ranges) <- names(mes_ranges)
  mes_introns <- transcript_introns[paste(names(mes_ranges),"mes",sep="_")]
  names(mes_introns) <- gsub("_mes","",names(mes_introns))
  annotation$mes$introns <- get_introns(mes_merged_ranges,mes_introns)
  
  #mee
  mee_ranges <- annotation$mee$genomic_ranges
  mee_introns <-  transcript_introns[names(mee_ranges)]
  names(mee_ranges) <- gsub("_mee","_place",names(mee_ranges))
  names(mee_ranges) <- gsub("_template","_mee",names(mee_ranges))
  names(mee_ranges) <- gsub("_place","_template",names(mee_ranges))
  
  mee_introns <- get_introns(mee_ranges,mee_introns)
  names(mee_introns) <- gsub("_mee","",names(mee_introns))
  names(mee_introns) <- gsub("_template","",names(mee_introns))
  mee_intron_list <- split(mee_introns,names(mee_introns))
  annotation$mee$introns <- IRanges(start=min(start(mee_intron_list)),
                                    end=max(end(mee_intron_list)),
                                    names=names(mee_intron_list))
  
  #a3
  a3_ranges <- annotation$a3$genomic_ranges
  a3_introns <- transcript_introns[paste(names(a3_ranges),"a3",sep="_")]
  names(a3_introns) <- gsub("_a3","",names(a3_introns))
  annotation$a3$introns <- get_introns(a3_ranges,a3_introns)
  
  #a5
  a5_ranges <- annotation$a5$genomic_ranges
  a5_introns <- transcript_introns[paste(names(a5_ranges),"a5",sep="_")]
  names(a5_introns) <- gsub("_a5","",names(a5_introns))
  annotation$a5$introns <- get_introns(a5_ranges,a5_introns)
  
  #afe
  afe_ids <- unique(gsub("_.*","",names(annotation$afe$genomic_ranges)))
  transcripts <- exonsBy(txdb,by="tx",use.names=T)
  afe_transcripts <- transcripts[paste(afe_ids,"afe",sep="_")]
  template_transcripts <- transcripts[paste(afe_ids,"template",sep="_")]
  afe_points <- sapply(afe_transcripts,function(x){
    if(as.character(strand(x))[1]=="+"){
      return(start(x)[1:2])
    }
    else{
      return(end(x)[c(2,1)])
    }
    
  })
  
  template_points <- sapply(template_transcripts,function(x){
    if(as.character(strand(x))[1]=="+"){
      return(start(x)[1:2])
    }
    else{
      return(end(x)[c(2,1)])
    }
    
  })
  points <- unname(rbind(afe_points,template_points))
  starts <- apply(points,2,min)
  ends <- apply(points,2,max)
  annotation$afe$introns <- IRanges(start=starts,
                     end = ends,
                     names = afe_ids)
  
  #ale
  ale_ids <- unique(gsub("_.*","",names(annotation$ale$genomic_ranges)))
  transcripts <- exonsBy(txdb,by="tx",use.names=T)
  ale_transcripts <- transcripts[paste(ale_ids,"ale",sep="_")]
  template_transcripts <- transcripts[paste(ale_ids,"template",sep="_")]
  ale_points <- sapply(ale_transcripts,function(x){
    if(as.character(strand(x))[1]=="+"){
      return(end(x)[c(length(x),length(x)-1)])
    }
    else{
      return(start(x)[c(length(x),length(x)-1)])
    }
    
  })
  
  template_points <- sapply(template_transcripts,function(x){
    if(as.character(strand(x))[1]=="+"){
      return(end(x)[c(length(x),length(x)-1)])
    }
    else{
      return(start(x)[c(length(x),length(x)-1)])
    }
    
  })
  points <- unname(rbind(ale_points,template_points))
  starts <- apply(points,2,min)
  ends <- apply(points,2,max)
  annotation$ale$introns <- IRanges(start=starts,
                                    end = ends,
                                    names = ale_ids)
  
  return(annotation)
}

get_introns <- function(exon_regions,transcript_introns){
  transcript_introns <- ranges(transcript_introns)
  exon_regions_list <- split(exon_regions,names(exon_regions))
  intron_regions <- unlist(subsetByOverlaps(transcript_introns,exon_regions_list))
  return(intron_regions)
}

pre_filter_type <- function(data,true_data){
  
  #type list 
  type_list <- as.list(true_data$type)
  names(type_list) <- true_data$gene_id
  
  #real type
  data[,real_type := unlist(type_list[gene_id],use.names=F)]
  
  #include?
  includes <- unique(data$type)
  print(paste0("The known event types are: ",paste(includes,collapse = ", ")))
  nrow_data <- nrow(data)
  data <- data[real_type %in% includes]
  print(paste0(nrow_data-nrow(data)," detected events were dropped"))
  
  return(data)
}




########### tools ##################
####################################

read_ds <- function(path,method){
  if(method=="cash"){
    return(read_cash(path))
  }
  if(method=="eventpointer"){
    return(read_eventpointer(path))
  }
  if(method=="aspli"){
    return(read_aspli(path))
  }
  if(method=="dexseq"){
    return(read_dexseq(path))
  }
  if(method=="majiq"){
    return(read_majiq(path))
  }
  if(method=="spladder"){
    return(read_spladder(path))
  }
  if(method == "jsplice"){
    return(read_jsplice(path))
  }
  stop(paste0("Method '",method,"' is not supported!"))
}

count_ds <- function(data,method,...){
  if(method=="cash"){
    return(count_cash(data,...))
  }
  if(method=="eventpointer"){
    return(count_eventpointer(data,...))
  }
  if(method=="majiq"){
    return(count_majiq(data,...))
  }
  if(method=="aspli"){
    return(count_aspli(data,...))
  }
  if(method=="spladder"){
    return(count_spladder(data,...))
  }
  if(method == "jsplice"){
    return(count_jsplice(data...))
  }
  stop(paste0("Method '",method,"' is not supported!"))
}

############ cash ######################
read_cash <- function(path){
  
  #read data
  file <- list.files(path,pattern = ".alldiff.txt$")
  if(length(file)!=1){
      stop("No input file found!")
  }
  path <-file.path(path,file)
  print(paste("Reading cash output:",path))
  data <- fread(path)
  
  #select relevant columns and split location
  data <- data[,.(gene_id=AccID,p_value=data$"P-Value",type=SplicingType,location=Location)]
  data[,location:=gsub(".*:","",data$location)]
  sep <- separate(data[,.(location)],location,c("start","end"),sep="-")
  data[,location:=NULL]
  data[,start:=as.integer(sep$start)]
  data[,end:=as.integer(sep$end)]
  
  #map for type parsing
  map <- list("AltEnd"="ale","A3SS"="a3","MXE"="mee","IR"="ir","A5SS"="a5","Cassette"="es","AltStart"="afe","Cassette_multi"="mes")
  data[,type:=as.factor(unlist(map[type],use.names = F))]
  
  return(data)
}

#count 
count_cash <- function(data,true_data,annotation,txdb){

  
  #lookup lists
  diff_spliced_list <- as.list(true_data$diff_spliced)
  names(diff_spliced_list) <- true_data$gene_id
  
  type_list <- as.list(true_data$type)
  names(type_list) <- true_data$gene_id
  
  #add strand
  data[,strand:=as.character(strand(genes(txdb)[data$gene_id]))]
  
  #gene_level
  counts <- data.table()
  gene_level <-unlist(diff_spliced_list[data$gene_id],use.names = F) 
  data[,gene_level:=gene_level]
  unique_data <- unique(data[,.(gene_id,gene_level)])
  counts$level <- "gene"
  counts$tp <- sum(unique_data$gene_level)
  counts$fp <- sum(unique_data$gene_level==F)
  counts$fn  <- sum(true_data$diff_spliced)-counts$tp
  
  #type_level
  counts_type <- data.table() 
  data[,type_level := gene_level&(unlist(type_list[gene_id],use.names = F)==type)]
  data[,real_type := unlist(type_list[gene_id],use.names=F)]
  unique_data <- unique(data[,.(gene_id,gene_level,type_level,type)])
  counts_type$level <- "type"
  counts_type$tp <- sum(unique_data$type_level)
  counts_type$fp <- sum(unique_data$type_level==F)
  counts_type$fn  <- sum(true_data$diff_spliced)-counts_type$tp
  counts <- rbind(counts,counts_type)
  
  #location_level
  counts_location <- data.table()
  data[type_level==F,location_level:=F]
  
  #es
  es_data <- data[type_level==T&type=="es"]
  es_ranges <- annotation$es$genomic_ranges[es_data$gene_id]
  es_data[,location_level:=((start(es_ranges)==start)&(end(es_ranges)==end))]
  data[type_level==T&type=="es",location_level:=es_data$location_level]
  es_data <- NULL
  es_ranges <- NULL
  
  #ir
  ir_data <- data[type_level==T&type=="ir"]
  ir_ranges <- annotation$ir$genomic_ranges[ir_data$gene_id]
  ir_data[,location_level:=(start(ir_ranges)==(start+1))&(end(ir_ranges)==(end-1))]
  data[type_level==T&type=="ir",location_level:=ir_data$location_level]
  ir_data <- NULL
  ir_ranges <- NULL
  
  #mes 
  mes_data <- data[type_level==T&type=="mes"]
  mes_ranges <- annotation$mes$genomic_ranges[mes_data$gene_id]
  mes_data[,location_level:=(min(start(mes_ranges))==start)&(max(end(mes_ranges))==end)]
  data[type_level==T&type=="mes",location_level:=mes_data$location_level]
  mes_data <- NULL
  mes_ranges <- NULL
  
  #afe
  afe_data <- data[type_level==T&type=="afe"] 
  afe_ranges_template <- annotation$afe$genomic_ranges[paste(afe_data$gene_id,"template",sep="_")]
  afe_ranges_afe <- annotation$afe$genomic_ranges[paste(afe_data$gene_id,afe_data$type,sep="_")]
  afe_data[,location_level:=
             ((start(afe_ranges_template)==start)&(end(afe_ranges_template)==end))
            |
             ((start(afe_ranges_afe)==start)&(end(afe_ranges_afe)==end))
  ]
  data[type_level==T&type=="afe",location_level:=afe_data$location_level]
  afe_data <- NULL
  afe_ranges_template <- NULL
  afe_ranges_afe <- NULL
  
  #ale
  ale_data <- data[type_level==T&type=="ale"] 
  ale_ranges_template <- annotation$ale$genomic_ranges[paste(ale_data$gene_id,"template",sep="_")]
  ale_ranges_ale <- annotation$ale$genomic_ranges[paste(ale_data$gene_id,ale_data$type,sep="_")]
  ale_data[,location_level:=
             ((start(ale_ranges_template)==start)&(end(ale_ranges_template)==end))
           |
             ((start(ale_ranges_ale)==start)&(end(ale_ranges_ale)==end))
           ]
  data[type_level==T&type=="ale",location_level:=ale_data$location_level]
  ale_data <- NULL
  ale_ranges_template <- NULL
  ale_ranges_afe <- NULL
  
  #mee
  mee_data <- data[type_level==T&type=="mee"] 
  mee_ranges_template <- annotation$mee$genomic_ranges[paste(mee_data$gene_id,"template",sep="_")]
  mee_ranges_mee <- annotation$mee$genomic_ranges[paste(mee_data$gene_id,mee_data$type,sep="_")]
  mee_data[,location_level:=
             ((start(mee_ranges_template)==start)&(end(mee_ranges_template)==end))
           |
             ((start(mee_ranges_mee)==start)&(end(mee_ranges_mee)==end))
           ]
  data[type_level==T&type=="mee",location_level:=mee_data$location_level]
  mee_data <- NULL
  mee_ranges_template <- NULL
  mee_ranges_afe <- NULL
  
  #a3
  a3_data <- data[type_level==T&type=="a3"]
  a3_ranges <- annotation$a3$genomic_ranges[a3_data$gene_id]
  strand_map_end = list("+"=-1,"-"=0)
  strand_map_start = list("+"=0,"-"=1)
  strand_start <- unlist(strand_map_start[a3_data$strand],use.names = F)
  strand_end <- unlist(strand_map_end[a3_data$strand],use.names = F)

  a3_data[,location_level:=((start(a3_ranges)==(start+strand_start))
                            &(end(a3_ranges)==(end+strand_end)))]
  data[type_level==T&type=="a3",location_level:=a3_data$location_level]
  a3_data <- NULL
  a3_ranges <- NULL
  
  #a5
  a5_data <- data[type_level==T&type=="a5"]
  a5_ranges <- annotation$a5$genomic_ranges[a5_data$gene_id]
  strand_map_start = list("+"=1,"-"=0)
  strand_map_end = list("+"=0,"-"=-1)
  strand_start <- unlist(strand_map_start[a5_data$strand],use.names = F)
  strand_end <- unlist(strand_map_end[a5_data$strand],use.names = F)
  
  a5_data[,location_level:=((start(a5_ranges)==(start+strand_start))
                            &(end(a5_ranges)==(end+strand_end)))]
  data[type_level==T&type=="a5",location_level:=a5_data$location_level]
  a5_data <- NULL
  a5_ranges <- NULL
  
  
  counts_location$level <- "location"
  counts_location$tp <- sum(data$location_level)
  counts_location$fp <- sum(data$location_level==F)
  counts_location$fn  <- sum(true_data$diff_spliced)-counts_location$tp
  counts <- rbind(counts,counts_location)
  
  return (counts)
}

########### eventpointer ####################

#read data
read_eventpointer <- function(path){
  file <- list.files(path,pattern = "DAS.txt$")
  if(length(file)!=1){
    stop("No input file found!")
  }
  path <-file.path(path,file)
  print(paste("Reading eventpointer output:",path))
  data <- fread(path,header=F,skip = 1)
  
  #formatting
  data <- data[,.(gene_id = V2, type = V3, location = V4, p_value =V5 )]
  data[,location:=gsub(".*:","",data$location)]
  sep <- separate(data[,.(location)],location,c("start","end"),sep="-")
  data[,location:=NULL]
  data[,start:=as.integer(sep$start)]
  data[,end:=as.integer(sep$end)]
  
  data <- data[startsWith(gene_id, "ENSG")]
  
  #map for type parsing
  map <- list("Alternative Last Exon"="ale",
              "Alternative 3' Splice Site"="a3",
              "Mutually Exclusive Exons"="mee",
              "Retained Intron"="ir",
              "Alternative 5' Splice Site"="a5",
              "Cassette Exon"="es",
              "Alternative First Exon"="afe",
              "Complex Event"="mes")
  data[,type:=as.factor(unlist(map[type],use.names = F))]
  
  
  return(data)
}

#count
count_eventpointer <- function(data,true_data,annotation,txdb){
  
  #lookup lists
  diff_spliced_list <- as.list(true_data$diff_spliced)
  names(diff_spliced_list) <- true_data$gene_id
  
  type_list <- as.list(true_data$type)
  names(type_list) <- true_data$gene_id
  
  #add strand
  data[,strand:=as.character(strand(genes(txdb)[data$gene_id]))]
  
  #gene_level
  counts <- data.table()
  gene_level <-unlist(diff_spliced_list[data$gene_id],use.names = F) 
  data[,gene_level:=gene_level]
  unique_data <- unique(data[,.(gene_id,gene_level)])
  counts$level <- "gene"
  counts$tp <- sum(unique_data$gene_level)
  counts$fp <- sum(unique_data$gene_level==F)
  counts$fn  <- sum(true_data$diff_spliced)-counts$tp
  
  
  #type_level
  counts_type <- data.table() 
  data[,type_level := gene_level&(unlist(type_list[gene_id],use.names = F)==type)]
  data[,real_type := unlist(type_list[gene_id],use.names=F)]
  unique_data <- unique(data[,.(gene_id,gene_level,type_level,type)])
  counts_type$level <- "type"
  counts_type$tp <- sum(unique_data$type_level)
  counts_type$fp <- sum(unique_data$type_level==F)
  counts_type$fn  <- sum(true_data$diff_spliced)-counts_type$tp
  counts <- rbind(counts,counts_type)
  
  #location_level
  counts_location <- data.table()
  data[type_level==F,location_level:=F]
  
  #es
  es_data <- data[type_level==T&type=="es"]
  es_ranges <- annotation$es$introns
  es_ranges <- es_ranges[es_data$gene_id]
  es_data[,location_level:=((start(es_ranges)==(start+1))&(end(es_ranges)==(end-1)))]
  data[type_level==T&type=="es",location_level:=es_data$location_level]
  es_data <- NULL
  es_ranges <- NULL

  
  #ir
  ir_data <- data[type_level==T&type=="ir"]
  ir_ranges <- annotation$ir$genomic_ranges[ir_data$gene_id]
  ir_data[,location_level:=(start(ir_ranges)==(start+1))&(end(ir_ranges)==(end-1))]
  data[type_level==T&type=="ir",location_level:=ir_data$location_level]
  ir_data <- NULL
  ir_ranges <- NULL
  
  #mes 
  mes_data <- data[type_level==T&type=="mes"]
  mes_ranges <- annotation$mes$introns
  mes_ranges <- mes_ranges[mes_data$gene_id]
  mes_data[,location_level:=((start(mes_ranges)==(start+1))&(end(mes_ranges)==(end-1)))]
  data[type_level==T&type=="mes",location_level:=mes_data$location_level]
  mes_data <- NULL
  mes_ranges <- NULL


  #mee
  mee_data <- data[type_level==T&type=="mee"]
  mee_ranges <- annotation$mee$introns[mee_data$gene_id]
  mee_data[,location_level:=(start(mee_ranges)==(start+1))&(end(mee_ranges)==(end-1))]
  data[type_level==T&type=="mee",location_level:=mee_data$location_level]
  mee_data <- NULL
  mee_ranges <- NULL

  #a3
  a3_data <- data[type_level==T&type=="a3"]
  a3_ranges <- annotation$a3$introns[a3_data$gene_id]
  a3_data[,location_level:=(start(a3_ranges)==(start+1))&(end(a3_ranges)==(end-1))]
  data[type_level==T&type=="a3",location_level:=a3_data$location_level]
  a3_data <- NULL
  a3_ranges <- NULL
  
  #a5
  a5_data <- data[type_level==T&type=="a5"]
  a5_ranges <- annotation$a5$introns[a5_data$gene_id]
  a5_data[,location_level:=(start(a5_ranges)==(start+1))&(end(a5_ranges)==(end-1))]
  data[type_level==T&type=="a5",location_level:=a5_data$location_level]
  a5_data <- NULL
  a5_ranges <- NULL
  
  #afe
  afe_data <- data[type_level==T&type=="afe"]
  afe_ranges <- annotation$afe$introns[afe_data$gene_id]
  afe_data[,location_level:=(start(afe_ranges)==start)|(end(afe_ranges)==end)]
  data[type_level==T&type=="afe",location_level:=afe_data$location_level]
  
  #ale
  ale_data <- data[type_level==T&type=="ale"]
  ale_ranges <- annotation$ale$introns[ale_data$gene_id]
  ale_data[,location_level:=(start(ale_ranges)==start)|(end(ale_ranges)==end)]
  data[type_level==T&type=="ale",location_level:=ale_data$location_level]
  
  counts_location$level <- "location"
  counts_location$tp <- sum(data$location_level)
  counts_location$fp <- sum(data$location_level==F)
  counts_location$fn  <- sum(true_data$diff_spliced)-counts_location$tp
  counts <- rbind(counts,counts_location)
  
  return(counts)
  
}



############ aspli ####################
read_aspli <- function(path){
  dir <- file.path(path,"exons")
  file <- list.files(dir,pattern = ".du.tab$")
  if(length(file)!=1){
    stop("No input file found!")
  }
  exon_path <-file.path(dir,file)
  print(paste("Reading aspli output:",exon_path))
  exon_data <- fread(exon_path)
  
  dir <- file.path(path,"introns")
  file <- list.files(dir,pattern = ".du.tab$")
  if(length(file)!=1){
    stop("No input file found!")
  }
  intron_path <-file.path(dir,file)
  print(paste("Reading aspli output:",intron_path))
  intron_data <- fread(intron_path)
  data <- rbind(exon_data,intron_data)
  data <- data[event!="-"]
  
  #rename events
  map <- list("Alt3ss"="a3", "Alt3ss*"="a3",  "Alt5ss"="a5" ,"Alt5ss*"="a5", "ES"="es","ES*"="es","IR"="ir","IR*"="ir" )
  data[,type:=unlist(map[event],use.names = F)]
  return(data[,.(gene_id=symbol,p_value=pvalue,start,end,type)])
}

count_aspli <- function(data,true_data,annotation,txdb){
  #lookup lists
  diff_spliced_list <- as.list(true_data$diff_spliced)
  names(diff_spliced_list) <- true_data$gene_id
  
  type_list <- as.list(true_data$type)
  names(type_list) <- true_data$gene_id
  
  #add strand
  data[,strand:=as.character(strand(genes(txdb)[data$gene_id]))]
  
  #include in true_data?
  includes <- unique(data$type)
  nrow_true_data <- nrow(true_data)
  true_data <- true_data[type %in% includes]
  print(paste0(nrow_true_data-nrow(true_data)," genes were dropped"))
  
  #gene level
  counts <- data.table()
  gene_level <-unlist(diff_spliced_list[data$gene_id],use.names = F) 
  data[,gene_level:=gene_level]
  unique_data <- unique(data[,.(gene_id,gene_level)])
  counts$level <- "gene"
  counts$tp <- sum(data$gene_level==T)
  counts$fp <- sum(data$gene_level==F)
  counts$fn  <- sum(true_data$diff_spliced)-counts$tp
  
  #type_level
  counts_type <- data.table() 
  #data[,real_type := unlist(type_list[gene_id],use.names=F)]
  data[,type_level := gene_level&real_type==type]
  unique_data <- unique(data[,.(gene_id,gene_level,type_level,type)])
  counts_type$level <- "type"
  counts_type$tp <- sum(unique_data$type_level)
  counts_type$fp <- sum(unique_data$type_level==F)
  counts_type$fn  <- sum(true_data$diff_spliced)-counts_type$tp
  counts <- rbind(counts,counts_type)
  
  #location_level
  counts_location <- data.table()
  data[type_level==F,location_level:=F]
  
  #es
  es_data <- data[type_level==T&type=="es"]
  es_ranges <- annotation$es$genomic_ranges[es_data$gene_id]
  es_data[,location_level:=((start(es_ranges)==start)&(end(es_ranges)==end))]
  data[type_level==T&type=="es",location_level:=es_data$location_level]
  es_data <- NULL
  es_ranges <- NULL
  
  #ir
  ir_data <- data[type_level==T&type=="ir"]
  ir_ranges <- annotation$ir$genomic_ranges[ir_data$gene_id]
  ir_data[,location_level:=((start(ir_ranges)==start)&(end(ir_ranges)==end))]
  data[type_level==T&type=="ir",location_level:=ir_data$location_level]
  ir_data <- NULL
  ir_ranges <- NULL
  
  #a3
  a3_data <- data[type_level==T&type=="a3"]
  a3_ranges <- annotation$a3$genomic_ranges[a3_data$gene_id]
  a3_data[,location_level:=((start(a3_ranges)==start)&(end(a3_ranges)==end))]
  data[type_level==T&type=="a3",location_level:=a3_data$location_level]
  a3_data <- NULL
  a3_ranges <- NULL
  
  #a5
  a5_data <- data[type_level==T&type=="a5"]
  a5_ranges <- annotation$a5$genomic_ranges[a5_data$gene_id]
  a5_data[,location_level:=((start(a5_ranges)==start)&(end(a5_ranges)==end))]
  data[type_level==T&type=="a5",location_level:=a5_data$location_level]
  a5_data <- NULL
  a5_ranges <- NULL
  
  counts_location$level <- "location"
  counts_location$tp <- sum(data$location_level==T)
  counts_location$fp <- sum(data$location_level==F)
  counts_location$fn  <- sum(true_data$diff_spliced)-counts_location$tp
  counts <- rbind(counts,counts_location)

  
  return(counts)
}

############# dexseq ################
read_dexseq <- function(path){
  print(paste("Reading dexseq output:",path))
  data <- fread(path,select = "groupID")
  return(data)
}

######### majiq ######################
read_majiq <- function(path){
  file <- list.files(path,pattern = ".tsv$")
  if(length(file)!=1){
    stop("No input file found!")
  }
  path <-file.path(path,file)
  print(paste("Reading aspli output:",path))
  data <- fread(path)

  setnames(data,"Gene ID","gene_id")
  
  #prepare es
  es <- data[ES==T]
  
  #prepare a3
  a3 <- data[A3SS == T]
  
  #prepare a5
  a5 <- data[A5SS==T]
  
  #prepare ir
  ir <- data[`IR coords`!=""]
  sep <- separate(ir[,.(`IR coords`)],`IR coords`,c("start","end"),sep="-")
  ir[,start:=as.integer(sep$start)]
  ir[,end:=as.integer(sep$end)]

  
  #prepare a5
  a5_junctions <- strsplit(a5$`Junctions coords`,";")
  locations <- sapply(a5_junctions,function(x){
    splitted <- strsplit(x,"-")
    junction_1 <- as.integer(splitted[[1]])
    junction_2 <- as.integer(splitted[[2]])
    point_1 <- junction_1[!(junction_1%in%junction_2)]
    point_2 <- junction_2[!(junction_2%in%junction_1)]
    start <- min(point_1,point_2)
    end <- max(point_1,point_2)
    return(c("start"=start,"end"=end))
  })
  a5[,start:=locations["start",]]
  a5[,end:=locations["end",]]
  
  #prepare a3
  a3_junctions <- strsplit(a3$`Junctions coords`,";")
  locations <- sapply(1:length(a3_junctions),function(x){
    shift <- 0
    if(a3[x]$A5SS){
      shift <- shift+2
      print("shift")
    }
    splitted <- strsplit(a3_junctions[[x]],"-")
    junction_1 <- as.integer(splitted[[1+shift]])
    junction_2 <- as.integer(splitted[[2+shift]])
    point_1 <- junction_1[!(junction_1%in%junction_2)]
    point_2 <- junction_2[!(junction_2%in%junction_1)]
    start <- min(point_1,point_2)
    end <- max(point_1,point_2)
    return(c("start"=start,"end"=end))
  })
  a3[,start:=locations["start",]]
  a3[,end:=locations["end",]]
  
  #prepare es
  es_junctions <- strsplit(es$`Junctions coords`,";")
  locations <- sapply(1:length(es_junctions),function(x){
    shift <- 0
    if(es[x]$A5SS){
      shift <- shift+2
      print("shift_1")
    }
    if(es[x]$A3SS){
      shift <- shift+2
      print("shift_2")
    }
    splitted <- strsplit(es_junctions[[x]],"-")
    junction_1 <- as.integer(splitted[[1+shift]])
    junction_2 <- as.integer(splitted[[2+shift]])
    start <- min(junction_1,junction_2)
    end <- max(junction_1,junction_2)
    return(c("start"=start,"end"=end))
  })
  es[,start:=locations["start",]]
  es[,end:=locations["end",]]

  
  data <- rbindlist(list(es=es,a3=a3,a5=a5,ir=ir),idcol="type",fill=T)
  return(data)
}

count_majiq <- function(data,true_data,annotation,txdb){
  #lookup lists
  diff_spliced_list <- as.list(true_data$diff_spliced)
  names(diff_spliced_list) <- true_data$gene_id
  
  type_list <- as.list(true_data$type)
  names(type_list) <- true_data$gene_id
  
  #add strand
  data[,strand:=as.character(strand(genes(txdb)[data$gene_id]))]
  
  #include in true_data?
  includes <- unique(data$type)
  nrow_true_data <- nrow(true_data)
  true_data <- true_data[type %in% includes]
  print(paste0(nrow_true_data-nrow(true_data)," genes were dropped"))
  
  #gene_level
  counts <- data.table()
  gene_level <-unlist(diff_spliced_list[data$gene_id],use.names = F) 
  data[,gene_level:=gene_level]
  unique_data <- unique(data[,.(gene_id,gene_level)])
  counts$level <- "gene"
  counts$tp <- sum(unique_data$gene_level)
  counts$fp <- sum(unique_data$gene_level==F)
  counts$fn  <- sum(true_data$diff_spliced)-counts$tp
  
  #type_level
  counts_type <- data.table() 
  data[,type_level := gene_level&real_type==type]
  unique_data <- unique(data[,.(gene_id,gene_level,type_level,type)])
  counts_type$level <- "type"
  counts_type$tp <- sum(unique_data$type_level==T)
  counts_type$fp <- sum(unique_data$type_level==F)
  counts_type$fn  <- sum(true_data$diff_spliced)-counts_type$tp
  counts <- rbind(counts,counts_type)
  
  #location_level
  counts_location <- data.table()
  data[type_level==F,location_level:=F]
  
  #ir
  ir_data <- data[type_level==T&type=="ir"]
  ir_ranges <- annotation$ir$genomic_ranges[ir_data$gene_id]
  ir_data[,location_level:=((start(ir_ranges)==(start+1))&(end(ir_ranges)==(end-1)))]
  data[type_level==T&type=="ir",location_level:=ir_data$location_level]
  ir_data <- NULL
  ir_ranges <- NULL
  
  #a5
  a5_data <- data[type_level==T&type=="a5"]
  a5_ranges <- annotation$a5$genomic_ranges[a5_data$gene_id]
  strand_map_start = list("+"=1,"-"=0)
  strand_map_end = list("+"=0,"-"=-1)
  strand_start <- unlist(strand_map_start[a5_data$strand],use.names = F)
  strand_end <- unlist(strand_map_end[a5_data$strand],use.names = F)
  
  a5_data[,location_level:=((start(a5_ranges)==(start+strand_start))
                            &(end(a5_ranges)==(end+strand_end)))]
  data[type_level==T&type=="a5",location_level:=a5_data$location_level]
  a5_data <- NULL
  a5_ranges <- NULL
  
  #a3
  a3_data <- data[type_level==T&type=="a3"]
  a3_ranges <- annotation$a3$genomic_ranges[a3_data$gene_id]
  strand_map_end = list("+"=-1,"-"=0)
  strand_map_start = list("+"=0,"-"=1)
  strand_start <- unlist(strand_map_start[a3_data$strand],use.names = F)
  strand_end <- unlist(strand_map_end[a3_data$strand],use.names = F)
  
  a3_data[,location_level:=((start(a3_ranges)==(start+strand_start))
                            &(end(a3_ranges)==(end+strand_end)))]
  data[type_level==T&type=="a3",location_level:=a3_data$location_level]
  a3_data <- NULL
  a3_ranges <- NULL
  
  #es
  es_data <- data[type_level==T&type=="es"]
  es_ranges <- annotation$es$introns
  es_ranges <- es_ranges[es_data$gene_id]
  es_data[,location_level:=((start(es_ranges)==(start+1))&(end(es_ranges)==(end-1)))]
  data[type_level==T&type=="es",location_level:=es_data$location_level]
  es_data <- NULL
  es_ranges <- NULL
  
  unique_data <-   unique_data <- unique(data[,.(gene_id,gene_level,type_level,type,location_level)])
  counts_location$level <- "location"
  counts_location$tp <- sum(unique_data$location_level==T)
  counts_location$fp <- sum(unique_data$location_level==F)
  counts_location$fn  <- sum(true_data$diff_spliced)-counts_location$tp
  counts <- rbind(counts,counts_location)

  return(counts)
}

#################### spladder #####################

read_spladder <- function(path){
  map <- list("mutex_exons"="mee","exon_skip"="es","intron_retention"="ir","mult_exon_skip"="mes","alt_3prime"="a3","alt_5prime"="a5")
  data <- rbindlist(lapply(names(map),function(type){
    file <- list.files(path,pattern = paste0("C3_",type,".gene_unique.tsv$"))
    if(length(file)!=1){
      stop(paste0("No input file found for type: ",type))
    }
    print(paste("Reading spladder output:",file.path(path,file)))
    sub_data <- fread(file.path(path,file))
    sub_data[,type:=map[[type]]]
    return(sub_data[,.(gene_id=gene,p_value=p_val,type)])
  }))
  return(data)
}

count_spladder <- function(data,true_data,annotation,txdb){
  #lookup lists
  diff_spliced_list <- as.list(true_data$diff_spliced)
  names(diff_spliced_list) <- true_data$gene_id
  
  type_list <- as.list(true_data$type)
  names(type_list) <- true_data$gene_id
  
  #add strand
  data[,strand:=as.character(strand(genes(txdb)[data$gene_id]))]
  
  #include in true_data?
  includes <- unique(data$type)
  nrow_true_data <- nrow(true_data)
  true_data <- true_data[type %in% includes]
  print(paste0(nrow_true_data-nrow(true_data)," genes were dropped"))
  
  
  #gene_level
  counts <- data.table()
  gene_level <-unlist(diff_spliced_list[data$gene_id],use.names = F) 
  data[,gene_level:=gene_level]
  unique_data <- unique(data[,.(gene_id,gene_level)])
  counts$level <- "gene"
  counts$tp <- sum(unique_data$gene_level==T)
  counts$fp <- sum(unique_data$gene_level==F)
  counts$fn  <- sum(true_data$diff_spliced)-counts$tp
  
  #type_level
  counts_type <- data.table() 
  data[,type_level := gene_level&real_type==type]
  unique_data <- unique(data[,.(gene_id,gene_level,type_level,type)])
  counts_type$level <- "type"
  counts_type$tp <- sum(unique_data$type_level==T)
  counts_type$fp <- sum(unique_data$type_level==F)
  counts_type$fn  <- sum(true_data$diff_spliced)-counts_type$tp
  counts <- rbind(counts,counts_type)
  
  return(counts)
}

################# jsplice ###########################

read_jsplice <- function(path){
  file <- list.files(path,pattern = ".txt$")
  if(length(file)!=1){
    stop("No/multiple input file/s found!")
  }
  path <-file.path(path,file)
  print(paste("Reading jsplice output:",path))
  data <- fread(path)
  
  #rename columns
  setnames(data,"Gene_name|ID","gene_id")
  setnames(data,"ASM_type","type")
  
  #remove without id
  nrow <- nrow(data)
  data <- data[gene_id!=""]
  print(paste0(nrow-nrow(data)," lines had no gene_id"))
  
  nrow <- nrow(data)
  data <- data[!grep("Unknown",data$type)]
  print(paste0(nrow-nrow(data)," lines had unknown as type"))
  
  

  return(data)
}
