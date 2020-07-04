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
  stop(paste0("Method '",method,"' is not supported!"))
}

count_ds <- function(data,method,...){
  if(method=="cash"){
    return(count_cash(data,...))
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
  print(a5_ranges)
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
read_eventpointer <- function(path){
  file <- list.files(path,pattern = "DAS.txt$")
  if(length(file)!=1){
    stop("No input file found!")
  }
  path <-file.path(path,file)
  print(paste("Reading eventpointer output:",path))
  data <- fread(path,header=F,skip = 1)
  return(data[,.(gene_id=V2,p_value=V5)])
}

############ aspli ####################
read_aspli <- function(path){
  path <- file.path(path,"genes")
  file <- list.files(path,pattern = ".de.tab$")
  if(length(file)!=1){
    stop("No input file found!")
  }
  path <-file.path(path,file)
  print(paste("Reading aspli output:",path))
  data <- fread(path)
  return(data[,.(gene_id=symbol,p_value=pvalue)])
}

############# dexseq ################
read_dexseq <- function(path){
  print(paste("Reading dexseq output:",path))
  data <- fread(path,select = "groupID")
  return(data)
}

