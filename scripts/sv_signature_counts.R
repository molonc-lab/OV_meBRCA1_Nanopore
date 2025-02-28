# Written by Felicity Newell
# v1.0
# Date: 09/02/2022

library(copynumber)
library(Palimpsest)
library(R.matlab)
library(bedr)

args=commandArgs(trailingOnly = T)

checkArgs<-function(args) {
  if(length(args) != 2)  {
    print("Needs 2 arguments:")
    print("Command is: Rscript --vanilla ../sv_signature_counts.R sampleFilePath outputFilePath")
    return(FALSE)
  } else if (!file.exists(args[1]))  {
    print(paste("Samplefile ",args[1], "does not exist"))
    return(FALSE)
  }
  return(TRUE)
}

process.file<-function(qsvData,sampleName) {
  qsv$distance=abs(qsv$chr_to_bkpt - qsv$chr_from_bkpt)
  qsv$distance[qsv$annotation == "CTX"]<-1
  qsv<-qsv[!grepl("GL",qsv$chr_from) & !grepl("GL",qsv$chr_to),]
  qsv$length_cat<-apply(qsv,1,function(x) {
    distance=as.numeric(x["distance"])
    type=x["annotation"]
    cat=""
    if (type == "CTX") {
      cat="-"
    } else if (distance < 1000) {
      cat="remove"
    } else if (distance >= 1000 & distance < 10000) {
      cat="1kb-10kb"
    } else if (distance >=10000 & distance < 100000) {
      cat="10kb-100kb"
    } else if(distance >=100000 & distance < 1000000) {
      cat="100kb-1mb"
    } else if(distance >=1000000 & distance < 10000000) {
      cat="1mb-10mb"
    } else {
      if (distance >=10000000) {
        cat=">10mb"
      }
    }
    return(cat)
  })
  
  qsv<-qsv[qsv$length_cat != "remove",]
  qsv<-qsv[qsv$annotation != "ITX",]
  
  if (nrow(qsv) > 0) {
    chr.from<-qsv[,c("chr_from","chr_from_bkpt","sv_id","annotation","length_cat")]
    chr.to<-qsv[,c("chr_to","chr_to_bkpt","sv_id","annotation","length_cat")]
    colnames(chr.to)<-c("chr_from","chr_from_bkpt","sv_id","annotation","length_cat")
    
    chr.all<-rbind(chr.from,chr.to)
    chr.all$chr_from<-as.numeric(chr.all$chr_from)
    chr.all<-chr.all[order(chr.all$chr_from,chr.all$chr_from_bkpt),]
    
    #Palimpsest - find clustered breakpoints
    bkp <- data.frame(chr=chr.all$chr_from,pos=chr.all$chr_from_bkpt,dist=NA,stringsAsFactors=F)
    tmp <- sub("X",23,sub("Y",24,bkp$chr));
    bkp <- bkp[order(bkp$chr,bkp$pos),]
    for(c in unique(bkp$chr)){
      ind <- which(bkp$chr==c)
      if(length(ind) > 1){
        bkp[ind,] <- bkp[ind[order(bkp[ind,"pos"])],]
        bkp[ind[-1],"dist"] <- bkp[ind[-1],"pos"]-bkp[ind[-length(ind)],"pos"]
      }
    }
    bed <- bkp[,c(1,2,2,3)]
    bed$chr<-paste("chr",bed$chr,sep="")
    indi <- bed2index(bed, sort = TRUE)
    tmp <- cluster.region(x = indi,distance = "1000000",verbose = T)
    tmp[,"chr"] <- data.frame(do.call('rbind', strsplit(as.character(tmp$index),':',fixed=TRUE)))[1]
    tmp[,c("start","end")] <- data.frame(do.call('rbind', strsplit(as.character(tmp$index),'-',fixed=TRUE)))[2]
    bkp$cluster <- tmp[,"regionIndex"]
    bkp$clustered <- "non_clustered"
    tt <- table(bkp$cluster)
    
    for(clst in names(which(tt >= 10))){
      bkp[which(bkp$cluster==clst),"clustered"] <-"clustered"
    }
    chr.all$palimpsest_cluster_type<-bkp$clustered
    
    #Set sv type
    chr.all$sv_type<-apply(chr.all,1,function(x) {
      a=x["annotation"]
      t=""
      if (a == "CTX") {
        t="translocation"
      } else if (a == "INV/ITX") {
        t="inversion"
      } else if(a == "DEL/ITX") {
        t="deletion"
      } else if (a == "DUP/INS/ITX") {
        t="duplication"
      }
      return (t)
    })
    
    #Set other end of clustered SVs to clustered
    clustered<-chr.all$sv_id[chr.all$palimpsest_cluster_type == "clustered"]
    chr.all$palimpsest_cluster_type[chr.all$sv_id %in% clustered]<-"clustered"
    
    chr.all$sv_sig_type<-paste(chr.all$sv_type,chr.all$palimpsest_cluster_type,chr.all$length_cat,sep=":")
    return(chr.all)
  } else {
    return(NULL)
  }
  
}

if(checkArgs(args)) {
  sampleFile=args[1]
  outFileName=args[2]

  #List of samples with qsv HighConfidence.dccfile path (column1) and sample name (column2). No header in file
  samples<-read.delim(sampleFile,stringsAsFactors = F,sep="\t",check.names=F,comment.char = "",fill=T,header=F)
  
  colnames(samples)=c("file","sampleName")
  
  subtypes.vec<-c("deletion:clustered:1kb-10kb","deletion:clustered:10kb-100kb","deletion:clustered:100kb-1mb","deletion:clustered:1mb-10mb","deletion:clustered:>10mb",
                  "duplication:clustered:1kb-10kb","duplication:clustered:10kb-100kb","duplication:clustered:100kb-1mb","duplication:clustered:1mb-10mb","duplication:clustered:>10mb",
                  "inversion:clustered:1kb-10kb","inversion:clustered:10kb-100kb","inversion:clustered:100kb-1mb","inversion:clustered:1mb-10mb","inversion:clustered:>10mb",
                  "translocation:clustered:-",
                  "deletion:non_clustered:1kb-10kb","deletion:non_clustered:10kb-100kb","deletion:non_clustered:100kb-1mb","deletion:non_clustered:1mb-10mb","deletion:non_clustered:>10mb",
                  "duplication:non_clustered:1kb-10kb","duplication:non_clustered:10kb-100kb","duplication:non_clustered:100kb-1mb","duplication:non_clustered:1mb-10mb","duplication:non_clustered:>10mb",
                  "inversion:non_clustered:1kb-10kb","inversion:non_clustered:10kb-100kb","inversion:non_clustered:100kb-1mb","inversion:non_clustered:1mb-10mb","inversion:non_clustered:>10mb",
                  "translocation:non_clustered:-")
  types.vec<-c(rep("clustered",16),rep("non_clustered",16))
  allcounts.df=data.frame("Type"=types.vec,"SubType"=subtypes.vec)
  
  if (length(unique(samples$sampleName)) != nrow(samples)) {
    stop("Sample IDs must be unique")
  } else {
    
    for (i in 1:nrow(samples)) {
      print(i)
      #Read in qsv files
      file<-samples$file[i]
      if (file.exists(file)) {
        qsv<-read.table(file,stringsAsFactors = F,sep="\t",check.names=F,comment.char = "#",fill=T,header=T)
      } else {
        stop(paste("File does not exist: ",file,sep=" "))
      }
      #Get sample counts
      if(nrow(qsv) > 0) {
        sample.counts=process.file(qsv,samples$sampleName[i])
        if (!is.null(sample.counts)) {
          counts<-as.data.frame(table(sample.counts$sv_sig_type))
          allcounts.df[,samples$sampleName[i]]<-counts$Freq[match(allcounts.df$SubType,counts$Var1)]
        } else {
          print(paste0(samples$sampleName[i], " does not have any SVs"))
          allcounts.df[,samples$sampleName[i]]=0
        }
      } else {
        print(paste0(samples$sampleName[i], " does not have any SVs"))
        allcounts.df[,samples$sampleName[i]]=0
      }
    }
    allcounts.df[is.na(allcounts.df)]=0
    
    print(paste("Writing to file: " ,outFileName))
    write.table(allcounts.df,outFileName,sep="\t",quote=F,row.names=F)
  }
}


