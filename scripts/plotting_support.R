library(ggtranscript)


plot_icgc_changed_meth <- function(icgc_meth_plot, gene_name, palette, icgc_mRNA_plot, transcript_annot){
  meth_plot <-icgc_meth_plot %>%
    filter(name == gene_name) %>%
    ggplot( aes(x = start, y= beta, group = Sample.ID, color =Donor))+
    geom_point()+
    geom_line(aes(linetype = Sample_Group))+
    scale_linetype_manual(values= c("Primary Tumour" = 1, "Relapse" = 3), "Sample Type")+
    scale_color_manual(values= palette,"Sample")+
    scale_y_continuous(limits = c(0, 1), name = "Methylation (beta)")+
    theme_bw()+
    theme(legend.position = "bottom", axis.title.x = element_blank(),legend.box="vertical", legend.margin = margin(), legend.spacing.x = unit(0, "mm"),  legend.spacing.y = unit(0.1, "mm"))

  
  gg_gene_annotation<-  transcript_annot %>%
    filter(gene_name %in% icgc_meth_plot$name)%>%
    ggplot(aes(
      xstart = start,
      xend = end,
      y = gene_name
    )) +
    geom_range() +
    geom_intron(
      data = to_intron( transcript_annot %>%
                          filter(gene_name %in% icgc_meth_plot$name) , "transcript_id"),
      aes(strand = strand)
    )+
    theme_bw()+
    coord_cartesian(xlim =  c(min(icgc_meth_plot$start), max(icgc_meth_plot$start)) )+
    theme(panel.grid = element_blank(),panel.border = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          axis.title.y = element_blank())
  
  
  icgc_mRNA_plot  <- icgc_mRNA_plot  %>%
    filter(`Associated Gene Name` == gene_name)
  
  if(dim(icgc_mRNA_plot)[1] >0){
    mRNA_plot <- icgc_mRNA_plot  %>%
      ggplot(aes(y = log_tpm, x =Sample_Group, label = donorLabel, group = donorLabel ))+
      geom_line(color  = "grey")+
      geom_point(aes(color = donorLabel))+
      scale_color_manual(values= palette,"Sample")+
      facet_wrap(~`Associated Gene Name`)+
      theme_bw()+
      ylab("mRNA (logTPM)")+
      theme(legend.position = "none", axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            strip.background = element_rect(fill = "white"))
    
    plot <- ggarrange(ggarrange(meth_plot+
                                  theme(legend.position = "none"),gg_gene_annotation, nrow = 2, heights = c(0.8, 0.2)),  mRNA_plot, widths = c(0.7,0.3))
  } else{
    plot <- ggarrange(meth_plot ,gg_gene_annotation, nrow = 2, heights = c(0.8, 0.2), common.legend = T, legend = "none")
  }
  
  return(plot)
}



draw_annotated_meth <- function(changed_cpg,name,palette, mRNA_pdx_plot, transcript_annot){
  cpg_methylation <- changed_cpg %>%
    filter(name== name)
  
  gene_name_within <- transcript_annot %>%
    filter(seqnames %in% unique(cpg_methylation$chrom), start <= max(cpg_methylation$chromEnd), end >= min(cpg_methylation$chromStart)) %>%
    pull(gene_name) %>%
    unique(.)
  
  if(length(gene_name_within) >0){
    
    gg_gene_annotation<-  transcript_annot %>%
      filter(gene_name %in% gene_name_within)%>%
      ggplot(aes(
        xstart = start,
        xend = end,
        y = gene_name
      )) +
      geom_range() +
      geom_intron(
        data = to_intron( transcript_annot %>%
                            filter(gene_name %in% gene_name_within) , "transcript_id"),
        aes(strand = strand)
      )+
      theme_bw()+
      coord_cartesian(xlim =  c(min(cpg_methylation$chromStart), max(cpg_methylation$chromEnd)) )+
      theme(panel.grid = element_blank(),panel.border = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks = element_blank(),
            axis.title.y = element_blank())
    
    
    gene_methylation <- cpg_methylation %>%
      ggplot(aes(x =start, y = fraction_mod/100, color = HR, group = alt_name))+
      geom_point(alpha = 0.3)+
      geom_line(alpha = 0.7)+
      theme_bw()+
      scale_color_manual(values = palette)+
      scale_y_continuous(labels = scales::percent, name = "%5mC")+
      coord_cartesian(xlim =  c(min(cpg_methylation$start), max(cpg_methylation$end))  )
    
    mRNA_pdx_plot <- mRNA_pdx_plot %>%
      filter(gene_name == name) 
    
    if(dim(mRNA_pdx_plot)[1] >0){
      mRNA_plot <- mRNA_pdx_plot %>%
        ggplot(aes(x = HR, y = log_tpm, color = HR, label = alt_name))+
        geom_point()+
        facet_wrap(~gene_name, scales = "free_y")+
        scale_colour_manual(values = palette)+
        theme_bw()+
        geom_text_repel(nudge_x = 0.4,segment.color = 'grey50')+
        theme(legend.position = "none", axis.title.x = element_blank(),strip.background = element_rect(fill = "white"),
              axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))+
        ylab("mRNA (logTPM)")
      plot <- ggarrange(ggarrange(gene_methylation+
                                    theme(axis.title.x = element_blank(), legend.position = "none"),gg_gene_annotation, nrow = 2, heights = c(0.8, 0.2)),  mRNA_plot, widths = c(0.7,0.3))
    }else{
      plot <- ggarrange(gene_methylation+
                          theme(axis.title.x = element_blank()), gg_gene_annotation, nrow = 2, align = "v", common.legend = T, heights = c(0.8, 0.2),
                        legend = "none")
    }
    return(plot)}
}

match_savana_bp_meth <- function(T_read_mapping_RDS, T_modkit_extract, N_modkit_extract, type, T_read_ref_mapping_RDS = NA){
  select_col_name <- c("chrom" , "SV_start", "SV_end", "REF","ALT", "ID","SV_type", "SV_len" )
  
  if(type == "DEL"){
  read_mapping <- readRDS(T_read_mapping_RDS) %>%
    filter(str_detect(ALT, "DEL")) %>%
    dplyr::rename(chrom = `#CHROM`)
  } else if (type == "BND") {
    select_col_name <- c("chrom" , "SV_start", "SV_end", "ALT", "ID","SV_type", "SV_len", "chr_end", "BND_type")
    read_mapping <- rbind(readRDS(T_read_mapping_RDS), readRDS(T_read_ref_mapping_RDS)) %>%
      mutate(SV_type = "BND", SV_len = NA, REF = NA)
  }
  else{
    read_mapping <- readRDS(T_read_mapping_RDS)
  }
  
  alt_reads <- read_mapping %>%
    dplyr::select(all_of(c(select_col_name, "read_id",  "read_type"))) %>%
    distinct() %>%
    left_join(fread(T_modkit_extract) %>%
                filter(mod_code == "m") %>%
                dplyr::select(read_id, forward_read_position, ref_position,chrom, mod_strand, ref_strand, read_length, mod_qual, base_qual), 
              by = join_by(chrom == chrom, read_id),
              relationship = "many-to-many") %>%
    mutate(mod = case_when(mod_qual >= 0.6 ~ "Me",
                           mod_qual <= 0.4 ~ "Un",
                           T ~ NA)) %>%
    filter(base_qual >=20, !is.na(mod),  ref_position != -1, (abs(ref_position -SV_start )<= 5000 | abs(ref_position -SV_end )<= 5000)) 
  if(type == "BND"){
    alt_reads <- alt_reads %>%
      filter( abs(SV_start- ref_position) <abs(SV_end- ref_position))
  }
  
  alt_reads <- alt_reads%>%
    group_by(ID, read_type) %>%
    filter(n()>2) %>%
    ungroup() %>%
    group_by(ID, read_type) %>%
    filter(n_distinct(read_id)>2) %>%
    ungroup() %>%
    group_by(ID) %>%
    filter(n_distinct(read_type)== 2) %>%
    ungroup()
  
  alt_P <- alt_reads%>%
    dplyr::select(all_of(c(select_col_name ,"ref_position"))) %>%
    distinct() %>%
    left_join(fread(N_modkit_extract) %>%
                filter(mod_code == "m") %>%
                dplyr::select(read_id, forward_read_position, ref_position,chrom, mod_strand, ref_strand, read_length, mod_qual, base_qual),
              by = join_by(chrom, ref_position),
              relationship = "many-to-many") %>%
    mutate(read_type = "parent") %>%
    mutate(mod = case_when(mod_qual >= 0.6 ~ "Me",
                           mod_qual <= 0.4 ~ "Un",
                           T ~ NA)) %>%
    filter(base_qual >=20, !is.na(mod),  ref_position != -1, (abs(ref_position -SV_start )<= 5000 | abs(ref_position -SV_end )<= 5000)) 
  
  if(type == "BND"){
    alt_P <- alt_P%>%
      filter(abs(SV_start- ref_position) <abs(SV_end- ref_position))
  }
  
  alt_P <- alt_P%>%
    group_by(ID, read_type) %>%
    filter(n()>2) %>%
    ungroup() %>%
    group_by(ID, read_type) %>%
    filter(n_distinct(read_id)>2) %>%
    ungroup() 
  
  total_reads <- rbind(alt_reads, alt_P) %>%
    mutate(read_type = factor(read_type, levels = c("alt", "ref", "parent"))) 
  
  return(total_reads)
}

isogenic_overlap_svs <-function(h453_reads, a630_reads, type){
  if (type == "INS") {
    selected_column = c("chrom", "SV_start", "REF", "SV_type")
  }else if (type == "DUP"){
    selected_column = c("chrom", "SV_start", "ALT", "SV_type")
  }else if (type == "BND"){
    selected_column = c("chrom", "SV_start", "ALT", "SV_type", "BND_type" , "chr_end")
  }
  else{
    selected_column = c("chrom", "SV_start", "REF", "ALT", "SV_type")
  }
  
  overlapped_reads <- h453_reads %>%
    dplyr::select(all_of(selected_column), SV_end_h453 = SV_end, ID_h453 =ID, SV_len_h453 =SV_len ) %>%
    distinct() %>%
    full_join(a630_reads %>%
                dplyr::select(all_of(selected_column), SV_end_a630 = SV_end, ID_a630 =ID, SV_len_a630 = SV_len )%>%
                distinct())%>%
    mutate(commonID = case_when(!is.na(ID_h453) &!is.na(ID_a630) ~ paste0(ID_h453, "_common"), 
                                !is.na(ID_h453) ~ paste0(ID_h453, "_h453"),
                                !is.na(ID_a630) ~ paste0(ID_a630, "_a630")),
           SV_len =pmax(SV_len_h453, SV_len_a630, na.rm =T),
           SV_end = pmax(SV_end_h453 ,SV_end_a630, na.rm = T))
  
  overlapped_reads_mod <-rbind(h453_reads %>%
                                 dplyr::select(-SV_end, -SV_len) %>%
                                 left_join(overlapped_reads %>%
                                       filter(!is.na(ID_h453)) %>%
                                       dplyr::select(ID = ID_h453, commonID, SV_len, SV_end)) %>%
                           dplyr::select(-ID) %>%
                           dplyr::rename(ID = commonID)%>%
                           dplyr::select(ID, SV_start, read_id, read_type, forward_read_position, ref_position, chrom, mod_strand, ref_strand, read_length, mod_qual, base_qual, mod, SV_len, SV_end) %>%
                           mutate(Sample = "OVCAR8 LP H4-53"),
                         a630_reads %>%
                           dplyr::select(-SV_end, -SV_len) %>%
                           right_join(overlapped_reads %>%
                                        filter(!is.na(ID_a630)) %>%
                                        dplyr::select(ID = ID_a630, commonID, SV_len, SV_end)) %>%
                           dplyr::select(-ID)%>%
                           dplyr::select(ID= commonID, SV_start, read_id, read_type, forward_read_position, ref_position, chrom, mod_strand, ref_strand, read_length, mod_qual, base_qual, mod, SV_len, SV_end)%>%
                           mutate(Sample = "OVCAR8 LP A6-30")) %>%
    filter(str_detect(ID, pattern = "_common"))
  
  overlapped_reads_pileup <- overlapped_reads_mod %>%
    mutate(ref_position = if_else(ref_strand== "-", ref_position -1, ref_position)) %>%
    dplyr::select(ref_position, chrom, mod, read_type, Sample, ID, SV_start, SV_end, SV_len)%>%
    group_by(Sample, ref_position, chrom,read_type, ID, SV_start, SV_end, SV_len) %>%
    filter(n()>2) %>%
    summarise(meth_perc = sum(mod == "Me")/n()) %>%
    ungroup() 
  
  if(type != "DUP"){
  overlapped_reads_pileup_filt <- overlapped_reads_pileup %>%
    filter(read_type == "parent", (meth_perc >=0.8| meth_perc <= 0.2) ) %>%
    dplyr::select(ref_position, chrom,  p_meth_perc = meth_perc, ID, SV_start, SV_end, SV_len) %>%
    distinct() %>%
    inner_join(overlapped_reads_pileup %>%
                 filter(read_type == "ref")) %>%
    filter(!is.na(ref_position)) %>%
    filter(abs(p_meth_perc - meth_perc)<=0.2) %>%
    dplyr::select(ref_position, chrom, Sample, p_meth_perc, ID, SV_start, SV_end, SV_len) %>%
    distinct() %>%
    inner_join(overlapped_reads_pileup %>%
                 filter(read_type != "parent")) %>%
    mutate(relat_pos = case_when(ref_position <= SV_start ~ ref_position -SV_start,
                                 ref_position >= SV_end ~ ref_position -SV_end,
                                 TRUE ~NA))    %>%
    filter(!is.na(relat_pos))%>%
    group_by(SV_len, ID, ref_position) %>%
    filter(n_distinct(read_type)==2) %>%
    ungroup() %>%
    filter(!is.na(ref_position)) %>%
    pivot_wider(values_from = meth_perc, names_from = Sample) %>%
    mutate(read_type = if_else(read_type == "ref", "REF", "ALT"),
           read_type = factor(read_type, levels = c("REF", "ALT"))) %>%
    filter(!is.na(`OVCAR8 LP A6-30`) & !is.na(`OVCAR8 LP H4-53`))
  } else{
    overlapped_reads_pileup_filt <- overlapped_reads_pileup %>%
      filter(read_type == "parent", (meth_perc >=0.8| meth_perc <= 0.2) ) %>%
      dplyr::select(ref_position, chrom,  p_meth_perc = meth_perc, ID, SV_start, SV_end, SV_len) %>%
      distinct() %>%
      inner_join(overlapped_reads_pileup %>%
                   filter(read_type == "ref")) %>%
      filter(!is.na(meth_perc), abs(p_meth_perc - meth_perc)<=0.2) %>%
      dplyr::select(ref_position, chrom, Sample, p_meth_perc, ID, SV_start, SV_end, SV_len) %>%
      distinct() %>%
      inner_join(overlapped_reads_pileup %>%
                   filter(read_type != "parent")) %>%
      mutate(relat_pos = case_when(ref_position <= SV_start ~ ref_position -SV_start,
                                   ref_position > SV_start &ref_position < SV_end ~ (ref_position-SV_start)/SV_len *500,
                                   ref_position >= SV_end ~ ref_position -SV_end +500,
                                   TRUE ~NA))    %>%
      group_by(SV_len, ID, ref_position) %>%
      filter(n_distinct(read_type)==2) %>%
      ungroup() %>%
      drop_na() %>%
      pivot_wider(values_from = meth_perc, names_from = Sample) %>%
      mutate(read_type = if_else(read_type == "ref", "REF", "ALT"),
             read_type = factor(read_type, levels = c("REF", "ALT"))) %>%
      drop_na() 
    }
  
  return(overlapped_reads_pileup_filt)
}

require(VennDiagram)
isogenic_overlap_svs_venn <- function(h453_reads, a630_reads ,type, output_file_name){
  
  if (type == "INS") {
    selected_column = c("chrom", "SV_start", "REF", "SV_type")
  }else if (type == "DUP"){
    selected_column = c("chrom", "SV_start", "ALT", "SV_type")
  }else if (type == "BND"){
    selected_column = c("chrom", "SV_start", "ALT", "SV_type", "BND_type", "chr_end")
  }
  else{
    selected_column = c("chrom", "SV_start", "REF", "ALT", "SV_type")
  }
  
  overlapped_reads <- h453_reads %>%
    dplyr::select(all_of(selected_column), SV_end_453 = SV_end, ID_h453 =ID, SV_len_h453 =SV_len ) %>%
    distinct() %>%
    full_join(a630_reads %>%
                dplyr::select(all_of(selected_column), SV_end_630 = SV_end, ID_a630 =ID, SV_len_a630 = SV_len )%>%
                distinct())%>%
    mutate(commonID = case_when(!is.na(ID_h453) &!is.na(ID_a630) ~ paste0(ID_h453, "_common"), 
                                !is.na(ID_h453) ~ paste0(ID_h453, "_h453"),
                                !is.na(ID_a630) ~ paste0(ID_a630, "_a630"))) 
  
  venn.diagram(list("OVCAR8 A6-meth" =  overlapped_reads %>%
                      filter(!is.na(ID_a630)) %>%
                      pull(commonID),
                    "OVCAR8 H4-KO" =   overlapped_reads%>%
                      filter(!is.na(ID_h453)) %>%
                      pull(commonID))
               ,fill = c("#007e5d", "#f8c828"),
               alpha = c(0.5, 0.5), cex =3, cat.dist = c(0.03, 0.03),cat.pos
               = c(-20, 20), main = paste0(type," (SV)"), main.pos = c(0.5, 0.92), main.cex = 2, cat.cex = 1.5,
               filename = output_file_name)
}


isogenic_overlap_indels_venn <- function(h453_reads, a630_reads ,type, output_file_name){
  
  selected_column = c("chrom", "start", "REF", "ALT", "type")
  
  overlapped_reads <- h453_reads %>%
    dplyr::select(all_of(selected_column), end_453 = end, ID_h453 =ID, indel_size ) %>%
    distinct() %>%
    full_join(a630_reads %>%
                dplyr::select(all_of(selected_column), end_630 = end, ID_a630 =ID, indel_size )%>%
                distinct())%>%
    mutate(commonID = case_when(!is.na(ID_h453) &!is.na(ID_a630) ~ paste0(ID_h453, "_common"), 
                                !is.na(ID_h453) ~ paste0(ID_h453, "_h453"),
                                !is.na(ID_a630) ~ paste0(ID_a630, "_a630"))) 
  
  venn.diagram(list("OVCAR8 A6-meth" =  overlapped_reads %>%
                      filter(!is.na(ID_a630)) %>%
                      pull(commonID),
                    "OVCAR8 H4-KO" =   overlapped_reads%>%
                      filter(!is.na(ID_h453)) %>%
                      pull(commonID))
               ,fill = c("#007e5d", "#f8c828"),
               alpha = c(0.5, 0.5), cex =3, cat.dist = c(0.03, 0.03),cat.pos
               = c(-20, 20), main = paste0(type," (INDEL)"), main.pos = c(0.5, 0.92), main.cex = 2, cat.cex = 1.5,
               filename = output_file_name)
}



require(MASS)
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


detect_indels_in_region <- function(reads, region_start, region_end, SampleID) {
  results <- lapply(seq_along(reads), function(i) {
    read <- reads[i]
    
    cigar_str <- as.character(cigar(read))[[1]]
    ops <- explodeCigarOps(cigar_str)[[1]]
    lens <- explodeCigarOpLengths(cigar_str)[[1]]
    
    ref_cursor <- start(read)
    query_cursor <- 1
    
    # Get read name or fallback
    read_id <- names(read)
    if (is.null(read_id) || length(read_id) == 0 || is.na(read_id) || read_id == "") {
      read_id <- paste0("read_", i)
    }
    
    read_indels <- list()
    region_covered <- FALSE
    
    for (j in seq_along(ops)) {
      op <- ops[j]
      len <- lens[j]
      
      if (op %in% c("M", "=", "X")) {
        match_start <- ref_cursor
        match_end <- ref_cursor + len - 1
        
        # Check if match region covers the entire input region
        if (match_start <= region_start && match_end >= region_end) {
          region_covered <- TRUE
        }
        
        ref_cursor <- ref_cursor + len
        query_cursor <- query_cursor + len
        
      } else if (op == "D") {
        del_start <- ref_cursor
        del_end <- ref_cursor + len - 1
        if (region_start <= del_end && region_end >= del_start) {
          read_indels[[length(read_indels) + 1]] <- data.frame(
            read_id = read_id,
            read_alt_type = "Deletion",
            IndelStart = del_start,
            IndelEnd = del_end,
            IndelLength = len,
            SampleID = SampleID,
            region_start = region_start,
            stringsAsFactors = FALSE
          )
        }
        ref_cursor <- ref_cursor + len
        
      } else if (op == "I") {
        if (region_start <= ref_cursor && region_end >= ref_cursor) {
          read_indels[[length(read_indels) + 1]] <- data.frame(
            read_id = read_id,
            read_alt_type = "Insertion",
            IndelStart = ref_cursor,
            IndelEnd = ref_cursor,
            IndelLength = len,
            SampleID = SampleID,
            region_start = region_start,
            stringsAsFactors = FALSE
          )
        }
        query_cursor <- query_cursor + len
      }
    }
    
    # If no indels found and region is covered, report as match
    if (length(read_indels) == 0 && region_covered) {
      return(data.frame(
        read_id = read_id,
        read_alt_type = "Match",
        IndelStart = region_start,
        IndelEnd = region_end,
        IndelLength = 0,
        SampleID = SampleID,
        region_start = region_start,
        stringsAsFactors = FALSE
      ))
    } else if (length(read_indels) > 0) {
      return(do.call(rbind, read_indels))
    } else {
      return(NULL)
    }
  })
  result_df <- do.call(rbind, results)
  
  if (is.null(result_df)) {
    result_df <- data.frame(
      read_id = character(0),
      read_alt_type = character(0),
      IndelStart = integer(0),
      IndelEnd = integer(0),
      IndelLength = integer(0),
      SampleID = SampleID,
      region_start = region_start,
      stringsAsFactors = FALSE
    )
  }
  
  return(result_df)
}

detect_snvs_in_region <- function(reads, region_start, region_end, SampleID) {
  
  results <- lapply(seq_along(reads), function(i) {
    read <- reads[i]
    
    cigar_str <- as.character(cigar(read))[[1]]
    ops  <- explodeCigarOps(cigar_str)[[1]]
    lens <- explodeCigarOpLengths(cigar_str)[[1]]
    
    ref_cursor <- start(read)
    query_cursor <- 1
    
    read_id <- names(read)
    if (is.null(read_id) || length(read_id) == 0 || is.na(read_id) || read_id == "") {
      read_id <- paste0("read_", i)
    }
    
    read_snvs <- list()
    region_covered <- FALSE
    
    for (j in seq_along(ops)) {
      op  <- ops[j]
      len <- lens[j]
      
      if (op %in% c("M", "=")) {
        match_start <- ref_cursor
        match_end   <- ref_cursor + len - 1
        
        if (match_start <= region_start && match_end >= region_end) {
          region_covered <- TRUE
        }
        
        ref_cursor   <- ref_cursor + len
        query_cursor <- query_cursor + len
        
      } else if (op == "X") {
        # Each base is a separate SNV
        for (k in seq_len(len)) {
          snv_pos <- ref_cursor
          
          if (snv_pos >= region_start && snv_pos <= region_end) {
            read_snvs[[length(read_snvs) + 1]] <- data.frame(
              read_id = read_id,
              read_alt_type = "SNV",
              SNV_pos = snv_pos,
              SampleID = SampleID,
              region_start = region_start,
              stringsAsFactors = FALSE
            )
          }
          
          ref_cursor   <- ref_cursor + 1
          query_cursor <- query_cursor + 1
        }
        
      } else if (op == "D") {
        ref_cursor <- ref_cursor + len
        
      } else if (op == "I") {
        query_cursor <- query_cursor + len
      }
    }
    
    # If no SNVs and region fully covered
    if (length(read_snvs) == 0 && region_covered) {
      return(data.frame(
        read_id = read_id,
        read_alt_type = "Match",
        SNV_pos = NA_integer_,
        SampleID = SampleID,
        region_start = region_start,
        stringsAsFactors = FALSE
      ))
    } else if (length(read_snvs) > 0) {
      return(do.call(rbind, read_snvs))
    } else {
      return(NULL)
    }
  })
  
  result_df <- do.call(rbind, results)
  
  if (is.null(result_df)) {
    result_df <- data.frame(
      read_id = character(0),
      read_alt_type = character(0),
      SNV_pos = integer(0),
      SampleID = SampleID,
      region_start = region_start,
      stringsAsFactors = FALSE
    )
  }
  
  return(result_df)
}




match_clair3_meth <- function(T_read_mapping_RDS, T_modkit_extract, N_modkit_extract){
  select_col_name <- c("chrom" , "start", "end", "REF","ALT", "ID","type", "indel_size" )
  
  read_mapping <- readRDS(T_read_mapping_RDS) %>%
    filter(read_type != "FAIL") %>%
    dplyr::rename(chrom = `#CHROM`) %>%
    mutate(start = POS, end = if_else(type  == "del", POS + indel_size, POS),
           read_type = case_when(read_type == "ALT" ~ "alt", 
                                 read_type == "REF" ~ "ref"))
  
  alt_reads <- read_mapping %>%
    dplyr::select(all_of(c(select_col_name, "read_id",  "read_type"))) %>%
    distinct() %>%
    inner_join(fread(T_modkit_extract) %>%
                filter(mod_code == "m") %>%
                dplyr::select(read_id, forward_read_position, ref_position,chrom, mod_strand, ref_strand, read_length, mod_qual, base_qual), 
              by = join_by(chrom, read_id),
              relationship = "many-to-many") %>%
    mutate(mod = case_when(mod_qual >= 0.6 ~ "Me",
                           mod_qual <= 0.4 ~ "Un",
                           T ~ NA)) %>%
    filter(base_qual >=20, !is.na(mod),  ref_position != -1, (abs(ref_position -start )<= 5000 | abs(ref_position -end )<= 5000)) %>%
    group_by(ID, read_type) %>%
    filter(n()>2) %>%
    ungroup() %>%
    group_by(ID, read_type) %>%
    filter(n_distinct(read_id)>2) %>%
    ungroup() %>%
    group_by(ID) %>%
    filter(n_distinct(read_type)== 2) %>%
    ungroup()
  
  alt_P <- alt_reads%>%
    dplyr::select(all_of(c(select_col_name ,"ref_position"))) %>%
    distinct() %>%
    left_join(fread(N_modkit_extract) %>%
                filter(mod_code == "m") %>%
                dplyr::select(read_id, forward_read_position, ref_position,chrom, mod_strand, ref_strand, read_length, mod_qual, base_qual),
              by = join_by(chrom, ref_position),
              relationship = "many-to-many") %>%
    mutate(read_type = "parent") %>%
    mutate(mod = case_when(mod_qual >= 0.6 ~ "Me",
                           mod_qual <= 0.4 ~ "Un",
                           T ~ NA)) %>%
    filter(base_qual >=20, !is.na(mod),  ref_position != -1, (abs(ref_position -start )<= 5000 | abs(ref_position -end )<= 5000)) %>%
    group_by(ID, read_type) %>%
    filter(n()>2) %>%
    ungroup() %>%
    group_by(ID, read_type) %>%
    filter(n_distinct(read_id)>2) %>%
    ungroup() 
  
  total_reads <- rbind(alt_reads, alt_P) %>%
    mutate(read_type = factor(read_type, levels = c("alt", "ref", "parent"))) 
  
  return(total_reads)
}


create_windows <- function(ID, chrom, SV_start, SV_end, SV_type) {
  # Left-side windows
  left_starts <- seq(from = SV_start - num_windows * window_size, 
                     to = SV_start - window_size, 
                     by = window_size)
  left_ends <- left_starts + window_size - 1
  
  left_starts[left_starts <= 0] <- NA
  left_ends[left_starts <= 0] <- NA
  
  # Right-side windows
  right_starts <- seq(from = SV_end + 1, 
                      to = SV_end + num_windows * window_size, 
                      by = window_size)
  right_ends <- right_starts + window_size - 1
  
  # Combine and return as data.frame
  data.frame(
    variant_id = ID,
    chrom = chrom,
    SV_start = SV_start,
    SV_end = SV_end,
    SV_type = SV_type,
    window_start = c(left_starts, right_starts),
    window_end = c(left_ends, right_ends)
  )
}

keep_bigger_overlaps <- function(gr1, gr2) {
  # Find overlaps between gr1 and gr2
  hits <- findOverlaps(gr1, gr2)
  
  # Extract indices of overlapping ranges
  gr1_idx <- queryHits(hits)
  gr2_idx <- subjectHits(hits)
  
  # Get widths of overlapping ranges
  gr1_widths <- width(gr1[gr1_idx])
  gr2_widths <- width(gr2[gr2_idx])
  
  # Decide which ranges to keep
  keep_gr1 <- gr1_idx[gr1_widths >= gr2_widths]  # keep gr1 if width equal or larger
  keep_gr2 <- gr2_idx[gr2_widths > gr1_widths]   # keep gr2 if strictly larger
  
  # Unique indices
  keep_gr1 <- unique(keep_gr1)
  keep_gr2 <- unique(keep_gr2)
  
  # Subset GRanges
  gr1_keep <- gr1[keep_gr1]
  gr2_keep <- gr2[keep_gr2]
  
  # Combine and return
  c(gr1_keep, gr2_keep)
}