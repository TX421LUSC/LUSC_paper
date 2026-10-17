## Script for Tx421 LUSC paper
#Fig 1d

#################################
### libraries and directories ###
#################################

library(tidyverse)
library(data.table)
library(fst)
library(GenomicRanges)
library(cowplot)
library(ggnewscale)


##########################
### path to input data ###
##########################

# cytoband info
cytoband.path <- "../data/hg19_cytobands.txt"
cytoband_df_indices_path <- "../data/cytoband_df_indices.RData"

# scna data
mphase_sample_path <- "../data/20221109_TRACERx421_scna_table.fst"
mphase_clonality_tumour_path <- "../data/20221108_mphase_clonality_tumour_cyto_cruk.RDS"

######################################################################
########LUSC input#############

squamous_data <- read.csv("../df-alvin/500_entire_tumour_with_grading.csv")

squamous_data %>% mutate(subtype_grade = case_when(keratinization %in% c("have_keratin") ~ "kerat_group",
                                                   keratinization %in% c("no_keratin") ~ "nonkerat_group")) %>% 
  mutate(tumour_id_muttable_cruk = cruk_tumourid) ->
  squamous_data1

nonkerat_group_vec <- unique(squamous_data1$tumour_id_muttable_cruk[squamous_data1$subtype_grade == "nonkerat_group"])
kerat_group_vec <- unique(squamous_data1$tumour_id_muttable_cruk[squamous_data1$subtype_grade == "kerat_group"])

########LUSC input
######################################################################

########################
### helper functions ###
########################

make_cohort_df_for_across_genome_plot_helper.new = function(cn_patient_df_dc, total_cohort_size, 
                                                            has_any_cn_of_choice = "loss.LOH_clonality") {
  print("Called make_cohort_df_for_across_genome_plot_helper")
  cn_patient_df_dc$has_any_cn_of_choice <- cn_patient_df_dc[[has_any_cn_of_choice]]
  cn_patient_df_dc_cn_of_choice_df = cn_patient_df_dc %>%
    filter(!is.na(has_any_cn_of_choice)) %>%
    mutate(CN_Subclonal=case_when(has_any_cn_of_choice == "clonal" ~ "Clonal",
                                  has_any_cn_of_choice == "subclonal" ~ "Subclonal")) %>%
    mutate(CN_Subclonal=factor(CN_Subclonal,levels=c("Clonal","Subclonal"))) %>%
    dplyr::select(tumour_id,cytoband_name,CN_Subclonal) %>%
    distinct() %>%
    group_by(cytoband_name,CN_Subclonal) %>%
    tally() %>%
    ungroup() %>%
    mutate(n=n/total_cohort_size) %>%
    mutate(cytoband_name=factor(cytoband_name,levels=c(cytoband_df$cytoband_name))) %>%
    filter(!is.na(cytoband_name)) %>%
    complete(cytoband_name,nesting(CN_Subclonal)) %>%
    mutate(n=ifelse(is.na(n),0,n)) %>%
    tidyr::spread(CN_Subclonal,n,fill=0) %>%
    mutate(Total_CN=Clonal+Subclonal,
           Clonal_Color="X")
  
  cn_patient_df_dc_cn_of_choice_df
  
}

createSimpleAcrossGenomePlot_new = function(cn_patient_df_dc,
                                            events_included="ArmAndFocal",
                                            LOHAsLoss=TRUE,
                                            profile_only=FALSE,
                                            fig_dir="",
                                            plot_title="",
                                            clonal_gain_colour ="grey50", subclonal_gain_colour="#f1a340",
                                            clonal_loss_colour ="grey50", subclonal_loss_colour="#f1a340", 
                                            total_cn_line_colour = "black",
                                            total_cn_line_type ="solid"){
  
  total_cohort_size = n_distinct(cn_patient_df_dc$tumour_id)
  print(total_cohort_size)
  
  if(LOHAsLoss==TRUE){
    events_included = paste(events_included,"LOHasLoss",sep="_")
  }
  
  
  ### GAINS PLOT
  cn_patient_df_dc_gains_df = make_cohort_df_for_across_genome_plot_helper.new(cn_patient_df_dc, total_cohort_size = total_cohort_size,
                                                                               has_any_cn_of_choice = "gain.amp_clonality")
  
  cn_patient_df_dc_gains = make_across_clonal_subclonal_genome_plot_itself_helper(cn_patient_df_dc_df = cn_patient_df_dc_gains_df,
                                                                                  gain_or_loss = "gain",
                                                                                  #LOHAsLoss = LOHAsLoss, 
                                                                                  total_cohort_size = total_cohort_size,
                                                                                  plot_title = plot_title,
                                                                                  clonal_gain_colour =clonal_gain_colour,
                                                                                  subclonal_gain_colour=subclonal_gain_colour,
                                                                                  clonal_loss_colour =clonal_loss_colour,
                                                                                  subclonal_loss_colour=subclonal_loss_colour,
                                                                                  total_cn_line_colour = "black",
                                                                                  total_cn_line_type ="solid")
  ### LOSS PLOT
  cn_patient_df_dc_losses_df = make_cohort_df_for_across_genome_plot_helper.new(cn_patient_df_dc, total_cohort_size = total_cohort_size,
                                                                                has_any_cn_of_choice = "loss.LOH_clonality")
  cn_patient_df_dc_losses = make_across_clonal_subclonal_genome_plot_itself_helper(cn_patient_df_dc_df = cn_patient_df_dc_losses_df,
                                                                                   gain_or_loss = "loss",
                                                                                   LOHAsLoss = LOHAsLoss, 
                                                                                   total_cohort_size = total_cohort_size, 
                                                                                   clonal_gain_colour =clonal_gain_colour,
                                                                                   subclonal_gain_colour=subclonal_gain_colour,
                                                                                   clonal_loss_colour =clonal_loss_colour,
                                                                                   subclonal_loss_colour=subclonal_loss_colour,
                                                                                   total_cn_line_colour = "black",
                                                                                   total_cn_line_type ="solid")
  
  cn_cohort_summary_df = full_join(cn_patient_df_dc_losses_df,cn_patient_df_dc_gains_df,by="cytoband_name",suffix=c("_loss","_gain"))
#  write.table(cn_cohort_summary_df,file=paste0(fig_dir,gsub(" ","",paste(plot_title,collapse="_")),"_cn_cohort_summary_df.txt"),sep="\t",quote=F,row.names=F)
  
  rel_heights_vec <-c(1,1)
  
  cn_hist_plots = plot_grid(cn_patient_df_dc_gains,
                            cn_patient_df_dc_losses,
                            align="v",
                            rel_heights=rel_heights_vec,
                            ncol=1)
  
  cn_hist_plots
  
}#end 


make_across_clonal_subclonal_genome_plot_itself_helper = function(cn_patient_df_dc_df,  gain_or_loss, LOHAsLoss=TRUE, 
                                                                  total_cohort_size, 
                                                                  plot_title="",
                                                                  clonal_gain_colour ="grey50", subclonal_gain_colour="#f1a340",
                                                                  clonal_loss_colour ="grey50", subclonal_loss_colour="#f1a340", 
                                                                  total_cn_line_colour = "black",
                                                                  total_cn_line_type ="solid") {
  print("Called make_across_clonal_subclonal_genome_plot_itself_helper")
  
  if(gain_or_loss %in% "gain") {
    clonal_colour <-  clonal_gain_colour
    subclonal_colour <- subclonal_gain_colour
  } else if (gain_or_loss %in% "loss") {
    clonal_colour <-  clonal_loss_colour
    subclonal_colour <- subclonal_loss_colour
  }
  
  cn_patient_df_dc_df_plot = cn_patient_df_dc_df %>%
    ggplot(aes(x=cytoband_name)) +
    geom_vline(xintercept = as.numeric(as.character(cytoband_df_indices$end_index)),linetype = "dashed", colour="grey") +
    geom_bar(aes(y=Clonal,fill=clonal_colour),stat="identity",alpha=0.7,width=1) +
    scale_fill_manual(values=c(clonal_colour),labels=c("Truncal")) +
    new_scale("fill") +
    geom_bar(aes(y=Subclonal,fill=subclonal_colour),stat="identity",alpha=0.7,width=1) +
    scale_fill_manual(values=c(subclonal_colour),labels=c("Subclonal")) +
    geom_path(aes(y=Clonal),color=clonal_colour,group="Clonal") +
    geom_path(aes(y=Subclonal),color=subclonal_colour,group="Subclonal") +
    geom_path(aes(y=Total_CN),color=total_cn_line_colour,group="Total_CN", linetype=total_cn_line_type) +
    theme_cowplot() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.line.x = element_blank(),
      legend.position = "right",
      legend.title = element_blank())
  
  if(gain_or_loss %in% "gain") {
    
    cn_patient_df_dc_df_plot_final = cn_patient_df_dc_df_plot +
      scale_x_discrete(drop=FALSE) +
      annotate("text", x=as.numeric(as.character(cytoband_df_indices$mean_index)),
               y = rep(c(1,0.9),nrow(cytoband_df_indices)/2), label = cytoband_df_indices$chromosome) +
      ylab("Prop. of Patients\nWith Gain/Amp") +
      ggtitle(gsub("\n", "", paste0(plot_title, " ",paste0(" - n=", total_cohort_size)))) +
      ylim(0,1)
    
  } else if (gain_or_loss %in% "loss") {
    
    cn_patient_df_dc_df_plot_final = cn_patient_df_dc_df_plot +
      scale_x_discrete(drop=FALSE) +
      ylab(ifelse(LOHAsLoss==TRUE,"Prop. of Patients\nWith Loss/LOH","Prop. of Patients\nWith Loss")) +
      scale_y_reverse(limits=c(1,0)) +
      ggtitle("")
    
  }
  
  cn_patient_df_dc_df_plot_final
  
}#end 


#################
### load data ###
#################

cytoband_df  <- read.table(cytoband.path, stringsAsFactors=FALSE, sep="\t", header=TRUE)

colnames(cytoband_df)  <-  c("chr","start","end","cytoband","X")
cytoband_df$cytoband_chr  <-  gsub("chr","",cytoband_df$chr)
cytoband_df  <-  cytoband_df[order(as.numeric(as.character(cytoband_df$cytoband_chr)),as.numeric(as.character(cytoband_df$start))),]
genome_size  <-  sum(as.numeric(cytoband_df$end-cytoband_df$start))
cytoband_df  <-  cytoband_df %>%
  mutate(cytoband_name = paste(cytoband_chr,cytoband,sep="")) %>%
  filter(!cytoband_chr%in%c("X","Y")) %>%
  select(cytoband_chr,start,end,cytoband_name)
cytoband_df$cytoband_arm = gsub(" ","",substr(gsub("q","q ",gsub("p","p ",cytoband_df$cytoband_name)),1,3))

load(cytoband_df_indices_path)

glimpse(cytoband_df_indices)

# load clonal/subclonal updated mphase data 
mphase_clonality_tumour_cyto_cruk <- readRDS(mphase_clonality_tumour_path)
mphase_tumours <- unique(mphase_clonality_tumour_cyto_cruk$tumour_id_mphase_cruk) # 403

### load SCNA per sample  ###
mphase_sample_df <- read.fst(mphase_sample_path, as.data.table=TRUE)

# number of tumours (including CN failed cases)
length(nonkerat_group_vec)

length(kerat_group_vec)


# keep only CN passed tumours 
nonkerat_group_vec <- intersect(mphase_tumours, nonkerat_group_vec )
kerat_group_vec <- intersect(mphase_tumours, kerat_group_vec )

# number of tumours (not including CN failed cases)
length(nonkerat_group_vec)

length(kerat_group_vec)


add_samp_vec_list <- list(nonkerat_group_vec,kerat_group_vec )
names(add_samp_vec_list)  <- c("nonkerat_group_vec","kerat_group_vec")

################################### 
## across genome plot (unpaired) ##
###################################

# tumours with multi-region scna data
multi_region_tumours <- mphase_sample_df %>%
  filter(!grepl('LN', sample)) %>%
  distinct(tumour_id, sample) %>%
  group_by(tumour_id) %>%
  summarise(count = n()) %>%
  filter(count >= 2)


# load existing patient summary of mphase with cytoband (max overlapped cytoband already selected)

all_prim_summaries.ori <- mphase_clonality_tumour_cyto_cruk %>%
  dplyr::rename(tumour_id = tumour_id_mphase_cruk ) %>%
  filter(tumour_id %in%  multi_region_tumours$tumour_id)


group_list <- names(add_samp_vec_list)  
group_list 
# [1] "nonkerat_group_vec"    "kerat_group_vec"


### Fig 1d Non-keratinized LUSC ###
subtype <- group_list[1]
  
  LUSC_tum_ids_vec_subtype <- unique(add_samp_vec_list[[subtype]])
  
  all_LUSC_prim_summaries <- all_prim_summaries.ori[all_prim_summaries.ori$tumour_id %in% LUSC_tum_ids_vec_subtype,] 
  
  all_LUSC_prim_summaries_cyto_df_list <-   all_LUSC_prim_summaries
  
  num_tum <- length(unique(all_LUSC_prim_summaries_cyto_df_list$tumour_id))
  print(paste("num_tum =", num_tum) )
  
  all_LUSC_plot_list <- createSimpleAcrossGenomePlot_new(cn_patient_df_dc = all_LUSC_prim_summaries_cyto_df_list,
                                                         events_included="ArmAndFocal",
                                                         LOHAsLoss=TRUE,
                                                         profile_only=FALSE,
                                                         plot_title=paste(subtype,"n=",num_tum),
                                                         fig_dir="")

  print(plot_grid(plotlist = list(all_LUSC_plot_list), align="v",ncol=1))
  
  
  
### Fig 1d Keratinized LUSC ###
subtype <- group_list[2]
  
  LUSC_tum_ids_vec_subtype <- unique(add_samp_vec_list[[subtype]])
  
  all_LUSC_prim_summaries <- all_prim_summaries.ori[all_prim_summaries.ori$tumour_id %in% LUSC_tum_ids_vec_subtype,] 
  
  all_LUSC_prim_summaries_cyto_df_list <-   all_LUSC_prim_summaries
  
  num_tum <- length(unique(all_LUSC_prim_summaries_cyto_df_list$tumour_id))
  print(paste("num_tum =", num_tum) )
  
  all_LUSC_plot_list <- createSimpleAcrossGenomePlot_new(cn_patient_df_dc = all_LUSC_prim_summaries_cyto_df_list,
                                                         events_included="ArmAndFocal",
                                                         LOHAsLoss=TRUE,
                                                         profile_only=FALSE,
                                                         plot_title=paste(subtype,"n=",num_tum),
                                                         fig_dir="")
  
  print(plot_grid(plotlist = list(all_LUSC_plot_list), align="v",ncol=1))
  