#Script for Tx421 LUSC paper
#Fig 2a,3a, Ext fig 4b

library(tidyverse)
library(ggpubr)
library(reshape2)


alltumour_dataset <- read.csv("../df-alvin/All_regions_500_regions-SCC_reviewed_With_mean_keratinMitosis.csv")
scc_dataset <- alltumour_dataset %>% filter(!is.na(Histology_per_tumour_id_muttable)) %>%
  filter(Histology_per_tumour_id_muttable == "Squamous cell carcinoma") %>%
  filter(included_500_regions == 1)
entiretumour_dataset <- read.csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv")


#keratinization Fig 2a
scc_dataset[, c("tumour_id_muttable_cruk", "region", "keratinization")] %>% 
  filter(!is.na(keratinization)) ->
  scc_selectforplot
#scc_dataset has 136 patients
#scc_selectforplot has 126 patients, because those regions without HE slides available for assessment are excluded.

scc_selectforplot %>% group_by(tumour_id_muttable_cruk, keratinization) %>% 
  summarize(no_of_regions = n()) ->
  scc_plot

scc_selectforplot %>% group_by(tumour_id_muttable_cruk) %>% 
  summarize(total_no_of_regions = n()) ->
  scc_regionnumberpercase

left_join(scc_plot, scc_regionnumberpercase, by = c("tumour_id_muttable_cruk" = "tumour_id_muttable_cruk")) ->
  join_scc

join_scc$proportion_regions <- join_scc$no_of_regions/join_scc$total_no_of_regions

join_scc$keratinization <- factor(join_scc$keratinization, levels = c("no_keratin", "have_keratin"))

join_scc %>% filter(!is.na(keratinization)) %>% 
  arrange(keratinization, proportion_regions) %>% 
  ungroup()->   #without ungroup row_number() will not work
  scc_plot1


#index number according to arrangement of poor_dx proportion
scc_plot1 %>% mutate(indexnumber = ifelse(keratinization == "no_keratin", row_number(), NA)) ->
  scc_plot2

scc_plot2 %>% dplyr::filter(!is.na(indexnumber)) %>% 
  dplyr::select(tumour_id_muttable_cruk, indexnumber) %>% 
  left_join(scc_plot2, ., by = c("tumour_id_muttable_cruk" = "tumour_id_muttable_cruk")) %>%
  arrange(indexnumber.y) ->
  scc_plot3

scc_plot3 %>%
  mutate(indexnumber.z = ifelse(is.na(indexnumber.y), row_number() - min(which(is.na(scc_plot3$indexnumber.y))) + 1 , indexnumber.y + length(scc_plot3$indexnumber.y[is.na(scc_plot3$indexnumber.y)]))) -> 
  scc_plot4

saveRDS(scc_plot4, "../df-alvin/tumour_regionalkeratinization.RDS") #this is saved for use in another PCA script. 

total_no_of_cases <- as.numeric(length(unique(scc_plot4$tumour_id_muttable_cruk)))

ggplot(scc_plot4, aes(fill=factor(keratinization, levels= c("have_keratin", "no_keratin")) , y=no_of_regions, x=indexnumber.z)) + 
  geom_bar(position="fill", stat="identity", width = 1) +
  ggtitle(paste0("Regional annotations of LUSC", ", keratinization, n = ", nrow(scc_regionnumberpercase))) +
  scale_fill_manual(values = c(  "have_keratin" = "#a6cee2", 
                                 "no_keratin" = "#f6999a"), drop = T) +
  scale_x_continuous(breaks = c(0.25*total_no_of_cases, 
                                0.5*total_no_of_cases,
                                0.75*total_no_of_cases,
                                total_no_of_cases),
                     labels = c("25%", "50%", "75%", "100%"),
                     name = "Percentage of tumours") +
  scale_y_continuous(name = "Proportion of tumour regions") +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_text(size = 13, family = "ArialMT"),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Fig 2a


#also plot classification by entire tumour
scc_plot5 <- scc_plot4 %>% rename(keratinization_regional = keratinization) %>% 
  left_join(entiretumour_dataset[, c("cruk_tumourid", "keratinization") ], by = c("tumour_id_muttable_cruk" = "cruk_tumourid"))

ggplot(scc_plot5, aes(fill=factor(keratinization, levels= c("have_keratin", "no_keratin")) , y=1, x=indexnumber.z)) + 
  geom_bar(position="fill", stat="identity", width = 1) +
  ggtitle(paste0("Regional annotations of LUSC", ", keratinization, n = ", nrow(scc_regionnumberpercase))) +
  scale_fill_manual(values = c(  "have_keratin" = "#a6cee2", 
                                 "no_keratin" = "#f6999a"), drop = T) +
  scale_x_continuous(breaks = c(0.25*total_no_of_cases, 
                                0.5*total_no_of_cases,
                                0.75*total_no_of_cases,
                                total_no_of_cases),
                     labels = c("25%", "50%", "75%", "100%"),
                     name = "Percentage of tumours") +
  scale_y_continuous(name = "Proportion of tumour regions") +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_text(size = 13, family = "ArialMT"),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Fig 2a top panel


#Fig 3a
scc_dataset[, c("tumour_id_muttable_cruk", "region", "tumour_budding")] %>% 
  filter(!is.na(tumour_budding)) ->
  scc_selectforplot

scc_selectforplot %>% group_by(tumour_id_muttable_cruk, tumour_budding) %>% 
  summarize(no_of_regions = n()) ->
  scc_plot

scc_selectforplot %>% group_by(tumour_id_muttable_cruk) %>% 
  summarize(total_no_of_regions = n()) ->
  scc_regionnumberpercase

left_join(scc_plot, scc_regionnumberpercase, by = c("tumour_id_muttable_cruk" = "tumour_id_muttable_cruk")) ->
  join_scc

join_scc$proportion_regions <- join_scc$no_of_regions/join_scc$total_no_of_regions

join_scc$tumour_budding <- factor(join_scc$tumour_budding, levels = c("have_budding", "no_budding"))

join_scc %>% filter(!is.na(tumour_budding)) %>% 
  arrange(tumour_budding, proportion_regions) %>% 
  ungroup()->   #without ungroup row_number() will not work
  scc_plot1


#index number according to arrangement of poor_dx proportion
scc_plot1 %>% mutate(indexnumber = ifelse(tumour_budding == "have_budding", row_number(), NA)) ->
  scc_plot2

scc_plot2 %>% dplyr::filter(!is.na(indexnumber)) %>% 
  dplyr::select(tumour_id_muttable_cruk, indexnumber) %>% 
  left_join(scc_plot2, ., by = c("tumour_id_muttable_cruk" = "tumour_id_muttable_cruk")) %>%
  arrange(indexnumber.y) ->
  scc_plot3

scc_plot3 %>%
  mutate(indexnumber.z = ifelse(is.na(indexnumber.y), row_number() - min(which(is.na(scc_plot3$indexnumber.y))) + 1 , indexnumber.y + length(scc_plot3$indexnumber.y[is.na(scc_plot3$indexnumber.y)]))) -> 
  scc_plot4

total_no_of_cases <- as.numeric(length(unique(scc_plot4$tumour_id_muttable_cruk)))

ggplot(scc_plot4, aes(fill=factor(tumour_budding, levels= c("no_budding", "have_budding")) , y=no_of_regions, x=indexnumber.z)) + 
  geom_bar(position="fill", stat="identity", width = 1) +
  ggtitle(paste0("Regional annotations of LUSC", ", tumour_budding, n = ", nrow(scc_regionnumberpercase))) +
  scale_fill_manual(values = c(  "no_budding" = "#fbb040", 
                                 "have_budding" = "#754c29"), drop = T) +
  scale_x_continuous(breaks = c(0.25*total_no_of_cases, 
                                0.5*total_no_of_cases,
                                0.75*total_no_of_cases,
                                total_no_of_cases),
                     labels = c("25%", "50%", "75%", "100%"),
                     name = "Percentage of tumours") +
  scale_y_continuous(name = "Proportion of tumour regions") +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_text(size = 13, family = "ArialMT"),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Fig 3a

#also plot classification by entire tumour
scc_plot5 <- scc_plot4 %>% rename(tumour_budding_regional = tumour_budding) %>% 
  left_join(entiretumour_dataset[, c("cruk_tumourid", "tumour_budding") ], by = c("tumour_id_muttable_cruk" = "cruk_tumourid"))

ggplot(scc_plot5, aes(fill=factor(tumour_budding, levels= c("no_budding", "have_budding")) , y=1, x=indexnumber.z)) + 
  geom_bar(position="fill", stat="identity", width = 1) +
  ggtitle(paste0("Regional annotations of LUSC", ", keratinization, n = ", nrow(scc_regionnumberpercase))) +
  scale_fill_manual(values = c(  "no_budding" = "#fbb040", 
                                 "have_budding" = "#754c29"), drop = T) +
  scale_x_continuous(breaks = c(0.25*total_no_of_cases, 
                                0.5*total_no_of_cases,
                                0.75*total_no_of_cases,
                                total_no_of_cases),
                     labels = c("25%", "50%", "75%", "100%"),
                     name = "Percentage of tumours") +
  scale_y_continuous(name = "Proportion of tumour regions") +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_text(size = 13, family = "ArialMT"),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Fig 3a top panel

#Ext fig 4b
scc_dataset[, c("tumour_id_muttable_cruk", "region", "single_cell_invasion")] %>% 
  filter(!is.na(single_cell_invasion)) ->
  scc_selectforplot

scc_selectforplot %>% group_by(tumour_id_muttable_cruk, single_cell_invasion) %>% 
  summarize(no_of_regions = n()) ->
  scc_plot

scc_selectforplot %>% group_by(tumour_id_muttable_cruk) %>% 
  summarize(total_no_of_regions = n()) ->
  scc_regionnumberpercase

left_join(scc_plot, scc_regionnumberpercase, by = c("tumour_id_muttable_cruk" = "tumour_id_muttable_cruk")) ->
  join_scc

join_scc$proportion_regions <- join_scc$no_of_regions/join_scc$total_no_of_regions

join_scc$single_cell_invasion <- factor(join_scc$single_cell_invasion, levels = c("have_sci", "no_sci"))

join_scc %>% filter(!is.na(single_cell_invasion)) %>% 
  arrange(single_cell_invasion, proportion_regions) %>% 
  ungroup()->   #without ungroup row_number() will not work
  scc_plot1


#index number according to arrangement of poor_dx proportion
scc_plot1 %>% mutate(indexnumber = ifelse(single_cell_invasion == "have_sci", row_number(), NA)) ->
  scc_plot2

scc_plot2 %>% dplyr::filter(!is.na(indexnumber)) %>% 
  dplyr::select(tumour_id_muttable_cruk, indexnumber) %>% 
  left_join(scc_plot2, ., by = c("tumour_id_muttable_cruk" = "tumour_id_muttable_cruk")) %>%
  arrange(indexnumber.y) ->
  scc_plot3

scc_plot3 %>%
  mutate(indexnumber.z = ifelse(is.na(indexnumber.y), row_number() - min(which(is.na(scc_plot3$indexnumber.y))) + 1 , indexnumber.y + length(scc_plot3$indexnumber.y[is.na(scc_plot3$indexnumber.y)]))) -> 
  scc_plot4

total_no_of_cases <- as.numeric(length(unique(scc_plot4$tumour_id_muttable_cruk)))

ggplot(scc_plot4, aes(fill=factor(single_cell_invasion, levels= c("no_sci", "have_sci")) , y=no_of_regions, x=indexnumber.z)) + 
  geom_bar(position="fill", stat="identity", width = 1) +
  ggtitle(paste0("Regional annotations of LUSC", ", single_cell_invasion, n = ", nrow(scc_regionnumberpercase))) +
  scale_fill_manual(values = c(  "no_sci" = "#ffde17", 
                                 "have_sci" = "#7f3f98"), drop = T) +
  scale_x_continuous(breaks = c(0.25*total_no_of_cases, 
                                0.5*total_no_of_cases,
                                0.75*total_no_of_cases,
                                total_no_of_cases),
                     labels = c("25%", "50%", "75%", "100%"),
                     name = "Percentage of tumours") +
  scale_y_continuous(name = "Proportion of tumour regions") +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_text(size = 13, family = "ArialMT"),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Ext fig 4b

#also plot classification by entire tumour
scc_plot5 <- scc_plot4 %>% rename(single_cell_invasion_regional = single_cell_invasion) %>% 
  left_join(entiretumour_dataset[, c("cruk_tumourid", "single_cell_invasion") ], by = c("tumour_id_muttable_cruk" = "cruk_tumourid"))

ggplot(scc_plot5, aes(fill=factor(single_cell_invasion, levels= c("no_sci", "have_sci")) , y=1, x=indexnumber.z)) + 
  geom_bar(position="fill", stat="identity", width = 1) +
  ggtitle(paste0("Regional annotations of LUSC", ", keratinization, n = ", nrow(scc_regionnumberpercase))) +
  scale_fill_manual(values = c( "no_sci" = "#ffde17", 
                                "have_sci" = "#7f3f98"), drop = T) +
  scale_x_continuous(breaks = c(0.25*total_no_of_cases, 
                                0.5*total_no_of_cases,
                                0.75*total_no_of_cases,
                                total_no_of_cases),
                     labels = c("25%", "50%", "75%", "100%"),
                     name = "Percentage of tumours") +
  scale_y_continuous(name = "Proportion of tumour regions") +
  theme(axis.text.x = element_text(size = 13, family = "ArialMT"),
        axis.title.x = element_text(size = 13, family = "ArialMT"),
        axis.text.y = element_text(size = 13, family = "ArialMT"),
        axis.title.y = element_text(size = 13, family = "ArialMT"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 13, family = "ArialMT"),
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title = element_blank())
#Ext fig 4b
