## Script for Tx421 LUSC paper
#Fig 3b, Ext 4c

library(tidyverse)
library(data.table)

input_plot_df_all <- read_csv("../df-alvin/500_entire_tumour_with_grading_with_wGIITMB_info.csv")


input_plot_df_all %>% dplyr::select(keratinization, tumour_budding) %>% table
# tumour_budding
# keratinization have_budding no_budding
# have_keratin           53         22
# no_keratin             30         31

input_plot_df_all %>% dplyr::select(keratinization, single_cell_invasion) %>% table
# single_cell_invasion
# keratinization have_sci no_sci
# have_keratin       32     43
# no_keratin         18     43


#invert x,y 
table(input_plot_df_all$tumour_budding, input_plot_df_all$keratinization) -> table1

df_gathered1 <- data.frame(var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2] ),
                           var2 = c(colnames(table1)[1],colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
                           count = c(table1[1], table1[3], table1[2], table1[4]))
df_gathered1$var2 <- factor( df_gathered1$var2 , levels=c("have_keratin", "no_keratin"))

df_gathered1$var1 <- factor( df_gathered1$var1 , levels=c("no_budding", "have_budding"))
fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

color_list <-  c("no_budding" = "#A6CEE3",
                 "greencolour" = "#B2DF8A",
                 "have_budding"= "#FB9A99")
g1 <- ggplot(data=df_gathered1, aes(x=var2, y=count, fill=var1, label = count)) +
  geom_bar(stat="identity", position = "fill", color = "black") +
  geom_bar_text( position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F)+
  theme_cowplot() +
  ylab("proportion") +
  # labs(x = var2, y = "proportion", fill = var1)   +
  ggtitle(paste0( " Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall=4))) +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=18))
print(g1)
#Fig 3b


table(input_plot_df_all$single_cell_invasion, input_plot_df_all$keratinization) -> table1

df_gathered1 <- data.frame(var1 = c(rownames(table1)[1], rownames(table1)[1], rownames(table1)[2], rownames(table1)[2] ),
                           var2 = c(colnames(table1)[1],colnames(table1)[2], colnames(table1)[1], colnames(table1)[2]),
                           count = c(table1[1], table1[3], table1[2], table1[4]))
df_gathered1$var2 <- factor( df_gathered1$var2 , levels=c("have_keratin", "no_keratin"))

df_gathered1$var1 <- factor( df_gathered1$var1 , levels=c("no_sci", "have_sci"))

fisher_p <- fisher.test(matrix(c(table1[1], table1[2], table1[3], table1[4]), nrow = 2))

color_list <-  c("no_sci" = "#A6CEE3",
                 "greencolour" = "#B2DF8A",
                 "have_sci"= "#FB9A99")
g1 <- ggplot(data=df_gathered1, aes(x=var2, y=count, fill=var1, label = count)) +
  geom_bar(stat="identity", position = "fill", color = "black") +
  geom_bar_text( position = "fill", reflow = TRUE) + 
  scale_fill_manual(values = color_list, na.value = '#bababa', drop = F)+
  theme_cowplot() +
  ylab("proportion") +
  # labs(x = var2, y = "proportion", fill = var1)   +
  ggtitle(paste0( " Fisher-exact test p-value: ", format(round(as.numeric(fisher_p$p.value), 4), nsmall=4))) +
  theme(axis.title.x = element_blank(),  
        axis.title.y = element_text(size=24),
        axis.text.x = element_text(size=20, angle = 45,  hjust=1), 
        axis.text.y = element_text(size=16),
        legend.position = "right",
        legend.title = element_blank(),
        legend.text = element_text(size=18))
print(g1)
#Ext fig 4c