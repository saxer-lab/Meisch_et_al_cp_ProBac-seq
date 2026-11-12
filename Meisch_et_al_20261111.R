# This script is for plotting the data for Meisch et al. 
# Written by Gerda Saxer 2026
# 

# Loading libraries
library(ggpubr)
library(rstatix)
library(tidyverse)
library(readxl)
library(broom)
library(patchwork)
library(cowplot)
library(rstatix)
library(purrr)
library(dbplyr)


######################################################################################################
####### Metabolites #############
######################################################################################################


############ FIGURE 1 ############
sup <- read_xlsx('Meisch_et_al_Data_Tables.xlsx', sheet = "Organic Acids", na = c("", "NA", "N/A") )
glimpse(sup)

sup %>%
  ggplot(aes(Compound, Value, fill=Sample))+
  stat_summary(fun='mean', geom = 'bar', position='dodge', width = 0.6)+
  stat_summary(geom = 'errorbar', fun.data=mean_cl_boot, position = position_dodge(.6), size=.3, width = .5)+
  theme_pubr(base_size = 10)+
  theme(legend.position = "right" )+
  xlab("Compound")+
  ylab("Concentration (g/L)")+
  scale_x_discrete(labels=c("AG01" = "L. acidophilus", "AG02"="B. lactis"))




######################################################################################################
####### Co-culture #############
######################################################################################################


########   FIGURE 2 ######################
coculture <- read_xlsx('Meisch_et_al_Data_Tables.xlsx', sheet = "Co-Culture ", na = c("", "NA", "N/A") )
glimpse(coculture)


culture.labs <- c('BC'='Ba+Cp','B' = 'Ba',  'C'='Cp', 'L'='La', 'LB'='La+Ba', 'LBC'='All three', 'LC'='La+Cp','CP' = 'C. perfringens CP 25037', 'DFM'='L. acidophilus AG01')

options(repr.plot.width =9, repr.plot.height =3.5) 
cocu <- coculture %>%
  ggplot(aes(fill=Fill, x = as.factor(Time), y = logCFU)) +
  facet_grid(~factor(Culture, levels = c('L', 'B', 'LB', 'C', 'LC', 'BC', 'LBC')), labeller = as_labeller(culture.labs))+
  stat_summary(fun = "mean", geom ="bar", position="dodge")+
  stat_summary(geom = "errorbar", fun.data=mean_cl_boot, position = position_dodge(.8), width=.5)+
  scale_y_continuous(limits = c(0, 10), breaks = c(0, 2, 4, 6, 8, 10))+
  geom_hline(yintercept = 6, linetype="dotted")+
  ylab("log10 (CFU/ml)") + 
  xlab("")+
  theme_pubr(base_size = 10)+
  theme(legend.position = "none")
#cocu

# ANOVA Viability among Cp containing cultures at T = 4
coculture %>% 
  filter(Time == 4) %>% 
  filter(Media =="CP") %>% 
  filter(Culture %in% c("C", "BC", "LC", "LBC")) %>% 
  aov(logCFU ~ Culture, .) %>% summary.aov

# ANOVA Viability among Cp containing cultures at T = 24
coculture %>% 
  filter(Time == 24) %>% 
  filter(Media =="CP") %>% 
  filter(Culture %in% c("C", "BC", "LC", "LBC")) %>% 
  aov(logCFU ~ Culture, .) %>% summary.aov

# Tukey HSD for Viability among Cp containing cultures at T = 24
coculture %>%   
  filter(Time == 24) %>% 
  filter(Media =="CP") %>% 
  filter(Culture %in% c("C", "BC", "LC", "LBC")) %>% 
  tukey_hsd(logCFU ~ Culture) %>% view()


## pH of co-cultures

ph <- coculture %>%
  ggplot(aes(x = as.factor(Time), y = pH)) +
  facet_grid(~factor(Culture, levels = c('L', 'B', 'LB', 'C', 'LC', 'BC', 'LBC')), labeller = as_labeller(culture.labs))+
  stat_summary(fun = "mean", geom ="point", position="dodge", shape='diamond', size=5)+
  stat_summary(geom = "errorbar", fun.data=mean_cl_boot, position = position_dodge(.8), size=0.5, width=.8)+
  ylim(3,7)+
  ylab("pH") + 
  xlab("Time")+
  theme_pubr(base_size = 10)
#ph

(cocu )/(ph +plot_layout(height = c(4,1))) + plot_annotation(tag_levels="A")

coculture %>% 
  filter(Time == 24) %>% 
  aov(pH ~ Culture, .) %>% summary.aov
pairwise_pH <- coculture %>%   
  filter(Time == 24) %>%
  tukey_hsd(pH ~ Culture) %>% view()


######################################################################################################
####### netB expression #############
######################################################################################################

netB <- read_xlsx('Meisch_et_al_Data_Tables.xlsx', sheet = "NetB expression", na = c("", "NA", "N/A") )
glimpse(netB)

control <- filter(netB, batch =='Control')
AG02_2.5 <- filter(netB, batch =='2.5% Ba')
AG02_5 <- filter(netB, batch =='5% Ba')
La <- filter(netB, batch == "3% La")

summary(netB)

netB$batch <- factor(netB$batch, levels = c('Control', '2.5% Ba', '5% Ba', '3% La'))
ggdensity(netB, x = 'AD4700001_02732', add="mean", rug=TRUE, color = 'batch') 
ggdensity(netB, x = 'AD4700001_02732', add='mean', rug =TRUE, color = 'batch', fill='batch', xlab="netB expression (log normalized)", palette = c('#0075CF', '#FF7F30','#90C715', '#93358D'))



############ Figure 3B ############

netB %>%
  ggplot(aes(x = factor(batch, 
                        levels = c("Control", "2.5% Ba", "5% Ba", "3% La")))) +
  geom_bar(aes(fill = batch)) +
  facet_wrap(~louvain, ncol = 4) +
  theme_pubr() +
  theme(axis.text.x = element_blank(), legend.title = element_blank()) +
  ylab("Number of Cells") +
  xlab("Treatment")


############ Fisher' exact test Treatment - Control - Cells expessing netB with expression > 0
threshold <- 0.0

# Summary table
summary_tbl <- netB %>%
  mutate(NetB = AD4700001_02732 > threshold) %>%
  group_by(batch) %>%
  summarise(Total_cells = n(),
    NetB_expressing = sum(NetB),
    Non_expressing = sum(!NetB),
    Percent_expressing = round(100 * NetB_expressing / Total_cells, 1),
    .groups = "drop")


control <- summary_tbl %>%
  filter(batch == "Control")

fisher_results <- summary_tbl %>%
  filter(batch != "Control") %>%
  rowwise() %>%
  mutate(p.value = fisher.test( matrix(
    c(control$NetB_expressing,
      control$Non_expressing,
      NetB_expressing,
      Non_expressing), 
    nrow = 2, byrow = TRUE))$p.value) %>%
  ungroup() %>%
  select(batch, p.value)

final_table <- summary_tbl %>%
  left_join(fisher_results, by = "batch")

final_table <- final_table %>%
  mutate( batch = factor( batch, levels = c("Control", "2.5% Ba", "5% Ba", "3% La"))) %>%
  arrange(batch)

final_table



############ Figure 4 and cut-off analyses ############

cutoffs <- seq(0, 0.2, by = 0.005)
samples <- list(Control = control, `2.5% Ba` = AG02_2.5, `5% Ba` = AG02_5,`3% La` = La)

final_df <- map_dfr(cutoffs, function(cutoff) {
  summary_tbl <- map_dfr(names(samples), function(trt) {
    df <- samples[[trt]]
    expressing <- sum(df$AD4700001_02732 > cutoff)
    total <- nrow(df)
    tibble(Sample = trt,
           cutoff = cutoff,
           `no expression` = total - expressing,
           `netB expression` = expressing,
           netB = expressing / total)
  })
  
  control_row <- summary_tbl %>%
    filter(Sample == "Control")
  
  summary_tbl %>%
    rowwise() %>%
    mutate(p.value = case_when(
      Sample == "Control" ~ NA_real_, TRUE ~ fisher.test(matrix(
        c(control_row$`netB expression`,
          control_row$`no expression`,
          `netB expression`,
          `no expression`),
        nrow = 2,byrow = TRUE))$p.value)) %>%
    ungroup()
  
})



nb <- ggplot(final_df, aes(cutoff, netB*100, color = Sample))+
  geom_point()+
  theme_pubr()+
  ylim(0,100)+
  theme(legend.position="none")+
  labs(x = "Cutoff expression", 
       y = expression(atop("Percent of cells expressing", 
                           italic("netB"))))
nb

f <- ggplot(final_df, aes(cutoff, p.value, color=Sample))+
  geom_point()+
  geom_line()+
  theme_pubr()+
  theme(legend.position="none")+
  geom_hline(yintercept=0.05)+
  ylab("F-Test p-value")+
  xlab("Cutoff expression")

nb + f


############ Figure 5B ############

##. Kruskal-Wallis test with Dunn post-hoc pairwise comparisons with Bonferroni correction

# Overall test
kruskal_test(netB, AD4700001_02732 ~ batch)
stat.test <- netB %>%
  dunn_test(AD4700001_02732 ~ batch, p.adjust.method = "bonferroni")
stat.test

stat.test <- netB %>%
  dunn_test(AD4700001_02732 ~ batch, p.adjust.method = "bonferroni") %>%
  filter(group1 == "Control" | group2 == "Control") %>%
  mutate( p.adj.label = sprintf("%.3f", p.adj))

netB %>%
  mutate(rank = rank(AD4700001_02732, ties.method = "average")) %>%
  group_by(batch) %>%
  summarise(
    mean_expr = mean(AD4700001_02732),
    median_expr = median(AD4700001_02732),
    mean_rank = mean(rank),
    .groups = "drop") %>%
  arrange(desc(mean_rank))

ggviolin(netB, x = "batch", y = "AD4700001_02732",
         ylab = "NetB expression (log transformed)", 
         xlab = "Treatment", 
         ggtheme = theme_pubr(), 
         fill = "batch") +
  stat_summary(fun = mean, geom = "point", shape = 95, size = 12, color = "white")+
  # stat_summary(fun = mean,geom = "point", shape = 18, size = 4, color = "white") +
  stat_pvalue_manual(stat.test, label = "p.adj.label", y.position = c(1.0, 1.1, 1.2), size = 4) +
  theme(legend.title = element_blank())



######################################################################################################
####### Cytotoxicity #############
######################################################################################################


############ Figure 7E ############

cytotox <- read_xlsx('Meisch_et_al_Data_Tables.xlsx', sheet = "Cytotoxicity", na = c("", "NA", "N/A") )
glimpse(cytotox)

## ANOVA
fit <- aov(Viability ~ Sample, data = cytotox)
summary(fit)

## Post-hoc Tukey HSD TEST
tukey_results <- TukeyHSD(fit) %>% broom::tidy()

tukey_results

stat.test2 <- tukey_results %>%
  filter(adj.p.value < 0.05) %>%
  mutate(adj.p.value.label = sprintf("%.3f", adj.p.value))

## Selecting significant pairwise comparisons
stat.test2 <- stat.test2 %>%
  separate(contrast, into = c("group1", "group2"), sep = "-")
stat.test2

stat.test2$y.position <- c(130, 140, 120)
stat_pvalue_manual(stat.test2,label = "adj.p.value.label")

# Order samples
cytotox$Sample <- factor(cytotox$Sample, levels = c("BHI", "CP", "3% AG01", "2.5% AG02", "5% AG02"))

# Mean and 95% CI
cytotox_ci <- cytotox %>%
  group_by(Sample) %>%
  summarise(n = sum(!is.na(Viability)), mean = mean(Viability, na.rm = TRUE),
            sd = sd(Viability, na.rm = TRUE), ci = qt(0.975, df = n - 1) * sd / sqrt(n),
            .groups = "drop")
cytotox_ci

# Plotting with 95% CI and significance
ggplot(cytotox_ci) +
  geom_col(aes(Sample, mean, fill=Sample), width = 0.5, show.legend = FALSE)+
  geom_errorbar(aes(Sample, mean, ymin=mean-ci, ymax=mean+ci), width = 0.2)+
  stat_pvalue_manual(stat.test2, label = "adj.p.value.label", size = 4)+
  scale_x_discrete(limits = c('BHI', 'CP', '3% AG01', '2.5% AG02', '5% AG02'), 
                   labels=c("BHI"="Media","CP"="Cp","3% AG01" = "Cp+3% La", "2.5% AG02" = "Cp+2.5% Ba", "5% AG02" = "Cp+5% Ba"))+
  theme_pubr()+
  ylim(0,150)+
  ylab("Epithelial Cell Viability (%)") + 
  xlab("")




######################################################################################################
####### Proof of Concept Trail #############
######################################################################################################


########## Figure 8: POC Lesion scores ################

lesions <- read_xlsx('Meisch_et_al_Data_Tables.xlsx', sheet = "POC_Lesion_Summary", na = c("", "NA", "N/A") )

lesions %>%
  mutate(Score = factor(Score), Treatment = recode(Treatment,
                                                   "NC" = "NC", "CC" = "CC", "AD-1E8" = "AD", "ED-1E8" = "ED"),
         Treatment = factor(Treatment, levels = c("NC", "CC", "AD", "ED"))) %>%
  ggplot(aes(fill = Score, x = Treatment, y = Birds)) +
  geom_bar(position = "fill", stat = "identity", width = 0.6) +
  scale_y_continuous(labels = scales::percent) +
  ylab("Percent of birds with lesions")+
  theme_pubr(base_size = 12) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")

