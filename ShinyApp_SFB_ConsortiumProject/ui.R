## Shiny App for visualizing SFB "HIT" consortium datasets
## Moritz Madern
## 08.08.2022



## Load packages and define relevant variables --------------------------------

# Load required packages
library(shiny)
library(shinythemes)
library(tidyverse)
library(cowplot)
library(scales)

# specify groups color schemes
colors_groups = c(WT_naive = "#DDF2D1", 
                  KO_naive = "#9FD4A3",
                  WT_Th1 = "#D6ECFF", 
                  KO_Th1 = "#9DCAEB", 
                  WT_Th2 = "#E6D5FF",
                  KO_Th2 = "#BE94E6",
                  WT_Th17 = "#fbd9d3", 
                  KO_Th17 = "#D16B75",
                  WT_Treg = "#D9C0B3",
                  KO_Treg = "#856B61",
                  REF= "black")
colors_time = set_names(x = c("#DDF2D1", "#68C5C7", "#369ABE", "#2A6CAF", "#253494"),
                        nm = "naive", "30m", "2h", "8h", "24h")         







## Read in and manipulate relevant data ----------------------------------------

# Read in Tcell activation proteomics data, prepare selection choices
df_proteins_TcellAct <- read.delim("data/aggregated_proteinGroups_TcellActivation.txt", header=TRUE, sep="\t")   # protein data wTcellActivation
geneNames_dup <- names(table(df_proteins_TcellAct$Gene.names))[table(df_proteins_TcellAct$Gene.names) > 1]
df_proteins_TcellAct$Gene.names[df_proteins_TcellAct$Gene.names %in% geneNames_dup] <- paste0(df_proteins_TcellAct$Gene.names, "_", df_proteins_TcellAct$Majority.protein.IDs)[df_proteins_TcellAct$Gene.names %in% geneNames_dup]   #handle duplicated gene names (by adding majority protein ID)
protein_choices_TcellAct<- df_proteins_TcellAct$Gene.names
bool_protein_reporter <- grepl(names(df_proteins_TcellAct), pattern="norm__aggregated")
samplenames_proteinTable_TcellAct = c("naive.1", "30m.3", "naive.2", "8h.1",   # specify sample names of input table
                                      "naive.3", "8h.2", "30m.1", "8h.3",
                                      "30m.2", "24h.1", "2h.1", "24h.2", 
                                      "2h.2", "24h.3", "2h.3", "boost")
names(df_proteins_TcellAct)[bool_protein_reporter] <- samplenames_proteinTable_TcellAct
groups_proteinTabe_TcellAct = c("naive", "30m", "naive", "8h", # specify group names of input table
                                "naive", "8h", "30m", "8h",
                                "30m", "24h", "2h", "24h", 
                                "2h", "24h", "2h", "boost")
samplenames_proteinTable_TcellAct_subselect = c("naive.1","naive.2","naive.3", # specify subselection/reordering for plotting
                                                "30m.1", "30m.2", "30m.3",
                                                "2h.1", "2h.2", "2h.3",
                                                "8h.1", "8h.2", "8h.3",
                                                "24h.1", "24h.2", "24h.3")
groups_proteinTable_TcellAct_subselect <- set_names(groups_proteinTabe_TcellAct, nm=samplenames_proteinTable_TcellAct)[samplenames_proteinTable_TcellAct_subselect]
df_phospho_TcellAct <- read.delim("data/SiteToProtein_Phospho (STY)Sites_TcellActivation.txt", header=TRUE, sep="\t", check.names = FALSE) # siteToProtein data woNaive
phosphoSite_choices_TcellAct <- df_phospho_TcellAct$unique_site_identifier
phosphoSite_choices_TcellAct_proteins <- df_phospho_TcellAct$Protein.group.IDs
df_acet_TcellAct <- read.delim("data/SiteToProtein_Acetyl (K)Sites_TcellActivation.txt", header=TRUE, sep="\t", check.names = FALSE) # siteToProtein data woNaive
acetSite_choices_TcellAct <- df_acet_TcellAct$unique_site_identifier
acetSite_choices_TcellAct_proteins <- df_acet_TcellAct$Protein.group.IDs


# Read in RNAseq data allSamples
df_RNA_all <- read.delim("data/Transcriptome_bulkRNAseq_allSamples.txt", header=TRUE, sep="\t")
df_RNA_all$Gene.names <- paste0(df_RNA_all$gene_name,"__", df_RNA_all$gene_id)  # create unique entries
gene_choices_all <- df_RNA_all$Gene.names
bool_RNA_counts <- grepl(names(df_RNA_all), pattern="normalized")
m_RNA_all <- df_RNA_all[,bool_RNA_counts] %>% as.matrix()
samplenames_RNAtable_all = c("WT_naive.1", "WT_naive.2", "WT_naive.3", # specify sample names of input table
                             "KO_naive.1", "KO_naive.2", "KO_naive.3", 
                             "WT_Th1.1", "WT_Th1.2", "WT_Th1.3", 
                             "KO_Th1.1", "KO_Th1.2", "KO_Th1.3", 
                             "WT_Th2.1", "WT_Th2.2", "WT_Th2.3", 
                             "KO_Th2.1", "KO_Th2.2", "KO_Th2.3", 
                             "WT_Th17.1", "WT_Th17.2", "WT_Th17.3", 
                             "KO_Th17.1", "KO_Th17.2", "KO_Th17.3", 
                             "WT_Treg.1", "WT_Treg.2", "WT_Treg.3", 
                             "KO_Treg.1", "KO_Treg.2", "KO_Treg.3")
names(df_RNA_all)[bool_RNA_counts] <- samplenames_RNAtable_all
groups_RNAtable_all = c("WT_naive", "WT_naive", "WT_naive", # specify sample names of input table
                        "KO_naive", "KO_naive", "KO_naive", 
                        "WT_Th1", "WT_Th1", "WT_Th1", 
                        "KO_Th1", "KO_Th1", "KO_Th1", 
                        "WT_Th2", "WT_Th2", "WT_Th2", 
                        "KO_Th2", "KO_Th2", "KO_Th2", 
                        "WT_Th17", "WT_Th17", "WT_Th17", 
                        "KO_Th17", "KO_Th17", "KO_Th17", 
                        "WT_Treg", "WT_Treg", "WT_Treg", 
                        "KO_Treg", "KO_Treg", "KO_Treg")


# Read in RNAseq data woNaive
df_RNA_woNaive <- read.delim("data/Transcriptome_bulkRNAseq_woNaive.txt", header=TRUE, sep="\t")
df_RNA_woNaive$Gene.names <- paste0(df_RNA_woNaive$gene_name,"__", df_RNA_woNaive$gene_id)  # create unique entries
gene_choices_woNaive <- df_RNA_woNaive$Gene.names
bool_RNA_counts <- grepl(names(df_RNA_woNaive), pattern="normalized")
m_RNA_woNaive <- df_RNA_woNaive[,bool_RNA_counts] %>% as.matrix()
samplenames_RNAtable_woNaive = c("WT_Th1.1", "WT_Th1.2", "WT_Th1.3", # specify sample names of input table
                                 "KO_Th1.1", "KO_Th1.2", "KO_Th1.3", 
                                 "WT_Th2.1", "WT_Th2.2", "WT_Th2.3", 
                                 "KO_Th2.1", "KO_Th2.2", "KO_Th2.3", 
                                 "WT_Th17.1", "WT_Th17.2", "WT_Th17.3", 
                                 "KO_Th17.1", "KO_Th17.2", "KO_Th17.3", 
                                 "WT_Treg.1", "WT_Treg.2", "WT_Treg.3", 
                                 "KO_Treg.1", "KO_Treg.2", "KO_Treg.3")
names(df_RNA_woNaive)[bool_RNA_counts] <- samplenames_RNAtable_woNaive
groups_RNAtable_woNaive = c("WT_Th1", "WT_Th1", "WT_Th1", # specify sample names of input table
                            "KO_Th1", "KO_Th1", "KO_Th1", 
                            "WT_Th2", "WT_Th2", "WT_Th2", 
                            "KO_Th2", "KO_Th2", "KO_Th2", 
                            "WT_Th17", "WT_Th17", "WT_Th17", 
                            "KO_Th17", "KO_Th17", "KO_Th17", 
                            "WT_Treg", "WT_Treg", "WT_Treg", 
                            "KO_Treg", "KO_Treg", "KO_Treg")


# Read in TMTset1 data woNaive, prepare selection choices
df_proteins_TMTset1 <- read.delim("data/ProteinGroups_filtered_normalized_TMTset1.txt", header=TRUE, sep="\t")   # protein data woNaive
protein_choices_TMTset1_woNaive <- df_proteins_TMTset1$Gene.names
bool_protein_reporter <- grepl(names(df_proteins_TMTset1), pattern="Reporter[.]intensity")
samplenames_proteinTable_TMTset1 = c("WT_naive.2", "WT_naive.1", "WT_naive.3", # specify sample names of input table
                                     "WT_Th17.1", "WT_Th1.2", "WT_Th17.2", 
                                     "WT_Th1.1", "WT_Th17.3", "WT_Th1.3", 
                                     "WT_Treg.1", "WT_Th2.1", "WT_Treg.2", 
                                     "WT_Th2.2", "WT_Treg.3", "WT_Th2.3", "REF")
names(df_proteins_TMTset1)[bool_protein_reporter] <- samplenames_proteinTable_TMTset1
groups_proteinTabe_TMTset1 = c("WT_naive", "WT_naive", "WT_naive", # specify group names of input table
                               "WT_Th17", "WT_Th1", "WT_Th17", 
                               "WT_Th1", "WT_Th17", "WT_Th1", 
                               "WT_Treg", "WT_Th2", "WT_Treg", 
                               "WT_Th2", "WT_Treg", "WT_Th2", "REF")
samplenames_proteinTable_TMTset1_woNaive = c("WT_Th1.1","WT_Th1.2","WT_Th1.3", # specify subselection for w/o naive
                                             "WT_Th2.1", "WT_Th2.2", "WT_Th2.3",
                                             "WT_Th17.1", "WT_Th17.2", "WT_Th17.3",
                                             "WT_Treg.1", "WT_Treg.2", "WT_Treg.3")
groups_proteinTable_TMTset1_woNaive <- set_names(groups_proteinTabe_TMTset1, nm=samplenames_proteinTable_TMTset1)[samplenames_proteinTable_TMTset1_woNaive]
df_phospho_TMTset1_woNaive <- read.delim("data/SiteToProtein_Phospho (STY)Sites_Set1_woNaive.txt", header=TRUE, sep="\t") # siteToProtein data woNaive
phosphoSite_choices_TMTset1_woNaive <- df_phospho_TMTset1_woNaive$unique_site_identifier
phosphoSite_choices_TMTset1_woNaive_proteins <- df_phospho_TMTset1_woNaive$Protein.group.IDs
df_acet_TMTset1_woNaive <- read.delim("data/SiteToProtein_Acetyl (K)Sites_Set1_woNaive.txt", header=TRUE, sep="\t") # siteToProtein data woNaive
acetSite_choices_TMTset1_woNaive <- df_acet_TMTset1_woNaive$unique_site_identifier
acetSite_choices_TMTset1_woNaive_proteins <- df_acet_TMTset1_woNaive$Protein.group.IDs








## Define Shiny UI --------------------------------------------------------
ui <- navbarPage("SFB HIT datasets",theme = shinytheme("cerulean"), collapsible = TRUE,
                 
                 
                 ## UI T-cell activation -------------------------------------
                 tabPanel("Proteomics T-cell Activation",
                          tags$h2("Dataset Overview"),
                          
                          # T-cell activation images 
                          fluidRow(
                            column(3,
                                   img(src="TcellActivation_experiment_design.png", align = "left", height="140%", width="140%")
                            ),
                            column(11
                            )    
                          ),
                          
                          
                          # Protein level T-cell activation  ---------------------
                          br(),
                          br(),
                          br(),
                          br(),
                          tags$h2("Protein level"),
                          
                          # T-cell Activation protein selection
                          selectInput(inputId = "selected_protein_TcellAct",
                                      label = "Specify protein's gene name:",
                                      choices = protein_choices_TcellAct,
                                      multiple = FALSE, width=500),
                          
                          # T-cell Activation protein ANOVA and plotting
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups on protein level via ANOVA:"),
                                   helpText("model: log2 Intensity ~ group "),
                                   textOutput(outputId = "text_protein_pval_ANOVA_TcellAct"),
                                   textOutput(outputId = "text_protein_adj.pval_ANOVA_TcellAct"),
                                   helpText(""),
                                   helpText("Calculated purity metrics for this protein:"),           
                                   textOutput(outputId = "text_protein_PPF_TcellAct"),
                                   textOutput(outputId = "text_protein_EIL_TcellAct"),
                                   helpText("Number of quantified PTMs for this protein:"),           
                                   textOutput(outputId = "text_protein_nrPhosphoSites_TcellAct"),
                                   textOutput(outputId = "text_protein_nrAcetSites_TcellAct"),
                            ),
                            column(8,
                                   plotOutput(outputId = "plot_protein_TcellAct", height="350px")
                            )      
                          ),
                          
                          
                          # Phospho Site level T-cell activation   -------------------
                          br(),
                          br(),
                          br(),
                          br(),
                          tags$h2("Phosphorylation site level"),
                          
                          # TT-cell Activation phosphoSite selection
                          selectInput(inputId = "selected_phosphoSite_TcellAct",
                                      label = "Specify phospho site of protein selected above:",
                                      choices = NULL, multiple = FALSE, width = 500, selected = ""),
                          helpText("Note: If no PTM site was quantified for this protein, the plots below will appear empty per default."),
                          
                          # T-cell Activation phosphoSite statistics and plotting
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups on phospho site level via ANOVA:"),
                                   helpText("model: log2 Intensity ~ group "),
                                   textOutput(outputId = "text_phospho_pval_ANOVA_TcellAct"),
                                   textOutput(outputId = "text_phospho_adj.pval_ANOVA_TcellAct"),
                                   helpText(""),
                                   helpText("Calculated purity metrics for this site:"),           
                                   textOutput(outputId = "text_phospho_PPF_TcellAct"),
                                   textOutput(outputId = "text_phospho_EIL_TcellAct"),
                            ),
                            column(8,
                                   plotOutput(outputId = "plot_phospho_TcellAct", height="350px")
                            )
                          ),
                          
                          # T-cell Activation grid plot site/protein normalization (non-batch corrected)
                          tags$h4("Phospho - site/protein normalization:"),
                          helpText("The plots below display results of site-to-protein normalization for the selected PTM site. To counteract false-positives resulting from differences in ion interference, intensities were additionally adjusted to equal interference (IF) levels prior to site-to-protein normalization (right-most column):"),
                          fluidRow(
                            column(12,
                                   plotOutput(outputId = "plot_grid_phospho_siteToProtein_TcellAct", height="200px"),
                            ),
                          ),
                          helpText("Testing for overall differences between IF-adjusted site/protein ratios via ANOVA:"),
                          helpText("model: ratio ~ group "),
                          textOutput(outputId = "text_phospho_pval_siteToProtein_IFadjust_TcellAct"),
                          textOutput(outputId = "text_phospho_adj.pval_siteToProtein_IFadjust_TcellAct"),
                          br(),
                          
                          
                          
                          # Acetyl Site level T-cell activation  -------------------
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          tags$h2("Acetylation site level"),
                          
                          # T-cell Activation acetSite selection
                          selectInput(inputId = "selected_acetSite_TcellAct",
                                      label = "Specify acetyl site of protein selected above:",
                                      choices = NULL, multiple = FALSE, width = 500, selected = ""),
                          helpText("Note: If no PTM site was quantified for this protein, the plots below will appear empty per default."),
                          
                          # T-cell Activation acetSite statistics and plotting
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups on acetyl site level via ANOVA:"),
                                   helpText("model: log2 Intensity ~ group "),
                                   textOutput(outputId = "text_acet_pval_ANOVA_TcellAct"),
                                   textOutput(outputId = "text_acet_adj.pval_ANOVA_TcellAct"),
                                   helpText(""),
                                   helpText("Calculated purity metrics for this site:"),           
                                   textOutput(outputId = "text_acet_PPF_TcellAct"),
                                   textOutput(outputId = "text_acet_EIL_TcellAct"),
                            ),
                            column(8,
                                   plotOutput(outputId = "plot_acet_TcellAct", height="350px")
                            )
                          ),
                          
                          # T-cell Activation grid plot site/protein normalization (non-batch corrected)
                          tags$h4("Acetyl - site/protein normalization:"),
                          helpText("The plots below display results of site-to-protein normalization for the selected modified peptide. To counteract false-positives due to differences in ion interference, intensities were additionally adjusted to equal interference (IF) levels prior to site-to-protein normalization (right-most column):"),
                          fluidRow(
                            column(12,
                                   plotOutput(outputId = "plot_grid_acet_siteToProtein_TcellAct", height="200px"),
                            ),
                          ),
                          helpText("Testing for overall differences between IF-adjusted site/protein ratios via ANOVA:"),
                          helpText("model: ratio ~ group "),
                          textOutput(outputId = "text_acet_pval_siteToProtein_IFadjust_TcellAct"),
                          textOutput(outputId = "text_acet_adj.pval_siteToProtein_IFadjust_TcellAct"),
                          br(),
                 ),
                 
                 
                 ## UI bulk RNAseq allSamples -------------------------------------
                 tabPanel("RNAseq (all samples)",
                          tags$h2("Dataset Overview"),
                         
                          # RNA all Samples images 
                          fluidRow(
                            column(4,
                                   img(src='bulkRNAseq_experiment_design.png', align = "left", height="105%", width="105%")
                            ),
                            column(8
                            )    
                          ),
                          
                          br(),
                          br(),
                          tags$h2("RNA level"),
                          
                          # RNA all gene selection
                          selectInput(inputId = "selected_gene_all",
                                      label = "Specify gene name:",
                                      choices = gene_choices_all,
                                      multiple = FALSE, width=500),
                          
                          # RNA plotting and ANOVA
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups via ANOVA:"),
                                   helpText("model: log2 counts ~ group + batch "),
                                   textOutput(outputId = "text_rna_pval_ANOVA_all"),
                                   textOutput(outputId = "text_rna_adj.pval_ANOVA_all"),
                                   ),
                            column(8,
                                   plotOutput(outputId = "plot_rna_all", height="350px")
                            )      
                          ),
 
                          # RNA - testing for differences between WT and KO
                          tags$h4("Addressing WT vs KO for selected gene"),
                          fluidRow(
                            column(11,
                                   helpText("Testing for differences between WT and KO across all groups via linear regression model F-test:"),
                                   helpText("model: log2 counts ~ celltype + batch + celltype:genotype"),
                                   textOutput(outputId = "text_rna_pval_WTvsKO_all"),
                                   textOutput(outputId = "text_rna_adj.pval_WTvsKO_all"),
                                   helpText("Note: This F-test tests agaist the null hypothesis H0 that there are no differences between WT and KO in any of the various cell types. The alternative hypothesis H1 accordingly corresponds to the hypothesis that, for at least one cell type, WT and KO are different."),
                                   ),
                            column(1)
                            )      
                          ),
                          
                     
                      
 
                 
                 ## UI bulk RNAseq woNaive -------------------------------------
                 tabPanel("RNAseq (w/o naive)",
                          tags$h2("Dataset Overview"),
                          
                          # RNA all Samples images 
                          fluidRow(
                            column(4,
                                   img(src='bulkRNAseq_experiment_design.png', align = "left", height="105%", width="105%")
                            ),
                            column(8
                            )    
                          ),
                          
                          br(),
                          br(),
                          tags$h2("RNA level"),
                          
                          # RNA all gene selection
                          selectInput(inputId = "selected_gene_woNaive",
                                      label = "Specify gene name:",
                                      choices = gene_choices_woNaive,
                                      multiple = FALSE, width=500),
                          
                          # RNA plotting and ANOVA
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups via ANOVA:"),
                                   helpText("model: log2 counts ~ group + batch "),
                                   textOutput(outputId = "text_rna_pval_ANOVA_woNaive"),
                                   textOutput(outputId = "text_rna_adj.pval_ANOVA_woNaive"),
                            ),
                            column(8,
                                   plotOutput(outputId = "plot_rna_woNaive", height="350px")
                            )      
                          ),
                          
                          # RNA - testing for differences between WT and KO
                          tags$h4("Addressing WT vs KO for selected gene"),
                          fluidRow(
                            column(11,
                                   helpText("Testing for differences between WT and KO across all groups via linear regression model F-test:"),
                                   helpText("model: log2 counts ~ celltype + batch + celltype:genotype"),
                                   textOutput(outputId = "text_rna_pval_WTvsKO_woNaive"),
                                   textOutput(outputId = "text_rna_adj.pval_WTvsKO_woNaive"),
                                   helpText("Note: This F-test tests agaist the null hypothesis H0 that there are no differences between WT and KO in any of the various cell types. The alternative hypothesis H1 accordingly corresponds to the hypothesis that for at least one cell type, WT and KO are different."),
                            ),
                            column(1)
                          )      
                 ),
                 
                 
                 
                 
                  
 
                 
                 
                 ## UI TMTset1 wo Naive -------------------------------------
                 tabPanel("Proteomics TMTset1 (w/o naive)",
                          tags$h2("Dataset Overview"),
                            
                          # TMTset1  images 
                          fluidRow(
                            column(3,
                                   img(src='TMTset1_experiment_design.png', align = "left", height="135%", width="135%")
                            ),
                            column(11
                            )    
                          ),
                          
                          
                          # Protein level TMTset1 woNaive  ---------------------
                          br(),
                          br(),
                          br(),
                          br(),
                          tags$h2("Protein level"),
                          
                          # TMTset1 woNaive protein selection
                          selectInput(inputId = "selected_protein_TMTset1_woNaive",
                                      label = "Specify protein's gene name:",
                                      choices = protein_choices_TMTset1_woNaive,
                                      multiple = FALSE, width=500),
                            
                          # TMTset1 woNaive protein ANOVA and plotting
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups on protein level via ANOVA:"),
                                   helpText("model: log2 Intensity ~ group "),
                                   textOutput(outputId = "text_protein_pval_ANOVA_TMTset1_woNaive"),
                                   textOutput(outputId = "text_protein_adj.pval_ANOVA_TMTset1_woNaive"),
                                   helpText("Number of quantified PTMs for this protein:"),           
                                   textOutput(outputId = "text_protein_nrPhosphoSites_TMTset1_woNaive"),
                                   textOutput(outputId = "text_protein_nrAcetSites_TMTset1_woNaive"),
                                   ),
                            column(8,
                                     plotOutput(outputId = "plot_protein_TMTset1_woNaive", height="350px")
                              )      
                            ),
                          helpText("Note: Due to MS3-quantification, ion interference at protein level is 0 (i.e. EIL = 0)."),
                            
                            
                          # Phospho Site level TMTset1 woNaive  -------------------
                          br(),
                          br(),
                          br(),
                          br(),
                          tags$h2("Phosphorylation site level"),
                            
                          # TMTset1 woNaive phosphoSite selection
                          selectInput(inputId = "selected_phosphoSite_TMTset1_woNaive",
                                      label = "Specify phospho site of protein selected above:",
                                      choices = NULL, multiple = FALSE, width = 500, selected = ""),
                          helpText("Note: If no PTM site was quantified for this protein, the plots below will appear empty per default."),
                          
                          # TMTset1 woNaive phosphoSite statistics and plotting
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups on phospho site level via ANOVA:"),
                                   helpText("model: log2 Intensity ~ group "),
                                   textOutput(outputId = "text_phospho_pval_ANOVA_TMTset1_woNaive"),
                                   textOutput(outputId = "text_phospho_adj.pval_ANOVA_TMTset1_woNaive"),
                                   helpText(""),
                                   helpText("Calculated purity metrics for this site:"),           
                                   textOutput(outputId = "text_phospho_PPF"),
                                   textOutput(outputId = "text_phospho_EIL"),
                            ),
                            column(8,
                                   plotOutput(outputId = "plot_phospho_TMTset1_woNaive", height="350px")
                            )
                          ),
                          
                          # TMTset1 phospho site woNaive grid plot site/protein normalization (non-batch corrected)
                          tags$h4("Phospho - site/protein normalization:"),
                          helpText("The plots below display results of site-to-protein normalization for the selected PTM site. To counteract false-positives resulting from differences in ion interference, intensities were additionally adjusted to equal interference (IF) levels prior to site-to-protein normalization (right-most column):"),
                          fluidRow(
                            column(12,
                                   plotOutput(outputId = "plot_grid_phospho_siteToProtein_TMTset1_woNaive", height="200px"),
                            ),
                          ),
                          helpText("Testing for overall differences between IF-adjusted site/protein ratios via ANOVA:"),
                          helpText("model: ratio ~ group "),
                          textOutput(outputId = "text_phospho_pval_siteToProtein_IFadjust_TMTset1_woNaive"),
                          textOutput(outputId = "text_phospho_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive"),
                          br(),
                            
                          # TMTset1 phospho site woNaive grid plot site/protein normalization (batch corrected)
                          tags$h4("Phospho - site/protein normalization (batch-corrected):"),
                          helpText("The same data after batch-correction via the comBat algorithm:"),
                          fluidRow(
                            column(12,
                                   plotOutput(outputId = "plot_grid_phospho_siteToProtein_TMTset1_woNaive_batchCorr", height="200px"),
                            ),
                          ),
                          helpText("Testing for overall differences between IF-adjusted batch-corrected site/protein ratios via ANOVA:"),
                          helpText("model: ratio ~ group "),
                          textOutput(outputId = "text_phospho_pval_siteToProtein_IFadjust_TMTset1_woNaive__batchCorr"),
                          textOutput(outputId = "text_phospho_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive_batchCorr"),
                            
                        
                          # Acetyl Site level woNaive  -------------------
                          br(),
                          br(),
                          br(),
                          br(),
                          br(),
                          tags$h2("Acetylation site level"),
                            
                          # TMTset1 woNaive acetSite selection
                          selectInput(inputId = "selected_acetSite_TMTset1_woNaive",
                                      label = "Specify acetyl site of protein selected above:",
                                      choices = NULL, multiple = FALSE, width = 500, selected = ""),
                          helpText("Note: If no PTM site was quantified for this protein, the plots below will appear empty per default."),
                            
                          # TMTset1 woNaive acetSite statistics and plotting
                          fluidRow(
                            column(4,
                                   helpText("Testing for overall differences between groups on acetyl site level via ANOVA:"),
                                   helpText("model: log2 Intensity ~ group "),
                                   textOutput(outputId = "text_acet_pval_ANOVA_TMTset1_woNaive"),
                                   textOutput(outputId = "text_acet_adj.pval_ANOVA_TMTset1_woNaive"),
                                   helpText(""),
                                   helpText("Calculated purity metrics for this site:"),           
                                   textOutput(outputId = "text_acet_PPF"),
                                   textOutput(outputId = "text_acet_EIL"),
                            ),
                            column(8,
                                   plotOutput(outputId = "plot_acet_TMTset1_woNaive", height="350px")
                            )
                          ),
                          
                          # TMTset1 acetyl site woNaive grid plot site/protein normalization (non-batch corrected)
                          tags$h4("Acetyl - site/protein normalization:"),
                          helpText("The plots below display results of site-to-protein normalization for the selected modified peptide. To counteract false-positives due to differences in ion interference, intensities were additionally adjusted to equal interference (IF) levels prior to site-to-protein normalization (right-most column):"),
                          fluidRow(
                            column(12,
                                    plotOutput(outputId = "plot_grid_acet_siteToProtein_TMTset1_woNaive", height="200px"),
                            ),
                          ),
                          helpText("Testing for overall differences between IF-adjusted site/protein ratios via ANOVA:"),
                          helpText("model: ratio ~ group "),
                          textOutput(outputId = "text_acet_pval_siteToProtein_IFadjust_TMTset1_woNaive"),
                          textOutput(outputId = "text_acet_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive"),
                          br(),
                            
                          # TMTset1 acetyl site woNaive grid plot site/protein normalization (batch corrected)
                          tags$h4("Acetyl - site/protein normalization (batch-corrected):"),
                          helpText("The same data after batch-correction via the comBat algorithm:"),
                          fluidRow(
                            column(12,
                                     plotOutput(outputId = "plot_grid_acet_siteToProtein_TMTset1_woNaive_batchCorr", height="200px"),
                            ),
                          ),
                          helpText("Testing for overall differences between IF-adjusted batch-corrected site/protein ratios via ANOVA:"),
                          helpText("model: ratio ~ group "),
                          textOutput(outputId = "text_acet_pval_siteToProtein_IFadjust_TMTset1_woNaive_batchCorr"),
                          textOutput(outputId = "text_acet_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive_batchCorr"),
                        )
)






