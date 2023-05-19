

# Load packages and define relevant variables --------------------------------

# Load required packages
library(shiny)
library(shinythemes)
library(tidyverse)
library(cowplot)
library(scales)


# specify group color schemes
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
                  REF="black")
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













## Define server ------------------------------------------------------------------
server <- function(input, output, session) {
  
  
  


  ## Server Protein level TcellAct  -------------------------------

  # Prepare protein pattern plot of selected protein i
  output$plot_protein_TcellAct <- renderPlot({

    # Specify selected protein's gene name
    i = input$selected_protein_TcellAct

    # Filter for protein i information
    df_proteins_i <- df_proteins_TcellAct %>% filter(Gene.names == i)

    # Select intensity values of protein i. Log2-transform
    v_protein_i <- df_proteins_i[samplenames_proteinTable_TcellAct_subselect] %>% as.numeric() 

    
    # replace 0 with NA. Then log2 transform
    v_protein_i[v_protein_i == 0] <- NA
    v_protein_i <- log2(v_protein_i)
    
    
    # Create dataframe for plotting
    df_gg_protein_i <- data.frame(x = factor(samplenames_proteinTable_TcellAct_subselect,levels=samplenames_proteinTable_TcellAct_subselect),
                                  groups=factor(groups_proteinTable_TcellAct_subselect), levels=unique(groups_proteinTable_TcellAct_subselect),
                                  y = v_protein_i)

    # Specify y-range, then prepare plot
    y <- df_gg_protein_i$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_prot <- ggplot(data=df_gg_protein_i) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_time[levels(df_gg_protein_i$groups)]) +
      theme_bw() +
      ggtitle(i) +
      ylab("log2 Intensity") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none")
    return(gg_prot)
  })

  # Prepare protein level ANOVA p-value Tcell Act
  output$text_protein_pval_ANOVA_TcellAct <- renderText({
    i = input$selected_protein_TcellAct
    df_proteins_i <- df_proteins_TcellAct %>% filter(Gene.names == i)
    return(paste0("ANOVA p-value: ", round(df_proteins_i$pval_ANOVA, digits=4)))
  })

  # Prepare protein level ANOVA adj.p-value Tcell Act
  output$text_protein_adj.pval_ANOVA_TcellAct <- renderText({
    i = input$selected_protein_TcellAct
    df_proteins_i <- df_proteins_TcellAct %>% filter(Gene.names == i)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_proteins_i$adj.pval_ANOVA, digits=4)))
  })

  # Prepare protein level EIL value Tcell Act
  output$text_protein_PPF_TcellAct <- renderText({
    i = input$selected_protein_TcellAct
    df_proteins_i <- df_proteins_TcellAct %>% filter(Gene.names == i)
    return(paste0("Precursor Purity Fraction (PPF): ", round(df_proteins_i$PPF, digits=2)))
  })
  
  # Prepare protein level EIL value Tcell Act
  output$text_protein_EIL_TcellAct <- renderText({
    i = input$selected_protein_TcellAct
    df_proteins_i <- df_proteins_TcellAct %>% filter(Gene.names == i)
    return(paste0("Estimated Interference Level (EIL): ", round(df_proteins_i$EIL, digits=2)))
  })
  
  # Prepare nr acetyl sites Tcell Act
  output$text_protein_nrAcetSites_TcellAct <- renderText({
    i = input$selected_protein_TcellAct
    proteinID_i <- df_proteins_TcellAct %>% filter(Gene.names == i) %>% select(id) %>% as.numeric()
    nrAcetSites <- sum(acetSite_choices_TcellAct_proteins == proteinID_i)
    return(paste0("Acetylation sites: ", round(nrAcetSites, digits=0)))
  })
  
  # Prepare number phospho sites TcellAct
  output$text_protein_nrPhosphoSites_TcellAct <- renderText({
    i = input$selected_protein_TcellAct
    proteinID_i <- df_proteins_TcellAct %>% filter(Gene.names == i) %>% select(id) %>% as.numeric()
    nrPhosphoSites <- sum(phosphoSite_choices_TcellAct_proteins == proteinID_i)
    return(paste0("Phosphorylation sites: ", round(nrPhosphoSites, digits=0)))
  })
  

  

## Server Phospho level TcellAct  -------------------------------

# Update phospho site choice based on protein selection
observeEvent(input$selected_protein_TcellAct,
             {
               protein_id <- df_proteins_TcellAct %>% filter(Gene.names == input$selected_protein_TcellAct) %>% select(id)
               if (length(phosphoSite_choices_TcellAct[phosphoSite_choices_TcellAct_proteins %in% protein_id]) > 0){
                 choices = phosphoSite_choices_TcellAct[phosphoSite_choices_TcellAct_proteins %in% protein_id]
               } else {
                 choices = " "
               }
               updateSelectizeInput(session, input = "selected_phosphoSite_TcellAct", choices = choices)
             })



  # Prepare phospho site pattern plot of selected site a
  output$plot_phospho_TcellAct <- renderPlot({

    # Specify selected site's unique name
    a = input$selected_phosphoSite_TcellAct

    # Filter for site's a information
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)

    # Select site intensities of site a
    v_phospho_a <- df_phospho_a[,paste0(samplenames_proteinTable_TcellAct_subselect)] %>%  as.numeric()

    # Prepare dataframe for plotting
    df_gg_phospho_a <- data.frame(x=factor(samplenames_proteinTable_TcellAct_subselect, levels=samplenames_proteinTable_TcellAct_subselect),
                                  groups=factor(groups_proteinTable_TcellAct_subselect), levels=unique(groups_proteinTable_TcellAct_subselect),
                                  y=log2(v_phospho_a))

    # Specify y-range, then prepare plot
    y <- df_gg_phospho_a$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_phospho <- ggplot(data=df_gg_phospho_a) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_time[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(a) +
      ylab("log2 Intensity") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none")
    return(gg_phospho)
  })


  # Prepare phospho level ANOVA p-value wo Naive
  output$text_phospho_pval_ANOVA_TcellAct <- renderText({
    a = input$selected_phosphoSite_TcellAct
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value:  ", round(df_phospho_a$pval_ANOVA, digits=4)))
  })


  # Prepare phospho level ANOVA adj.p-value wo Naive
  output$text_phospho_adj.pval_ANOVA_TcellAct <- renderText({
    a = input$selected_phosphoSite_TcellAct
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value:  ", round(df_phospho_a$adj.pval_ANOVA, digits=4)))
  })


  # Prepare phospho level ANOVA adj.p-value wo Naive
  output$text_phospho_PPF_TcellAct <- renderText({
    a = input$selected_phosphoSite_TcellAct
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("Precursor Purity Fraction (PPF)   :   ", round(df_phospho_a$PPF, digits=2)))
  })


  # Prepare phospho level ANOVA adj.p-value wo Naive
  output$text_phospho_EIL_TcellAct <- renderText({
    a = input$selected_phosphoSite_TcellAct
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("Estimated Interference Level (EIL)   :   ", round(df_phospho_a$EIL, digits=2)))
  })
  
  
  
  
  
  


  # Prepare phospho site grid plot
  output$plot_grid_phospho_siteToProtein_TcellAct <- renderPlot({

    # Specify selected site's unique name
    a = input$selected_phosphoSite_TcellAct

    # Filter for site's a information
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)

    # Select site intensities of site a
    v_phospho_a <- df_phospho_a[,paste0(samplenames_proteinTable_TcellAct_subselect)] %>%  as.numeric()
    v_protein_a <- df_phospho_a[,paste0(samplenames_proteinTable_TcellAct_subselect, "__underlyingProtein")] %>%  as.numeric()
    v_siteToProtein_a <- df_phospho_a[,paste0(samplenames_proteinTable_TcellAct_subselect, "__siteToProtein")] %>%  as.numeric()
    v_siteToProtein_IFadjust_a <- df_phospho_a[,paste0(samplenames_proteinTable_TcellAct_subselect, "__siteToProtein_IFadjust")] %>%  as.numeric()

    # Prepare dataframe for plotting
    df_gg_phospho_a <- data.frame(x=factor(samplenames_proteinTable_TcellAct_subselect, levels=samplenames_proteinTable_TcellAct_subselect),
                                  groups=factor(groups_proteinTable_TcellAct_subselect, levels=unique(groups_proteinTable_TcellAct_subselect)),
                                  v_phospho_a = v_phospho_a,
                                  v_protein_a = v_protein_a,
                                  v_siteToProtein_a = v_siteToProtein_a,
                                  v_siteToProtein_IFadjust_a = v_siteToProtein_IFadjust_a)
    df_gg_phospho_a$v_phospho_a <- df_gg_phospho_a$v_phospho_a/median(df_gg_phospho_a$v_phospho_a, na.rm=TRUE)
    df_gg_phospho_a$v_protein_a <- df_gg_phospho_a$v_protein_a/median(df_gg_phospho_a$v_protein_a, na.rm=TRUE)
    y_max <- max(c(2.5,
                   df_gg_phospho_a$v_phospho_a + 0.35,
                   df_gg_phospho_a$v_protein_a + 0.35,
                   df_gg_phospho_a$v_siteToProtein_a + 0.35,
                   df_gg_phospho_a$v_siteToProtein_IFadjust_a + 0.35),
                 na.rm=TRUE)

    # prepare ggplot of site level
    gg_acet <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_phospho_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_phospho_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # prepare ggplot of protein level
    gg_prot <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_protein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_protein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("protein level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # prepare ggplot of siteToProtein level
    gg_siteToProt <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # prepare ggplot of siteToProtein level IF-adjusted
    gg_siteToprotein_IFadjust <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_IFadjust_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_IFadjust_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein IF-adjusted ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # build grid and return
    gg_grid <- cowplot::plot_grid(gg_acet, gg_prot, gg_siteToProt, gg_siteToprotein_IFadjust,
                                  nrow=1, ncol=4)
    return(gg_grid)
  })


  # Prepare phospho level IF-adjusted site/protein ratio ANOVA adj.p-value wo Naive
  output$text_phospho_adj.pval_siteToProtein_IFadjust_TcellAct <- renderText({
    a = input$selected_phosphoSite_TcellAct
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_phospho_a$adj.pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })


  # Prepare phospho level IF-adjusted site/protein ratio ANOVA p-value wo Naive
  output$text_phospho_pval_siteToProtein_IFadjust_TcellAct <- renderText({
    a = input$selected_phosphoSite_TcellAct
    df_phospho_a <- df_phospho_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value: ", round(df_phospho_a$pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })



  
  ## Server Acetyl level TcellAct -------------------------------------------------
  
  # Update phospho site choice based on protein selection
  observeEvent(input$selected_protein_TcellAct,
               {
                 protein_id <- df_proteins_TcellAct %>% filter(Gene.names == input$selected_protein_TcellAct) %>% select(id)
                 if (length(acetSite_choices_TcellAct[acetSite_choices_TcellAct_proteins %in% protein_id]) > 0){
                   choices = acetSite_choices_TcellAct[acetSite_choices_TcellAct_proteins %in% protein_id]
                 } else {
                   choices = " "
                 }
                 updateSelectizeInput(session, input = "selected_acetSite_TcellAct", choices = choices)
               })
  
  
  
  # Prepare acet site pattern plot of selected site a
  output$plot_acet_TcellAct <- renderPlot({
    
    # Specify selected site's unique name
    a = input$selected_acetSite_TcellAct
    
    # Filter for site's a information
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)
    
    # Select site intensities of site a
    v_acet_a <- df_acet_a[,paste0(samplenames_proteinTable_TcellAct_subselect)] %>%  as.numeric()
    
    # Prepare dataframe for plotting
    df_gg_acet_a <- data.frame(x=factor(samplenames_proteinTable_TcellAct_subselect, levels=samplenames_proteinTable_TcellAct_subselect),
                                  groups=factor(groups_proteinTable_TcellAct_subselect), levels=unique(groups_proteinTable_TcellAct_subselect),
                                  y=log2(v_acet_a))
    
    # Specify y-range, then prepare plot
    y <- df_gg_acet_a$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_acet <- ggplot(data=df_gg_acet_a) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_time[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(a) +
      ylab("log2 Intensity") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none")
    return(gg_acet)
  })
  
  
  # Prepare acet level ANOVA p-value wo Naive
  output$text_acet_pval_ANOVA_TcellAct <- renderText({
    a = input$selected_acetSite_TcellAct
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value:  ", round(df_acet_a$pval_ANOVA, digits=4)))
  })
  
  
  # Prepare acet level ANOVA adj.p-value wo Naive
  output$text_acet_adj.pval_ANOVA_TcellAct <- renderText({
    a = input$selected_acetSite_TcellAct
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value:  ", round(df_acet_a$adj.pval_ANOVA, digits=4)))
  })
  
  
  # Prepare acet level ANOVA adj.p-value wo Naive
  output$text_acet_PPF_TcellAct <- renderText({
    a = input$selected_acetSite_TcellAct
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("Precursor Purity Fraction (PPF)   :   ", round(df_acet_a$PPF, digits=2)))
  })
  
  
  # Prepare acet level ANOVA adj.p-value wo Naive
  output$text_acet_EIL_TcellAct <- renderText({
    a = input$selected_acetSite_TcellAct
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("Estimated Interference Level (EIL)   :   ", round(df_acet_a$EIL, digits=2)))
  })
  
  
  # Prepare acet site grid plot
  output$plot_grid_acet_siteToProtein_TcellAct <- renderPlot({

    # Specify selected site's unique name
    a = input$selected_acetSite_TcellAct

    # Filter for site's a information
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)

    # Select site intensities of site a
    v_acet_a <- df_acet_a[,paste0(samplenames_proteinTable_TcellAct_subselect)] %>%  as.numeric()
    v_protein_a <- df_acet_a[,paste0(samplenames_proteinTable_TcellAct_subselect, "__underlyingProtein")] %>%  as.numeric()
    v_siteToProtein_a <- df_acet_a[,paste0(samplenames_proteinTable_TcellAct_subselect, "__siteToProtein")] %>%  as.numeric()
    v_siteToProtein_IFadjust_a <- df_acet_a[,paste0(samplenames_proteinTable_TcellAct_subselect, "__siteToProtein_IFadjust")] %>%  as.numeric()

    # Prepare dataframe for plotting
    df_gg_acet_a <- data.frame(x=factor(samplenames_proteinTable_TcellAct_subselect, levels=samplenames_proteinTable_TcellAct_subselect),
                               groups=factor(groups_proteinTable_TcellAct_subselect, levels=unique(groups_proteinTable_TcellAct_subselect)),
                               v_acet_a = v_acet_a,
                               v_protein_a = v_protein_a,
                               v_siteToProtein_a = v_siteToProtein_a,
                               v_siteToProtein_IFadjust_a = v_siteToProtein_IFadjust_a)
    df_gg_acet_a$v_acet_a <- df_gg_acet_a$v_acet_a/median(df_gg_acet_a$v_acet_a, na.rm=TRUE)
    df_gg_acet_a$v_protein_a <- df_gg_acet_a$v_protein_a/median(df_gg_acet_a$v_protein_a, na.rm=TRUE)
    y_max <- max(c(2.5,
                   df_gg_acet_a$v_acet_a + 0.35,
                   df_gg_acet_a$v_protein_a + 0.35,
                   df_gg_acet_a$v_siteToProtein_a + 0.35,
                   df_gg_acet_a$v_siteToProtein_IFadjust_a + 0.35),
                 na.rm=TRUE)

    # prepare ggplot of site level
    gg_acet <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_acet_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_acet_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # prepare ggplot of protein level
    gg_prot <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_protein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_protein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("protein level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # prepare ggplot of siteToProtein level
    gg_siteToProt <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # prepare ggplot of siteToProtein level IF-adjusted
    gg_siteToprotein_IFadjust <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_IFadjust_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_IFadjust_a, col=groups), cex=5) +
      scale_color_manual(values=colors_time[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein IF-adjusted ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")

    # build grid and return
    gg_grid <- cowplot::plot_grid(gg_acet, gg_prot, gg_siteToProt, gg_siteToprotein_IFadjust, nrow=1, ncol=4)
    return(gg_grid)
  })

  # Prepare acet level IF-adjusted site/protein ratio ANOVA adj.p-value wo Naive
  output$text_acet_adj.pval_siteToProtein_IFadjust_TcellAct <- renderText({
    a = input$selected_acetSite_TcellAct
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_acet_a$adj.pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })

  # Prepare acet level IF-adjusted site/protein ratio ANOVA p-value wo Naive
  output$text_acet_pval_siteToProtein_IFadjust_TcellAct <- renderText({
    a = input$selected_acetSite_TcellAct
    df_acet_a <- df_acet_TcellAct %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value: ", round(df_acet_a$pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })

  
  
  ## Server RNA level allSamples  --------------------------------
  output$plot_rna_all <- renderPlot({
    
    # Specify selected gene                                         
    r = input$selected_gene_all
    
    # Filter for gene r information
    df_RNA_all_r <- df_RNA_all %>% filter(Gene.names == r)
    
    # Extract count data of gene r. Log2 transform
    v_RNA_all_r <- df_RNA_all_r[samplenames_RNAtable_all] %>%  as.numeric() 
    v_RNA_all_r <- log2(v_RNA_all_r + 1)
    
    # Create dataframe for plotting
    df_gg_RNA_all_r <- data.frame(x = factor(samplenames_RNAtable_all, levels=samplenames_RNAtable_all),
                                  groups = factor(groups_RNAtable_all, levels=names(colors_groups)),
                                  y = v_RNA_all_r)
    
    # Specify y-range, then prepare plot
    y <- df_gg_RNA_all_r$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_rna_all <- ggplot(data=df_gg_RNA_all_r) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_groups[levels(df_gg_RNA_all_r$groups)]) +
      theme_bw() +
      ggtitle(r) +
      ylab("log2 counts") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none") 
    return(gg_rna_all)
  })
  
  # Prepare RNA level ANOVA p-value all samples
  output$text_rna_pval_ANOVA_all <- renderText({
    r = input$selected_gene_all
    df_rna_r <- df_RNA_all %>% filter(Gene.names == r)
    return(paste0("ANOVA p-value: ", round(df_rna_r$ANOVA_pval, digits=4)))
  })  
  
  # Prepare RNA level ANOVA adj.p-value all samples
  output$text_rna_adj.pval_ANOVA_all <- renderText({
    r = input$selected_gene_all
    df_rna_r <- df_RNA_all %>% filter(Gene.names == r)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_rna_r$ANOVA_adj.pval, digits=4)))
  })  
  
  # Prepare RNA level WT vs KO adj.p-value all samples
  output$text_rna_pval_WTvsKO_all <- renderText({
    r = input$selected_gene_all
    df_rna_r <- df_RNA_all %>% filter(Gene.names == r)
    return(paste0("F-test p-value: ", round(df_rna_r$pval_ANOVA_interaction, digits=4)))
  })  
  
  # Prepare RNA level WT vs KO adj.p-value all samples
  output$text_rna_adj.pval_WTvsKO_all <- renderText({
    r = input$selected_gene_all
    df_rna_r <- df_RNA_all %>% filter(Gene.names == r)
    return(paste0("F-test FDR-adj. p-value: ", round(df_rna_r$adj.pval_ANOVA_interaction, digits=4)))
  })  
  
  
  
  
  
  ## Server RNA level woNaive --------------------------------
  output$plot_rna_woNaive <- renderPlot({
    
    # Specify selected gene                                         
    r = input$selected_gene_woNaive
    
    # Filter for gene r information
    df_RNA_woNaive_r <- df_RNA_woNaive %>% filter(Gene.names == r)
    
    # Extract count data of gene r. Log2 transform
    v_RNA_woNaive_r <- df_RNA_woNaive_r[samplenames_RNAtable_woNaive] %>%  as.numeric() 
    v_RNA_woNaive_r <- log2(v_RNA_woNaive_r + 1)
    
    # Create dataframe for plotting
    df_gg_RNA_woNaive_r <- data.frame(x = factor(samplenames_RNAtable_woNaive, levels=samplenames_RNAtable_woNaive),
                                      groups = factor(groups_RNAtable_woNaive, levels=names(colors_groups)),
                                      y = v_RNA_woNaive_r)
    
    # Specify y-range, then prepare plot
    y <- df_gg_RNA_woNaive_r$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_rna_woNaive <- ggplot(data=df_gg_RNA_woNaive_r) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_groups[levels(df_gg_RNA_woNaive_r$groups)]) +
      theme_bw() +
      ggtitle(r) +
      ylab("log2 counts") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none") 
    return(gg_rna_woNaive)
  })
  
  # Prepare RNA level ANOVA p-value all samples
  output$text_rna_pval_ANOVA_woNaive <- renderText({
    r = input$selected_gene_woNaive
    df_rna_r <- df_RNA_woNaive %>% filter(Gene.names == r)
    return(paste0("ANOVA p-value: ", round(df_rna_r$ANOVA_pval, digits=4)))
  })  
  
  # Prepare RNA level ANOVA adj.p-value all samples
  output$text_rna_adj.pval_ANOVA_woNaive <- renderText({
    r = input$selected_gene_woNaive
    df_rna_r <- df_RNA_woNaive %>% filter(Gene.names == r)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_rna_r$ANOVA_adj.pval, digits=4)))
  })  
  
  # Prepare RNA level WT vs KO adj.p-value all samples
  output$text_rna_pval_WTvsKO_woNaive <- renderText({
    r = input$selected_gene_woNaive
    df_rna_r <- df_RNA_woNaive %>% filter(Gene.names == r)
    return(paste0("F-test p-value: ", round(df_rna_r$pval_ANOVA_interaction, digits=4)))
  })  
  
  # Prepare RNA level WT vs KO adj.p-value all samples
  output$text_rna_adj.pval_WTvsKO_woNaive <- renderText({
    r = input$selected_gene_woNaive
    df_rna_r <- df_RNA_woNaive %>% filter(Gene.names == r)
    return(paste0("F-test FDR-adj. p-value: ", round(df_rna_r$adj.pval_ANOVA_interaction, digits=4)))
  })  
  
  
  
  
  
  
  
  ## Server Protein level TMTset1 woNaive  -------------------------------
  
  # Prepare protein pattern plot of selected protein i
  output$plot_protein_TMTset1_woNaive <- renderPlot({
    
    # Specify selected protein's gene name
    i = input$selected_protein_TMTset1_woNaive
    
    # Filter for protein i information
    df_proteins_i <- df_proteins_TMTset1 %>% filter(Gene.names == i)
    
    # Selevt intensity values of protein i
    v_protein_i <- df_proteins_i[samplenames_proteinTable_TMTset1_woNaive] %>% as.numeric()
    
    # Create dataframe for plotting
    df_gg_protein_i <- data.frame(x = factor(samplenames_proteinTable_TMTset1_woNaive, levels=samplenames_proteinTable_TMTset1_woNaive),
                                          groups=factor(groups_proteinTable_TMTset1_woNaive, levels=unique(groups_proteinTable_TMTset1_woNaive)),
                                          y = v_protein_i)
    
    # Specify y-range, then prepare plot
    y <- df_gg_protein_i$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_prot <- ggplot(data=df_gg_protein_i) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_groups[levels(df_gg_protein_i$groups)]) +
      theme_bw() +
      ggtitle(i) +
      ylab("log2 Intensity") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none") 
    return(gg_prot)
  })
  
  # Prepare protein level ANOVA p-value wo Naive
  output$text_protein_pval_ANOVA_TMTset1_woNaive <- renderText({
    i = input$selected_protein_TMTset1_woNaive
    df_proteins_i <- df_proteins_TMTset1 %>% filter(Gene.names == i)
    return(paste0("ANOVA p-value: ", round(df_proteins_i$pval_ANOVA_woNaiveGroup, digits=4)))
  })
  
  # Prepare protein level ANOVA adj.p-value wo Naive
  output$text_protein_adj.pval_ANOVA_TMTset1_woNaive <- renderText({
    i = input$selected_protein_TMTset1_woNaive
    df_proteins_i <- df_proteins_TMTset1 %>% filter(Gene.names == i)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_proteins_i$adj.pval_ANOVA_woNaiveGroup, digits=4)))
  })
  
  # Prepare nr acetyl sites Tcell Act
  output$text_protein_nrAcetSites_TMTset1_woNaive <- renderText({
    i = input$selected_protein_TMTset1_woNaive
    proteinID_i <- df_proteins_TMTset1 %>% filter(Gene.names == i) %>% select(id) %>% as.numeric()
    nrAcetSites <- sum(acetSite_choices_TMTset1_woNaive_proteins == proteinID_i)
    return(paste0("Acetylation sites: ", round(nrAcetSites, digits=0)))
  })
  
  # Prepare number phospho sites TcellAct
  output$text_protein_nrPhosphoSites_TMTset1_woNaive <- renderText({
    i = input$selected_protein_TMTset1_woNaive
    proteinID_i <- df_proteins_TMTset1 %>% filter(Gene.names == i) %>% select(id) %>% as.numeric()
    nrPhosphoSites <- sum(phosphoSite_choices_TMTset1_woNaive_proteins == proteinID_i)
    return(paste0("Phosphorylation sites: ", round(nrPhosphoSites, digits=0)))
  })
  
  
  
  
  
  
  ## Server Phospho level TMTset1 woNaive  -------------------------------
  
  # Update phospho site choice based on protein selection
  observeEvent(input$selected_protein_TMTset1_woNaive, 
               {
                 protein_id <- df_proteins_TMTset1   %>% filter(Gene.names == input$selected_protein_TMTset1_woNaive) %>% select(id)
                 if (length(phosphoSite_choices_TMTset1_woNaive [phosphoSite_choices_TMTset1_woNaive_proteins %in% protein_id]) > 0){
                   choices = phosphoSite_choices_TMTset1_woNaive [phosphoSite_choices_TMTset1_woNaive_proteins %in% protein_id]
                 } else {
                   choices = " "
                 }
                 updateSelectizeInput(session, input = "selected_phosphoSite_TMTset1_woNaive", choices = choices)
               })
  
  
  # Prepare phospho site pattern plot of selected site a
  output$plot_phospho_TMTset1_woNaive <- renderPlot({
    
    # Specify selected site's unique name
    a = input$selected_phosphoSite_TMTset1_woNaive
    
    # Filter for site's a information
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    
    # Select site intensities of site a
    v_phospho_a <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive)] %>%  as.numeric()
    
    # Prepare dataframe for plotting
    df_gg_phospho_a <- data.frame(x=factor(samplenames_proteinTable_TMTset1_woNaive, levels=samplenames_proteinTable_TMTset1_woNaive),
                                  groups=factor(groups_proteinTable_TMTset1_woNaive, levels=unique(groups_proteinTable_TMTset1_woNaive)),
                                  y=log2(v_phospho_a))
    
    # Specify y-range, then prepare plot
    y <- df_gg_phospho_a$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_phospho <- ggplot(data=df_gg_phospho_a) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(a) +
      ylab("log2 Intensity") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none") 
    return(gg_phospho)
  })
  
  
  # Prepare phospho level ANOVA p-value wo Naive
  output$text_phospho_pval_ANOVA_TMTset1_woNaive <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value:  ", round(df_phospho_a$pval_ANOVA, digits=4)))
  })
  
  
  # Prepare phospho level ANOVA adj.p-value wo Naive
  output$text_phospho_adj.pval_ANOVA_TMTset1_woNaive <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value:  ", round(df_phospho_a$adj.pval_ANOVA, digits=4)))
  })
  
  
  # Prepare phospho level ANOVA adj.p-value wo Naive
  output$text_phospho_PPF <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("Precursor Purity Fraction (PPF)   :   ", round(df_phospho_a$PPF, digits=2)))
  })
  
  
  # Prepare phospho level ANOVA adj.p-value wo Naive
  output$text_phospho_EIL <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("Estimated Interference Level (EIL)   :   ", round(df_phospho_a$EIL, digits=2)))
  })
  
  
  # Prepare phospho site grid plot
  output$plot_grid_phospho_siteToProtein_TMTset1_woNaive <- renderPlot({
    
    # Specify selected site's unique name
    a = input$selected_phosphoSite_TMTset1_woNaive
    
    # Filter for site's a information
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    
    # Select site intensities of site a
    v_phospho_a <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive)] %>%  as.numeric()
    v_protein_a <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__underlyingProtein")] %>%  as.numeric()
    v_siteToProtein_a <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein")] %>%  as.numeric()
    v_siteToProtein_IFadjust_a <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein_IFadjust")] %>%  as.numeric()
    
    # Prepare dataframe for plotting
    df_gg_phospho_a <- data.frame(x=factor(samplenames_proteinTable_TMTset1_woNaive, levels=samplenames_proteinTable_TMTset1_woNaive),
                                  groups=factor(groups_proteinTable_TMTset1_woNaive, levels=unique(groups_proteinTable_TMTset1_woNaive)),
                                  v_phospho_a = v_phospho_a,
                                  v_protein_a = v_protein_a,
                                  v_siteToProtein_a = v_siteToProtein_a,
                                  v_siteToProtein_IFadjust_a = v_siteToProtein_IFadjust_a)
    df_gg_phospho_a$v_phospho_a <- df_gg_phospho_a$v_phospho_a/median(df_gg_phospho_a$v_phospho_a, na.rm=TRUE)
    df_gg_phospho_a$v_protein_a <- df_gg_phospho_a$v_protein_a/median(df_gg_phospho_a$v_protein_a, na.rm=TRUE)
    y_max <- max(c(2.5,  
                   df_gg_phospho_a$v_phospho_a + 0.35,
                   df_gg_phospho_a$v_protein_a + 0.35,
                   df_gg_phospho_a$v_siteToProtein_a + 0.35,
                   df_gg_phospho_a$v_siteToProtein_IFadjust_a + 0.35),
                 na.rm=TRUE)
    
    # prepare ggplot of site level
    gg_acet <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_phospho_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_phospho_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of protein level
    gg_prot <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_protein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_protein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("protein level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level
    gg_siteToProt <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level IF-adjusted
    gg_siteToprotein_IFadjust <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_IFadjust_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_IFadjust_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein IF-adjusted ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # build grid and return
    gg_grid <- cowplot::plot_grid(gg_acet, gg_prot, gg_siteToProt, gg_siteToprotein_IFadjust,
                                  nrow=1, ncol=4)
    return(gg_grid)
  })
  
  
  # Prepare phospho level IF-adjusted site/protein ratio ANOVA adj.p-value wo Naive
  output$text_phospho_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_phospho_a$adj.pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })
  
  
  # Prepare phospho level IF-adjusted site/protein ratio ANOVA p-value wo Naive
  output$text_phospho_pval_siteToProtein_IFadjust_TMTset1_woNaive <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value: ", round(df_phospho_a$pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })
  
  
  # Prepare phospho site-to-protein grid plot (batch-corrected)
  output$plot_grid_phospho_siteToProtein_TMTset1_woNaive_batchCorr <- renderPlot({
    
    # Specify selected site's unique name
    a = input$selected_phosphoSite_TMTset1_woNaive
    
    # Filter for site's a information
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    
    # Select site intensities of site a
    v_phospho_a_batchCorr <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__batchCorr")] %>%  as.numeric()
    v_protein_a_batchCorr <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__underlyingProtein__batchCorr")] %>%  as.numeric()
    v_siteToProtein_a_batchCorr <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein__batchCorr")] %>%  as.numeric()
    v_siteToProtein_IFadjust_a_batchCorr <- df_phospho_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein_IFadjust__batchCorr")] %>%  as.numeric()
    
    # Prepare dataframe for plotting
    df_gg_phospho_a <- data.frame(x=factor(samplenames_proteinTable_TMTset1_woNaive, levels=samplenames_proteinTable_TMTset1_woNaive),
                                  groups=factor(groups_proteinTable_TMTset1_woNaive, levels=unique(groups_proteinTable_TMTset1_woNaive)),
                                  v_phospho_a_batchCorr = v_phospho_a_batchCorr,
                                  v_protein_a_batchCorr = v_protein_a_batchCorr,
                                  v_siteToProtein_a_batchCorr = v_siteToProtein_a_batchCorr,
                                  v_siteToProtein_IFadjust_a_batchCor = v_siteToProtein_IFadjust_a_batchCorr)
    df_gg_phospho_a$v_phospho_a_batchCorr <- df_gg_phospho_a$v_phospho_a_batchCorr/median(df_gg_phospho_a$v_phospho_a_batchCorr, na.rm = TRUE)
    df_gg_phospho_a$v_protein_a_batchCorr <- df_gg_phospho_a$v_protein_a_batchCorr/median(df_gg_phospho_a$v_protein_a_batchCorr, na.rm = TRUE)
    y_max <- max(c(2.5, 
                   df_gg_phospho_a$v_phospho_a_batchCorr + 0.35,
                   df_gg_phospho_a$v_protein_a_batchCorr + 0.35,
                   df_gg_phospho_a$v_siteToProtein_a_batchCorr + 0.35,
                   df_gg_phospho_a$v_siteToProtein_IFadjust_a_batchCor + 0.35),
                 na.rm=TRUE)
    
    # prepare ggplot of site level (batch-corrected)
    gg_phospho_batchCorr <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x, y=v_phospho_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x, y=v_phospho_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site level")) +
      ylab("relative intensity (batch-corrected)") +
      xlab("") +
      ylim(-0.1,y_max) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=9, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of protein level (batch-corrected)
    gg_protein_batchCorr <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x, y=v_protein_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x, y=v_protein_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("protein level")) +
      ylab("relative intensity (batch-corrected)") +
      xlab("") +
      ylim(-0.1,y_max) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=9, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level (batch-corrected)
    gg_siteToprotein_batchCorr <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein ratio")) +
      ylab("ratio (batch-corrected)") +
      ylim(-0.1,y_max) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=7.5),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level IF-adjusted (batch-corrected)
    gg_siteToprotein_IFadjust_batchCorr <- ggplot(data=df_gg_phospho_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_IFadjust_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_IFadjust_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_phospho_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein IF-adjusted ratio")) +
      ylab("ratio (batch-corrected)") +
      ylim(-0.1,y_max) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=7.5),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # build grid and return
    gg_grid <- cowplot::plot_grid(gg_phospho_batchCorr, gg_protein_batchCorr, gg_siteToprotein_batchCorr, gg_siteToprotein_IFadjust_batchCorr,
                                  nrow=1, ncol=4)
    return(gg_grid)
  })
  
  
  # Prepare phospho level batch-corrected IF-adjusted site/protein ratio ANOVA adj.p-value wo Naive
  output$text_phospho_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive_batchCorr <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_phospho_a$adj.pval_ANOVA_siteToProtein_IFadjust__batchCorr, digits=4)))
  })
  
  
  # Prepare phospho level batch-correctedIF-adjusted site/protein ratio ANOVA p-value wo Naive
  output$text_phospho_pval_siteToProtein_IFadjust_TMTset1_woNaive_batchCorr <- renderText({
    a = input$selected_phosphoSite_TMTset1_woNaive
    df_phospho_a <- df_phospho_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value: ", round(df_phospho_a$pval_ANOVA_siteToProtein_IFadjust__batchCorr, digits=4)))
  })
  
  
  
  
  
  
  ## Server Acetyl level TMTset1 woNaive -------------------------------------------------
  
  # Update acetyl site choice based on protein selection
  observeEvent(input$selected_protein_TMTset1_woNaive, 
               {
                 protein_id <- df_proteins_TMTset1   %>% filter(Gene.names == input$selected_protein_TMTset1_woNaive) %>% select(id)
                 if (length(acetSite_choices_TMTset1_woNaive [acetSite_choices_TMTset1_woNaive_proteins %in% protein_id]) > 0){
                   choices = acetSite_choices_TMTset1_woNaive [acetSite_choices_TMTset1_woNaive_proteins %in% protein_id]
                 } else {
                   choices = " "
                 }
                 updateSelectizeInput(session, input = "selected_acetSite_TMTset1_woNaive", choices = choices)
               })
  
  # Prepare acet site pattern plot of selected site a
  output$plot_acet_TMTset1_woNaive <- renderPlot({
    
    # Specify selected site's unique name
    a = input$selected_acetSite_TMTset1_woNaive
    
    # Filter for site's a information
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    
    # Select site intensities of site a
    v_acet_a <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive)] %>%  as.numeric()
    
    # Prepare dataframe for plotting
    df_gg_acet_a <- data.frame(x=factor(samplenames_proteinTable_TMTset1_woNaive, levels=samplenames_proteinTable_TMTset1_woNaive),
                               groups=factor(groups_proteinTable_TMTset1_woNaive, levels=unique(groups_proteinTable_TMTset1_woNaive)),
                               y=log2(v_acet_a))
    
    # Specify y-range, then prepare plot
    y <- df_gg_acet_a$y
    span_original <- abs(diff(range(y, na.rm=TRUE)))
    if (span_original > 1){
      y_lim <- range(y, na.rm = TRUE) + c(-0.25, 0.25)
    } else {
      y_lim <- median(y, na.rm = TRUE) + c(-0.75, 0.75)
    }
    gg_acet <- ggplot(data=df_gg_acet_a) +
      geom_line(aes(x=x,y=y, group=1), color="black") +
      geom_point(aes(x=x,y=y, col=groups), cex=10) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(a) +
      ylab("log2 Intensity") +
      xlab("") +
      ylim(y_lim[1], y_lim[2]) +
      theme(plot.title = element_text(color="black", size=16.5, face="bold"),
            axis.title.y = element_text(color="black", size=15, face="bold"),
            axis.text.y =element_text(size=13),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=13),
            legend.position="none") 
    return(gg_acet)
  })
  
  
  # Prepare acet level ANOVA p-value wo Naive
  output$text_acet_pval_ANOVA_TMTset1_woNaive <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value:  ", round(df_acet_a$pval_ANOVA, digits=4)))
  })
  
  
  # Prepare acet level ANOVA adj.p-value wo Naive
  output$text_acet_adj.pval_ANOVA_TMTset1_woNaive <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value:  ", round(df_acet_a$adj.pval_ANOVA, digits=4)))
  })
  
  
  # Prepare acet level ANOVA adj.p-value wo Naive
  output$text_acet_PPF <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("Precursor Purity Fraction (PPF):   ", round(df_acet_a$PPF, digits=2)))
  })
  
  
  # Prepare acet level ANOVA adj.p-value wo Naive
  output$text_acet_EIL <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("Estimated Interference Level (EIL):   ", round(df_acet_a$EIL, digits=2)))
  })
  
  
  # Prepare acet site grid plot
  output$plot_grid_acet_siteToProtein_TMTset1_woNaive <- renderPlot({
    
    # Specify selected site's unique name
    a = input$selected_acetSite_TMTset1_woNaive
    
    # Filter for site's a information
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    
    # Select site intensities of site a
    v_acet_a <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive)] %>%  as.numeric()
    v_protein_a <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__underlyingProtein")] %>%  as.numeric()
    v_siteToProtein_a <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein")] %>%  as.numeric()
    v_siteToProtein_IFadjust_a <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein_IFadjust")] %>%  as.numeric()
    
    # Prepare dataframe for plotting
    df_gg_acet_a <- data.frame(x=factor(samplenames_proteinTable_TMTset1_woNaive, levels=samplenames_proteinTable_TMTset1_woNaive),
                               groups=factor(groups_proteinTable_TMTset1_woNaive, levels=unique(groups_proteinTable_TMTset1_woNaive)),
                               v_acet_a = v_acet_a,
                               v_protein_a = v_protein_a,
                               v_siteToProtein_a = v_siteToProtein_a,
                               v_siteToProtein_IFadjust_a = v_siteToProtein_IFadjust_a)
    df_gg_acet_a$v_acet_a <- df_gg_acet_a$v_acet_a/median(df_gg_acet_a$v_acet_a, na.rm=TRUE)
    df_gg_acet_a$v_protein_a <- df_gg_acet_a$v_protein_a/median(df_gg_acet_a$v_protein_a, na.rm=TRUE)
    y_max <- max(c(2.5,  
                   df_gg_acet_a$v_acet_a + 0.35,
                   df_gg_acet_a$v_protein_a + 0.35,
                   df_gg_acet_a$v_siteToProtein_a + 0.35,
                   df_gg_acet_a$v_siteToProtein_IFadjust_a + 0.35),
                 na.rm=TRUE)
    
    # prepare ggplot of site level
    gg_acet <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_acet_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_acet_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of protein level
    gg_prot <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_protein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_protein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("protein level")) +
      ylab("relative intensity") +
      xlab("") +
      ylim(-0.1,y_max ) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level
    gg_siteToProt <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level IF-adjusted
    gg_siteToprotein_IFadjust <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_IFadjust_a, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_IFadjust_a, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein IF-adjusted ratio")) +
      ylab("ratio") +
      ylim(-0.1,y_max ) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # build grid and return
    gg_grid <- cowplot::plot_grid(gg_acet, gg_prot, gg_siteToProt, gg_siteToprotein_IFadjust, nrow=1, ncol=4)
    return(gg_grid)
  })
  
  
  # Prepare acet level IF-adjusted site/protein ratio ANOVA adj.p-value wo Naive
  output$text_acet_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_acet_a$adj.pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })
  
  
  # Prepare acet level IF-adjusted site/protein ratio ANOVA p-value wo Naive
  output$text_acet_pval_siteToProtein_IFadjust_TMTset1_woNaive <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value: ", round(df_acet_a$pval_ANOVA_siteToProtein_IFadjust, digits=4)))
  })
  
  
  # Prepare acet site-to-protein grid plot (batch-corrected)
  output$plot_grid_acet_siteToProtein_TMTset1_woNaive_batchCorr <- renderPlot({
    
    # Specify selected site's unique name
    a = input$selected_acetSite_TMTset1_woNaive
    
    # Filter for site's a information
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    
    # Select site intensities of site a
    v_acet_a_batchCorr <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__batchCorr")] %>%  as.numeric()
    v_protein_a_batchCorr <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__underlyingProtein__batchCorr")] %>%  as.numeric()
    v_siteToProtein_a_batchCorr <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein__batchCorr")] %>%  as.numeric()
    v_siteToProtein_IFadjust_a_batchCorr <- df_acet_a[,paste0(samplenames_proteinTable_TMTset1_woNaive, "__siteToProtein_IFadjust__batchCorr")] %>%  as.numeric()
    
    # Prepare dataframe for plotting
    df_gg_acet_a <- data.frame(x=factor(samplenames_proteinTable_TMTset1_woNaive, levels=samplenames_proteinTable_TMTset1_woNaive),
                               groups=factor(groups_proteinTable_TMTset1_woNaive, levels=unique(groups_proteinTable_TMTset1_woNaive)),
                               v_acet_a_batchCorr = v_acet_a_batchCorr,
                               v_protein_a_batchCorr = v_protein_a_batchCorr,
                               v_siteToProtein_a_batchCorr = v_siteToProtein_a_batchCorr,
                               v_siteToProtein_IFadjust_a_batchCor = v_siteToProtein_IFadjust_a_batchCorr)
    df_gg_acet_a$v_acet_a_batchCorr <- df_gg_acet_a$v_acet_a_batchCorr/median(df_gg_acet_a$v_acet_a_batchCorr, na.rm = TRUE)
    df_gg_acet_a$v_protein_a_batchCorr <- df_gg_acet_a$v_protein_a_batchCorr/median(df_gg_acet_a$v_protein_a_batchCorr, na.rm = TRUE)
    y_max <- max(c(2.5, 
                   df_gg_acet_a$v_acet_a_batchCorr + 0.35,
                   df_gg_acet_a$v_protein_a_batchCorr + 0.35,
                   df_gg_acet_a$v_siteToProtein_a_batchCorr + 0.35,
                   df_gg_acet_a$v_siteToProtein_IFadjust_a_batchCor + 0.35),
                 na.rm=TRUE)
    
    # prepare ggplot of site level (batch-corrected)
    gg_acet_batchCorr <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x, y=v_acet_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x, y=v_acet_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site level")) +
      ylab("relative intensity (batch-corrected)") +
      xlab("") +
      ylim(-0.1,y_max) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=9, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of protein level (batch-corrected)
    gg_protein_batchCorr <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x, y=v_protein_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x, y=v_protein_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("protein level")) +
      ylab("relative intensity (batch-corrected)") +
      xlab("") +
      ylim(-0.1,y_max) +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=9, face="bold"),
            axis.text.y =element_text(size=8),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level (batch-corrected)
    gg_siteToprotein_batchCorr <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein ratio")) +
      ylab("ratio (batch-corrected)") +
      ylim(-0.1,y_max) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=7.5),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # prepare ggplot of siteToProtein level IF-adjusted (batch-corrected)
    gg_siteToprotein_IFadjust_batchCorr <- ggplot(data=df_gg_acet_a) +
      geom_hline(yintercept = 1, linetype="dashed", alpha=0.3) +
      geom_line(aes(x=x,y=v_siteToProtein_IFadjust_a_batchCorr, group=1), color="black") +
      geom_point(aes(x=x,y=v_siteToProtein_IFadjust_a_batchCorr, col=groups), cex=5) +
      scale_color_manual(values=colors_groups[levels(df_gg_acet_a$groups)]) +
      theme_bw() +
      ggtitle(paste0("site/protein IF-adjusted ratio")) +
      ylab("ratio (batch-corrected)") +
      ylim(-0.1,y_max) +
      xlab("") +
      theme(plot.title = element_text(color="black", size=12, face="bold"),
            axis.title.y = element_text(color="black", size=10, face="bold"),
            axis.text.y =element_text(size=7.5),
            axis.text.x = element_text(angle = 90, vjust = 0.6, hjust=0.6,size=8),
            legend.position="none")
    
    # build grid and return
    gg_grid <- cowplot::plot_grid(gg_acet_batchCorr, gg_protein_batchCorr, gg_siteToprotein_batchCorr, gg_siteToprotein_IFadjust_batchCorr, nrow=1, ncol=4)
    return(gg_grid)
  })
  
  
  # Prepare acet level batch-corrected IF-adjusted site/protein ratio ANOVA adj.p-value wo Naive
  output$text_acet_adj.pval_siteToProtein_IFadjust_TMTset1_woNaive_batchCorr <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA FDR-adj. p-value: ", round(df_acet_a$adj.pval_ANOVA_siteToProtein_IFadjust__batchCorr, digits=4)))
  })
  
  
  # Prepare acet level batch-correctedIF-adjusted site/protein ratio ANOVA p-value wo Naive
  output$text_acet_pval_siteToProtein_IFadjust_TMTset1_woNaive_batchCorr <- renderText({
    a = input$selected_acetSite_TMTset1_woNaive
    df_acet_a <- df_acet_TMTset1_woNaive %>% filter(unique_site_identifier == a)
    return(paste0("ANOVA p-value: ", round(df_acet_a$pval_ANOVA_siteToProtein_IFadjust__batchCorr, digits=4)))
  })
  
  
}

