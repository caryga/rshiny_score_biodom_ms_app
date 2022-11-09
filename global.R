#
# This application enables biological domain characterization by users and
# plots information about specific scored targets using the data collected by 
# the TREAT-AD Emory-Sage-SGC center, including:
#   - Overall, Genetics, and Multi-Omics scores
#   - Genetic association score evidence
#   - Meta-analysis predictions of RNA and Protein expression changes
#   - AD Biodomain annotations 
#

# packages ----------------------------------------------------------------

library(synapser)
library(network)
library(sna)
library(shiny)
library(shinythemes)
library(DT)
library(plotly)
library('org.Hs.eg.db')
library(clusterProfiler)
## TODO: enable other enrichment tools, e.g. gprofiler2
# library(gprofiler2) 
library(tidyverse)

# settings ----------------------------------------------------------------

# graphics theme
theme_set(theme_bw())
theme_update( axis.text = element_text(size = 12) )

# data --------------------------------------------------------------------

### load data from synapse ----
# synLogin()
# 
# ##
# # Example Datasets for Enrichment
# ##
# 
# # TMT Proteomics Modules
# tmt.mods <- readxl::read_xlsx(synGet('syn28551982')$path, sheet = 4, skip = 2) %>%
#   mutate(gene = str_split_fixed(UniqueID, '\\|',2)[,1],
#          mod = str_split_fixed(kMEtableSortVector, '\\|', 3) %>% as_tibble()) %>%
#   select(gene, mod) %>%
#   unnest(mod) %>%
#   rename_with(~c('gene','mod','col','eigen'), everything())
# 
# consens.clust <- tibble(
#   consens = c('ConsensusClusterA',
#               'ConsensusClusterB',
#               'ConsensusClusterC',
#               'ConsensusClusterD',
#               'ConsensusClusterE' ),
#   mod = list(c('IFGyellow', 'PHGyellow', 'TCXblue'),
#              c('CBEturquoise','DLPFCblue','FPturquoise','IFGturquoise','PHGturquoise','STGblue','TCXturquoise'),
#              c('CBEyellow','DLPFCyellow','FPyellow','IFGbrown','PHGbrown','STGbrown','TCXgreen'),
#              c('CBEbrown','DLPFCbrown','FPblue','IFGblue','PHGgreen','STGyellow','TCXyellow'),
#              c('CBEblue','DLPFCturquoise','FPbrown','PHGblue','STGturquoise','TCXbrown'))) %>%
#   unnest(mod)
# 
# # AMP-AD transcriptomic modules
# tx.mods <- read_csv( synapser::synTableQuery('SELECT * FROM syn11932957')$filepath, show_col_types = F) %>%
#   select(-ROW_ID, -ROW_VERSION) %>% 
#   left_join(., consens.clust, by=c('Module'='mod')) %>%
#   relocate(consens, .after = Module) %>% 
#   select(gene = external_gene_name, mod = Module, consens, ensg = GeneID) %>% distinct()
# 
# # AMP-AD transcriptomic submodules
# tx.submods <- bind_rows(
#   readRDS( synapser::synGet('syn23660898')$path), #Mayo
#   readRDS( synapser::synGet('syn23660900')$path), #MSSM
#   readRDS( synapser::synGet('syn23660902')$path)  #ROSMAP
# ) %>%
#   mutate(mod = str_split_fixed(Submodule, '_', 2) %>% as.data.frame() %>% pull(V1) ) %>%
#   relocate(mod, .before= Submodule) %>%
#   inner_join(consens.clust, ., by = c('mod' = 'mod')) %>% 
#   left_join(., tx.mods %>% select(ensg, gene), by = c('Gene' = 'ensg')) %>% 
#   select(gene, mod, submod = Submodule, consens, ensg = Gene) %>% distinct()
#   
# ##
# # Biodomain definitions
# ##
# 
# biodom <- full_join(
#   # biodomains
#   readRDS(synGet('syn25428992')$path),
#   # domain labels
#   read_csv(synGet('syn26856828')$path, show_col_types = F),
#   by = c('Biodomain'='domain')
#   ) %>%
#   mutate(Biodomain = case_when(Biodomain == 'none' ~ NA_character_, T ~ Biodomain))
# 
# ## TODO: Mmus defs? syn26592124
# 
# ##
# # TAD scores
# ##
# scores <-  synapser::synTableQuery('select * from syn25575156')$filepath %>% read_csv() %>%
#     select(ENSG, GeneName, Overall, rank= Overall_rank,
#            Omics = OmicsScore, Genetics = GeneticsScore
#            # , Literature = LiteratureScore, Neuropath = NeuropathScore
#            ) %>%
#     pivot_longer(cols = c(Overall, Omics, Genetics
#                           # , Literature, Neuropath
#                           )) %>%
#     mutate(name = fct_relevel(name, c(
#       # 'Literature',
#       'Omics', 'Genetics','Overall')))
# 
# ##
# # Gene-Biodomain mappings
# ##
# biodom_genes <- readRDS(  synapser::synGet('syn25428992')$path  ) %>%
#     select(Biodomain, GO_ID, symbol = hgnc_symbol) %>% unnest_longer(symbol) %>% filter(!is.na(symbol), symbol != '') %>%
#     group_by(symbol, Biodomain) %>% summarise(n_term = length(unique(GO_ID))) %>% ungroup() %>%
#     left_join(.,
#               readRDS(  synapser::synGet('syn25428992')$path  ) %>%
#                   select(Biodomain, GO_ID, symbol=hgnc_symbol) %>%
#                   group_by(Biodomain) %>% summarise(bd_terms = length(unique(GO_ID))) %>% ungroup(),
#               by = 'Biodomain') %>%
#     mutate(pct = 100*(n_term / bd_terms) ) %>%
#     left_join(.,
#               readRDS(  synapser::synGet('syn25428992')$path  ) %>%
#                 select(Biodomain, symbol = hgnc_symbol, n_symbol = n_hgncSymbol, GO_ID, GOterm_Name ) %>%
#                 unnest_longer(symbol),
#               by = c('symbol','Biodomain')
#               )
# dom.cols <- read_csv( synapser::synGet('syn26856828')$path )
# 
# ##
# # TAD DE meta-analyses
# ##
# omics <- read_csv( synapser::synTableQuery('select * from syn22758536')$filepath ) %>%
#     select(GName, RNA_TE, RNA_fdr_CorPVal, Pro_TE, Pro_fdr_CorPVal) %>%
#     pivot_longer(
#         c(RNA_TE, Pro_TE), names_sep = '_', names_to = c('x', 'y'), values_to = 'TE'
#     ) %>%
#     pivot_longer(
#         c(RNA_fdr_CorPVal, Pro_fdr_CorPVal), names_sep = '_', names_to = c('a', 'b'), values_to = 'fdr'
#     ) %>%
#     distinct() %>% filter(x == a) %>% select(-y, -a, -b) %>%
#     mutate(x = str_replace_all(x, 'Pro','Protein'))
# 
# ##
# # TAD Genetic Evidence
# ##
# gen <- synapser::synTableQuery('SELECT * FROM syn26844312')$filepath %>% read_csv() %>%
#     filter(!is.na(GeneName), !duplicated(GeneName)) %>%
#     mutate(min_gwasP = case_when(min_gwasP == 0 ~ 1e-100,
#                                  min_gwasP == 1 ~ NA_real_,
#                                  T~min_gwasP),
#            min_gwasP = -log10(min_gwasP),
#            min_gwasP_rank = rank( min_gwasP , na.last = 'keep')/length(which(.$min_gwasP != 1)),
#            min_qtlFDR = case_when(min_qtlFDR == 0 ~ 1e-100,
#                                   min_qtlFDR == 1 ~ NA_real_,
#                                   T~min_qtlFDR),
#            min_qtlFDR = -log10(min_qtlFDR),
#            min_qtlFDR_rank = rank( min_qtlFDR , na.last = 'keep') /length(which(!is.na(min_qtlFDR)))
#            )%>%
#     select(GeneName, score_rank,
#            min_gwasP_rank, #min_gwasP, meanRank_gwasP,
#            min_qtlFDR_rank, #min_qtlFDR, meanRank_qtlFDR,
#            coding_variant_summary, noncoding_variant_summary,
#            Hsap_pheno_score, Ortholog_pheno_score) %>%
#     pivot_longer(-c(GeneName,score_rank)) %>%
#     mutate( name = fct_relevel(name,
#                                c('Ortholog_pheno_score', 'Hsap_pheno_score',
#                                  'noncoding_variant_summary','coding_variant_summary',
#                                        'min_qtlFDR_rank','min_gwasP_rank')),
#             value = case_when( value == 0 ~ NA_real_, T ~ value)) %>%
#     arrange( name )
# 
# ##
# # mutate scores to change plotting order and switch in NA values for 0's
# ##
# scores <- scores %>%
#     mutate( name = fct_relevel(name, c(
#       # 'Neuropath', 'Literature',
#       'Omics','Genetics','Overall')),
#             value = case_when( value == 0 ~ NA_real_, T ~ value)) %>%
#     arrange( name )
# 
# # save.image('appData.RData')
# saveRDS(biodom, 'data/biodom.rds')
# saveRDS(biodom_genes, 'data/biodom_genes.rds')
# saveRDS(dom.cols, 'data/dom_cols.rds')
# saveRDS(scores, 'data/scores.rds')
# saveRDS(gen, 'data/gen.rds')
# saveRDS(omics, 'data/omics.rds')
# saveRDS(consens.clust, 'data/consens_clust.rds')
# saveRDS(tx.mods, 'data/tx_mods.rds')
# saveRDS(tx.submods, 'data/tx_submods.rds')
# saveRDS(tmt.mods, 'data/tmt_mods.rds')

### load pre-compiled data ----
# load('appData.RData')

biodom <- readRDS('data/biodom.rds')
biodom_genes <- readRDS('data/biodom_genes.rds')
dom.cols <- readRDS('data/dom_cols.rds')
scores <- readRDS('data/scores.rds')
gen <- readRDS( 'data/gen.rds')
omics <- readRDS('data/omics.rds')
consens.clust <- readRDS('data/consens_clust.rds')
tx.mods <- readRDS('data/tx_mods.rds')
tx.submods <- readRDS('data/tx_submods.rds')
tmt.mods <- readRDS('data/tmt_mods.rds')


# functions ---------------------------------------------------------------

# global to handle keyboard <enter> as a mouse click of "submit"
jscode <- '
$(function() {
  var $els = $("[data-proxy-click]");
  $.each(
    $els,
    function(idx, el) {
      var $el = $(el);
      var $proxy = $("#" + $el.data("proxyClick"));
      $el.keydown(function (e) {
        if (e.keyCode == 13) {
          $proxy.click();
        }
      });
    }
  );
});
'

# this function 
#	(1) performs GSEA for against all GO terms,
# (2) annotates results with biodomain-specific information

bd.gse <- function(gene_list, 
                   key_type = 'SYMBOL', 
                   score_type = 'pos'){
  
  # run GSEA analysis
  enr <- clusterProfiler::gseGO(
    geneList = gene_list, 
    ont = 'all',
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
    keyType = key_type,
    # minGSSize = 10,
    # maxGSSize = 500,
    eps = 0,
    pvalueCutoff = 1,
    # by = 'fgsea', nproc = 8, 
    scoreType = score_type
  )
  
  enr@result <- enr@result %>% 
    mutate(pathway = Description, 
           pval = pvalue, 
           padj = p.adjust, 
           ES = enrichmentScore, 
           size = setSize, 
           leadingEdge = core_enrichment) %>% 
    mutate(term = pathway) %>% 
    left_join(., biodom %>% select(pathway=GOterm_Name, Biodomain), by='pathway') %>% 
    full_join(., biodom %>% select(Biodomain, abbr, label, color) %>% distinct() , by = 'Biodomain' ) %>%
    relocate(Biodomain, .after = pathway) %>% 
    # mutate(
    #   Biodomain = case_when(is.na(Biodomain)~'none', T~ Biodomain),
    #   color = case_when(is.na(color)~'#7f7f7f', T~color)
    #   # , Biodomain = fct_relevel(Biodomain, as.character(bdt$domain) )
    # ) %>%
    # arrange(Biodomain) %>%
    rowwise() %>% mutate(leadingEdge = str_split(leadingEdge, '/')) 
  
  return(enr)
  
}

# this function 
#	(1) calculates hypergeometric GO term enrichment for all GO terms,
# (2) annotates results with biodomain-specific information
bd.enr <- function(gene_list, key_type = 'SYMBOL'){
  
  # run GO enricment analysis
  enr <- clusterProfiler::enrichGO(
    gene = gene_list, 
    OrgDb = org.Hs.eg.db::org.Hs.eg.db, 
    keyType = key_type, 
    pvalueCutoff = 1)
  
  enr@result <- enr@result %>% 
    left_join(., biodom %>% select(ID = GO_ID, Biodomain), by = 'ID') %>%
    full_join(., biodom %>% select(Biodomain, abbr, label, color) %>% distinct() , by = 'Biodomain' ) %>% 
    mutate(term = Description,
           padj = p.adjust)
  
  return(enr)
  
}


# This function quantifies the number of significantly enriched GO terms
# by AD biological domain

bd.tally <- function( enrRes, biodomDefTbl = biodom ){
  bdt <- bind_cols(
    biodomDefTbl %>% select(Biodomain, abbr, label, color) %>% distinct(),
    n_term = map_dbl( 
      unique(biodomDefTbl$Biodomain),
      ~ biodomDefTbl %>% filter(Biodomain == .x) %>% pull(GOterm_Name) %>% length()),
    n_sig_term = map_dbl( 
      unique(biodomDefTbl$Biodomain),
      ~ enrRes %>% filter( p.adjust < 0.05, Biodomain == .x) %>% pull(term) %>% unique() %>% length())) %>% 
    mutate(n_sig_term = ifelse(label == 'none', 
                               length(setdiff(enrRes$term[ enrRes$p.adjust < 0.05 ], biodom$GOterm_Name)), n_sig_term )) %>% 
    mutate(Biodomain = fct_reorder(Biodomain, n_sig_term, .desc = F)) %>% 
    arrange(Biodomain)
  
  bdt$prop <- bdt$n_sig_term / bdt$n_term
  
  return(bdt)
}

# this function 
# (2) plots violin plot of all biodom terms

bd.violin <- function(enr_result, plot_type , biodomDefTbl = biodom, 
                    padj_threshold = 0.05, size_cutoff = 5, include_none_bd = T){
  
  bdt <- bd.tally(enr_result, biodomDefTbl)
  # %>% arrange(prop)
  
  e <- enr_result %>% ungroup() %>% 
    filter( case_when( include_none_bd == F ~ label != 'none', T ~ !is.na(label)) ) %>%
    mutate(label = fct_relevel(
      label, 
      as.character(bdt$label) %>% str_subset(.,'none', negate = T) %>% c( 'none', .) ) #
    ) %>% 
    arrange(label)
  
  if(plot_type == 'enr'){
    
    p2 <- ggplot(e, aes( label, -log10(p.adjust) )) + theme_minimal()+
        scale_color_identity(guide = 'none')+ scale_fill_identity(guide = 'none')+
        scale_x_discrete(drop = F)+
        scale_size(
          limits = e$Count %>% range(na.rm = T),
          range = c(0.05,5))+
        geom_jitter( data = subset(e, p.adjust > 0.05 ),
                     color = 'grey85',  size = .3 )+
        geom_violin(data = subset(e, p.adjust <= 0.05 ),
                    scale='width', aes(col = color), lwd = .5 )+
        geom_jitter( data = subset(e, p.adjust <= 0.05 ),
                     aes(fill = color, size = Count, text = term),
                     color = 'grey20', alpha = .4, shape = 21 ) +
        labs(x='') + coord_flip() +
        theme(legend.position = 'none')
    
    p3 <- plotly::ggplotly(p2, tooltip = c('text','size'))
    
  } else if (plot_type == 'sig'){
    
    p2 <- e %>% arrange(desc(NES)) %>% 
      ggplot(aes( label, -log10(padj) ))+ theme_minimal()+
      scale_color_identity(guide = 'none')+ scale_fill_identity(guide = 'none')+
      scale_x_discrete(drop = F)+
      scale_size_continuous(
        limits = e$NES[e$padj <= padj_threshold] %>% range(na.rm = T),
        range = c(0.05,5), trans = 'exp')+
      geom_jitter( data = subset(e, padj > padj_threshold ),
                   color = 'grey85',  size = .3)+
      geom_violin(data = subset(e, padj <= padj_threshold ), 
                  scale='width', aes(col = color), lwd = .5)+ #color = 'grey50',
      geom_jitter( data = subset(e, padj <= padj_threshold ),
                   aes(fill = color, size = NES, text = term), 
                  color = 'grey20', alpha = .4, shape = 21) +
      geom_hline(yintercept = -log10(0.05), lty = 2)+
      labs(x='')+ coord_flip()+
      theme(legend.position = 'none')
    
    p3 <- plotly::ggplotly(p2, tooltip = c('text','size'))
    
  } else if (plot_type == 'nes'){
    
    p2 <- e %>% arrange(padj) %>% 
      ggplot(aes( label, NES ))+ theme_minimal()+
      scale_color_identity(guide = 'none')+ scale_fill_identity('none')+
      scale_x_discrete(drop = F)+
      scale_size(
        limits = -log10(e$padj[e$padj < padj_threshold]) %>% range(na.rm = T),
        range = c(0.05,5))+
      geom_jitter( data = subset(e, padj > padj_threshold ),
                   color = 'grey85',  size = .3)+
      geom_violin(data = subset(e, padj <= padj_threshold & NES < 0 ),
                  scale='width', aes(col = color), lwd = .5)+ 
      geom_violin(data = subset(e, padj <= padj_threshold & NES > 0 ),
                  scale='width', aes(col = color), lwd = .5)+ 
      geom_jitter( data = subset(e, padj <= padj_threshold ),
                   aes(fill = color, size = -log10(padj), text = term), 
                   color = 'grey20', alpha = .4, shape = 21)+ 
      geom_hline(yintercept = 0, lty = 2)+
      labs(x='')+ coord_flip()+
      theme(legend.position = 'none')
    
    p3 <- plotly::ggplotly(p2, tooltip = c('text','size'))
    
  }
  
  return(p3)
  
}


# # this function 
# # plots gsea results in table format
# 
# gse.table <- function(gene_list, enr_result, biodomDefTbl = biodom, 
#                       padj_threshold = 0.05, size_cutoff = 5,
#                       score_type = 'pos', vp = 'sig',
#                       include_none_bd = F, top_term_highlight = F){
#   
#   biodom.annotated <- biodomDefTbl %>% 
#     select(GOterm_Name, hgnc_symbol, n_hgncSymbol) %>% 
#     filter(!is.na(n_hgncSymbol)) %>% distinct() %>% 
#     pull(hgnc_symbol, name = GOterm_Name)
#   
#   bdt <- bd.tally(enr_result, biodomDefTbl)
#   # %>% arrange(prop)
#   
#   if(score_type == 'pos'){
#     
#     idx <- enr_result %>% 
#       filter(padj < padj_threshold, size > size_cutoff,
#              if( include_none_bd == F ) Biodomain != 'none' else !is.na(Biodomain) ) %>% 
#       arrange(desc(NES)) %>% select(pathway) %>% distinct() %>% 
#       pull(pathway) %>% .[1:20]
#     
#   } else if( score_type == 'std' ) {
#     
#     idx <- c(
#       enr_result %>% 
#         filter(padj < padj_threshold, size > size_cutoff, NES > 0,
#                if( include_none_bd == F ) Biodomain != 'none' else !is.na(Biodomain)) %>%  
#         arrange(desc(NES)) %>% select(pathway) %>% distinct() %>% 
#         pull(pathway) %>% .[1:10],
#       enr_result %>% 
#         filter(padj < padj_threshold, size > size_cutoff, NES < 0,
#                if( include_none_bd == F ) Biodomain != 'none' else !is.na(Biodomain)) %>%  
#         arrange(NES) %>% select(pathway) %>% distinct() %>% 
#         pull(pathway) %>% .[10:1]
#       
#     )
#   }
#   
#   pathways <- map(
#     idx,
#     ~ enr_result %>% select(pathway, leadingEdge) %>% distinct() %>% 
#       filter(pathway == .x) %>% pull(leadingEdge) %>% unlist() %>% 
#       c(., biodom.annotated[[.x]]) %>% unique()
#   )
#   names(pathways) = idx
#   
#   p <- fgsea::plotGseaTable( pathways , gene_list, enr_result, colwidths = c(5, 3, 0.8, 0, 1.2), render = F )
#   
#   return(p)
#   
# }
#   
# bd.dotplot <- function(enr_result){
#   
#   p <- dotplot(enr_result, showCategory = 20, orderBy = 'p.adjust')+ 
#     theme(axis.text.y = element_text(size = 9),
#           axis.text.x = element_text(size = 9))
#   
#   return(p)
#   
# }


  #     theme(axis.text.y = element_text(size = 9),
  #           axis.text.x = element_text(size = 9)),