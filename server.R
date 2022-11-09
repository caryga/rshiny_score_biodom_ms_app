
server <- function(input, output, session) {
  
  rv <- reactiveValues()
  
  
  ## enrichment tests ----
  rv$type = NULL
  rv$gene_list = NULL
  rv$enr_res = NULL
  rv$violin_plot = NULL
  rv$term_plot = NULL
  
  
  observe({
    
  if(input$enr_type == 'enr' | input$tmt_prepop == T ){  
      rv$type <- 'enr' 
  } else if( input$enr_type == 'pos' | (input$tad_prepop == T & input$tad_prepop_ %in% c('Overall','Genetics','Omics')) ){
      rv$type <- 'pos'
  } else if( input$enr_type == 'std' | (input$tad_prepop == T & input$tad_prepop_ %in% c('RNA','Protein')) ){
      rv$type <- 'std'
  }
    
  }) %>% bindEvent(input$submit_genelist, ignoreInit = T, ignoreNULL = T )
  
  ### genelist ----
  observe({
    if(input$tad_prepop == T & rv$type == 'pos' ){
        rv$gene_list <- scores %>% 
          filter(name == input$tad_prepop_) %>% 
          arrange(desc(value)) %>% pull(value, name = GeneName) %>% .[1:input$tad_slider] 
    } else if(input$tad_prepop == T & rv$type == 'std'){
        rv$gene_list <- omics %>% 
          filter(x == input$tad_prepop_, fdr < 0.05) %>% 
          arrange(desc(TE)) %>% 
          pull(TE, name = GName)
    } else if(input$tmt_prepop == T){
        rv$gene_list <- tmt.mods %>% 
          filter(mod == input$tmt_mod) %>% 
          pull(gene)
    } else if(input$txmod_prepop == T){
      rv$gene_list <- tx.mods %>% 
        filter(mod == input$tx_mod) %>% 
        pull(gene)
    } else if(input$txsubmod_prepop == T){
      rv$gene_list <- tx.submods %>% 
        filter(submod == input$txsub_mod) %>% 
        pull(gene)
    } else if(!is.null(input$file1)){
      if(rv$type %in% c('enr')){
        rv$gene_list <- read_csv(input$file1$datapath, col_names = input$header) %>%
          rename_with(., ~ c('gene'), everything()) %>% pull(gene)
      } else if(rv$type %in% c('pos','std')){
        rv$gene_list <- read_csv(input$file1$datapath, col_names = input$header) %>%
          rename_with(., ~ c('gene','score'), everything()) %>% pull(score, name = gene)
      }
    } else if(input$txt_input != ''){
      if(rv$type %in% c('enr')){
        rv$gene_list <- tibble(input$txt_input) %>% 
          rename_with(., ~ c('gene'), everything()) %>% pull(gene)
      } else if(rv$type %in% c('pos','std')){
        rv$gene_list <- tibble(input$txt_input) %>% 
          rename_with(., ~ c('gene','score'), everything()) %>% pull(score, name = gene)
      }
      }
  }) %>% bindEvent(input$submit_genelist, ignoreInit = T, ignoreNULL = T )
  
  ### enrichment ----
  observe({
    
    req(rv$gene_list)
    
    if( rv$type == 'enr' ) {
      withProgress(message = 'enriching...', #detail = "(please wait)",
                   style = 'notification', value = NULL, {
                     rv$enr_res <- bd.enr(rv$gene_list, key_type = 'SYMBOL')
                   })
    } else if( rv$type %in% c('pos','std') ) {
      withProgress(message = 'enriching...', #detail = "(please wait)",
                   style = 'notification', value = NULL, {
                     rv$enr_res <- bd.gse(rv$gene_list, key_type = 'SYMBOL', score_type = rv$type )
                   })
    } 
    
  }) %>% bindEvent(input$submit_genelist, ignoreInit = T, ignoreNULL = T )
  
  ### violin plot ----
  observe({
    
    req(rv$enr_res)
    
    if (rv$type == 'enr') {
      # output$bd_violins <- renderPlot({ bd.violin(rv$enr_res@result, plot_type = 'enr', include_none_bd = F) })
      output$bd_violins <- plotly::renderPlotly( bd.violin(rv$enr_res@result, plot_type = 'enr', include_none_bd = T) )
      
    } else if (rv$type  == 'std') {
      # output$bd_violins <- renderPlot({ bd.violin(rv$enr_res@result, plot_type = 'nes', include_none_bd = F) })
      output$bd_violins <- plotly::renderPlotly( bd.violin(rv$enr_res@result, plot_type = 'nes', include_none_bd = T) )
    
    } else if (rv$type  == 'pos') {
      # output$bd_violins <- renderPlot({ bd.violin(rv$enr_res@result, plot_type = 'nes', include_none_bd = F) })
      output$bd_violins <- plotly::renderPlotly( bd.violin(rv$enr_res@result, plot_type = 'sig', include_none_bd = T) )
    }
    
  }) %>% bindEvent(input$submit_genelist, ignoreInit = T, ignoreNULL = T )

  
  # ### term plot ----
  # observe({
  #   
  #   req(rv$gene_list, rv$enr_res)
  #   
  #   if (rv$type == 'enr') {
  #     output$term_plot <- renderPlot({ bd.dotplot(rv$enr_res) })
  #   } else if (rv$type %in% c('pos','std')) {
  #     output$term_plot <- renderPlot({ gse.table(rv$gene_list, rv$enr_res@result, score_type = rv$type, include_none_bd = F) })
  #   }
  #   
  # }) %>% bindEvent(input$submit_genelist, ignoreInit = T, ignoreNULL = T )
      
        
  ### enr result data table ----
  output$enr_table <- DT::renderDataTable( server = FALSE, {
    
    req(rv$enr_res)
    
    rv$enr_res@result %>%
      select(Biodomain, term, p.adjust,
             contains('NES'), contains('GeneRatio'),
             contains('leadingEdge'), contains('geneID')) %>%
      mutate(Biodomain = ifelse(is.na(Biodomain), 'none', Biodomain)) %>%
      DT::datatable(.,
                    filter = 'top',
                    
                    extensions = 'Buttons',
                    options = list(
                      paging = TRUE,
                      searching = TRUE,
                      fixedColumns = F,
                      autoWidth = TRUE,
                      ordering = TRUE,
                      dom = 'Bfrtip',
                      buttons = c('copy', 'csv', 'excel')
                    ),
                    
                    class = "display"
      
                    ) %>%
      DT::formatSignif(., columns = c('p.adjust'), digits = 3)
    
    })
  
  output$mytable1  <- DT::renderDataTable(
    DT::datatable(
      { plots.dfs()[[1]] }))
  
  
  ## target info ----
  
  rv$generate_score_expr_plot = NULL
  rv$generate_biodom_plot = NULL
  
  ### scores & expression plot ----
  
  observe({
    
    tep.tg = input$symbol %>% str_to_upper()
    
    if( tep.tg %in% scores$GeneName ){
      
      score.plot <- ggplot(scores, aes( value, name )) + 
        geom_violin(scale = 'width', aes(fill = name), alpha = .3, 
                    draw_quantiles = c(.5,.9), trim = F)+ 
        geom_point(data = subset(scores, GeneName == tep.tg), size = 3, alpha = .6) + 
        viridis::scale_fill_viridis(discrete = T, guide ='none', option = 'D', direction = -1) + 
        labs(y='', x= 'score', 
             title = paste0(tep.tg,'\nTREAT-AD scores')
             , subtitle = paste0('Rank #',
                                 scores$rank[which(scores$GeneName==tep.tg)] %>% unique(),
                                 '\nlines: 50th & 90th quantiles'))
      
    } else {
      
      score.plot <- ggplot(scores, aes( value, name)) + 
        geom_violin(scale = 'width', aes(fill = name), alpha = .3, 
                    draw_quantiles = c(.5,.9), trim = T)+ 
        geom_point(data = subset(scores, GeneName == tep.tg), size = 3, alpha = .6) + 
        viridis::scale_fill_viridis(discrete = T, guide ='none', option = 'D', direction = -1) + 
        labs(y='', x= 'score', 
             title = paste0(tep.tg,'\nTREAT-AD scores')
             , subtitle = paste0('**GENE NOT SCORED**'))
      
    }
    
    if( tep.tg %in% omics$GName ){
      
      expr.plot <- ggplot(omics, aes( TE, x )) + 
        geom_violin(fill = 'grey60', alpha = .3) +
        geom_point(
          data = subset(omics, GName == tep.tg & fdr < 0.05),
          size = 3, color = 'red', alpha = .6) +
        geom_point(
          data = subset(omics, GName == tep.tg & fdr > 0.05),
          size = 3, color = 'grey65', alpha = .6) +
        geom_vline(xintercept = 0, lty = 2) +
        labs(y = '', 
             x = 'meta-analysis log fold change\nAD vs Control',
             title = paste0(tep.tg,'\nDEG meta-analysis'),
             subtitle = 'red: FDR < 0.05') 
      
    } else {
      
      expr.plot <- ggplot(omics, aes( TE, x )) + 
        geom_violin(fill = 'grey60', alpha = .3) +
        geom_point(
          data = subset(omics, GName == tep.tg & fdr < 0.05),
          size = 3, color = 'red', alpha = .6) +
        geom_point(
          data = subset(omics, GName == tep.tg & fdr > 0.05),
          size = 3, color = 'grey65', alpha = .6) +
        geom_vline(xintercept = 0, lty = 2) +
        labs(y = '', 
             x = 'meta-analysis log fold change\nAD vs Control',
             title = paste0(tep.tg,'\nDEG meta-analysis'),
             subtitle = '**EXPRESSION NOT DETECTED**') 
    }
    
    if( tep.tg %in% gen$GeneName ){
      
      gen.plot <- ggplot(gen, aes(value, name)) + 
        geom_violin(scale = 'width', aes(fill = name), alpha = .3,
                    draw_quantiles = c(.5,.95), trim = F)+
        geom_point(data = subset(gen, GeneName == tep.tg),
                   size = 3, alpha = .6 )+
        viridis::scale_fill_viridis(discrete = T, guide ='none', option = 'D', direction = -1) +
        labs(y='', x= 'score', 
             title = paste0(tep.tg,'\nGenetic Evidence'),
             subtitle = paste0('Rank #',
                               gen$score_rank[which(gen$GeneName==tep.tg)] %>% unique(),
                               '\nlines: 50th & 95th quantiles\n95th quantile min_gwasP ~ p = 5e-8')
        )
      
    } else {
      
      gen.plot <- ggplot(gen, aes(value, name)) + 
        geom_violin(scale = 'width', aes(fill = name), alpha = .3,
                    draw_quantiles = c(.5,.95), trim = F)+
        geom_point(data = subset(gen, GeneName == tep.tg),
                   size = 3, alpha = .6 )+
        viridis::scale_fill_viridis(discrete = T, guide ='none', option = 'D', direction = -1) +
        labs(y='', x= 'score', 
             title = paste0(tep.tg,'\nGenetic Evidence' ),
             subtitle = '**GENE NOT SCORED**') 
    }
    
    rv$generate_score_expr_plot <- cowplot::plot_grid(
      score.plot, 
      cowplot::plot_grid(
        gen.plot,expr.plot, ncol = 1,
        rel_heights = c(1,.7), align = 'v' ), 
      ncol=2 )
    
  }) %>% bindEvent(input$update, ignoreInit = F)
  
  ### biodom plot ----    
  observe({
    
    tep.tg = input$symbol %>% str_to_upper()
    
    if(tep.tg %in% biodom_genes$symbol){
      
      bd.plot <- biodom_genes %>% 
        filter(symbol == tep.tg, Biodomain != 'none') %>% 
        select(-GO_ID,-GOterm_Name, -n_symbol) %>%  distinct() %>% 
        full_join(., dom.cols, by = c('Biodomain'='domain')) %>% 
        mutate( pct = case_when(is.na(pct)~0, T~pct),
                Biodomain = fct_reorder(Biodomain, pct) ) %>% 
        ggplot(., aes( pct, Biodomain )) + 
        geom_segment( aes(yend=Biodomain, xend=0), color = 'grey50') + 
        geom_point(size = 3, alpha = .9, aes(fill = color), shape = 21) + 
        scale_fill_identity() + expand_limits(x=0) + 
        labs(y='', x = '% Biodom\nGO terms', 
             title = paste0(tep.tg,'\nBiodomains')) 
      
      net <- biodom_genes %>% 
        filter(symbol == tep.tg) %>% 
        select(Biodomain, GOterm_Name) %>% 
        network::network(., directed = F)
      
      v.attr = tibble( v = net %v% 'vertex.names' )
      v.attr = bind_cols(
        v.attr,
        sizes = biodom_genes %>% filter(symbol == tep.tg) %>% 
          .[match( net %v% 'vertex.names',.$GOterm_Name),] %>% pull(n_symbol),
        class = if_else(v.attr$v %in% unique(biodom_genes$Biodomain), 'bd','term')
      ) %>% 
        left_join(., dom.cols, by = c('v'='domain')) %>% 
        mutate(sizes = case_when( is.na(sizes) ~ 2000, T ~ as.double(sizes)),
               sizes = sizes / max(sizes),
               # sizes = rank(sizes, ties.method = 'max'),
               color = case_when( is.na(color) ~ 'grey80', T ~ color)) 
      
      net %v% 'size' = v.attr$sizes
      net %v% 'color' = v.attr$color
      net %v% 'class' = v.attr$class
      
      bd.nw <- suppressWarnings(
        GGally::ggnet2( 
          net, edge.size = .5, shape.palette = c(15,19), 
          node.size = 'size', node.color = 'color', node.shape = 'class',
          label = unique(biodom_genes$Biodomain) , label.size = 4) + 
          theme(legend.position = 'none') +
          geom_point(aes(text = v.attr$v ), color = 'grey80', alpha = 0, size = 10) +
          ggtitle(paste0(tep.tg,'- Biodomain term network'))
      )
      
    } else { 
      
      bd.plot  <-  biodom_genes %>% 
        filter(symbol == 'SMOC1') %>% 
        select(-GO_ID,-GOterm_Name, -n_symbol) %>%  distinct() %>% 
        full_join(., dom.cols, by = c('Biodomain'='domain')) %>% 
        mutate( pct = NA_real_,
                Biodomain = fct_reorder(Biodomain, pct) ) %>% 
        ggplot(., aes( pct, Biodomain )) + 
        geom_segment( aes(yend=Biodomain, xend=0), color = 'grey50') + 
        geom_point(size = 3, alpha = .9, aes(fill = color), shape = 21) + 
        scale_fill_identity() + expand_limits(x=0) + 
        labs(y='', x = '% Biodom\nGO terms', 
             title = paste0(tep.tg, '\nBiodomains',tep.tg),
             subtitle = '**NO BIODOM ANNOTATED**')  
      
      bd.nw = NULL
      
    }
    
    rv$generate_biodom_plot <- bd.plot
    
    rv$generate_bd_nw = bd.nw 
    
  }) %>%  bindEvent(input$update, ignoreInit = F)
  
  # plot rendering ---- 
  output$score_expr_plot <- renderPlot( rv$generate_score_expr_plot, width = 600  )
  output$biodom_plot <- renderPlot( rv$generate_biodom_plot )
  output$biodom_nw <- plotly::renderPlotly( plotly::ggplotly(rv$generate_bd_nw, tooltip = 'text')  ) # %>% layout(width = 350)
  # output$cell_type_plot <- renderPlot( rv$generate_cell_type_plot, width = 600 )
  # output$phenotypes_plot <- renderPlot( rv$generate_phenotypes_plot, width = 600 )
  
  output$linkOut <- renderUI({ 
    list(
      HTML('<strong>Link Outs: </strong>'),
      HTML(paste0('<a href = "https://www.genecards.org/cgi-bin/carddisp.pl?gene=', input$symbol, '" target="_blank">GeneCards</a>')),
      ' | ',
      HTML(paste0('<a href = "https://pubmed.ncbi.nlm.nih.gov/?term=alzheimer%20', input$symbol, '" target="_blank">PubMed</a>'))
    )
  }) 
  
  
}

