#
# User interface ####
#
ui <- fluidPage(
  theme = shinytheme("spacelab"), # pick a theme later
  navbarPage(
    "TREAT-AD: risk scores & biological domains",
    
    #
    ## Enrichment Test ####
    #
      tabPanel(
        "Biological Domain Enrichment",
        
        ### Input ----
        fluidRow(
          
          column( 
            width = 12,
            
            wellPanel(
              
              h4(strong("(1) Select Enrichment Type:")),
              radioButtons("enr_type", NULL, selected = 'enr', inline = T,
                           choiceNames = c('Hypergeometric', 'GSEA (positive values)', 'GSEA (standard)'),
                           choiceValues = c('enr', 'pos', 'std')),
              h5(em("ensure that enrichment type selection matches the format of the genelist provided")),
              tags$hr(),
              # Input: use pre-defined gene lists
              h4(strong("(2a) Use pre-populated gene list:")),
              # checkboxInput('prepop', '?', F),
              
              fluidRow(
                
                column(6,
                       checkboxInput('tad_prepop', strong('TREAT-AD target scores (GSEA)'), F),
                       radioButtons('tad_prepop_', '',
                                    choiceNames = c('Target Risk Scores (positive)',
                                                    'Genetic Risk Scores (positive)',
                                                    'Multi-Omics Risk Scores (positive)',
                                                    'Transcriptomic Meta-analysis, FDR < 0.05 (standard)',
                                                    'Proteomic Meta-analysis, FDR < 0.05 (standard)'),
                                    choiceValues = c('Overall','Genetics','Omics','RNA','Protein')
                       ),
                       h3(""),
                       sliderInput('tad_slider', 'Use top N scored targets (for scores):', 
                                   max = scores %>% filter(name == 'Overall') %>% nrow(), 
                                   min = 1, round = T,
                                   value = 5000, step = 100, width = '75%')),
                
                column(6,
                       checkboxInput('tmt_prepop', strong('AMP-AD TMT Proteomics Modules (Hypergeometric)'), F),
                       selectInput('tmt_mod', NULL , unique(tmt.mods$mod), selected = "M42 Matrisome" ),
                       strong('From:'),
                       em(a(href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8825285/','Johnson et al 2022'))
                       
                       ,
                       checkboxInput('txmod_prepop', strong('AMP-AD Transcriptomics Modules (Hypergeometric)'), F),
                       selectInput('tx_mod', NULL , unique(tx.mods$mod), selected = "STGblue" ),
                       strong('From:'),
                       em(a(href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7428328/','Wan et al 2020'))
                       
                       ,
                       checkboxInput('txsubmod_prepop', strong('AMP-AD Transcriptomics Sub-Modules (Hypergeometric)'), F),
                       selectInput('txsub_mod', NULL , unique(tx.submods$submod), selected = "PHGturquoise_2" ),
                       strong('From:'),
                       em(a(href='https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7295244/','Milind et al 2020'))
                       
                       )
                
                
              ),
              
              h3( "" ),
              
              tags$hr(),
              h4(strong("(2b) Supply your own gene list:")),
              fluidRow(
                column(6,
                  fluidRow(
                    column(8,
                           # Input: Select a file 
                           fileInput("file1", "Upload File (.csv)",
                                     multiple = FALSE,
                                     accept = c("text/csv","text/plain",
                                                "text/comma-separated-values,text/plain",".csv"))
                           ),
                    column(2,
                           # Input: Checkbox if file has header 
                           checkboxInput("header", "Includes Header Row?", F)
                           )
                    )
                ),
                column(6,
                       textAreaInput('txt_input','Paste Input Gene List:',
                                     resize = 'vertical', height = '150px', #width = '500px',
                                     placeholder = 'enter a list of genes like (hypergeometric):\nAPP\nAPOE...\n\nor genes with values in decreasing order like (GSEA):\nAPP\t4.71\nAPOE\t4.68...' 
                       ))
              ),
              # h3( "" ),
              tags$hr(),
              actionButton(inputId = 'submit_genelist' , label = 'Submit', class="btn-primary")
              
            )
          )
        ), 
        
        ### Output ----
        fluidRow(
          tabsetPanel( 
            # tabPanel('Biological Domain Enrichments', plotOutput("bd_violins", width = '50%', height = '600px' )),
            tabPanel('Enrichment Plot', plotly::plotlyOutput("bd_violins", width = '75%', height = '500px' )),
            # tabPanel('GO Term Enrichments', plotOutput("term_plot", width = '50%', height = '600px' )),
            tabPanel('Enrichment Results Table',DT::dataTableOutput("enr_table"))
          )
        )
      ),
      
    #
    ## Target-specific information ####
    #
    tabPanel(
      "Target-specific information",
      
      ### Input ----
      fluidRow( 
        column( width = 12,
                wellPanel( 
                  div(tags$head(tags$script(HTML(jscode))),
                      
                      div(style="display: inline-block; width: 300px ;", 
                          tagAppendAttributes(
                            textInput( inputId = "symbol", #width = 5,
                                       label = "Gene symbol:",
                                       placeholder = 'e.g. APP or PLEC or SMOC1'),
                            `data-proxy-click` = "update")),
                      div(style="display: inline-block; width: 95px ;", 
                          actionButton("update", "Plot", class="btn-primary"))) , #, class = "btn-success"
                  br(),
                  strong("notes:"), 
                  # br(),
                  # em("\u2022 Overall is Genetics + Omics only"),
                  br(),
                  em("\u2022 No point indicates gene not scored or measured"),
                  br(),
                  em("\u2022 If a gene isn't displayed, try clicking on the 'GeneCards' link below to ensure the correct symbol")
                )
        )
      ),
      
      ### Output ----
      fluidRow( 
        tabsetPanel(
          type = "tabs", 
          tabPanel( "Scores & Expression", plotOutput("score_expr_plot")), 
          tabPanel( "Biodomains",
            fluidRow(
              column(4, plotOutput("biodom_plot"))
              , column(6,  plotly::plotlyOutput('biodom_nw'))
              ))
          )
      )
    ),
  
  # For more information links
  fluidRow(
    column(width = 12, 
           hr(),
           strong('For more information:'),
           # br(),
           em(a(href='https://www.synapse.org/#!Synapse:syn25575156/tables/','Overall Scores'), ' | ', 
              a(href='https://www.synapse.org/#!Synapse:syn26844312/tables/','Genetics Scores'), ' | ', 
              a(href='https://www.synapse.org/#!Synapse:syn22758536/tables/','Multi-omics Scores')),
           # hr(),
           uiOutput('linkOut')
    )
  )
)
)

