x = scores %>% filter(name == 'Overall', !is.na(GeneName)) %>% 
  arrange(desc(value)) %>% slice(1:5000) %>% 
  pull(Overall, name = GeneName) %>%
  # bd.gse(gene_list = ., key_type = 'SYMBOL', score_type = 'pos')
  pull(GeneName) %>% 
  gprofiler2::gost(query = ., organism = 'hsapiens', ordered_query = T, significant = F, sources = 'GO' )

x$result <- x$result %>% 
  left_join(., biodom %>% select(term_id = GO_ID, Biodomain), by = 'term_id') %>%
  full_join(., biodom %>% select(Biodomain, abbr, label, color) %>% distinct() , by = 'Biodomain' ) %>% 
  mutate(term = term_name,
         padj = p_value) %>% 
  filter( label != 'none') %>% 
  mutate(label = fct_relevel(label, bd.tally(x$result %>% filter(significant == 'TRUE') %>% pull(term), biodom) %>% 
                               pull(label) %>% as.character() %>% str_subset(.,'none', negate = T))) %>% 
  arrange(label)

ggplot(x$result, aes( label, -log10(padj) )) + theme_minimal()+
  scale_color_identity(guide = 'none')+ scale_fill_identity(guide = 'none')+
  scale_x_discrete(drop = F)+
  scale_size(
    limits = x$result$intersection_size %>% range(na.rm = T),
    range = c(0.05,5))+
  geom_jitter( data = subset(x$result, padj > 0.05 ),
               color = 'grey85',  size = .3 )+
  geom_violin(data = subset(x$result, padj <= 0.05 ),
              scale='width', aes(col = color), lwd = .5 )+
  geom_jitter( data = subset(x$result, padj <= 0.05 ),
               aes(fill = color, size = intersection_size, text = term),
               color = 'grey20', alpha = .4, shape = 21 ) +
  labs(x='') + coord_flip() +
  theme(legend.position = 'none')


bd.violin(x@result, type = 'nes')
tad.scores %>% slice(1:5000) %>% filter(!is.na(GeneName)) %>% pull(Overall, name = GeneName) %>% gse.table(., x@result) %>% gridExtra::grid.arrange(.)

y = tmt.mods %>% filter(grepl('M42',mod)) %>% pull(gene) %>% 
  # bd.enr(gene_list = ., key_type = 'SYMBOL')
  gprofiler2::gost(query = ., organism = 'hsapiens', ordered_query = F, significant = F, sources = 'GO' )

y$result <- y$result %>% 
  left_join(., biodom %>% select(term_id = GO_ID, Biodomain), by = 'term_id') %>%
  full_join(., biodom %>% select(Biodomain, abbr, label, color) %>% distinct() , by = 'Biodomain' ) %>% 
  mutate(term = term_name,
         padj = p_value)

bd.violin(y@result, type = 'enr')
dotplot(y)


withProgress(message = 'fitting full linear model...', #detail = "(please wait)",
             style = 'notification', value = NULL, {
               
               newfit <- linear_fit(
                 cohort.data(), 
                 cohort.names(), 
                 mmus.betamatrix(), 
                 mmus.vars())
               
               lin({newfit})
               
             })

