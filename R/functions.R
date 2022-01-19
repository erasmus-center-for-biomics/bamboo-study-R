

#' Read metaphlan data
#' 
#' @param path
#'
#' 
#' @export
#' 
read_metaphlan <- function(path) {
  read_tsv(
    path,
    col_names = c('clade_name', 'ncbi_taxid_trace', 'relative_abundance', 'additional_species'),
    col_types = 'ccdc', 
    comment = '#'
  ) %>%
    mutate(
      ID = str_match(
        basename(path), 
        'I[0-9]{2}-[0-9]{4}-[A-Z]?[0-9]{2}')[,1],
      clade_id = str_remove(clade_name, '^.*[|]'),
      clade_type = str_remove(clade_id, '__.*$'),
      clade_id = str_remove(clade_id, '^.__'),
      ncbi_taxid = str_remove(ncbi_taxid_trace , '^.*[|]'),
      parent_taxid = str_match(ncbi_taxid_trace , '([^|]*)[|][^|]*$')[,2]
    )
}


#' Create a phyloseq taxonomy object
#'
#' @param data the metaphlan data
#'
#' @export
to_taxtable <- function(data) {
  data %>%
    filter(clade_type == 's') %>%
    group_by(clade_id, clade_name) %>%
    summarise() %>%
    ungroup() %>%
    separate(
      clade_name,
      c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
      sep = '[|]',
      remove = TRUE
    ) %>%
    mutate(
      Domain = str_remove(Domain, '^.__'),
      Phylum = str_remove(Phylum, '^.__'),
      Class = str_remove(Class, '^.__'),
      Order = str_remove(Order, '^.__'),
      Family = str_remove(Family, '^.__'),
      Genus = str_remove(Genus, '^.__'),
      Species = str_remove(Species, '^.__')
    ) %>%
    arrange(clade_id) %>%
    (function(tbl){
      m <- as.matrix(tbl[,2:ncol(tbl)])
      rownames(m) <- tbl[[1]]
      m
    }) %>%
    tax_table()
}

#' Create a phyloseq OTU table
#' 
#' @param data the metaphlan data
#' @param clades clades to include
#'
#' @export
to_otutable <- function(data, clades) {
  data %>%
    filter(clade_id %in% clades) %>%
    pivot_wider(
      id_cols = clade_id,
      names_from = DNA_ID,
      values_from = relative_abundance,
      values_fill = list(
        'relative_abundance' = 0.0
      ) 
    ) %>%
    arrange(clade_id) %>%
    (function(tbl){
      m <- as.matrix(tbl[,2:ncol(tbl)])
      rownames(m) <- tbl[[1]]
      m
    }) %>%
    otu_table(taxa_are_rows = TRUE)
}

#' Prepare the phyloseq sample data
#'
#' @param metadata
#' @param samples samples to include
#' 
#' @export
to_sampledata <- function(metadata, samples) {
  metadata %>%
    filter(DNA_ID %in% samples) %>%
    arrange(DNA_ID) %>%
    as.data.frame() %>%
    (function(d){
      rownames(d) <- d$DNA_ID 
      d
    }) %>%
    sample_data()
}


#' Get the alpha diversity metrics
#'
#'
alpha_diversity <- function(pobj) {
  bind_cols(
    diversity(
      pobj,  
      c('inverse_simpson', 'shannon', 'gini_simpson', 'coverage')
    ),
    richness(pobj),
    evenness(pobj),
    dominance(pobj),
    rarity(pobj)
  ) %>% 
    (function(x){
      x %>%
        as_tibble() %>%
        mutate(DNA_ID = rownames(x))
    })
}



density_or_scatter <- function(d, crit) {
  if(n_distinct(d$mvalue) <= 3) {
    d %>%
      ggplot(aes(x = value, fill = as.factor(mvalue))) +
      geom_density(alpha = 0.5, adjust = 0.5) +
      geom_point(y = 0, pch = 21, size = 2) +
      labs(
        x = crit$axis,
        y = crit$criteria,
        caption = sprintf(
          "P-value: %s", 
          format(crit$p.value, digits = 2))
      ) +
      theme_minimal() +
      theme(
        panel.border = element_rect(colour = 'grey20', fill = NA)
      )
  } else {
    d %>%
      ggplot(aes(x = value, y = mvalue)) +
      geom_point() +
      labs(
        x = crit$axis,
        y = crit$criteria,
        caption = sprintf(
          "P-value: %s", 
          format(crit$p.value, digits = 2))
      ) +
      theme_minimal() +
      theme(
        panel.border = element_rect(colour = 'grey20', fill = NA)
      )
  }
}


#
# Targets helpers
#


tar_write_tsv <- function(x, f, ...) {
  write_tsv(x, f, ...)
  f
}

tar_ggsave <- function(x, f, ...) {
  ggsave(x, filename = f, ...)
  f
}

