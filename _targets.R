library(targets)

source('R/functions.R')

options(tidyverse.quiet = TRUE)
tar_option_set(packages = c(
  "tidyverse",
  "phyloseq",
  "microbiome",
  "wesanderson",
  "broom",
  "cowplot",
  "ggrepel",
  "ggforce"
 )) 


if(!dir.exists("output"))
  dir.create("output")

list(

  #
  # Specify the input files
  #

  tar_target(
    metadata_file,
    'data/bamboo_metadata.csv',
    format = 'file'
  ),
  
  tar_target(
    metaphlan_files,
    file.path('data', list.files('data', pattern = '.tsv$')),
    format = 'file'
  ),
  
  #
  # Read the data
  #
  tar_target(
    metadata,
    read_delim(metadata_file, ";")
  ),
  
  tar_target(
    metaphlan,
    left_join(
      metaphlan_files %>%
        map_df(read_metaphlan),
      metadata %>% select(ID, DNA_ID),
      by = 'ID')
  ),
  
  # filter viruses from the metaphlan and rescale the relative abundance
  tar_target(
    metaphlan_flt,
    metaphlan %>%
      mutate(
        `is virus` = str_detect(clade_name, '^k__Viruses')
      ) %>%
      filter(!`is virus`) %>%
      group_by(ID, clade_type) %>%
      mutate(
        relative_abundance = relative_abundance / sum(relative_abundance) * 100
      ) %>%
      ungroup()
  ),
  
  #
  # Create a Phyloseq object
  #
  tar_target(
    taxtable, 
    to_taxtable(metaphlan_flt)
  ),
  
  tar_target(
    otutable, 
    to_otutable(metaphlan_flt, rownames(taxtable))
  ),
  
  tar_target(
    sampledata,
    to_sampledata(metadata, colnames(otutable))
  ),
  
  tar_target(
    phyloobj,
    phyloseq(otutable, taxtable, sampledata)
  ),
  
  #
  # Get the alpha diversity metrics
  #
  
  tar_target(
    alphadiv,
    alpha_diversity(phyloobj)
  ),
  
  
  #
  # Plot the microbial composition
  #
  
  tar_target(
    topgenus,
    metaphlan %>%
      filter(clade_type == 'g') %>%
      group_by(clade_id, ncbi_taxid, clade_name) %>%
      summarise(ra = sum(relative_abundance)) %>%
      ungroup() %>%
      mutate(
        `is virus` = str_detect(clade_name, '^k__Viruses'),
        r = rank(-ra)) %>%
      arrange(r)
  ),
  
  tar_target(
    metagenus,
    metaphlan %>%
      filter(clade_type == 'g') %>%
      mutate(
        `is virus` = str_detect(clade_name, '^k__Viruses'),
        `in top` = clade_id %in% topgenus$clade_id[topgenus$r < 20],
        label = ifelse(`in top`, clade_id, "other")
      ) %>%
      group_by(DNA_ID, `is virus`, label) %>%
      summarise(relative_abundance = sum(relative_abundance)) %>%
      ungroup() %>%
      group_by(DNA_ID,  `is virus`) %>%
      mutate(RnA = relative_abundance / sum(relative_abundance)) %>%
      ungroup() %>%
      mutate(
        label = parse_factor(label, levels = c(topgenus$clade_id[topgenus$r < 20], "other"))
      ) %>%
      left_join(
        metadata %>%
          select(DNA_ID, PatientTitle, Intervention, Stage),
        by = 'DNA_ID'
      ) %>%
      mutate(
        Ilabel = case_when(
          Intervention == "0" ~ "Run-in VF",
          Intervention == "A" ~ "Intervention MF",
          Intervention == "C" ~ "Intervention VF",
          TRUE ~ "unknown"
        ),
        I = sprintf("(%d) %s", Stage, Ilabel),
        Subject = as.integer(str_remove(PatientTitle, '^ART-0+'))
      )
  ),
  
  # a colour map of 3 colours per base colour
  tar_target(
    tri_colour_map,
    c(
      # blue
      '#1E90FF', '#4668a5', '#0000FF', 
      # green
      '#90EE90', '#008000', '#6B8E23',
      # red
      '#FF0000', '#FA8072', '#CD5C5C',
      # teal
      '#66CDAA', '#20B2AA', '#008080',
      # orange
      '#FFD700', '#FFA500', '#FF8C00',
      # purple
      '#E6E6FA', '#800080', '#6A5ACD',
      # brown
      '#CD853F', '#A52A2A', '#D2691E',
      # pink 
      '#FFC0CB', '#DB7093', '#C71585',
      # end with black
      'black'
    )
  ),
  
  tar_target(
    virus_labels,
    metagenus %>%
      filter(`is virus`) %>%
      filter(label != "other") %>%
      mutate(x = as.character(label)) %>%
      pull(x)
  ),
  
  tar_target(
    genus_colour_map,
    metagenus %>%
      filter(!`is virus`) %>% 
      pull(label) %>%
      levels() %>%
      (function(l, rm){
        l <- l[!l %in% rm]
        x <- c(tri_colour_map[1:length(l)-1], 'black')
        names(x) <- l
        x
      })(rm=virus_labels)
  ),
  
  tar_target(
    sample_grid,
    metagenus %>% 
      group_by(Subject) %>% 
      summarise(
        T = case_when(
          any(str_detect(I, '\\(1\\).*MF')) ~ "MF",
          TRUE ~ "VF"
        )
      ) %>%
      ungroup() %>%
      group_by(T) %>%
      mutate(
        i = 1:n(),
        column = ceiling(i / 4),
        row = i  - (column * 4) + 4
      ) %>%
      ungroup() %>%
      mutate(column = sprintf("%s %d", T, column))
  ),
  
  tar_target(
    plt_metagenus,
    metagenus %>% 
      filter(!`is virus`) %>% 
      left_join(sample_grid, by = "Subject") %>%
      split(.$T) %>%
      map(function(d) {
        
        o <- d %>%
          group_by(Stage, Ilabel) %>%
          summarise() %>%
          ungroup() %>%
          arrange(Stage)
        
        txt <- d %>%
          group_by(Subject, column, row) %>%
          summarise() %>%
          mutate(
            RnA = 0.9, 
            Ilabel = o$Ilabel[2]
          )
        
        d %>%
          rename(genus = label) %>%
          mutate(Ilabel = parse_factor(Ilabel, levels = o$Ilabel)) %>%
          ggplot() + 
          geom_bar(aes(x = Ilabel, y = RnA, fill = genus), stat = 'identity') +
          geom_label(aes(x = Ilabel, y = RnA, label = Subject), data = txt) +
          labs(
            x = "",
            y = "",
          ) +
          scale_fill_manual(values = genus_colour_map) +
          facet_grid(row ~ column, scales = 'free_x') +
          theme_minimal() +
          theme(
            aspect.ratio = 1,
            panel.grid = element_blank(),
            strip.text = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust=0),
            axis.title.y = element_blank()
          ) 
      })
  ),
  
  tar_target(
    path_metagenus_1,
    tar_ggsave(plt_metagenus[[1]], "output/metagenus_MF.pdf", width = 12, height = 12),
    format = 'file'
  ),
  
  tar_target(
    path_metagenus_2,
    tar_ggsave(plt_metagenus[[2]], "output/metagenus_VF.pdf", width = 12, height = 12),
    format = 'file'
  ),
  
  #
  # Get the beta diversity metrics
  #
  
  tar_target(
    todo, 
    c("bray")
  ),
  
  tar_target(
    betadist,
    todo %>%
      set_names(nm = todo) %>%
      map(function(tp) {
        distance(phyloobj, tp)
      })
  ),
  
  #
  # Determine which genus differ between clusters 
  # based on Bray Curtis distance   
  #

  # cluster the samples
  tar_target(
    hcl,
    hclust(betadist[["bray"]], method = "ward.D2")
  ),
  
  tar_target(
    hk,
    cutree(hcl, k = 3)
  ),
  
  # convert the m
  tar_target(
    genus_m, 
    metaphlan %>%
      filter(clade_type == 'g') %>%
      mutate(
        `is virus` = str_detect(clade_name, '^k__Viruses')
      ) %>%
      filter(!`is virus`) %>%
      group_by(DNA_ID, clade_id) %>%
      summarise(relative_abundance = sum(relative_abundance)) %>%
      ungroup() %>%
      pivot_wider(
        id_cols = "clade_id",
        names_from = "DNA_ID",
        values_from = "relative_abundance",
        values_fill = list("relative_abundance" = 0.0)
      ) %>%
      (function(x){
        m <- as.matrix(x[,2:ncol(x)])
        rownames(m) <- x[[1]]
        m
      })
  ),
  
  # determine differential genuses between the clusters
  tar_target(
    genus_kruskal,
    apply(genus_m, 1, function(x, g = hk){
      kruskal.test(x, g) %>% tidy()
    }) %>%
      bind_rows() %>%
      mutate(
        padj = p.adjust(p.value, method = "holm"),
        clade_id = rownames(genus_m),
        `top clade` = clade_id %in% metagenus$label
      ) %>%
      arrange(p.value)
  ),
  
  # write the kruskal results to a file
  tar_target(
    path_genus_kruskal,
    tar_write_tsv(genus_kruskal, "output/genus_kruskal.xls"),
    format = "file"
  ),

  # plot the clusters
  tar_target(
    plt_cluster,
    genus_m %>%
      as_tibble() %>%
      mutate(
        clade_id = rownames(genus_m)
      ) %>%
      filter(clade_id %in% c("Aeriscardovia", "Klebsiella", "Olsenella", "Thermoleophilum", "Enorma")) %>%
      pivot_longer(
        !clade_id,
        names_to = "DNA_ID",
        values_to = "score"
      ) %>%
      left_join(
        tibble(
          DNA_ID = names(hk),
          Cluster = as.factor(hk)
        ),
        by = "DNA_ID"
      ) %>%
      ggplot(aes(x = DNA_ID, y = score)) +
      geom_bar(stat = "identity") +
      facet_grid(clade_id ~ Cluster, scales = "free", space = "free_x") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 90)
      )
  ),
  
  #
  # Determine whether and which categorical  
  #   parameters are over-represented in the clusters
  #
  
  # sanitize the input to a categorical list
  tar_target(
    catset,
    tibble(
      Cluster = as.factor(hk),
      DNA_ID = names(hk)
    ) |>
      left_join(
        metadata, 
        by = "DNA_ID"
      ) |> 
      mutate(
        Gender = ifelse(Gender_male == 1, "M", "F"),
        Breastfed = ifelse(Breastfed_yes == 1, "Y", "N"),
        Delivery = ifelse(Delivery_Normal == 1, "N", "O"),
        `Bowel frequency` = str_c(
          Bowel_more_3_day,
          Bowel_2_3_day,
          Bowel_once_day,
          Bowel_3_6_week,
          Bowel_less_2_week,
          sep = "/"
        ),
        Consistency = str_c(
          Consistency_Watery,
          Consistency_Soft,
          Consistency_Formed,          
          sep = "/"
        ),
        Color = str_c(
          Color_1,
          Color_2,
          Color_3,
          Color_4,
          Color_5,
          Color_6,
          sep = "/"
        ),
        Amount = str_c(
          Amount_smear,
          Amount_Upto25,
          Amount_25_50,
          Amount_50,
          sep = "/"
        )
      ) |>
      select(
        DNA_ID,
        Cluster,
        PatientID,
        Gender,
        Breastfed, 
        Delivery,
        # `Bowel per day`,
        `Bowel frequency`,
        Consistency,
        Color,
        Amount
      ) |> 
      pivot_longer(
        c(
          Gender, 
          Breastfed, 
          Delivery, 
          # `Bowel per day`,
          `Bowel frequency`,
          Consistency,
          # Consistency2,
          Color,
          Amount
        ),
        names_to = "P",
        values_to = "V"
      )
  ),
  
  # use a chisquare test to determine significant over 
  # or under representation
  tar_target(
    category_cluster_test,
    unique(catset[['P']]) |>
      map_df(function(p){
        catset |> 
          filter(`P` == p) |> 
          group_by(Cluster, V) |>
          summarise(n = n_distinct(DNA_ID)) |>
          ungroup() |> 
          pivot_wider(
            id_cols = Cluster, 
            names_from = "V", 
            values_from = "n", 
            values_fill = list("n" = 0)
          ) |>
          (function(x){
            m <- as.matrix(x[,2:ncol(x)])
            rownames(m) <- x$Cluster
            m
          })() |>
          chisq.test(simulate.p.value = TRUE, rescale.p = TRUE) |>
          broom::tidy() |>
          mutate(
            category = p
          )
      }) |>
      mutate(
        method = str_replace_all(method, "[\t,\n]", " "),
        padj = p.adjust(p.value, method = "fdr")
      )
  ),
  
  # write the statistical results
  tar_target(
    path_category_cluster_test,
    tar_write_tsv(category_cluster_test, "output/chisquare_cluster_to_parameters.xls"),
    format = "file"
  ),
  
  #
  # Test continuous clinical parameters
  #

  # sanitize input
  tar_target(
    contset,
    tibble(
      Cluster = as.factor(hk),
      DNA_ID = names(hk)
    ) |>
      left_join(
        metadata, 
        by = "DNA_ID"
      ) |>
      select(
        DNA_ID,
        Cluster,
        PatientID,
        MilkVolume = MilkVolume_week,
        StoolFrequency = Average_StoolFreq 
      )
  ),

  # test using a kruskal test  
  tar_target(
    continuous_cluster_test,
    c("MilkVolume", "StoolFrequency") |>
      map_df(function(v){
        kruskal.test(
          contset[["Cluster"]], 
          contset[[v]]) |>
          tidy() |>
          mutate(
            variable = v
          )
      })
  ),
  
  # sanitize the output and perform multiple testing correction
  tar_target(
    cluster_test,
    bind_rows(
      category_cluster_test |> 
        select(variable = category, p.value),
      continuous_cluster_test |>
        select(variable, p.value)
    ) |>
      mutate(padj = p.adjust(p.value, method = "fdr"))
  ),
  
  tar_target(
    path_cluster_test,
    tar_write_tsv(cluster_test, "output/cluster_to_parameters.xls"),
    format = "file"
  ),
  
  #
  # Perform multidimensional scaling and retrieve the coordinates
  #
  
  tar_target(
    betamds,
    betadist %>%
      set_names(nm = todo) %>%
      map(function(x) {
        ordinate(phyloobj, "MDS", x)
      })
  ),
  
  tar_target(
    betacoord,
    1:length(betamds) %>%
      map_df(function(i) {
        x <- betamds[[i]]
        x$vectors %>%
          (function(m){
            left_join(
              tibble(
                DNA_ID = rownames(m), 
                as_tibble(m)), 
              metadata, by = 'DNA_ID')
          }) %>%
          mutate(
            method = names(betamds)[i]
          )
      }) %>%
      split(.$method) %>%
      map_df(function(d) {
        d %>%
          pivot_longer(
            starts_with('Axis'),
            names_to = 'axis',
            values_to = 'value'
          )
      })
  ),
  
  # write the coordinates for future reference
  tar_target(
    path_betacoord,
    tar_write_tsv(betacoord, "output/betacoord.tsv"),
    format = 'file'
  ),
  
  # visualize the results
  tar_target(
    plt_mds,
    betacoord %>%
      filter(method == "bray") %>%
      filter(axis %in% c("Axis.1", "Axis.2")) %>%
      pivot_wider(names_from = c("axis"), values_from = "value") %>%
      mutate(
        Subject = as.integer(str_remove(PatientTitle, '^ART-0+'))
      ) %>%
      ggplot(aes(x = Axis.1, y = Axis.2)) +
      geom_point(aes(colour = as.factor(Breastfed_yes), pch = as.factor(Gender_male)), size = 3) +
      geom_text_repel(aes(label = Subject), max.overlaps = NA) +
      guides(
        shape = guide_legend(title = "Male"),
        colour = guide_legend(title = "Breastfed")) +
      scale_colour_manual(values = wes_palette('Zissou1', 2, 'continuous')) + 
      theme_minimal() + 
      theme(
        aspect.ratio = 1,
        panel.border = element_rect(fill = NA, colour = "gray50")
      )
  ),
  
  tar_target(
    path_mds,
    tar_ggsave(plt_mds, "output/bray_mds_no_virus.pdf", width = 10, height = 9),
    format = "file"
  ),

  #
  # end of script
  #
  
  NULL
)
