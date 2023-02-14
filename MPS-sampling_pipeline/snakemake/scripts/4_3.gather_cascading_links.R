# 4_3.gather_cascading_links.R
# Gather cascading genome clusterings into one single file.


# save.image("~/4_3.gather_cascading_links.R")
# stop()
# load("~/4_3.gather_cascading_links.R")


##### Input #####


library(magrittr)


gen_Lincomb_links <- snakemake@input[["gen_Lincomb_links"]] %>%
  data.table::fread()
Lincomb_precomp_links <- snakemake@input[["Lincomb_precomp_links"]] %>%
  data.table::fread()
Lincomb_MPS_cluster_links <- snakemake@input[["Lincomb_MPS_cluster_links"]] %>%
  data.table::fread()


##### genome_cascading_links ######


# Simply merge all <data.tables>.
genome_cascading_links <- data.table::merge.data.table(
  x = gen_Lincomb_links,
  y = Lincomb_precomp_links,
  by = "Lincombination_ID",
  all = TRUE,
  order = TRUE)
genome_cascading_links <- data.table::merge.data.table(
  x = genome_cascading_links,
  y = Lincomb_MPS_cluster_links,
  by = "Lincombination_ID",
  all = TRUE,
  order = TRUE)
# Re-order the columns.
genome_cascading_links <- genome_cascading_links[
  , .(genomeAccession,
      Lincombination_ID,
      preconnected_component_ID,
      MPS_cluster_ID)]
# Fill missing values among each genome clustering level.
#   Put the opposite of the combination ID if needed.
genome_cascading_links[is.na(preconnected_component_ID),
                       preconnected_component_ID := -Lincombination_ID]
genome_cascading_links[MPS_cluster_ID == "",
                       MPS_cluster_ID := -Lincombination_ID]


##### Output #####


data.table::fwrite(
  x = genome_cascading_links,
  file = snakemake@output[["genome_cascading_links"]],
  sep = ";")


