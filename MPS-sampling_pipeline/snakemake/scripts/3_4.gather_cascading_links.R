# 3_4.gather_cascading_links.R
# Build MPS-clusters through hierarchical clustering with complete-linkage.


# save.image("~/3_4.gather_cascading_links.RDS")
# stop()
# load("~/3_4.gather_cascading_links.RDS")


##### Input #####


library(magrittr)
print("sessionInfo()")
print(sessionInfo())
print(".libPaths()")
print(.libPaths())
print("loadedNamespeces()")
print(loadedNamespaces())
cat("\n")


cores <- snakemake@params[["cores"]]
cat("Parallelization on", cores, "cores.\n")


if (cores > 1) {
  doParallel::registerDoParallel(cores = cores)
  parallel <- TRUE
} else {
  parallel <- FALSE
}


gen_Lincomb_links <- snakemake@input[["gen_Lincomb_links"]] %>%
  data.table::fread()
Lincomb_MPS_cluster_links <- snakemake@input[["Lincomb_MPS_cluster_links"]] %>%
  data.table::fread()


##### Compute ######


# Simply merge all <data.tables>.
genome_cascading_links <- data.table::merge.data.table(
  x = gen_Lincomb_links,
  y = Lincomb_MPS_cluster_links,
  by = "Lincombination_ID",
  all = TRUE,
  order = TRUE)
# Re-order the columns.
genome_cascading_links <- genome_cascading_links[
  , .(genomeAccession,
      Lincombination_ID,
      MPS_cluster_ID)]


##### Output #####


data.table::fwrite(
  x = genome_cascading_links,
  file = snakemake@output[["genome_cascading_links"]],
  sep = ";")
