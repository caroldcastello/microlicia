# Validate input files and compare species names with the MCC tree.
# Run from the repository root: source("scripts/00_validate_inputs.R")

needed <- c("ape", "yaml")
missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Install required packages first: ", paste(missing, collapse = ", "))
}

cfg <- yaml::read_yaml("config/analysis.yml")
paths <- unlist(cfg$files, use.names = TRUE)
not_found <- paths[!file.exists(paths)]
if (length(not_found)) {
  stop("Input file(s) not found: ", paste(not_found, collapse = ", "),
       ". See data/input/README.md and config/analysis.yml.")
}

tree <- ape::read.nexus(paths[["mcc_tree"]])
if (inherits(tree, "multiPhylo")) stop("Expected one MCC tree; multiple trees were read.")
if (is.null(tree$edge.length)) stop("The MCC tree has no branch lengths.")
if (anyDuplicated(tree$tip.label)) stop("The MCC tree contains duplicated tip labels.")

occ <- read.csv(paths[["occurrences"]], check.names = FALSE)
traits <- read.csv(paths[["stamen_traits"]], check.names = FALSE)
taxon_col <- cfg$columns$taxon
for (nm in c("occurrences", "stamen_traits")) {
  dat <- if (nm == "occurrences") occ else traits
  if (is.null(taxon_col) || !taxon_col %in% names(dat)) {
    stop("Taxon column '", taxon_col, "' is missing from ", nm, ".")
  }
  tx <- trimws(as.character(dat[[taxon_col]]))
  if (anyNA(tx) || any(!nzchar(tx))) stop("Missing/blank taxon names found in ", nm, ".")
}

for (col_key in c("longitude", "latitude")) {
  col <- cfg$columns[[col_key]]
  if (is.null(col) || !col %in% names(occ)) {
    stop("Occurrence column '", col, "' is missing; check config/analysis.yml.")
  }
  if (!is.numeric(occ[[col]])) {
    stop("Occurrence column '", col, "' must be numeric: ", col)
  }
}

occ_taxa <- unique(trimws(as.character(occ[[taxon_col]])))
trait_taxa <- unique(trimws(as.character(traits[[taxon_col]])))
tree_taxa <- tree$tip.label
compare_taxa <- function(data_taxa, tree_taxa, label) {
  data.frame(dataset = label, n_taxa = length(data_taxa),
             matched_to_tree = sum(data_taxa %in% tree_taxa),
             absent_from_tree = sum(!data_taxa %in% tree_taxa))
}

out_dir <- "results/00_input_checks"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(rbind(compare_taxa(occ_taxa, tree_taxa, "occurrences"),
                compare_taxa(trait_taxa, tree_taxa, "stamen_traits")),
          file.path(out_dir, "taxon_match_summary.csv"), row.names = FALSE)
write.csv(data.frame(taxon = setdiff(occ_taxa, tree_taxa)),
          file.path(out_dir, "occurrence_taxa_absent_from_tree.csv"), row.names = FALSE)
write.csv(data.frame(taxon = setdiff(trait_taxa, tree_taxa)),
          file.path(out_dir, "trait_taxa_absent_from_tree.csv"), row.names = FALSE)
write.csv(data.frame(taxon = setdiff(tree_taxa, union(occ_taxa, trait_taxa))),
          file.path(out_dir, "tree_taxa_without_occurrence_or_trait_data.csv"), row.names = FALSE)

tip_depths <- ape::node.depth.edgelength(tree)[seq_along(tree_taxa)]
write.csv(data.frame(n_tips = length(tree_taxa), has_branch_lengths = TRUE,
                     is_ultrametric = ape::is.ultrametric(tree),
                     root_to_tip_range = diff(range(tip_depths))),
          file.path(out_dir, "mcc_tree_check.csv"), row.names = FALSE)
message("Input checks complete. Review files in ", out_dir,
        " and resolve name mismatches before analysis.")
