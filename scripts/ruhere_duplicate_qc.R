# Optional geographic QC: flag duplicate species-coordinate records without deleting them.
# Run from repository root after installing ape and RuHere.

required <- c("ape", "RuHere", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install packages first: ", paste(missing, collapse = ", "))

input_occ <- "data/input/ocorrencias.xlsx.csv"
input_tree <- "data/input/teste_datacao_MCC.tree"
out_csv <- "results/ruhere_duplicate_qc.csv"
out_json <- "results/ruhere_duplicate_qc_summary.json"
dir.create("results", showWarnings = FALSE, recursive = TRUE)

occ <- read.csv(input_occ, sep = ";", fileEncoding = "CP1252", check.names = FALSE,
                stringsAsFactors = FALSE)
tree <- ape::read.nexus(input_tree)
if (inherits(tree, "multiPhylo")) stop("Expected one MCC tree.")
if (!all(c("taxon_tree", "Latitude_decimal", "Longitude_decimal") %in% names(occ))) {
  stop("Occurrence table must include taxon_tree and decimal latitude/longitude columns.")
}
if (any(!occ$taxon_tree %in% tree$tip.label)) {
  stop("At least one taxon_tree label does not match a tree tip exactly.")
}
occ$Latitude_decimal <- suppressWarnings(as.numeric(occ$Latitude_decimal))
occ$Longitude_decimal <- suppressWarnings(as.numeric(occ$Longitude_decimal))
occ$coordinate_valid <- is.finite(occ$Latitude_decimal) & is.finite(occ$Longitude_decimal) &
  occ$Latitude_decimal >= -90 & occ$Latitude_decimal <= 90 &
  occ$Longitude_decimal >= -180 & occ$Longitude_decimal <= 180
occ$duplicate_coordinate_flag <- NA
occ$.source_row <- seq_len(nrow(occ))
valid <- occ[occ$coordinate_valid, , drop = FALSE]

if (nrow(valid)) {
  checked <- RuHere::flag_duplicates(
    occ = valid,
    species = "taxon_tree",
    long = "Longitude_decimal",
    lat = "Latitude_decimal"
  )
  # RuHere's duplicated_flag is TRUE for the retained representative and FALSE for
  # additional records with the same species and coordinates. Keep both classes.
  occ$duplicate_coordinate_flag[checked$.source_row] <- checked$duplicated_flag
  n_flagged <- sum(!checked$duplicated_flag, na.rm = TRUE)
} else {
  n_flagged <- 0L
}
occ$.source_row <- NULL
write.csv(occ, out_csv, row.names = FALSE, fileEncoding = "UTF-8")
summary <- list(
  package = "RuHere::flag_duplicates",
  n_occurrence_rows = nrow(occ),
  n_rows_with_numerically_valid_coordinates = sum(occ$coordinate_valid),
  n_duplicate_coordinate_rows_flagged = n_flagged,
  rows_deleted = 0L,
  note = "Flags are a screening aid; inspect collection metadata and vouchers before removing records."
)
jsonlite::write_json(summary, out_json, pretty = TRUE, auto_unbox = TRUE)
message("RuHere screening complete. No input rows were removed.")
