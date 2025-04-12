# GEO parser script to run in parallel on hundreds of GSE datasets
# Save this as: parse_geo_batch.R

suppressPackageStartupMessages({
  library(GEOquery)
  library(tidyverse)
  library(furrr)
  library(future)
})

plan(multisession, workers = as.integer(Sys.getenv("NSLOTS")))

# ---------- Setup ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript parse_geo_batch.R gse_list.txt output_dir")

gse_list_file <- args[1]
output_base_dir <- args[2]
if (!dir.exists(output_base_dir)) dir.create(output_base_dir, recursive = TRUE)

gse_ids <- read_lines(gse_list_file)

# ---------- Functions ----------
experimentData_to_df <- function(gse_id) {
  tryCatch({
    expData_obj <- get(paste0(gse_id, "_exp"))
    metadata_list <- list()
    for (slot_name in slotNames(expData_obj)) {
      value <- slot(expData_obj, slot_name)
      if (is.list(value)) value <- paste(unlist(value), collapse = "; ")
      metadata_list[[slot_name]] <- value
    }
    if (!is.null(expData_obj@other) && length(expData_obj@other) > 0) {
      other_flat <- lapply(expData_obj@other, function(x) if (is.list(x)) paste(unlist(x), collapse = "; ") else x)
      metadata_list <- c(metadata_list, other_flat)
    }
    metadata_df <- as.data.frame(metadata_list, stringsAsFactors = FALSE) %>%
      rename(gse_id = geo_accession) %>%
      relocate(gse_id)
  }, error = function(e) {
    message("Failed to extract dataset-level metadata for ", gse_id, ": ", conditionMessage(e))
    tibble(gse_id = gse_id)
  })
}

extract_sample_metadata <- function(pheno_df, gse_id) {
  tryCatch({
    base_fields <- c("geo_accession", "source_name_ch1", "organism_ch1")
    keywords <- c("tissue", "organ", "source", "age", "sex", "gender", "disease", "state", "condition", "brain", "region")
    all_colnames <- colnames(pheno_df)
    characteristics_cols <- all_colnames[str_detect(all_colnames, "characteristics")]
    first_row_values <- pheno_df %>% slice(1) %>% select(all_of(characteristics_cols)) %>% mutate(across(everything(), ~tolower(as.character(.)))) %>% unlist()
    matched_cols <- names(first_row_values)[
      map_lgl(seq_along(characteristics_cols), function(i) {
        str_detect(first_row_values[i], paste(keywords, collapse = "|"))
      })
    ]
    supplementary_cols <- all_colnames[str_detect(all_colnames, "supplementary")]
    selected_cols <- intersect(unique(c(base_fields, matched_cols, supplementary_cols)), colnames(pheno_df))
    
    pheno_df %>% select(all_of(selected_cols)) %>%
      mutate(sample_row = row_number()) %>%
      pivot_longer(cols = contains("characteristics_ch1"), names_to = "char_col", values_to = "char_value") %>%
      separate(char_value, into = c("key", "value"), sep = ":", extra = "merge", fill = "right") %>%
      mutate(key = str_trim(key) %>% str_to_lower() %>% str_replace_all(" ", "_"),
             value = str_trim(value)) %>%
      filter(!is.na(key)) %>%
      pivot_wider(id_cols = sample_row, names_from = key, values_from = value) %>%
      left_join(pheno_df %>% mutate(sample_row = row_number()) %>% select(sample_row, any_of(c(base_fields, supplementary_cols))), by = "sample_row") %>%
      select(any_of(c(base_fields, supplementary_cols)), everything(), -sample_row) %>%
      relocate(any_of(base_fields)) %>%
      mutate(gse_id = gse_id) %>%
      relocate(gse_id) %>%
      rename(gsm_id = geo_accession)
  }, error = function(e) {
    message("Failed to extract sample metadata for ", gse_id, ": ", conditionMessage(e))
    tibble(gse_id = gse_id)
  })
}

parse_single_gse <- function(gse_id, base_dir) {
  dest_dir <- file.path(base_dir, gse_id)
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)
  
  message("Processing ", gse_id)
  gse <- tryCatch({
    getGEO(GEO = gse_id, destdir = dest_dir, GSEMatrix = TRUE, getGPL = FALSE, parseCharacteristics = TRUE)
  }, error = function(e) {
    message("Failed to download ", gse_id, ": ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(gse)) return(NULL)
  
  eset <- gse[[1]]
  assign(paste0(gse_id, "_eset"), eset, envir = .GlobalEnv)
  assign(paste0(gse_id, "_exp"), experimentData(eset), envir = .GlobalEnv)
  assign(paste0(gse_id, "_pheno"), pData(eset), envir = .GlobalEnv)
  assign(paste0(gse_id, "_feature"), fData(eset), envir = .GlobalEnv)
  assign(paste0(gse_id, "_expr"), exprs(eset), envir = .GlobalEnv)
  
  assign(paste0(gse_id, "_exp_df"), experimentData_to_df(gse_id), envir = .GlobalEnv)
  assign(paste0(gse_id, "_exp_df_selected"), get(paste0(gse_id, "_exp_df")) %>% select(any_of(c("gse_id", "title", "abstract", "url", "pubMedIds", "overall_design", "sample_id", "status", "summary", "supplementary_file"))), envir = .GlobalEnv)
  
  assign(paste0(gse_id, "_pheno_filtered"), extract_sample_metadata(get(paste0(gse_id, "_pheno")), gse_id), envir = .GlobalEnv)
  
  save(list = c(
    paste0(gse_id, "_eset"), paste0(gse_id, "_exp"), paste0(gse_id, "_pheno"), paste0(gse_id, "_feature"), paste0(gse_id, "_expr"),
    paste0(gse_id, "_exp_df"), paste0(gse_id, "_exp_df_selected"), paste0(gse_id, "_pheno_filtered")),
    file = file.path(dest_dir, paste0(gse_id, "_objects.RData"))
  )
  
  write_csv(get(paste0(gse_id, "_exp_df_selected")), file.path(dest_dir, paste0(gse_id, "_dataset_metadata.csv")))
  write_csv(get(paste0(gse_id, "_pheno_filtered")), file.path(dest_dir, paste0(gse_id, "_sample_metadata.csv")))
  
  return(gse_id)
}

# ---------- Run Batch ----------

# Run and track both success and failure
results <- future_map(gse_ids, function(gse_id) {
  tryCatch({
    parse_single_gse(gse_id, base_dir = output_base_dir)
    list(gse_id = gse_id, status = "success", error = NA)
  }, error = function(e) {
    list(gse_id = gse_id, status = "failed", error = conditionMessage(e))
  })
})

# Convert list to tibble
results_tbl <- bind_rows(results)

# Write success list
results_tbl %>%
  filter(status == "success") %>%
  pull(gse_id) %>%
  write_lines(file.path(output_base_dir, "success_gse_ids.txt"))

# Write failure list
results_tbl %>%
  filter(status == "failed") %>%
  mutate(error = replace_na(error, "unknown error")) %>%
  write_csv(file.path(output_base_dir, "failed_gse_ids.csv"))


