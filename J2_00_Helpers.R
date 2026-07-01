# ============================================================
# J2_00_Helpers.R
# SHARED HELPERS FOR THE Journal_2_rev PIPELINE
# ============================================================
# Source this file at the top of EVERY J2_0X script in this pipeline:
#   source("J2_00_Helpers.R")
#
# Provides:
#   - Output directory setup under "Journal_2_rev/"
#   - theme_j2_en(): ggplot theme - English, NO title/subtitle, black
#     (not grey) italic footnote caption containing method notes only
#   - save_png_j2(): consistent PNG export
#   - wrap_caption(): wraps long caption strings onto multiple lines so
#     they are not clipped at the right edge of the saved PNG
#   - finalize_outputs(): writes EVERY table in a script to THREE formats
#     at once - .xlsx (multi-sheet), .tex (booktabs LaTeX), .docx (Word,
#     via officer::body_add_table + safe_df(), per project convention -
#     flextable is never used)
#
# Conventions enforced pipeline-wide:
#   - All plot text, table headers, and table content: English
#   - No plot titles/subtitles; only a bottom caption with the
#     statistical method, in black (#000000), never grey
#   - safe_df() always applied before officer::body_add_table()
# ============================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(openxlsx)
  library(officer)
})

# ---- Output directories ----------------------------------------------
DIR_ROOT  <- "Journal_2_rev"
DIR_XLSX  <- file.path(DIR_ROOT, "xlsx")
DIR_PLOTS <- file.path(DIR_ROOT, "plots")
DIR_RDS   <- file.path(DIR_ROOT, "rds")
DIR_TEX   <- file.path(DIR_ROOT, "tex")
DIR_WORD  <- file.path(DIR_ROOT, "word")
for (d in c(DIR_ROOT, DIR_XLSX, DIR_PLOTS, DIR_RDS, DIR_TEX, DIR_WORD)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---- Plot theme: English, no titles, black footnote -------------------
theme_j2_en <- function() {
  theme_minimal(base_size = 11) +
    theme(
      plot.title       = element_blank(),
      plot.subtitle    = element_blank(),
      plot.caption     = element_text(size = 8, color = "#000000",
                                       face = "italic", hjust = 0,
                                       margin = margin(t = 8)),
      axis.text        = element_text(color = "#000000"),
      axis.title       = element_text(color = "#000000"),
      legend.text      = element_text(color = "#000000"),
      legend.title     = element_text(color = "#000000"),
      strip.text       = element_text(color = "#000000"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3)
    )
}

save_png_j2 <- function(p, nm, w = 8, h = 6) {
  path <- file.path(DIR_PLOTS, paste0(nm, ".png"))
  tryCatch({
    ggsave(path, p, width = w, height = h, dpi = 300, bg = "white")
    cat(sprintf("  OK  %s\n", path))
  }, error = function(e) {
    while (dev.cur() > 1) dev.off()
    cat(sprintf("  FAILED %s: %s\n", nm, e$message))
  })
}

# ---- wrap_caption(): prevents long plot captions from being clipped at
# the right edge of the saved PNG. ggplot2 does not auto-wrap labs(caption=)
# text, so any caption longer than the plot's rendered width at the given
# font size will overflow past the panel and be cut off rather than wrapped.
# Use as: labs(..., caption = wrap_caption("long caption text here"))
wrap_caption <- function(text, width = 95) {
  paste(strwrap(text, width = width), collapse = "\n")
}

# ---- safe_df(): mandatory sanitizer before officer::body_add_table ----
safe_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(col) {
    if (is.numeric(col)) {
      out <- ifelse(is.na(col), "", formatC(col, digits = 4, format = "f"))
    } else {
      out <- ifelse(is.na(col), "", as.character(col))
    }
    out
  })
  df
}

# ---- LaTeX helpers ------------------------------------------------------
latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%$&_#{}])", "\\\\\\1", x)
  x
}

latex_table_block <- function(df, caption, label, digits = 4) {
  df_fmt <- df
  for (cl in names(df_fmt)) {
    if (is.numeric(df_fmt[[cl]])) {
      df_fmt[[cl]] <- formatC(df_fmt[[cl]], digits = digits, format = "f")
    }
  }
  ncol_t <- ncol(df_fmt)
  align  <- paste(rep("l", ncol_t), collapse = "")
  header <- paste(latex_escape(names(df_fmt)), collapse = " & ")
  body_lines <- apply(df_fmt, 1, function(r) paste(latex_escape(r), collapse = " & "))
  c(
    "\\begin{table}[ht]",
    "\\centering",
    sprintf("\\caption{%s}", caption),
    sprintf("\\label{%s}", label),
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(body_lines, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

# ---- Per-format writers --------------------------------------------------
write_xlsx_multi <- function(table_list, path) {
  wb <- createWorkbook()
  hs <- createStyle(textDecoration = "bold", fgFill = "#1F4E79",
                     fontColour = "white", halign = "center")
  for (nm in names(table_list)) {
    sheet_nm <- substr(gsub("[^A-Za-z0-9_]", "_", nm), 1, 31)
    addWorksheet(wb, sheet_nm)
    df <- table_list[[nm]]$df
    writeData(wb, sheet_nm, df)
    addStyle(wb, sheet_nm, hs, rows = 1, cols = 1:ncol(df), gridExpand = TRUE)
  }
  saveWorkbook(wb, path, overwrite = TRUE)
  cat(sprintf("  OK  %s\n", path))
}

write_tex_multi <- function(table_list, path) {
  all_lines <- c(
    "% Auto-generated by Journal_2_rev pipeline.",
    "% Requires \\usepackage{booktabs} in the LaTeX preamble.",
    ""
  )
  for (nm in names(table_list)) {
    tbl <- table_list[[nm]]
    all_lines <- c(all_lines,
                    latex_table_block(tbl$df, tbl$caption, paste0("tab:", nm)), "")
  }
  writeLines(all_lines, path)
  cat(sprintf("  OK  %s\n", path))
}

write_docx_multi <- function(table_list, path, title) {
  doc <- officer::read_docx()
  doc <- officer::body_add_par(doc, title, style = "heading 1")

  # "Table Grid" is not guaranteed to exist in every officer/Word locale's
  # blank template. Detect a usable table style once; fall back to no
  # style (officer's plain default table) if none is found, rather than
  # letting the whole export fail.
  available_styles <- tryCatch(officer::styles_info(doc)$style_name, error = function(e) character(0))
  table_style <- if ("Table Grid" %in% available_styles) "Table Grid" else NA_character_

  for (nm in names(table_list)) {
    tbl <- table_list[[nm]]
    doc <- officer::body_add_par(doc, tbl$caption, style = "heading 2")
    doc <- tryCatch({
      if (!is.na(table_style)) {
        officer::body_add_table(doc, safe_df(tbl$df), style = table_style)
      } else {
        officer::body_add_table(doc, safe_df(tbl$df))
      }
    }, error = function(e) {
      cat(sprintf("  NOTE: body_add_table styled insertion failed (%s); retrying without style.\n", e$message))
      officer::body_add_table(doc, safe_df(tbl$df))
    })
    doc <- officer::body_add_par(doc, "", style = "Normal")
  }
  print(doc, target = path)
  cat(sprintf("  OK  %s\n", path))
}

# ---- Master export: call ONCE at the end of every script ----------------
# table_list: named list, each element = list(df = <data.frame>, caption = "English caption")
finalize_outputs <- function(table_list, script_tag, title) {
  cat(sprintf("\n--- Exporting all tables for %s (.xlsx + .tex + .docx) ---\n", script_tag))
  write_xlsx_multi(table_list, file.path(DIR_XLSX, paste0(script_tag, ".xlsx")))
  write_tex_multi(table_list,  file.path(DIR_TEX,  paste0(script_tag, ".tex")))
  write_docx_multi(table_list, file.path(DIR_WORD, paste0(script_tag, ".docx")), title)
}

# ---- Shared closed-form TE estimator (Battese & Coelli, 1988) -----------
# Used as a robust fallback whenever frontier::efficiencies() fails or
# returns an object of unexpected length - this happens occasionally when
# the MLE covariance matrix is ill-behaved, and failing SILENTLY (as a
# bare tryCatch(..., error = function(e) NULL) would) is exactly the kind
# of bug this helper exists to prevent pipeline-wide.
compute_TE_BC88 <- function(e, sigmaSq, gamma) {
  sigma_u2    <- gamma * sigmaSq
  sigma_v2    <- (1 - gamma) * sigmaSq
  sigma_star2 <- sigma_u2 * sigma_v2 / sigmaSq
  sigma_star  <- sqrt(sigma_star2)
  mu_star     <- -e * sigma_u2 / sigmaSq
  z           <- mu_star / sigma_star
  exp(-mu_star + 0.5 * sigma_star2) * pnorm(z - sigma_star) / pnorm(z)
}

# ---- Misc shared statistical helpers -------------------------------------
sig_label <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  if (p < 0.10)  return(".")
  return("ns")
}

cat("J2_00_Helpers.R loaded. Output root:", DIR_ROOT, "\n")
cat("  wrap_caption() available:", exists("wrap_caption"), "\n")
