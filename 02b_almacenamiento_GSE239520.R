# 02b — Storage time before DNA extraction in the third cohort (GSE239520)
#
# GSE239520 deposits, for each sample, the number of days it was left at room temperature before DNA
# extraction. This script reads that field from the SOFT file, describes its distribution by group, tests
# whether long storage (> 300 days) is associated with group (two-sided Fisher exact test) and writes the
# table. The variable enters no model; its effect on the locus is examined by exclusion in 02c.
# Paper section: Discussion 4.1; Supplementary Table S8
# Input:  GSE239520 SOFT file (CM)
# Output: outputs/S4_almacenamiento_GSE239520.csv (Supplementary Table S8)

# ---- 1. Packages and paths --------------------------------------------------------------------------
suppressMessages(library(data.table))
CM  <- "data/evaluation"                                     # GEO files (edit to the local path)
OUT <- "outputs"
dir.create(OUT, showWarnings = FALSE)

# ---- 2. Phenotype from the SOFT file ----------------------------------------------------------------
li  <- gsub("\r", "", readLines(file.path(CM, "GSE239520.txt")))   # SOFT files carry Windows line endings
idx <- grep("^\\^SAMPLE", li)                                # first line of each sample block
gsm <- sub(".*= ", "", grep("^!Sample_geo_accession", li, value = TRUE))
blo <- cut(seq_along(li), breaks = c(idx, length(li) + 1), right = FALSE, labels = FALSE)   # sample of each line
ch  <- data.table(b = blo[grep("^!Sample_characteristics_ch1", li)],
                  t = sub("^!Sample_characteristics_ch1 = ", "",
                          grep("^!Sample_characteristics_ch1", li, value = TRUE)))
ch[, campo := trimws(sub(":.*", "", t))]                     # field name
ch[, valor := trimws(sub("^[^:]*:", "", t))]                 # field value
ch[, geo := gsm[b]]
ph <- as.data.table(dcast(ch, geo ~ campo, value.var = "valor"))   # one row per GSM
setnames(ph, make.names(names(ph)))                          # syntactically valid column names
ph[, malt := ifelse(grepl("maltreat", group, ignore.case = TRUE), 1,     # 1 = maltreatment
             ifelse(grepl("typically", group, ignore.case = TRUE), 0, NA))]   # 0 = comparison
cat(sprintf("GSE239520: %d muestras | %d maltrato / %d comparación\n",
            nrow(ph), sum(ph$malt == 1, na.rm = TRUE), sum(ph$malt == 0, na.rm = TRUE)))

# ---- 3. Storage time by group and Supplementary Table S8 --------------------------------------------
dias <- suppressWarnings(as.numeric(ph[["days.left.at.rt.before.the.dna.extraction"]]))
if (all(is.na(dias))) {                                      # the field name may vary across GEO versions
  cand <- grep("days", names(ph), value = TRUE, ignore.case = TRUE)
  stopifnot(length(cand) == 1)                               # it must be unambiguous
  dias <- suppressWarnings(as.numeric(ph[[cand[1]]]))
}
stopifnot(!all(is.na(dias)))
alm <- data.table(geo = ph$geo, malt = ph$malt, dias = dias)   # one row per sample
alm <- alm[!is.na(malt) & !is.na(dias)]                      # complete cases only
alm[, largo := dias > 300]                                   # long storage: no values fall between 8 and 386
alm[, gsmn := as.integer(sub("^GSM", "", geo))]              # accession number, to order the series
tab <- table(alm$malt, alm$largo)                            # group by long storage
fis <- fisher.test(tab)                                      # two-sided Fisher exact test
cat("\n=== días a temperatura ambiente antes de la extracción ===\n")
for (g in c(1, 0)) {                                         # maltreatment, then comparison
  d <- sort(alm[malt == g, dias])
  cat(sprintf("  %-12s n=%2d | rango %d-%d | almacenamiento largo (>300 d): %d\n",
              ifelse(g == 1, "maltrato", "comparación"), length(d),
              min(d), max(d), sum(d > 300)))
}
h <- sort(alm$dias)
cat(sprintf("  hueco sin observaciones: %d a %d días\n",
            h[sum(h < 300)], h[sum(h < 300) + 1]))           # edges of the gap
pos <- match(alm[largo == TRUE, gsmn], sort(alm$gsmn))       # position of the long-storage samples in the series
cat(sprintf("  los largos caen entre los accesos %d y %d de la serie (de %d)\n",
            min(pos), max(pos), nrow(alm)))
cat(sprintf("  máximo de días después de ese tramo: %d\n",
            max(alm[match(gsmn, sort(alm$gsmn)) > max(pos), dias])))   # no later sample exceeds 7 days
cat(sprintf("  Fisher exacto bilateral: p = %.3g\n", fis$p.value))
fwrite(alm[order(gsmn)], file.path(OUT, "S4_almacenamiento_GSE239520.csv"))
cat("\nEscrito:", file.path(OUT, "S4_almacenamiento_GSE239520.csv"), "\n")
