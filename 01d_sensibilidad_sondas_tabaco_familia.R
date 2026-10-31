# 01d — Three further sensitivity analyses of the discovery cohort
#
# Reproduces the loading of 01b (same probe filter, same 116 participants, same design matrix) and adds:
#   S7  the genome-wide analysis restricted to the probes with no missing value, so no value is imputed;
#   S8  a methylation-derived smoking score (EpiSmokEr; Bollepalli et al., 2019) and, in parallel, the
#       methylation of cg05575921 in AHRR as covariates, because smoking was not recorded;
#   S9  one participant per extended family (the member closest in age to the median of the replication
#       cohort), which removes the family structure.
# Each analysis reports the genome-wide count and the eight ZG16B promoter probes.
# Paper section: Methods 2.6; Results 3.7; Supplementary Tables S4, S5 and S7
# Input:  GSE226085 processed matrix and SOFT file (BU); GSE239521 SOFT file (CM);
#         outputs/T3_efectos_cohorte1.csv (the 22 probes of the primary analysis, written by 01b)
# Output: outputs/S7_sondas_completas_{resumen,locus}.csv (Supplementary Table S4),
#         S8_tabaco_{resumen,locus}.csv and S8_puntaje_tabaco_por_muestra.csv (Supplementary Table S5),
#         S9_una_por_familia_{resumen,locus,seleccion}.csv (Supplementary Table S7), sessionInfo_01d.txt

# ---- 1. Packages, paths and locus -------------------------------------------------------------------
suppressMessages({
  library(data.table); library(limma); library(EpiDISH); library(splines)
  library(impute)                                            # k-nearest-neighbor imputation
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})
set.seed(42)                                                 # same seed as 01 and 01b
BU  <- "data/GSE226085"                                      # GSE226085 (edit to the local path)
CM  <- "data/evaluation"                                     # SOFT file of the replication cohort (edit to the local path)
OUT <- "outputs"
dir.create(OUT, showWarnings = FALSE)

LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",   # TSS200 -> negative
          "cg19128189","cg08829024","cg12475088","cg04543077"),  # TSS1500 -> positive
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))
AHRR <- "cg05575921"                                         # smoking probe in AHRR (Hernández et al., 2023)

# ---- 2. Phenotypes from the GEO SOFT files (discovery and replication cohorts) ----------------------
fenotipos <- list()
for (g in c("GSE226085", "GSE239521")) {
  fichero <- if (g == "GSE226085") file.path(BU, "GSE226085.txt") else file.path(CM, "GSE239521.txt")
  li  <- gsub("\r", "", readLines(fichero))                  # SOFT files carry Windows line endings
  idx <- grep("^\\^SAMPLE", li)                              # first line of each sample block
  gsm <- sub(".*= ", "", grep("^!Sample_geo_accession", li, value = TRUE))
  blo <- cut(seq_along(li), breaks = c(idx, length(li) + 1), # sample of each line
             right = FALSE, labels = FALSE)
  ch  <- data.table(b = blo[grep("^!Sample_characteristics_ch1", li)],
                    t = sub("^!Sample_characteristics_ch1 = ", "",
                            grep("^!Sample_characteristics_ch1", li, value = TRUE)))
  ch[, campo := trimws(sub(":.*", "", t))]
  ch[, valor := trimws(sub("^[^:]*:", "", t))]
  ch[, geo   := gsm[b]]
  ph <- as.data.table(dcast(ch, geo ~ campo, value.var = "valor"))   # long to wide: one row per GSM
  setnames(ph, make.names(names(ph)))
  fenotipos[[g]] <- ph[]
}

# ---- 3. Discovery cohort, loaded as in 01b ----------------------------------------------------------
# The full raw matrix of the analytic sample (Braw_ok) is also kept, because EpiSmokEr needs probes that
# the missingness filter may have removed
ph <- fenotipos[["GSE226085"]]
ph[, `:=`(violencia = suppressWarnings(as.numeric(violence.trauma)), edad = suppressWarnings(as.numeric(age)),
          sexo = gender, familia = extended.family,
          lote = plate, col = individual, gen = generation)]
ph[, usar := grepl("^Y|^y|TRUE|1", passed.methylation.quality.control) &   # authors' quality control
             grepl("^Y|^y|TRUE|1", included.in.final.analytical.sample)]   # and their final analytic sample
cab <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"), nrows = 0)
B <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"),
           select = setdiff(names(cab), grep("Detection", names(cab), value = TRUE)),   # detection p columns not loaded
           showProgress = FALSE)
cg <- B[[1]]; B <- as.matrix(B[, -1]); rownames(B) <- cg
n_inicial <- nrow(B)                                         # probe count before filtering
com <- intersect(colnames(B), ph$col)
ph  <- ph[match(com, col)]; B <- B[, com, drop = FALSE]      # phenotype aligned with the matrix
stopifnot(identical(colnames(B), ph$col))
data(centEpiFibIC.m)                                         # EpiDISH epithelial/fibroblast/immune reference
Bf <- B[rowMeans(is.na(B)) < 0.05 & !grepl("^rs|^ch\\.", rownames(B)), , drop = FALSE]   # article filter
ph[, Epi := epidish(beta.m = Bf, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
ok <- ph$usar & !is.na(ph$violencia) & !is.na(ph$edad) & !is.na(ph$Epi)   # analytic sample
ok[is.na(ok)] <- FALSE
cat(sprintf("Encadenado de sondas: %d -> %d | muestra analítica n = %d | %d familias\n",
            n_inicial, nrow(Bf), sum(ok), uniqueN(ph$familia[ok])))
stopifnot(nrow(Bf) == 736065, sum(ok) == 116)                # figures reported in the article
Braw_ok <- B[, ok, drop = FALSE]                             # full raw matrix, analytic sample (for S8)
Bok     <- Bf[, ok, drop = FALSE]                            # filtered, with its missing values (for S7 and KNN)
rm(B, Bf); gc()
d <- data.frame(y    = ph$violencia[ok], edad = ph$edad[ok],
                sexo = droplevels(factor(ph$sexo[ok])), Epi = ph$Epi[ok],
                lote = droplevels(factor(ph$lote[ok])))
dis <- model.matrix(~ y + ns(edad, 3) + sexo + lote + ns(Epi, 4), data = d)   # design of the article
stopifnot(qr(dis)$rank == ncol(dis))
primarios <- fread(file.path(OUT, "T3_efectos_cohorte1.csv"))$cpg   # the 22 probes of the primary analysis
cat(sprintf("CpGs significativos del análisis primario cargados: %d\n", length(primarios)))

# ---- 4. S7: probes with complete data only (no imputation) ------------------------------------------
Bcc <- Bok[rowSums(is.na(Bok)) == 0, , drop = FALSE]         # probes without any missing value
cat(sprintf("S7 · sondas completas: %d de %d (%.1f%%)\n", nrow(Bcc), nrow(Bok), 100 * nrow(Bcc) / nrow(Bok)))
stopifnot(all(LOCUS$cpg %in% rownames(Bcc)))                 # the locus must be complete
Mcc <- pmin(pmax(Bcc, 1e-4), 1 - 1e-4); Mcc <- log2(Mcc / (1 - Mcc))   # M-values
fit7  <- eBayes(lmFit(Mcc, dis)); tt7 <- topTable(fit7, coef = 2, number = Inf, sort.by = "none")
fit7b <- eBayes(lmFit(Bcc, dis)); tb7 <- topTable(fit7b, coef = 2, number = Inf, sort.by = "none")   # same model on beta: delta-beta
hits7 <- rownames(tt7)[tt7$adj.P.Val < 0.05]
res7 <- data.table(analisis = "sondas completas", n_sondas = nrow(Bcc),
                   lambda = round(median(qchisq(tt7$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1), 3),
                   n_fdr05 = length(hits7),
                   primarios_presentes = sum(primarios %in% rownames(Bcc)),
                   primarios_recuperados = sum(primarios %in% hits7),
                   zg16b_fdr05 = sum(LOCUS$cpg %in% hits7))
cat("\n=== S7 · resumen ===\n"); print(res7)
loc7 <- merge(LOCUS, data.table(cpg = rownames(tt7), M_est = tt7$logFC, p = tt7$P.Value, fdr = tt7$adj.P.Val), by = "cpg")
loc7 <- merge(loc7, data.table(cpg = rownames(tb7), dBeta = tb7$logFC), by = "cpg")
loc7[, `:=`(dir_ok = sign(M_est) == dir_pred, analisis = "sondas completas")]; loc7 <- loc7[order(region, cpg)]
cat("\n=== S7 · el locus de ZG16B con solo sondas completas ===\n"); print(loc7)
fwrite(res7, file.path(OUT, "S7_sondas_completas_resumen.csv"))
fwrite(loc7, file.path(OUT, "S7_sondas_completas_locus.csv"))
rm(Bcc, Mcc, fit7, fit7b); gc()

# ---- 5. KNN imputation, as in 01b, for S8 and S9 ----------------------------------------------------
t0 <- Sys.time()
invisible(capture.output(Bknn <- impute.knn(Bok, k = 10, maxp = 1500, rng.seed = 42)$data))   # k = 10, fixed seed
cat(sprintf("KNN (k=10) sobre %d x %d en %.1f min | perdidos restantes: %d\n",
            nrow(Bknn), ncol(Bknn), as.numeric(difftime(Sys.time(), t0, units = "mins")), sum(is.na(Bknn))))
stopifnot(sum(is.na(Bknn)) == 0)
M <- pmin(pmax(Bknn, 1e-4), 1 - 1e-4); M <- log2(M / (1 - M))   # M-values for inference
rm(Bok); gc()
fitP <- eBayes(lmFit(M, dis)); ttP <- topTable(fitP, coef = 2, number = Inf, sort.by = "none")   # primary model, for comparison
lambdaP <- median(qchisq(ttP$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
cat(sprintf("Primario reproducido    : λ = %.3f | CpGs a FDR<0.05: %d  (debe dar 1.033 y 22)\n",
            lambdaP, sum(ttP$adj.P.Val < 0.05)))
rm(fitP); gc()

# ---- 6. S8: methylation-derived smoking score and AHRR as covariates --------------------------------
# EpiSmokEr on the full raw matrix, with missing values replaced by the probe median for this calculation
# only; sample sheet with sex coded 1 = male, 2 = female; the continuous score is used. The cohort did not
# record smoking, so two coherence checks are printed: score by generation and correlation with the beta
# of cg05575921 (AHRR). Model A: design + score; model B: design + M-value of cg05575921.
if (!requireNamespace("EpiSmokEr", quietly = TRUE)) {        # installation from GitHub
  cat("EpiSmokEr no está instalado; intentando remotes::install_github('sailalithabollepalli/EpiSmokEr')\n")
  tryCatch(remotes::install_github("sailalithabollepalli/EpiSmokEr", upgrade = "never", quiet = TRUE),
           error = function(e) cat("FALLO la instalación:", conditionMessage(e),
                                   "\nInstalar a mano en R:  remotes::install_github('sailalithabollepalli/EpiSmokEr')\n"))
}
stopifnot(requireNamespace("EpiSmokEr", quietly = TRUE))
suppressMessages(library(EpiSmokEr))
# EpiSmokEr exports a function named `all` that masks base::all; base::all is used explicitly from here on
all <- base::all
Bfill <- Braw_ok                                             # full raw matrix, only for the score
na_f  <- which(is.na(Bfill), arr.ind = TRUE)
if (nrow(na_f) > 0) {                                        # probe-median fill
  med <- apply(Bfill[unique(na_f[, 1]), , drop = FALSE], 1, median, na.rm = TRUE)
  Bfill[na_f] <- med[match(rownames(Bfill)[na_f[, 1]], names(med))]
}
Bfill <- Bfill[complete.cases(Bfill), , drop = FALSE]        # in case a probe was entirely missing
cat(sprintf("S8 · matriz para EpiSmokEr: %d sondas x %d muestras, NA restantes: %d\n",
            nrow(Bfill), ncol(Bfill), sum(is.na(Bfill))))
ss <- data.frame(row.names = colnames(Bfill),                # minimal sample sheet
                 sex = ifelse(grepl("^M", d$sexo, ignore.case = TRUE), 1L, 2L))   # 1 = male, 2 = female
stopifnot(base::all(table(ss$sex) > 0))
out <- epismoker(dataset = Bfill, samplesheet = ss, method = "all")
out <- as.data.frame(out)
cat("Columnas devueltas por EpiSmokEr:", paste(names(out), collapse = ", "), "\n")
col_sc <- grep("score", names(out), ignore.case = TRUE, value = TRUE)   # column holding the score
col_sc <- if (any(grepl("smokingScore", col_sc))) "smokingScore" else col_sc[1]
stopifnot(length(col_sc) == 1)
id_col <- if ("SampleName" %in% names(out)) out$SampleName else rownames(out)
smk <- out[[col_sc]][match(colnames(Bfill), id_col)]         # aligned with the column order
stopifnot(!anyNA(smk))
cat(sprintf("Puntaje continuo usado: '%s' | media %.2f, DE %.2f, rango %.2f a %.2f\n",
            col_sc, mean(smk), sd(smk), min(smk), max(smk)))
rm(Bfill); gc()

# Coherence checks, in place of the comparison with self-reported smoking the cohort does not allow
ahrr_beta <- if (AHRR %in% rownames(Braw_ok)) Braw_ok[AHRR, ] else rep(NA_real_, ncol(Braw_ok))
gen_ok <- ph$gen[ok]
cat("\n=== S8 · coherencia del puntaje ===\n")
for (g in sort(unique(gen_ok))) {                            # score by generation (higher in mothers expected)
  v <- smk[gen_ok == g]
  cat(sprintf("  %-25s n = %2d | media %.2f\n", g, length(v), mean(v)))
}
cat(sprintf("Correlación (Spearman) puntaje ~ edad: %.2f\n", cor(smk, d$edad, method = "spearman")))
if (!anyNA(ahrr_beta))
  cat(sprintf("Correlación (Spearman) puntaje ~ β de cg05575921 (AHRR): %.2f  <- se espera negativa\n",
              cor(smk, ahrr_beta, method = "spearman")))
fwrite(data.table(col = colnames(Braw_ok), generacion = gen_ok, edad = d$edad, sexo = as.character(d$sexo),
                  puntaje_tabaco = round(smk, 4), beta_AHRR_cg05575921 = round(ahrr_beta, 4)),
       file.path(OUT, "S8_puntaje_tabaco_por_muestra.csv"))
rm(Braw_ok); gc()
stopifnot(!anyNA(ahrr_beta))
mA <- pmin(pmax(ahrr_beta, 1e-4), 1 - 1e-4); mA <- log2(mA / (1 - mA))   # AHRR on the M-value scale

# Reference rows (primary model), then models A and B in a loop
res8 <- data.table(analisis = "primario (referencia)", n_sondas = nrow(M), lambda = round(lambdaP, 3),
                   n_fdr05 = sum(ttP$adj.P.Val < 0.05), primarios_recuperados = length(primarios),
                   zg16b_fdr05 = sum(LOCUS$cpg %in% rownames(ttP)[ttP$adj.P.Val < 0.05]))
loc8 <- merge(LOCUS, data.table(cpg = rownames(ttP), M_est = ttP$logFC, p = ttP$P.Value, fdr = ttP$adj.P.Val), by = "cpg")
loc8[, `:=`(dir_ok = sign(M_est) == dir_pred, analisis = "primario (referencia)")]; loc8 <- loc8[order(region, cpg)]
modelos <- list("manuscrito + puntaje EpiSmokEr" = cbind(dis, tabaco = as.numeric(scale(smk))),
                "manuscrito + AHRR cg05575921"   = cbind(dis, ahrr = as.numeric(scale(mA))))
for (etiqueta in names(modelos)) {
  disX <- modelos[[etiqueta]]
  stopifnot(qr(disX)$rank == ncol(disX))
  fitX  <- eBayes(lmFit(M, disX));    ttX <- topTable(fitX,  coef = 2, number = Inf, sort.by = "none")
  fitXb <- eBayes(lmFit(Bknn, disX)); tbX <- topTable(fitXb, coef = 2, number = Inf, sort.by = "none")   # delta-beta
  hitsX <- rownames(ttX)[ttX$adj.P.Val < 0.05]
  res8 <- rbind(res8, data.table(analisis = etiqueta, n_sondas = nrow(M),
                                 lambda = round(median(qchisq(ttX$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1), 3),
                                 n_fdr05 = length(hitsX), primarios_recuperados = sum(primarios %in% hitsX),
                                 zg16b_fdr05 = sum(LOCUS$cpg %in% hitsX)))
  locX <- merge(LOCUS, data.table(cpg = rownames(ttX), M_est = ttX$logFC, p = ttX$P.Value, fdr = ttX$adj.P.Val), by = "cpg")
  locX <- merge(locX, data.table(cpg = rownames(tbX), dBeta = tbX$logFC), by = "cpg")
  locX[, `:=`(dir_ok = sign(M_est) == dir_pred, analisis = etiqueta)]
  loc8 <- rbind(loc8, locX[order(region, cpg)], fill = TRUE)
  rm(fitX, fitXb); gc()
}
cat("\n=== S8 · resumen genómico ===\n"); print(res8)
cat("\n=== S8 · el locus de ZG16B bajo los tres modelos ===\n")
print(loc8[, .(analisis, cpg, region, dir_pred, M_est = round(M_est, 4), p = signif(p, 3), fdr = signif(fdr, 3), dir_ok)])
fwrite(res8, file.path(OUT, "S8_tabaco_resumen.csv"))
fwrite(loc8, file.path(OUT, "S8_tabaco_locus.csv"))

# ---- 7. S9: one participant per extended family -----------------------------------------------------
# Rule fixed before running: in each extended family, the participant whose age is closest to the median
# age of the replication cohort (ties toward the younger). The design of the article is kept if it leaves at
# least 25 residual degrees of freedom; otherwise it is simplified to linear terms, and the choice is printed.
ph2 <- fenotipos[["GSE239521"]]
ed2 <- suppressWarnings(as.numeric(ph2[[grep("^age", names(ph2), value = TRUE)[1]]]))
objetivo <- median(ed2, na.rm = TRUE)                        # about 14.1 years
cat(sprintf("S9 · mediana de edad de la cohorte de réplica: %.2f años (n = %d)\n", objetivo, sum(!is.na(ed2))))
cand <- data.table(idx = which(ok), familia = ph$familia[ok], edad = d$edad, dist = abs(d$edad - objetivo))
setorder(cand, familia, dist, edad)                          # per family: closest first, then younger
sel  <- cand[, .SD[1], by = familia]                         # one per family
sel9 <- match(sel$idx, which(ok))                            # positions within the analytic sample
d9   <- d[sel9, ]
cat(sprintf("S9 · n = %d participantes, uno por familia | edad: mediana %.1f, rango %.0f a %.0f | en 10-19 años: %d\n",
            nrow(d9), median(d9$edad), min(d9$edad), max(d9$edad), sum(d9$edad >= 10 & d9$edad <= 19)))
print(table(generacion = ph$gen[ok][sel9]))
d9$sexo <- droplevels(d9$sexo); d9$lote <- droplevels(d9$lote)
dis9 <- model.matrix(~ y + ns(edad, 3) + sexo + lote + ns(Epi, 4), data = d9)   # design of the article
if (nrow(dis9) - ncol(dis9) < 25 || qr(dis9)$rank < ncol(dis9)) {   # too few residual degrees of freedom
  dis9 <- model.matrix(~ y + edad + sexo + lote + Epi, data = d9)   # pre-specified linear fallback
  cat("S9 · diseño simplificado a términos lineales por grados de libertad\n")
} else cat("S9 · se conserva el diseño del manuscrito\n")
cat(sprintf("S9 · parámetros %d, grados de libertad residuales %d\n", ncol(dis9), nrow(dis9) - ncol(dis9)))
stopifnot(qr(dis9)$rank == ncol(dis9))
fit9  <- eBayes(lmFit(M[, sel9], dis9));    tt9 <- topTable(fit9,  coef = 2, number = Inf, sort.by = "none")
fit9b <- eBayes(lmFit(Bknn[, sel9], dis9)); tb9 <- topTable(fit9b, coef = 2, number = Inf, sort.by = "none")
hits9 <- rownames(tt9)[tt9$adj.P.Val < 0.05]
res9 <- data.table(analisis = "una por familia", n = nrow(d9), edad_mediana = median(d9$edad),
                   edad_min = min(d9$edad), edad_max = max(d9$edad),
                   lambda = round(median(qchisq(tt9$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1), 3),
                   n_fdr05 = length(hits9), primarios_recuperados = sum(primarios %in% hits9),
                   zg16b_fdr05 = sum(LOCUS$cpg %in% hits9),
                   zg16b_dir_ok = sum(sign(tt9[LOCUS$cpg, "logFC"]) == LOCUS$dir_pred),
                   zg16b_p_nominal = sum(tt9[LOCUS$cpg, "P.Value"] < 0.05))
cat("\n=== S9 · resumen ===\n"); print(res9)
loc9 <- merge(LOCUS, data.table(cpg = rownames(tt9), M_est = tt9$logFC, p = tt9$P.Value, fdr = tt9$adj.P.Val), by = "cpg")
loc9 <- merge(loc9, data.table(cpg = rownames(tb9), dBeta = tb9$logFC), by = "cpg")
loc9[, `:=`(dir_ok = sign(M_est) == dir_pred, analisis = "una por familia")]; loc9 <- loc9[order(region, cpg)]
cat("\n=== S9 · el locus de ZG16B con un participante por familia ===\n")
print(loc9[, .(cpg, region, dir_pred, M_est = round(M_est, 4), dBeta = round(dBeta, 4), p = signif(p, 3), fdr = signif(fdr, 3), dir_ok)])
fwrite(res9, file.path(OUT, "S9_una_por_familia_resumen.csv"))
fwrite(loc9, file.path(OUT, "S9_una_por_familia_locus.csv"))
fwrite(data.table(familia = sel$familia, col = ph$col[ok][sel9], edad = d9$edad, generacion = ph$gen[ok][sel9],
                  violencia = d9$y), file.path(OUT, "S9_una_por_familia_seleccion.csv"))

# ---- 8. Session information -------------------------------------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_01d.txt"))
cat("Salidas escritas en", normalizePath(OUT), ":\n")
print(list.files(OUT, pattern = "^S[789]_|sessionInfo_01d"))
