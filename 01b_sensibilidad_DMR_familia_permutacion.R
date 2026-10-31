# 01b — Primary analysis on the imputed matrix, regions and intersection, within-family model, descriptives
#
# Reproduces the loading of 01 and adds: the k-nearest-neighbor imputation of the isolated missing values
# and the primary genome-wide model on the imputed matrix; (A) delta-beta per unit of exposure (the same
# model fitted on beta); (B) a second model with the within-family correlation estimated by
# duplicateCorrelation; (C) differentially methylated regions with DMRcate and the DMP-DMR intersection;
# (D) family-level permutations, written out below but not run (see the note there); (E) descriptives,
# group comparability, the variance partition of the violence score, the inventory of ZG16B probes and
# the sensitivity analysis that reinstates the ten participants with technical replicates.
# Paper section: Methods 2.3-2.6; Results 3.1, 3.2, 3.4 and 3.7; Table 2; Supplementary Tables S1, S2, S6
# Input:  GSE226085 processed matrix and SOFT file (BU); GSE239520 and GSE239521 SOFT files (CM);
#         outputs/tronco_ZG16B.RDS (01) and outputs/cohorte_discordante.RDS (02)
# Output: outputs/T0_comparabilidad_grupos.csv, T3_efectos_cohorte1.csv, T4_ZG16B_correlacion_familiar.csv,
#         T5_DMR_cohorte1.csv, T6_interseccion_DMP_DMR.csv, S2_inventario_sondas_ZG16B.csv,
#         S3_sensibilidad_replicas_tecnicas.csv,
#         descriptivos_01b.RDS, sensibilidad_01b.RDS, sessionInfo_01b.txt

# ---- 1. Packages, paths and locus -------------------------------------------------------------------
suppressMessages({
  library(data.table); library(limma); library(EpiDISH); library(splines)
  library(DMRcate)                                           # region-level analysis
  library(statmod)                                           # required by duplicateCorrelation
  library(impute)                                            # k-nearest-neighbor imputation
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})
set.seed(42)                                                 # same seed as 01
BU  <- "data/GSE226085"                                      # GSE226085 (edit to the local path)
CM  <- "data/evaluation"                                     # SOFT files of the evaluation cohorts (edit to the local path)
OUT <- "outputs"
dir.create(OUT, showWarnings = FALSE)
N_PERM <- 30                                                 # permutations, only if section D is reactivated
LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",   # TSS200 -> negative
          "cg19128189","cg08829024","cg12475088","cg04543077"),  # TSS1500 -> positive
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))

# ---- 2. Phenotypes from the GEO SOFT files of the three cohorts -------------------------------------
fenotipos <- list()
for (g in c("GSE226085", "GSE239521", "GSE239520")) {
  fichero <- if (g == "GSE226085") file.path(BU, "GSE226085.txt") else file.path(CM, paste0(g, ".txt"))
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

# ---- 3. Discovery cohort, loaded as in 01 -----------------------------------------------------------
ph <- fenotipos[["GSE226085"]]
ph[, `:=`(violencia = suppressWarnings(as.numeric(violence.trauma)), edad = suppressWarnings(as.numeric(age)),
          sexo = gender, familia = extended.family,
          lote = plate, col = individual)]
ph[, usar := grepl("^Y|^y|TRUE|1", passed.methylation.quality.control) &   # authors' quality control
             grepl("^Y|^y|TRUE|1", included.in.final.analytical.sample)]   # and their final analytic sample
cab <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"), nrows = 0)   # header only
B <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"),
           select = setdiff(names(cab), grep("Detection", names(cab), value = TRUE)),   # detection p columns not loaded
           showProgress = FALSE)
cg <- B[[1]]                                                 # first column: probe identifiers
B  <- as.matrix(B[, -1])
rownames(B) <- cg
n_inicial <- nrow(B)                                         # probe count before filtering

# Kept aside for the technical-replicate sensitivity analysis (section E): the matrix has 160 columns,
# 140 named by individual and 20 with suffix _1/_2 (two profiles of ten participants), which the name
# matching below leaves out. Only the eight locus probes and the EpiDISH reference probes are needed.
data(centEpiFibIC.m)                                         # EpiDISH epithelial/fibroblast/immune reference
B_zg  <- B[intersect(LOCUS$cpg, rownames(B)), , drop = FALSE]
B_ref <- B[intersect(rownames(centEpiFibIC.m), rownames(B)), , drop = FALSE]
ph_todo <- copy(ph)                                          # full phenotype, replicates included
cat(sprintf("Apartado para la sensibilidad: %d sondas del locus y %d de referencia, x %d columnas\n",
            nrow(B_zg), nrow(B_ref), ncol(B_zg)))

com <- intersect(colnames(B), ph$col)
ph  <- ph[match(com, col)]                                   # phenotype aligned with the matrix
B   <- B[, com, drop = FALSE]
stopifnot(identical(colnames(B), ph$col))
B <- B[rowMeans(is.na(B)) < 0.05 &                           # 5% or more missing values out
       !grepl("^rs|^ch\\.", rownames(B)), , drop = FALSE]    # SNP-targeting and non-CpG probes out
ph[, Epi := epidish(beta.m = B, ref.m = centEpiFibIC.m,      # cell-type deconvolution
                    method = "RPC")$estF[, "Epi"]]           # epithelial fraction
ok <- ph$usar & !is.na(ph$violencia) & !is.na(ph$edad) & !is.na(ph$Epi)   # analytic sample
ok[is.na(ok)] <- FALSE
cat(sprintf("Encadenado de sondas: %d -> %d tras filtrar perdidos, SNP y no-CpG\n",
            n_inicial, nrow(B)))                             # 865918 -> 736065
cat(sprintf("Muestra analítica: n = %d de %d | %d familias extensas\n",
            sum(ok), ncol(B), uniqueN(ph$familia[ok])))      # 116 of 140, 42 families
stopifnot(nrow(B) == 736065, sum(ok) == 116)                 # figures reported in the article

# ---- 4. Missing values and primary model ------------------------------------------------------------
# Isolated missing values remain after the filter. They are imputed with KNN (impute.knn, k = 10, fixed
# seed); the earlier convention (missing -> beta 0.5) is computed alongside to document how much the count
# depends on it. The eight ZG16B probes have no missing value.
d <- data.frame(y    = ph$violencia[ok],                     # continuous exposure, 0-19
                edad = ph$edad[ok],
                sexo = droplevels(factor(ph$sexo[ok])),
                Epi  = ph$Epi[ok],
                lote = droplevels(factor(ph$lote[ok])))      # plate
dis <- model.matrix(~ y + ns(edad, 3) + sexo + lote + ns(Epi, 4), data = d)   # design of the article
stopifnot(qr(dis)$rank == ncol(dis))
Bok <- B[, ok, drop = FALSE]                                 # matrix of the analytic sample
rm(B); gc()

n_na  <- sum(is.na(Bok))
n_son <- sum(rowSums(is.na(Bok)) > 0)
zg_na <- sum(rowSums(is.na(Bok[intersect(LOCUS$cpg, rownames(Bok)), , drop = FALSE])) > 0)
cat(sprintf("Perdidos tras el filtro: %d valores en %d sondas (%.1f%% de las %d)\n",
            n_na, n_son, 100 * n_son / nrow(Bok), nrow(Bok)))
cat(sprintf("  Sondas de ZG16B con algún perdido: %d de 8   <- si es 0, el locus no depende\n",
            zg_na))
stopifnot(zg_na == 0)                                        # the locus must not depend on the imputation

# Earlier convention first (missing -> M = 0, beta = 0.5), freed immediately
Mc <- pmin(pmax(Bok, 1e-4), 1 - 1e-4); Mc <- log2(Mc / (1 - Mc)); Mc[!is.finite(Mc)] <- 0
ttC <- topTable(eBayes(lmFit(Mc, dis)), coef = 2, number = Inf, sort.by = "none")
cat(sprintf("\nConvención del tronco   : λ = %.3f | CpGs a FDR<0.05: %d\n",
            median(qchisq(ttC$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1),
            sum(ttC$adj.P.Val < 0.05)))
rm(Mc); gc()

# KNN imputation, the analysis reported
t0 <- Sys.time()
Bknn <- impute.knn(Bok, k = 10, maxp = 1500, rng.seed = 42)$data   # k = 10 neighbors, fixed seed
rm(Bok); gc()
cat(sprintf("KNN (k=10) sobre %d x %d en %.1f min | perdidos restantes: %d\n",
            nrow(Bknn), ncol(Bknn), as.numeric(difftime(Sys.time(), t0, units = "mins")),
            sum(is.na(Bknn))))
stopifnot(sum(is.na(Bknn)) == 0)

# Primary model on the imputed matrix (M-values), and the same model on beta for delta-beta
M <- pmin(pmax(Bknn, 1e-4), 1 - 1e-4); M <- log2(M / (1 - M))
fitM <- eBayes(lmFit(M, dis))
fitB <- eBayes(lmFit(Bknn, dis))
ttM  <- topTable(fitM, coef = 2, number = Inf, sort.by = "none")
lambdaM <- median(qchisq(ttM$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)   # genomic inflation factor
cat(sprintf("PRIMARIO (KNN)          : λ = %.3f | CpGs a FDR<0.05: %d\n",
            lambdaM, sum(ttM$adj.P.Val < 0.05)))
zg_cmp <- merge(LOCUS,                                       # the locus under both missing-value conventions
                data.table(cpg = rownames(ttM),
                           knn_beta = round(ttM$logFC, 4),
                           knn_fdr  = signif(ttM$adj.P.Val, 3),
                           tronco_beta = round(ttC$logFC, 4),
                           tronco_fdr  = signif(ttC$adj.P.Val, 3)), by = "cpg")
cat("\n=== El locus de ZG16B bajo las dos convenciones de perdidos ===\n")
print(zg_cmp[order(region, cpg)])

# ---- 5. A: delta-beta per unit of exposure (same design, fitted on beta) ----------------------------
ttB <- topTable(fitB, coef = 2, number = Inf, sort.by = "none")
stopifnot(identical(rownames(ttB), rownames(ttM)))           # same probe order
rm(fitB, Bknn); gc()                                         # DMRcate needs the memory
an <- as.data.table(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19),   # hg19 annotation
                    keep.rownames = "cpg")
efectos <- data.table(cpg   = rownames(ttM),
                      M_est = ttM$logFC,                     # coefficient on the M-value scale
                      dBeta = ttB$logFC,                     # delta-beta per unit of violence
                      p     = ttM$P.Value,
                      fdr   = ttM$adj.P.Val)                 # Benjamini-Hochberg FDR
efectos <- merge(efectos, an[, .(cpg, chr, pos,
                                 gene  = UCSC_RefGene_Name,
                                 grupo = UCSC_RefGene_Group,
                                 isla  = Relation_to_Island)], by = "cpg")
for (i in seq_len(nrow(efectos))) {                          # the manifest repeats gene and region once per transcript
  set(efectos, i, "gene",  paste(unique(strsplit(efectos$gene[i],  ";")[[1]]), collapse = "/"))
  set(efectos, i, "grupo", paste(unique(strsplit(efectos$grupo[i], ";")[[1]]), collapse = "/"))
}
T3 <- efectos[fdr < 0.05][order(p)]                          # the significant probes
cat("\n=== TABLA 3 · los 22 CpGs con su Δβ interpretable ===\n")
print(T3[, .(cpg, gene, grupo, chr, pos,
             M = round(M_est, 3), dBeta = round(dBeta, 4),
             p = signif(p, 3), fdr = signif(fdr, 3))])
fwrite(T3, file.path(OUT, "T3_efectos_cohorte1.csv"))
zg <- merge(LOCUS, efectos, by = "cpg")[order(pos)]          # the ZG16B locus
cat("\n=== Locus ZG16B · las 8 sondas con su Δβ ===\n")
print(zg[, .(cpg, region, pos, dir_pred,
             M = round(M_est, 3), dBeta = round(dBeta, 4),
             fdr = signif(fdr, 3),
             dir_ok = sign(M_est) == dir_pred)])

# Concordance test of the evaluation cohorts (one-sided exact binomial on the sign), from the saved artifacts
for (f in c("tronco_ZG16B.RDS", "cohorte_discordante.RDS")) {
  ruta <- file.path(OUT, f)
  if (!file.exists(ruta)) { cat("· falta", f, "\n"); next }
  r  <- readRDS(ruta)
  if (!is.null(r$cohorte2)) r <- r$cohorte2                  # 01 nests the replication cohort in $cohorte2
  k  <- if (!is.null(r$n_dir)) r$n_dir else sum(r$locus$dir_ok)   # concordant probes
  n  <- if (!is.null(r$locus)) nrow(r$locus) else 8L
  bt <- binom.test(k, n, 0.5, alternative = "greater")
  cat(sprintf("%-28s %d/%d en la dirección predicha | binomial p = %.4f\n",
              f, k, n, bt$p.value))
}

# ---- 6. B: within-family correlation (limma duplicateCorrelation) -----------------------------------
# Second model, identical to the primary one except for the correlation within extended families in the
# variance structure. Both FDRs are reported for the locus.
sub <- sample(nrow(M), min(50000, nrow(M)))                  # random subsample to estimate the correlation
dc  <- duplicateCorrelation(M[sub, , drop = FALSE], dis,
                            block = ph$familia[ok])          # the extended family is the block
rho <- dc$consensus.correlation
cat(sprintf("\nCorrelación intra-familia (consenso sobre %d sondas): ρ = %.4f\n",
            length(sub), rho))
fitF <- eBayes(lmFit(M, dis, block = ph$familia[ok], correlation = rho))   # refit with the correlation
ttF  <- topTable(fitF, coef = 2, number = Inf, sort.by = "none")
cat(sprintf("Modelo con correlación familiar: λ = %.3f | CpGs a FDR<0.05: %d (primario: %d)\n",
            median(qchisq(ttF$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1),
            sum(ttF$adj.P.Val < 0.05), sum(ttM$adj.P.Val < 0.05)))
comp <- merge(zg[, .(cpg, region, M_primario = M_est, dBeta, fdr_primario = fdr)],
              data.table(cpg = rownames(ttF), M_familia = ttF$logFC,
                         fdr_familia = ttF$adj.P.Val), by = "cpg")
setorder(comp, region, cpg)
cat("\n=== Locus ZG16B · modelo primario frente a modelo con correlación familiar ===\n")
print(comp[, .(cpg, region, dBeta = round(dBeta, 4),
               primario = sprintf("%.3f (FDR %.3g)", M_primario, fdr_primario),
               familiar = sprintf("%.3f (FDR %.3g)", M_familia, fdr_familia))])
fwrite(comp, file.path(OUT, "T4_ZG16B_correlacion_familiar.csv"))
cat(sprintf("\n  De las 8 sondas de ZG16B, mantienen FDR<0.05: %d en el primario, %d con familia\n",
            sum(comp$fdr_primario < 0.05), sum(comp$fdr_familia < 0.05)))

# ---- 7. C: regions (DMRcate) and DMP-DMR intersection -----------------------------------------------
# Same design matrix and coefficient as the position-level analysis; only a gene present in both
# analyses is carried forward. The arraytype name changes across DMRcate versions, so it is read from
# the installed function.
validos <- eval(formals(DMRcate::cpg.annotate)$arraytype)
ARR <- intersect(c("EPICv1", "EPIC"), validos)[1]
stopifnot(!is.na(ARR))
cat(sprintf("DMRcate %s admite arraytype = {%s} -> se usa \"%s\"\n",
            as.character(packageVersion("DMRcate")),
            paste(validos, collapse = ", "), ARR))
anot <- cpg.annotate(object = M, datatype = "array", what = "M",
                     analysis.type = "differential", design = dis,
                     coef = 2, arraytype = ARR, fdr = 0.05)  # seeding FDR 0.05
cat(sprintf("\nSondas que superan el semillado (FDR<0.05): %d\n",
            sum(anot@ranges$is.sig)))
dmr <- dmrcate(anot, lambda = 1000, C = 2)                   # Gaussian kernel bandwidth and scaling factor
reg <- if (length(dmr@coord) > 0) as.data.table(extractRanges(dmr, genome = "hg19")) else
       data.table()
cat(sprintf("Regiones diferencialmente metiladas: %d\n", nrow(reg)))
if (nrow(reg) > 0) {
  reg <- reg[order(min_smoothed_fdr)]
  print(reg[1:min(10, nrow(reg)),
            .(seqnames, start, end, width, no.cpgs,
              overlapping.genes, meandiff = round(meandiff, 4),
              fdr = signif(min_smoothed_fdr, 3))])
  fwrite(reg, file.path(OUT, "T5_DMR_cohorte1.csv"))
  gr_dmp <- an[cpg %in% T3$cpg, .(cpg, chr, pos)]            # significant DMPs with their position
  gr_dmp[, en_dmr := FALSE]
  for (i in seq_len(nrow(gr_dmp)))                           # does each DMP fall inside a region?
    set(gr_dmp, i, "en_dmr", any(reg$seqnames == gr_dmp$chr[i] & reg$start <= gr_dmp$pos[i] & reg$end >= gr_dmp$pos[i]))
  inter <- merge(gr_dmp[en_dmr == TRUE], T3[, .(cpg, gene, grupo, fdr)], by = "cpg")
  cat("\n=== INTERSECCIÓN DMP-DMR · los CpGs significativos que además caen en una región ===\n")
  print(inter[order(fdr)])
  fwrite(inter, file.path(OUT, "T6_interseccion_DMP_DMR.csv"))
  cat(sprintf("\n  %d de los %d DMP caen dentro de un DMR | ¿ZG16B entre ellos? %s\n",
              nrow(inter), nrow(T3),
              ifelse(any(grepl("ZG16B", inter$gene)), "SÍ", "NO")))
} else {                                                     # null outcome, reported as such
  cat("  DMRcate no formó ninguna región. Se reporta así: la capa de regiones no\n",
      "  respalda la de posiciones, y la intersección queda vacía.\n")
}

# ---- 8. D: family-level permutations (not run) ------------------------------------------------------
# Kept as comments so that the procedure stays visible; it is reactivated by removing the leading '# '.
# Permutation is not the convention in a probe-level EWAS (the convention is lambda + FDR), and the
# non-independence of families is handled by the within-family model of section B.
# fam <- as.character(ph$familia[ok])                        # family of each sample
# vio_fam <- tapply(d$y, fam, function(z) z[1])              # exposure per family
# stopifnot(!any(is.na(vio_fam)))
# nulo <- data.table(perm = seq_len(N_PERM), lambda = NA_real_, n_fdr = NA_integer_)
# for (i in seq_len(N_PERM)) {
#   barajado <- setNames(sample(as.vector(vio_fam)), names(vio_fam))   # shuffle the UNNAMED vector, relabel after
#   y_perm   <- barajado[fam]
#   dp  <- d; dp$y <- as.numeric(y_perm)
#   dsp <- model.matrix(~ y + ns(edad, 3) + sexo + lote + ns(Epi, 4), data = dp)
#   ttp <- topTable(eBayes(lmFit(M, dsp)), coef = 2, number = Inf, sort.by = "none")
#   nulo[i, `:=`(lambda = median(qchisq(ttp$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1),
#                n_fdr = sum(ttp$adj.P.Val < 0.05))]
#   cat(sprintf("  perm %2d/%d · λ = %.3f · CpGs FDR<0.05 = %d\n", i, N_PERM, nulo$lambda[i], nulo$n_fdr[i]))
# }
# stopifnot(length(unique(nulo$lambda)) > 1)                 # a null that does not vary means a broken permutation
# p_emp <- (1 + sum(nulo$n_fdr >= sum(ttM$adj.P.Val < 0.05))) / (1 + N_PERM)   # empirical p
# fwrite(nulo, file.path(OUT, "S2_permutaciones_familia.csv"))

# ---- 9. E: group comparability in the evaluation cohorts --------------------------------------------
filas <- list()
for (g in c("GSE239521", "GSE239520")) {
  p2 <- copy(fenotipos[[g]])
  p2[, malt := ifelse(grepl("maltreat", group, ignore.case = TRUE), 1L,   # 1 = maltreatment
              ifelse(grepl("typically", group, ignore.case = TRUE), 0L, NA_integer_))]
  p2[, ed := suppressWarnings(as.numeric(age))]
  p2 <- p2[!is.na(malt)]
  tt <- t.test(ed ~ malt, data = p2)                         # Welch t test for age
  ch <- suppressWarnings(chisq.test(table(p2$gender, p2$malt)))   # chi-squared test for sex
  filas[[g]] <- data.table(cohorte = g, n = nrow(p2),
             edad_expuestos = sprintf("%.2f ± %.2f", mean(p2[malt==1]$ed), sd(p2[malt==1]$ed)),
             edad_control   = sprintf("%.2f ± %.2f", mean(p2[malt==0]$ed), sd(p2[malt==0]$ed)),
             p_edad = signif(tt$p.value, 3), p_sexo = signif(ch$p.value, 3))
}
comparabilidad <- rbindlist(filas)
cat("\n=== Comparabilidad de los grupos (cohortes 2 y 3) ===\n"); print(comparabilidad)
fwrite(comparabilidad, file.path(OUT, "T0_comparabilidad_grupos.csv"))

# Observed against the expected background: promoter share among the analyzed and the significant probes
an_bg <- an[cpg %in% rownames(ttM)]                          # the background actually analyzed
pf <- 100 * mean(grepl("TSS200|TSS1500|5'UTR|1stExon", an_bg$UCSC_RefGene_Group))
ps <- 100 * mean(grepl("TSS200|TSS1500|5'UTR|1stExon", an[cpg %in% T3$cpg]$UCSC_RefGene_Group))
cat(sprintf("\nSondas en región promotora: fondo %.2f%% (n=%d) | significativas %.2f%% (n=%d)\n",
            pf, nrow(an_bg), ps, nrow(T3)))
n_bonf <- sum(ttM$P.Value < 0.05 / nrow(ttM))                # probes passing Bonferroni
cat(sprintf("CpGs que superan Bonferroni (%.2e): %d\n", 0.05/nrow(ttM), n_bonf))
saveRDS(list(comparabilidad = comparabilidad, prom_fondo = pf, prom_sig = ps,
             n_bonf = n_bonf, ttM_p = ttM$P.Value, ttM_cpg = rownames(ttM)),
        file.path(OUT, "descriptivos_01b.RDS"))

# ---- 10. E: sensitivity with the ten participants with technical replicates reinstated --------------
# The analytic sample excludes ten participants profiled twice. The locus model (eight probes, not the
# genome) is refitted keeping the first profile of each, to show that the exclusion did not favor the result.
ind    <- sub("_[12]$", "", colnames(B_zg))                  # suffix removed
primer <- !duplicated(ind)                                   # one profile per individual
ph_ext <- ph_todo[match(ind[primer], col)]
Bz     <- B_zg[, primer, drop = FALSE]
ph_ext[, Epi := epidish(beta.m = B_ref[, primer, drop = FALSE],
                        ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
ok_e <- ph_ext$usar & !is.na(ph_ext$violencia) & !is.na(ph_ext$edad) & !is.na(ph_ext$Epi)
ok_e[is.na(ok_e)] <- FALSE
cat(sprintf("\nSensibilidad de réplicas: n = %d (frente a %d del análisis principal)\n",
            sum(ok_e), sum(ok)))
d_e <- data.frame(y = ph_ext$violencia[ok_e], edad = ph_ext$edad[ok_e],
                  sexo = droplevels(factor(ph_ext$sexo[ok_e])),
                  Epi = ph_ext$Epi[ok_e],
                  lote = droplevels(factor(ph_ext$lote[ok_e])))
dis_e <- model.matrix(~ y + ns(edad, 3) + sexo + lote + ns(Epi, 4), data = d_e)   # same design
stopifnot(qr(dis_e)$rank == ncol(dis_e))
Mz <- pmin(pmax(Bz[, ok_e, drop = FALSE], 1e-4), 1 - 1e-4); Mz <- log2(Mz / (1 - Mz))
tt_e <- topTable(eBayes(lmFit(Mz, dis_e)), coef = 2, number = Inf, sort.by = "none")
sens_rep <- merge(LOCUS, data.table(cpg = rownames(tt_e),
                                    M_ampliado = round(tt_e$logFC, 4),
                                    p_ampliado = signif(tt_e$P.Value, 3)), by = "cpg")
sens_rep <- merge(sens_rep, zg[, .(cpg, M_principal = M_est, p_principal = p)], by = "cpg")
sens_rep[, dir_ok := sign(M_ampliado) == dir_pred]
setorder(sens_rep, region, cpg)
cat("=== El locus con y sin los participantes replicados ===\n")
print(sens_rep[, .(cpg, region, M_principal = round(M_principal, 4), M_ampliado,
                   p_principal = signif(p_principal, 3), p_ampliado, dir_ok)])
cat(sprintf("  Signo predicho conservado en %d de 8 sondas | correlación de los coeficientes r = %.4f\n",
            sum(sens_rep$dir_ok), cor(sens_rep$M_principal, sens_rep$M_ampliado)))
fwrite(sens_rep, file.path(OUT, "S3_sensibilidad_replicas_tecnicas.csv"))
rm(B_zg, B_ref); gc()

# ---- 11. E: descriptives of the analytic sample and variance partition of the exposure --------------
ph_ok <- ph[ok]
cat(sprintf("\nComposición: %d mujeres / %d hombres | edad M=%.2f DE=%.2f rango %.1f-%.1f\n",
            sum(ph_ok$sexo == "Female"), sum(ph_ok$sexo == "Male"),
            mean(ph_ok$edad), sd(ph_ok$edad), min(ph_ok$edad), max(ph_ok$edad)))
cat("  generaciones: "); print(table(ph_ok$generation))
cat(sprintf("  violencia: M=%.2f DE=%.2f rango %.0f-%.0f\n",
            mean(ph_ok$violencia), sd(ph_ok$violencia),
            min(ph_ok$violencia), max(ph_ok$violencia)))
# Is exposure an individual or a family characteristic? Sums of squares between and within families
gm  <- mean(ph_ok$violencia)
agg <- ph_ok[, .(m = mean(violencia), k = .N, ss = sum((violencia - mean(violencia))^2)),
             by = familia]
ss_b <- sum(agg$k * (agg$m - gm)^2)                          # between families
ss_w <- sum(agg$ss)                                          # within families
idn  <- ph_ok[, .(u = uniqueN(violencia), k = .N), by = familia]   # identical score within the family?
cat(sprintf("  varianza de la violencia: %.2f%% entre familias, %.2f%% dentro\n",
            100*ss_b/(ss_b+ss_w), 100*ss_w/(ss_b+ss_w)))
cat(sprintf("  familias con score idéntico: %d de %d (y %d de las %d con más de un miembro)\n",
            sum(idn$u == 1), nrow(idn), sum(idn$u == 1 & idn$k > 1), sum(idn$k > 1)))
# Inventory of the ZG16B probes on the array
zg_all <- an[grepl("ZG16B", UCSC_RefGene_Name)][order(pos)]
zg_all[, analizada := cpg %in% rownames(ttM)]                # did it survive the filter?
zg_all[, grupo := ""]
for (i in seq_len(nrow(zg_all)))
  set(zg_all, i, "grupo", paste(unique(strsplit(zg_all$UCSC_RefGene_Group[i], ";")[[1]]), collapse = "/"))
cat(sprintf("\nZG16B en el array: %d sondas anotadas | %d en promotor | %d analizadas\n",
            nrow(zg_all), sum(grepl("TSS", zg_all$grupo)), sum(zg_all$analizada)))
print(zg_all[, .(cpg, pos, grupo, analizada)])
fwrite(zg_all[, .(cpg, pos, grupo, analizada)], file.path(OUT, "S2_inventario_sondas_ZG16B.csv"))

# ---- 12. Saved artifacts and session information ----------------------------------------------------
saveRDS(list(rho_familia = rho,
             locus_primario = zg, locus_familia = comp,
             n_dmr = nrow(reg), efectos = T3),
        file.path(OUT, "sensibilidad_01b.RDS"))
writeLines(capture.output(sessionInfo()),
           file.path(OUT, "sessionInfo_01b.txt"))
cat("\n=== RESUMEN PARA EL MANUSCRITO ===\n")
cat(sprintf("A · Δβ de los 22 CpGs             -> T3_efectos_cohorte1.csv\n"))
cat(sprintf("B · correlación intra-familia ρ = %.4f; sondas de ZG16B con FDR<0.05: %d de 8\n",
            rho, sum(comp$fdr_familia < 0.05)))
cat(sprintf("C · regiones diferencialmente metiladas: %d\n", nrow(reg)))
cat(sprintf("D · permutación por familia: no ejecutada (ver la sección D)\n"))
