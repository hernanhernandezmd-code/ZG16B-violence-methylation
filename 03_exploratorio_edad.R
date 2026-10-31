# 03 — Exploratory analysis of age at the eight ZG16B promoter probes
#
# Post hoc analysis prompted by the discordance between the two evaluation cohorts: the prediction held in
# adolescents (saliva) and was inverted in young children (buccal epithelium). It asks whether the exposure
# effect at the eight promoter probes depends on age. It does not change the discovery analysis or the
# pre-specified test. The published coefficients are reproduced first as a control (stopifnot). Then, in
# each cohort: (a) sign concordance within age groups (discovery: <10, 10-19 and >=20 years; evaluation
# cohorts: split at the median age); (b) an exposure-by-age interaction (age centered at the median),
# counting the probes whose interaction moves the exposure effect toward the predicted sign. Per-probe
# ordinary least squares with lm(): eight probes, not a genome-wide scan, so no variance moderation.
# Paper section: Results 3.8 (Exploratory analysis of age); Supplementary Table S10
# Input:  GSE226085 processed matrix and SOFT file (BU); GSE239520 and GSE239521 IDATs and SOFT files (CM);
#         outputs/T4_ZG16B_correlacion_familiar.csv and outputs/TABLA4_locus_tres_cohortes.csv (controls)
# Output: outputs/E1_edad_resumen.csv, E1_edad_locus.csv, E1_interaccion_resumen.csv,
#         E1_interaccion_locus.csv, sessionInfo_03.txt

# ---- 1. Packages, paths and locus -------------------------------------------------------------------
suppressMessages({ library(data.table); library(minfi); library(wateRmelon)
                   library(EpiDISH); library(splines) })
set.seed(42)                                                 # same seed as the other scripts
BU  <- "data/GSE226085"                                      # GSE226085 (edit to the local path)
CM  <- "data/evaluation"                                     # GSE239520 and GSE239521 (edit to the local path)
OUT <- "outputs"

# The eight analyzed promoter probes and the sign fixed as prediction
LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",   # TSS200 -> negative
          "cg19128189","cg08829024","cg12475088","cg04543077"),  # TSS1500 -> positive
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))
data(centEpiFibIC.m)                                         # EpiDISH epithelial/fibroblast/immune reference

# ---- 2. Phenotypes and IDAT basenames from the GEO SOFT files ---------------------------------------
fenotipos <- list()
for (g in c("GSE226085", "GSE239521", "GSE239520")) {
  fichero <- if (g == "GSE226085") file.path(BU, "GSE226085.txt") else file.path(CM, paste0(g, ".txt"))
  li  <- gsub("\r", "", readLines(fichero))                  # SOFT files carry Windows line endings
  idx <- grep("^\\^SAMPLE", li)
  gsm <- sub(".*= ", "", grep("^!Sample_geo_accession", li, value = TRUE))
  blo <- cut(seq_along(li), breaks = c(idx, length(li) + 1), right = FALSE, labels = FALSE)   # sample of each line
  ch  <- data.table(b = blo[grep("^!Sample_characteristics_ch1", li)],
                    t = sub("^!Sample_characteristics_ch1 = ", "",
                            grep("^!Sample_characteristics_ch1", li, value = TRUE)))
  ch[, campo := trimws(sub(":.*", "", t))][, valor := trimws(sub("^[^:]*:", "", t))][, geo := gsm[b]]
  ph <- as.data.table(dcast(ch, geo ~ campo, value.var = "valor"))
  setnames(ph, make.names(names(ph)))
  sup <- grep("^!Sample_supplementary_file", li, value = TRUE)
  if (length(sup)) {                                         # IDAT basenames, when the series deposits IDATs
    b2 <- blo[grep("^!Sample_supplementary_file", li)]
    bs <- data.table(geo = gsm[b2], f = sub(".*= ", "", sup))[grepl("_Grn", f)]
    bs[, base := sub("_Grn.idat.gz$", "", basename(f))]
    ph <- merge(ph, bs[, .(geo, base)], by = "geo", all.x = TRUE)
  }
  fenotipos[[g]] <- ph[]
}

# ---- 3. Discovery cohort (GSE226085) ----------------------------------------------------------------
datos <- list()                                              # M-values and covariates of each cohort
ph1 <- fenotipos[["GSE226085"]]
ph1[, `:=`(y = suppressWarnings(as.numeric(violence.trauma)), edad = suppressWarnings(as.numeric(age)),
           sexo = gender, lote = plate, col = individual, gen = generation)]
ph1[, usar := grepl("^Y|^y|TRUE|1", passed.methylation.quality.control) &
              grepl("^Y|^y|TRUE|1", included.in.final.analytical.sample)]
cab <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"), nrows = 0)
B1  <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"),
             select = setdiff(names(cab), grep("Detection", names(cab), value = TRUE)), showProgress = FALSE)
cg <- B1[[1]]; B1 <- as.matrix(B1[, -1]); rownames(B1) <- cg
com <- intersect(colnames(B1), ph1$col); ph1 <- ph1[match(com, col)]; B1 <- B1[, com, drop = FALSE]
# Same probe filter as the article: <5% missing, no SNP or non-CpG probes
B1  <- B1[rowMeans(is.na(B1)) < 0.05 & !grepl("^rs|^ch\\.", rownames(B1)), , drop = FALSE]
ph1[, Epi := epidish(beta.m = B1, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
ok <- ph1$usar & !is.na(ph1$y) & !is.na(ph1$edad) & !is.na(ph1$Epi); ok[is.na(ok)] <- FALSE
stopifnot(sum(ok) == 116)                                    # analytic sample of the article
M1 <- pmin(pmax(B1[LOCUS$cpg, ok], 1e-4), 1 - 1e-4); M1 <- log2(M1 / (1 - M1))   # M-values
rm(B1); gc()
datos[["descubrimiento"]] <- list(M = M1, d = data.frame(y = ph1$y[ok], edad = ph1$edad[ok], sexo = factor(ph1$sexo[ok]),
                                                         Epi = ph1$Epi[ok], lote = factor(ph1$lote[ok]), gen = ph1$gen[ok]))
cat("\nEdades del descubrimiento por generación:\n")
for (gn in sort(unique(ph1$gen[ok]))) {
  v <- ph1$edad[ok][ph1$gen[ok] == gn]
  cat(sprintf("  %-25s n = %2d | mediana %.1f | rango %.1f-%.1f\n", gn, length(v), median(v), min(v), max(v)))
}

# ---- 4. Evaluation cohorts (GSE239521 saliva, GSE239520 buccal) -------------------------------------
# Same processing as 01 and 02: IDATs -> dasen -> EpiDISH; maltreatment coded 1, comparison 0
nombres <- c(GSE239521 = "réplica (adolescentes, saliva)", GSE239520 = "tercera (niños, bucal)")
for (g in names(nombres)) {
  ph <- fenotipos[[g]]
  bases <- file.path(CM, g, ph$base)
  okf <- file.exists(paste0(bases, "_Grn.idat.gz")) & file.exists(paste0(bases, "_Red.idat.gz"))
  ph <- ph[okf]; bases <- bases[okf]
  rg <- read.metharray(basenames = bases, force = TRUE); colnames(rg) <- ph$geo
  ms <- preprocessRaw(rg); rm(rg); gc()
  B  <- getBeta(dasen(ms)); rm(ms); gc()
  ph[, Epi := epidish(beta.m = B, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
  ph[, y := ifelse(grepl("maltreat", group, ignore.case = TRUE), 1,
            ifelse(grepl("typically", group, ignore.case = TRUE), 0, NA))]
  ph[, edad := suppressWarnings(as.numeric(age))]
  ok <- !is.na(ph$y) & !is.na(ph$Epi) & !is.na(ph$edad)
  M <- pmin(pmax(B[LOCUS$cpg, ok], 1e-4), 1 - 1e-4); M <- log2(M / (1 - M)); rm(B); gc()
  d <- data.frame(y = ph$y[ok], edad = ph$edad[ok], sexo = factor(ph$gender[ok]), Epi = ph$Epi[ok])
  datos[[nombres[[g]]]] <- list(M = M, d = d)
  cat(sprintf("%s · edad por grupo: maltrato mediana %.1f (%.1f-%.1f) | comparación mediana %.1f (%.1f-%.1f)\n", g,
              median(d$edad[d$y == 1]), min(d$edad[d$y == 1]), max(d$edad[d$y == 1]),
              median(d$edad[d$y == 0]), min(d$edad[d$y == 0]), max(d$edad[d$y == 0])))
}

# ---- 5. Analysis plan: every model fitted below -----------------------------------------------------
# One row per model: cohort, group label, age range or half, degrees of freedom of the age and epithelial
# splines (lowered in the subgroups) and whether the exposure-by-age interaction enters
plan <- list()
plan[[length(plan) + 1]] <- list(coh = "descubrimiento", grupo = "todos", tipo = "todos", gl = c(3, 4))
for (g in list(c("menores de 10 años", 0, 10), c("10 a 19 años", 10, 20), c("20 años o más", 20, 200)))
  plan[[length(plan) + 1]] <- list(coh = "descubrimiento", grupo = g[1], tipo = "rango",
                                   lim = as.numeric(g[2:3]), gl = c(1, 2))
plan[[length(plan) + 1]] <- list(coh = "descubrimiento", grupo = "interacción", tipo = "interaccion", gl = c(3, 4))
for (coh in nombres) {
  plan[[length(plan) + 1]] <- list(coh = coh, grupo = "todos", tipo = "todos", gl = c(3, 4))
  plan[[length(plan) + 1]] <- list(coh = coh, grupo = "menor", tipo = "mitad", gl = c(1, 2))
  plan[[length(plan) + 1]] <- list(coh = coh, grupo = "mayor", tipo = "mitad", gl = c(1, 2))
  plan[[length(plan) + 1]] <- list(coh = coh, grupo = "interacción", tipo = "interaccion", gl = c(3, 4))
}

# ---- 6. Fitting: one lm() per probe for each row of the plan ----------------------------------------
R <- list(); L <- list(); I <- list(); IL <- list(); todos <- list()
for (pl in plan) {
  M <- datos[[pl$coh]]$M; d <- datos[[pl$coh]]$d
  s <- rep(TRUE, nrow(d)); grupo <- pl$grupo
  if (pl$tipo == "rango") s <- d$edad >= pl$lim[1] & d$edad < pl$lim[2]
  if (pl$tipo == "mitad") {                                  # split at the median age of the cohort
    med <- median(d$edad)
    s <- if (pl$grupo == "menor") d$edad < med else d$edad >= med
    grupo <- sprintf("mitad %s (edad %s %.1f)", pl$grupo, if (pl$grupo == "menor") "<" else ">=", med)
  }
  ds <- d[s, ]; Ms <- M[, s, drop = FALSE]
  if (pl$tipo == "rango") { ds <- droplevels(ds); if (nrow(ds) < 15) next }
  interaccion <- pl$tipo == "interaccion"
  if (interaccion) ds$edad_c <- ds$edad - median(ds$edad)    # age centered at the cohort median
  # Exposure + age (spline) + epithelial fraction (spline) [+ sex] [+ plate] [+ exposure x age]
  f <- paste("m ~ y + ns(edad,", pl$gl[1], ") + ns(Epi,", pl$gl[2], ")")
  if (length(unique(ds$sexo)) > 1) f <- paste(f, "+ sexo")  # adults in the discovery cohort are all women
  if ("lote" %in% names(ds) && length(unique(ds$lote)) > 1) f <- paste(f, "+ lote")
  if (interaccion) f <- paste(f, "+ y:edad_c")
  z <- list()
  for (cg in LOCUS$cpg) {
    ds$m <- Ms[cg, ]
    co <- summary(lm(as.formula(f), data = ds))$coefficients
    z[[cg]] <- data.table(cpg = cg, coef = co["y", 1], p = co["y", 4],
                          int = if (interaccion) co["y:edad_c", 1] else NA_real_,
                          p_int = if (interaccion) co["y:edad_c", 4] else NA_real_)
  }
  z <- merge(LOCUS, rbindlist(z), by = "cpg")
  coh <- pl$coh
  if (interaccion) {
    # The interaction moves the effect toward the prediction when it has the predicted sign:
    # more negative with age at TSS200, more positive at TSS1500
    z[, hacia_pred := sign(int) == dir_pred]
    k <- sum(z$hacia_pred)
    IL[[length(IL) + 1]] <- z[, .(cohorte = coh, cpg, region, dir_pred, int = signif(int, 3), p_int = signif(p_int, 3), hacia_pred)]
    I[[length(I) + 1]] <- data.table(cohorte = coh, n = nrow(ds), interaccion_hacia_pred = sprintf("%d de 8", k),
                                     p_binomial = signif(binom.test(k, 8, 0.5, alternative = "greater")$p.value, 3),
                                     p_int_min = signif(min(z$p_int), 3))
    next
  }
  if (pl$tipo == "todos") todos[[coh]] <- z                  # kept for the controls below
  z[, dir_ok := sign(coef) == dir_pred]                      # sign concordance with the prediction
  k <- sum(z$dir_ok)
  L[[length(L) + 1]] <- z[, .(cohorte = coh, grupo, cpg, region, dir_pred, coef = round(coef, 4), p = signif(p, 3), dir_ok)]
  R[[length(R) + 1]] <- data.table(cohorte = coh, grupo, n = nrow(ds),
                                   expuestos = if (all(ds$y %in% 0:1)) sum(ds$y) else NA,
                                   edad_mediana = round(median(ds$edad), 1),
                                   edad_rango = sprintf("%.1f-%.1f", min(ds$edad), max(ds$edad)),
                                   tss200_ok = sum(z[region == "TSS200"]$dir_ok), tss1500_ok = sum(z[region == "TSS1500"]$dir_ok),
                                   concordancia = sprintf("%d de 8", k),
                                   p_binomial = signif(binom.test(k, 8, 0.5, alternative = "greater")$p.value, 3))
}

# ---- 7. Controls: the full-cohort coefficients must equal the published ones ------------------------
t4 <- fread(file.path(OUT, "T4_ZG16B_correlacion_familiar.csv"))   # discovery, primary model
control <- merge(todos[["descubrimiento"]][, .(cpg, coef)], t4[, .(cpg, M_primario)], by = "cpg")
cat("\nControl: coeficientes del descubrimiento frente a los publicados (T4)\n")
print(control[, .(cpg, coef = round(coef, 4), M_primario = round(M_primario, 4))])
stopifnot(max(abs(control$coef - control$M_primario)) < 1e-6)
publicado <- fread(file.path(OUT, "TABLA4_locus_tres_cohortes.csv"))   # coefficients of Table 3
for (g in names(nombres)) {
  ref <- publicado[, .(cpg = `Probe ID`, pub = if (g == "GSE239521") `Repl. coef` else `Third coef`)]
  cmp <- merge(todos[[nombres[[g]]]][, .(cpg, coef)], ref, by = "cpg")
  cat(sprintf("Control %s: diferencia máxima con lo publicado = %.4f\n", g, max(abs(cmp$coef - cmp$pub))))
}

# ---- 8. Results and outputs -------------------------------------------------------------------------
res <- rbindlist(R); loc <- rbindlist(L); ires <- rbindlist(I); iloc <- rbindlist(IL)
cat("\n=== CONCORDANCIA DE SIGNO POR SUBGRUPO DE EDAD ===\n"); print(res)
cat("\n=== INTERACCIÓN EXPOSICIÓN x EDAD (¿el efecto se acerca al signo predicho con la edad?) ===\n"); print(ires)
print(iloc)
fwrite(res,  file.path(OUT, "E1_edad_resumen.csv")); fwrite(loc,  file.path(OUT, "E1_edad_locus.csv"))
fwrite(ires, file.path(OUT, "E1_interaccion_resumen.csv")); fwrite(iloc, file.path(OUT, "E1_interaccion_locus.csv"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_03.txt"))
