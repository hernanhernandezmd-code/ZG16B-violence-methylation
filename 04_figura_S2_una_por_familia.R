# 04 — Supplementary Fig. S2: ZG16B promoter probes with all participants and with one per family
#
# Shows the main limitation and its check in one figure: the eight promoter coefficients in the full
# discovery sample (n = 116) and in one participant per extended family (n = 42, the participants of
# S9_una_por_familia_seleccion.csv written by 01d). Same design matrix as the article and 01d. The
# coefficients from lm() equal the limma estimates and are checked against Table 3 (T4) and
# Supplementary Table S7 (S9_una_por_familia_locus.csv) with stopifnot; the error bars are 95% ordinary
# least squares confidence intervals.
# Paper section: Results 3.7; Limitations; Supplementary Fig. S2
# Input:  GSE226085 processed matrix and SOFT (BU); outputs/S9_una_por_familia_seleccion.csv,
#         outputs/T4_ZG16B_correlacion_familiar.csv, outputs/S9_una_por_familia_locus.csv
# Output: outputs/HernandezFigS2.tif (600 dpi),
#         outputs/FigS2_una_por_familia.csv, outputs/sessionInfo_04.txt

# ---- 1. Packages, paths and locus -------------------------------------------------------------------
suppressMessages({ library(data.table); library(EpiDISH); library(splines)
                   library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) })
set.seed(42)
BU  <- "data/GSE226085"                                      # GSE226085 (edit to the local path)
OUT <- "outputs"
LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",   # TSS200 -> negative
          "cg19128189","cg08829024","cg12475088","cg04543077"),  # TSS1500 -> positive
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))

# ---- 2. Phenotype and methylation of the discovery cohort -------------------------------------------
# SOFT characteristics to columns (same parsing as 01)
li  <- gsub("\r", "", readLines(file.path(BU, "GSE226085.txt")))
idx <- grep("^\\^SAMPLE", li); gsm <- sub(".*= ", "", grep("^!Sample_geo_accession", li, value = TRUE))
blo <- cut(seq_along(li), breaks = c(idx, length(li) + 1), right = FALSE, labels = FALSE)
ch  <- data.table(b = blo[grep("^!Sample_characteristics_ch1", li)],
                  t = sub("^!Sample_characteristics_ch1 = ", "", grep("^!Sample_characteristics_ch1", li, value = TRUE)))
ch[, campo := trimws(sub(":.*", "", t))][, valor := trimws(sub("^[^:]*:", "", t))][, geo := gsm[b]]
ph <- as.data.table(dcast(ch, geo ~ campo, value.var = "valor")); setnames(ph, make.names(names(ph)))
ph[, `:=`(y = suppressWarnings(as.numeric(violence.trauma)), edad = suppressWarnings(as.numeric(age)),
          sexo = gender, lote = plate, col = individual)]
ph[, usar := grepl("^Y|^y|TRUE|1", passed.methylation.quality.control) &
             grepl("^Y|^y|TRUE|1", included.in.final.analytical.sample)]

cab <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"), nrows = 0)
B   <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"),
             select = setdiff(names(cab), grep("Detection", names(cab), value = TRUE)), showProgress = FALSE)
cg <- B[[1]]; B <- as.matrix(B[, -1]); rownames(B) <- cg
com <- intersect(colnames(B), ph$col); ph <- ph[match(com, col)]; B <- B[, com, drop = FALSE]
B   <- B[rowMeans(is.na(B)) < 0.05 & !grepl("^rs|^ch\\.", rownames(B)), , drop = FALSE]   # article filter
data(centEpiFibIC.m)
ph[, Epi := epidish(beta.m = B, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
ok <- ph$usar & !is.na(ph$y) & !is.na(ph$edad) & !is.na(ph$Epi); ok[is.na(ok)] <- FALSE
stopifnot(sum(ok) == 116)                                    # analytic sample of the article
M <- pmin(pmax(B[LOCUS$cpg, ok], 1e-4), 1 - 1e-4); M <- log2(M / (1 - M))   # M-values
rm(B); gc()
d <- data.frame(y = ph$y[ok], edad = ph$edad[ok], sexo = factor(ph$sexo[ok]),
                Epi = ph$Epi[ok], lote = factor(ph$lote[ok]), col = ph$col[ok])

# ---- 3. Coefficients and 95% intervals in both samples ----------------------------------------------
sel <- fread(file.path(OUT, "S9_una_por_familia_seleccion.csv"))$col   # one participant per family (01d)
stopifnot(length(sel) == 42, all(sel %in% d$col))
muestras <- list("completa" = rep(TRUE, nrow(d)), "una por familia" = d$col %in% sel)
filas <- list()
for (etiqueta in names(muestras)) {                          # full sample, then one per family
  s  <- muestras[[etiqueta]]
  ds <- droplevels(d[s, ])
  for (cg in LOCUS$cpg) {                                    # one model per probe
    ds$m <- M[cg, s]
    f <- lm(m ~ y + ns(edad, 3) + sexo + lote + ns(Epi, 4), data = ds)   # design of the article
    ci <- confint(f)["y", ]
    filas[[length(filas) + 1]] <- data.table(cpg = cg, muestra = etiqueta, n = nrow(ds), coef = coef(f)["y"],
                                             lo = ci[1], hi = ci[2], p = summary(f)$coefficients["y", 4])
  }
}
r <- rbindlist(filas)

# Control: coefficients must equal the published ones (T4 full sample, S9 one per family)
t4 <- fread(file.path(OUT, "T4_ZG16B_correlacion_familiar.csv")); s9 <- fread(file.path(OUT, "S9_una_por_familia_locus.csv"))
chk <- merge(merge(r[muestra == "completa", .(cpg, c1 = coef)], t4[, .(cpg, M_primario)]),
             merge(r[muestra == "una por familia", .(cpg, c2 = coef)], s9[, .(cpg, M_est)]))
print(chk[, .(cpg, c1 = round(c1, 4), M_primario = round(M_primario, 4), c2 = round(c2, 4), M_est = round(M_est, 4))])
stopifnot(max(abs(chk$c1 - chk$M_primario)) < 1e-6, max(abs(chk$c2 - chk$M_est)) < 1e-6)

an <- as.data.table(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19), keep.rownames = "cpg")[cpg %in% LOCUS$cpg, .(cpg, pos)]
r  <- merge(merge(r, LOCUS, by = "cpg"), an, by = "cpg")
r[, kb := pos / 1000 + ifelse(muestra == "completa", -0.012, 0.012)]   # small horizontal offset per sample
fwrite(r[, .(cpg, region, pos, muestra, n, coef = round(coef, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 3))],
       file.path(OUT, "FigS2_una_por_familia.csv"))

# ---- 4. Figure ---------------------------------------------------------------------------------------
# Same conventions as Fig. 3: circles TSS200, triangles TSS1500, vertical line at the region boundary
tiff(file.path(OUT, "HernandezFigS2.tif"), width = 6.5, height = 4, units = "in", res = 600,
     compression = "lzw")
par(mar = c(4.2, 4.6, 1.2, 1), family = "sans", cex.axis = .9, tcl = -0.35, mgp = c(2.6, .7, 0))
yl <- range(c(r$lo, r$hi)) * 1.1
plot(r$kb, r$coef, type = "n", ylim = yl, xlim = range(r$pos) / 1000 + c(-.06, .06),
     xlab = "Position on chromosome 16 (kb, hg19)", ylab = "Coefficient per violent event (M-value)")
abline(h = 0, col = "grey60"); abline(v = 2879.771, col = "grey75")
segments(r$kb, r$lo, r$kb, r$hi, col = ifelse(r$muestra == "completa", "black", "grey45"))
simbolo <- ifelse(r$region == "TSS200", ifelse(r$muestra == "completa", 19, 1),
                  ifelse(r$muestra == "completa", 17, 2))              # filled = all, open = one per family
points(r$kb, r$coef, pch = simbolo, col = ifelse(r$muestra == "completa", "black", "grey35"), cex = 1.3)
legend("bottomleft", c("All participants (n = 116)", "One per family (n = 42)"), pch = c(19, 1),
       col = c("black", "grey35"), bty = "n", cex = .85)
legend("bottomright", c("TSS1500", "TSS200"), pch = c(17, 19), bty = "n", cex = .85, horiz = TRUE)
invisible(dev.off())
cat("\nFigura S2 escrita\n")
print(r[order(muestra, pos), .(cpg, region, muestra, coef = round(coef, 4), lo = round(lo, 4), hi = round(hi, 4), p = signif(p, 3))])
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_04.txt"))
