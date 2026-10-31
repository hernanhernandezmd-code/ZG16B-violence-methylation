# 02c — Sensitivity to storage time in the third cohort (GSE239520)
#
# Repeats the evaluation of the fixed prediction after excluding the 16 comparison samples left at room
# temperature for more than 300 days before DNA extraction (described in 02b), and compares the result
# with the full run. Processing and model are those of 02.
# Paper section: Discussion 4.1; Supplementary Table S9
# Input:  GSE239520 IDATs and SOFT file (CM)
# Output: outputs/S5_sensibilidad_almacenamiento.csv (Supplementary Table S9)

# ---- 1. Packages, paths and locus -------------------------------------------------------------------
suppressMessages({
  library(data.table); library(minfi); library(wateRmelon); library(limma)
  library(EpiDISH); library(splines)
})
set.seed(42); CM <- "data/evaluation"; OUT <- "outputs"   # GEO files (edit CM to the local path)
dir.create(OUT, showWarnings = FALSE)

LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",   # TSS200 -> negative
          "cg19128189","cg08829024","cg12475088","cg04543077"),  # TSS1500 -> positive
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))

# ---- 2. Phenotype and IDAT basenames from the SOFT file ---------------------------------------------
li  <- gsub("\r", "", readLines(file.path(CM, "GSE239520.txt")))   # SOFT files carry Windows line endings
idx <- grep("^\\^SAMPLE", li)                                # first line of each sample block
gsm <- sub(".*= ", "", grep("^!Sample_geo_accession", li, value = TRUE))
blo <- cut(seq_along(li), breaks = c(idx, length(li) + 1),   # sample of each line
           right = FALSE, labels = FALSE)
ch  <- data.table(b = blo[grep("^!Sample_characteristics_ch1", li)],
                  t = sub("^!Sample_characteristics_ch1 = ", "",
                          grep("^!Sample_characteristics_ch1", li, value = TRUE)))
ch[, campo := trimws(sub(":.*", "", t))][, valor := trimws(sub("^[^:]*:", "", t))]
ch[, geo := gsm[b]]
ph <- as.data.table(dcast(ch, geo ~ campo, value.var = "valor"))
setnames(ph, make.names(names(ph)))
sup <- grep("^!Sample_supplementary_file", li, value = TRUE)   # IDAT file names
b2  <- blo[grep("^!Sample_supplementary_file", li)]
bs  <- data.table(geo = gsm[b2], f = sub(".*= ", "", sup))[grepl("_Grn", f)]
bs[, base := sub("_Grn.idat.gz$", "", basename(f))]
ph  <- merge(ph, bs[, .(geo, base)], by = "geo")              # the 85 samples

# ---- 3. Processing identical to 02 ------------------------------------------------------------------
bases <- file.path(CM, "GSE239520", ph$base)
okf <- file.exists(paste0(bases,"_Grn.idat.gz")) & file.exists(paste0(bases,"_Red.idat.gz"))
ph <- ph[okf]; bases <- bases[okf]                           # samples with both IDATs present
rg <- read.metharray(basenames = bases, force = TRUE); colnames(rg) <- ph$geo
ms <- preprocessRaw(rg); rm(rg); gc()
B <- getBeta(dasen(ms)); rm(ms); gc()                        # dasen normalization (wateRmelon)
data(centEpiFibIC.m)                                         # EpiDISH reference
ph[, Epi := epidish(beta.m = B, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
ph[, malt := ifelse(grepl("maltreat", group, ignore.case=TRUE), 1,
            ifelse(grepl("typically", group, ignore.case=TRUE), 0, NA))]
ph[, edad := suppressWarnings(as.numeric(age))]
dias <- suppressWarnings(as.numeric(ph[["days.left.at.rt.before.the.dna.extraction"]]))   # days at room temperature
if (all(is.na(dias))) {                                      # the field name may vary across GEO versions
  cand <- grep("days", names(ph), value = TRUE, ignore.case = TRUE)
  stopifnot(length(cand) == 1); dias <- suppressWarnings(as.numeric(ph[[cand[1]]]))
}
ph[, dias := dias][, largo := dias > 300]                    # long storage
pres <- LOCUS$cpg[LOCUS$cpg %in% rownames(B)]                # the eight probes
M <- B[pres, , drop = FALSE]
M[] <- log2(pmin(pmax(M,1e-4),1-1e-4)/(1-pmin(pmax(M,1e-4),1-1e-4)))   # M-values
M[!is.finite(M)] <- 0; rm(B); gc()

# ---- 4. The two runs: full cohort and without the 16 long-storage samples ---------------------------
corridas <- list("completa" = rep(TRUE, nrow(ph)),           # must give 0 of 8, as in 02
                 "sin almacenamiento largo" = !ph$largo)     # the question of this script
res <- list()
for (etiqueta in names(corridas)) {
  ok <- corridas[[etiqueta]] & !is.na(ph$malt) & !is.na(ph$Epi) & !is.na(ph$edad)
  d <- data.frame(y = ph$malt[ok], edad = ph$edad[ok],
                  sexo = droplevels(factor(ph$gender[ok])), Epi = ph$Epi[ok])
  for (v in names(d)) if (is.factor(d[[v]]) && nlevels(d[[v]]) < 2) d[[v]] <- NULL   # drop single-level factors
  dis <- model.matrix(~ y + ns(edad,3) + sexo + ns(Epi,4), data = d)   # same model as 02
  stopifnot(qr(dis)$rank == ncol(dis))
  tt <- limma::topTable(eBayes(lmFit(M[, ok, drop=FALSE], dis)), coef = 2,
                        number = Inf, sort.by = "none")
  zg <- merge(LOCUS, data.table(cpg = rownames(tt), beta = round(tt$logFC,4),
                                p = signif(tt$P.Value,3)), by = "cpg")
  zg[, dir_ok := sign(beta) == dir_pred]; setorder(zg, region, p)
  n_ok <- sum(zg$dir_ok)
  bt <- binom.test(n_ok, nrow(zg), 0.5, alternative = "greater")   # same test as 02
  cat(sprintf("\n=== %s (n=%d: %d maltrato / %d comparación) ===\n", etiqueta, sum(ok),
              sum(ph$malt[ok]==1), sum(ph$malt[ok]==0)))
  print(zg)
  cat(sprintf("  dirección predicha: %d de %d | binomial unilateral p = %.4f\n",
              n_ok, nrow(zg), bt$p.value))
  res[[etiqueta]] <- zg[, corrida := etiqueta][]
}
fwrite(rbindlist(res), file.path(OUT, "S5_sensibilidad_almacenamiento.csv"))
cat("\nEscrito:", file.path(OUT, "S5_sensibilidad_almacenamiento.csv"), "\n")
