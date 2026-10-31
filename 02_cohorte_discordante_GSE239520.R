# 02 — Third cohort (GSE239520): evaluation of the fixed directional prediction
#
# Evaluates, without modification, the directional prediction fixed in 01 at the eight ZG16B promoter
# probes, in buccal epithelium from 85 children assessed for maltreatment. Processing is identical to that
# of the replication cohort: IDATs -> preprocessRaw -> dasen -> EpiDISH, then limma with age, sex and
# epithelial fraction (natural splines). The statistic is the number of probes carrying the predicted
# sign, tested against one half with a one-sided exact binomial test.
# Paper section: Results 3.6; Supplementary Table S3
# Input:  GSE239520 IDATs and SOFT file (CM)
# Output: outputs/S1_cohorte_discordante_GSE239520.csv (Supplementary Table S3), cohorte_discordante.RDS

# ---- 1. Packages, paths and locus -------------------------------------------------------------------
suppressMessages({
  library(data.table); library(minfi); library(wateRmelon); library(limma)
  library(EpiDISH); library(splines)
})
set.seed(42); CM <- "data/evaluation"; OUT <- "outputs"   # GEO files (edit CM to the local path)
dir.create(OUT, showWarnings = FALSE)

# The eight analyzed promoter probes and the sign fixed as prediction in 01
LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",   # TSS200 -> negative
          "cg19128189","cg08829024","cg12475088","cg04543077"),  # TSS1500 -> positive
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))

# ---- 2. Phenotype and IDAT basenames from the SOFT file ---------------------------------------------
li  <- gsub("\r", "", readLines(file.path(CM, "GSE239520.txt")))   # SOFT files carry Windows line endings
idx <- grep("^\\^SAMPLE", li)                                # first line of each sample block
gsm <- sub(".*= ", "", grep("^!Sample_geo_accession", li, value = TRUE))
blo <- cut(seq_along(li), breaks = c(idx, length(li) + 1), right = FALSE, labels = FALSE)   # sample of each line
ch  <- data.table(b = blo[grep("^!Sample_characteristics_ch1", li)],
                  t = sub("^!Sample_characteristics_ch1 = ", "",
                          grep("^!Sample_characteristics_ch1", li, value = TRUE)))
ch[, campo := trimws(sub(":.*", "", t))][, valor := trimws(sub("^[^:]*:", "", t))][, geo := gsm[b]]
ph <- as.data.table(dcast(ch, geo ~ campo, value.var = "valor")); setnames(ph, make.names(names(ph)))
sup <- grep("^!Sample_supplementary_file", li, value = TRUE)
b2  <- blo[grep("^!Sample_supplementary_file", li)]
bs  <- data.table(geo = gsm[b2], f = sub(".*= ", "", sup))[grepl("_Grn", f)]   # one IDAT pair per sample
bs[, base := sub("_Grn.idat.gz$", "", basename(f))]
ph  <- merge(ph, bs[, .(geo, base)], by = "geo")

# ---- 3. Processing, identical to the replication cohort ---------------------------------------------
bases <- file.path(CM, "GSE239520", ph$base)
okf <- file.exists(paste0(bases,"_Grn.idat.gz")) & file.exists(paste0(bases,"_Red.idat.gz"))
ph <- ph[okf]; bases <- bases[okf]                           # samples with both IDATs present
rg <- read.metharray(basenames = bases, force = TRUE); colnames(rg) <- ph$geo
ms <- preprocessRaw(rg); rm(rg); gc()
B <- getBeta(dasen(ms)); rm(ms); gc()                        # dasen normalization (wateRmelon)
data(centEpiFibIC.m)                                         # EpiDISH epithelial/fibroblast/immune reference
ph[, Epi := epidish(beta.m = B, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
ph[, malt := ifelse(grepl("maltreat", group, ignore.case=TRUE), 1,
            ifelse(grepl("typically", group, ignore.case=TRUE), 0, NA))]   # 1 = maltreatment, 0 = comparison
ph[, edad := suppressWarnings(as.numeric(age))]
cat(sprintf("GSE239520 · %d sondas x %d muestras | %d maltrato / %d control\n",
            nrow(B), ncol(B), sum(ph$malt==1,na.rm=TRUE), sum(ph$malt==0,na.rm=TRUE)))
cat(sprintf("  fracción epitelial: mediana %.3f | rango %.3f-%.3f\n",
            median(ph$Epi), min(ph$Epi), max(ph$Epi)))

# M-values of the eight probes, with betas clipped away from 0 and 1
pres <- LOCUS$cpg[LOCUS$cpg %in% rownames(B)]
M <- B[pres, , drop = FALSE]
M[] <- log2(pmin(pmax(M,1e-4),1-1e-4)/(1-pmin(pmax(M,1e-4),1-1e-4))); M[!is.finite(M)] <- 0
rm(B); gc()
ok <- !is.na(ph$malt) & !is.na(ph$Epi) & !is.na(ph$edad)
d <- data.frame(y = ph$malt[ok], edad = ph$edad[ok],
                sexo = droplevels(factor(ph$gender[ok])), Epi = ph$Epi[ok])
for (v in names(d)) if (is.factor(d[[v]]) && nlevels(d[[v]]) < 2) d[[v]] <- NULL   # drop single-level factors
dis <- model.matrix(~ y + ns(edad,3) + sexo + ns(Epi,4), data = d)
stopifnot(qr(dis)$rank == ncol(dis))                         # full rank, as in 01 and 01b
tt <- limma::topTable(eBayes(lmFit(M[, ok, drop=FALSE], dis)), coef = 2, number = Inf, sort.by = "none")

# ---- 4. Result: sign concordance and binomial test --------------------------------------------------
zg <- merge(LOCUS, data.table(cpg = rownames(tt), beta = round(tt$logFC,4),
                              t = round(tt$t,2), p = signif(tt$P.Value,3)), by = "cpg")
zg[, dir_ok := sign(beta) == dir_pred]; setorder(zg, region, p)
cat(sprintf("\n=== GSE239520 · epitelio bucal (n=%d) ===\n", sum(ok))); print(zg)
n_ok <- sum(zg$dir_ok); bt <- binom.test(n_ok, nrow(zg), 0.5, alternative = "greater")
cat(sprintf("\n  TSS200: %d de 4 negativas | TSS1500: %d de 4 positivas\n",
            zg[region=="TSS200", sum(dir_ok)], zg[region=="TSS1500", sum(dir_ok)]))
cat(sprintf("  DIRECCIÓN PREDICHA: %d de %d | binomial p = %.4f  %s\n",
            n_ok, nrow(zg), bt$p.value, ifelse(bt$p.value < 0.05, "REPLICA", "NO replica")))
fwrite(zg, file.path(OUT, "S1_cohorte_discordante_GSE239520.csv"))
saveRDS(list(locus = zg, n = sum(ok), n_dir = n_ok, p_binom = bt$p.value,
             epi_mediana = median(ph$Epi)), file.path(OUT, "cohorte_discordante.RDS"))
# The storage variable deposited with this series is described in 02b and examined by exclusion in 02c.
