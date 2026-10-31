# 01 — Discovery (GSE226085) and targeted replication (GSE239521)
#
# Discovery cohort: buccal epithelium, with exposure to violence as a continuous count (0-19), limma model
# with age (natural splines), sex, plate and epithelial fraction (natural splines). Replication cohort:
# saliva, with binary child maltreatment. The directional prediction at the eight ZG16B promoter probes is
# fixed below, before the replication cohort is read, and is evaluated with a one-sided exact binomial
# test on the sign. The exposure of the discovery cohort is a family-level characteristic; the
# within-family correlation model and the primary analysis on the imputed matrix are in 01b.
# Paper section: Methods 2.2-2.5; Results 3.3 and 3.5; Table 3
# Input:  GSE226085 processed matrix and SOFT file (BU); GSE239521 IDATs and SOFT file (CM)
# Output: outputs/T1_hits_genomicos_cohorte1.csv, T2_locus_ZG16B_dos_cohortes.csv, tronco_ZG16B.RDS,
#         sessionInfo_tronco.txt

# ---- 1. Packages, paths and the fixed locus ---------------------------------------------------------
suppressMessages({
  library(data.table); library(minfi); library(wateRmelon); library(limma)
  library(EpiDISH); library(splines)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
})
set.seed(42)
BU  <- "data/GSE226085"            # GSE226085, processed matrix deposited by the authors (edit to the local path)
CM  <- "data/evaluation"           # GSE239520 and GSE239521, IDATs per sample (edit to the local path)
OUT <- "outputs"
dir.create(OUT, showWarnings = FALSE)

# The prediction: negative coefficients at the four TSS200 probes, positive at the four TSS1500 probes
LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",
          "cg19128189","cg08829024","cg12475088","cg04543077"),
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))

# ---- 2. Phenotypes of both cohorts from their GEO SOFT files ----------------------------------------
fenotipos <- list()
for (g in c("GSE226085", "GSE239521")) {
  fichero <- if (g == "GSE226085") file.path(BU, "GSE226085.txt") else file.path(CM, "GSE239521.txt")
  li  <- gsub("\r", "", readLines(fichero))                  # SOFT files carry Windows line endings
  idx <- grep("^\\^SAMPLE", li)                              # first line of each sample block
  gsm <- sub(".*= ", "", grep("^!Sample_geo_accession", li, value = TRUE))
  blo <- cut(seq_along(li), breaks = c(idx, length(li) + 1), right = FALSE, labels = FALSE)   # sample of each line
  ch  <- data.table(b = blo[grep("^!Sample_characteristics_ch1", li)],
                    t = sub("^!Sample_characteristics_ch1 = ", "",
                            grep("^!Sample_characteristics_ch1", li, value = TRUE)))
  ch[, campo := trimws(sub(":.*", "", t))][, valor := trimws(sub("^[^:]*:", "", t))][, geo := gsm[b]]
  ph <- as.data.table(dcast(ch, geo ~ campo, value.var = "valor"))   # one row per GSM
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

# ---- 3. Cohort 1: discovery (GSE226085) -------------------------------------------------------------
# Deposited processed matrix; probe filter (5% or more missing, rs and ch. probes); EpiDISH (centEpiFibIC);
# analytic sample = authors' quality control + final analytic sample + complete covariates
ph1 <- fenotipos[["GSE226085"]]
ph1[, `:=`(violencia = suppressWarnings(as.numeric(violence.trauma)), edad = suppressWarnings(as.numeric(age)),
           sexo = gender, familia = extended.family, lote = plate, col = individual)]
ph1[, usar := grepl("^Y|^y|TRUE|1", passed.methylation.quality.control) &
              grepl("^Y|^y|TRUE|1", included.in.final.analytical.sample)]

cab <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"), nrows = 0)
B1 <- fread(file.path(BU, "GSE226085_geo_igp_processed_matrix_v1.csv.gz"),
            select = setdiff(names(cab), grep("Detection", names(cab), value = TRUE)),   # detection p columns not loaded
            showProgress = FALSE)
cg <- B1[[1]]; B1 <- as.matrix(B1[, -1]); rownames(B1) <- cg
com <- intersect(colnames(B1), ph1$col); ph1 <- ph1[match(com, col)]; B1 <- B1[, com, drop = FALSE]
stopifnot(identical(colnames(B1), ph1$col))                  # phenotype aligned with the matrix
B1 <- B1[rowMeans(is.na(B1)) < 0.05 & !grepl("^rs|^ch\\.", rownames(B1)), , drop = FALSE]
data(centEpiFibIC.m)
ph1[, Epi := epidish(beta.m = B1, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
cat(sprintf("COHORTE 1 · %d sondas x %d muestras | fracción epitelial mediana %.3f\n",
            nrow(B1), ncol(B1), median(ph1$Epi)))

# M-values with betas clipped away from 0 and 1
M1 <- B1; M1[] <- log2(pmin(pmax(M1,1e-4),1-1e-4)/(1-pmin(pmax(M1,1e-4),1-1e-4))); M1[!is.finite(M1)] <- 0
rm(B1); gc()
ok1 <- ph1$usar & !is.na(ph1$violencia) & !is.na(ph1$edad) & !is.na(ph1$Epi); ok1[is.na(ok1)] <- FALSE

# Model: violence + age (spline) + sex + plate + epithelial fraction (spline); single-level factors dropped
d <- data.frame(y = ph1$violencia[ok1], edad = ph1$edad[ok1], sexo = droplevels(factor(ph1$sexo[ok1])),
                Epi = ph1$Epi[ok1], lote = droplevels(factor(ph1$lote[ok1])))
for (v in names(d)) if (is.factor(d[[v]]) && nlevels(d[[v]]) < 2) d[[v]] <- NULL
trm <- c("y", "ns(edad, 3)")
if ("sexo" %in% names(d)) trm <- c(trm, "sexo")
if ("lote" %in% names(d)) trm <- c(trm, "lote")
trm <- c(trm, "ns(Epi, 4)")
dis <- model.matrix(as.formula(paste("~", paste(trm, collapse = " + "))), data = d)
stopifnot(qr(dis)$rank == ncol(dis))
tt1 <- limma::topTable(eBayes(lmFit(M1[, ok1, drop = FALSE], dis)), coef = 2, number = Inf, sort.by = "none")
# Genomic inflation factor: median observed chi-square over its expectation
lambda1 <- median(qchisq(tt1$P.Value, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)
cat(sprintf("  EWAS genómico (n=%d): λ = %.2f | CpGs a FDR<0.05: %d | Bonferroni: %d\n",
            sum(ok1), lambda1, sum(tt1$adj.P.Val < 0.05),
            sum(tt1$P.Value < 0.05/nrow(tt1))))

# ---- 4. Genome-wide hits of cohort 1 (FDR < 0.05) ---------------------------------------------------
an <- as.data.table(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19), keep.rownames = "cpg")
sig <- rownames(tt1)[tt1$adj.P.Val < 0.05]
hits <- merge(data.table(cpg = sig, beta = round(tt1[sig,"logFC"],4),
                         p = signif(tt1[sig,"P.Value"],3), fdr = signif(tt1[sig,"adj.P.Val"],3)),
              an[, .(cpg, chr, pos, gene = UCSC_RefGene_Name, grupo = UCSC_RefGene_Group,
                     isla = Relation_to_Island)], by = "cpg")
for (i in seq_len(nrow(hits))) {                             # one entry per gene and per gene region
  set(hits, i, "gene",  paste(unique(strsplit(hits$gene[i],  ";")[[1]]), collapse = "/"))
  set(hits, i, "grupo", paste(unique(strsplit(hits$grupo[i], ";")[[1]]), collapse = "/"))
}
cat("\n=== TABLA 1 · CpGs a FDR<0.05 en el metiloma completo (cohorte 1) ===\n")
print(hits[order(p)])
fwrite(hits[order(p)], file.path(OUT, "T1_hits_genomicos_cohorte1.csv"))

# ---- 5. The ZG16B locus in cohort 1: the eight promoter probes and their sign -----------------------
zg1 <- merge(LOCUS, data.table(cpg = rownames(tt1), beta = round(tt1$logFC,4),
                               p = signif(tt1$P.Value,3), fdr = signif(tt1$adj.P.Val,3)), by = "cpg")
zg1 <- merge(zg1, an[, .(cpg, pos)], by = "cpg"); setorder(zg1, pos)
zg1[, dir_ok := sign(beta) == dir_pred]
cat("\n=== TABLA 2 · locus ZG16B en la cohorte 1 (descubrimiento) ===\n"); print(zg1)
cat(sprintf("\n  TSS200: %d de 4 negativas | TSS1500: %d de 4 positivas | FDR<0.05: %d sondas\n",
            zg1[region=="TSS200", sum(beta<0)], zg1[region=="TSS1500", sum(beta>0)],
            sum(zg1$fdr < 0.05)))
rm(M1); gc()

# ---- 6. Cohort 2: replication (GSE239521, saliva) ---------------------------------------------------
# IDATs -> preprocessRaw -> dasen; same model without plate; only the sign of the eight probes is evaluated
# against the prediction fixed above
ph2 <- fenotipos[["GSE239521"]]
bases <- file.path(CM, "GSE239521", ph2$base)
okf <- file.exists(paste0(bases,"_Grn.idat.gz")) & file.exists(paste0(bases,"_Red.idat.gz"))
ph2 <- ph2[okf]; bases <- bases[okf]                         # samples with both IDATs present
cat(sprintf("COHORTE 2 · IDATs: %d\n", length(bases)))
rg <- read.metharray(basenames = bases, force = TRUE); colnames(rg) <- ph2$geo
ms <- preprocessRaw(rg); rm(rg); gc()
B2 <- getBeta(dasen(ms)); rm(ms); gc()                       # dasen normalization (wateRmelon)
ph2[, Epi := epidish(beta.m = B2, ref.m = centEpiFibIC.m, method = "RPC")$estF[, "Epi"]]
ph2[, malt := ifelse(grepl("maltreat", group, ignore.case=TRUE), 1,
             ifelse(grepl("typically", group, ignore.case=TRUE), 0, NA))]   # 1 = maltreatment, 0 = comparison
ph2[, edad := suppressWarnings(as.numeric(age))]
cat(sprintf("  %d sondas x %d muestras | %d maltrato / %d control | fracción epitelial mediana %.3f\n",
            nrow(B2), ncol(B2), sum(ph2$malt==1,na.rm=TRUE), sum(ph2$malt==0,na.rm=TRUE),
            median(ph2$Epi)))

pres <- LOCUS$cpg[LOCUS$cpg %in% rownames(B2)]
M2 <- B2[pres, , drop = FALSE]
M2[] <- log2(pmin(pmax(M2,1e-4),1-1e-4)/(1-pmin(pmax(M2,1e-4),1-1e-4))); M2[!is.finite(M2)] <- 0
rm(B2); gc()
ok2 <- !is.na(ph2$malt) & !is.na(ph2$Epi) & !is.na(ph2$edad)
d <- data.frame(y = ph2$malt[ok2], edad = ph2$edad[ok2], sexo = droplevels(factor(ph2$gender[ok2])),
                Epi = ph2$Epi[ok2])
for (v in names(d)) if (is.factor(d[[v]]) && nlevels(d[[v]]) < 2) d[[v]] <- NULL
trm <- c("y", "ns(edad, 3)")
if ("sexo" %in% names(d)) trm <- c(trm, "sexo")
trm <- c(trm, "ns(Epi, 4)")
dis <- model.matrix(as.formula(paste("~", paste(trm, collapse = " + "))), data = d)
stopifnot(qr(dis)$rank == ncol(dis))
tt2 <- limma::topTable(eBayes(lmFit(M2[, ok2, drop = FALSE], dis)), coef = 2, number = Inf, sort.by = "none")
zg2 <- merge(LOCUS, data.table(cpg = rownames(tt2), beta = round(tt2$logFC,4),
                               t = round(tt2$t,2), p = signif(tt2$P.Value,3)), by = "cpg")
zg2[, dir_ok := sign(beta) == dir_pred]; setorder(zg2, region, p)
cat("\n=== TABLA 3 · réplica dirigida del locus ZG16B en saliva (cohorte 2) ===\n"); print(zg2)

n_ok <- sum(zg2$dir_ok)
bt <- binom.test(n_ok, nrow(zg2), 0.5, alternative = "greater")   # one-sided exact binomial test
cat(sprintf("\n  TSS200: %d de 4 negativas | TSS1500: %d de 4 positivas\n",
            zg2[region=="TSS200", sum(dir_ok)], zg2[region=="TSS1500", sum(dir_ok)]))
cat(sprintf("  DIRECCIÓN PREDICHA: %d de %d | binomial de una cola p = %.4f %s\n",
            n_ok, nrow(zg2), bt$p.value, ifelse(bt$p.value < 0.05, "<- REPLICA", "")))

# ---- 7. Joint result and saved artifacts ------------------------------------------------------------
conj <- rbind(zg1[, .(cohorte = "1 · GSE226085 (bucal, violencia)", cpg, region, beta, p, dir_ok)],
              zg2[, .(cohorte = "2 · GSE239521 (saliva, maltrato)", cpg, region, beta, p, dir_ok)])
fwrite(conj, file.path(OUT, "T2_locus_ZG16B_dos_cohortes.csv"))
saveRDS(list(cohorte1 = list(hits = hits, locus = zg1, lambda = lambda1, n = sum(ok1)),
             cohorte2 = list(locus = zg2, n = sum(ok2), n_dir = n_ok, p_binom = bt$p.value),
             sesion = sessionInfo()),
        file.path(OUT, "tronco_ZG16B.RDS"))
cat("\n=== RESULTADO DEL TRONCO ===\n")
cat(sprintf("  Cohorte 1 (descubrimiento, bucal, n=%d): λ=%.2f, %d CpGs a FDR<0.05,\n",
            sum(ok1), lambda1, sum(tt1$adj.P.Val < 0.05)))
cat(sprintf("     de los cuales %d en el promotor de ZG16B; patrón espacial TSS200 %d/4 - TSS1500 %d/4\n",
            sum(zg1$fdr < 0.05), zg1[region=="TSS200", sum(beta<0)], zg1[region=="TSS1500", sum(beta>0)]))
cat(sprintf("  Cohorte 2 (réplica, saliva, n=%d): %d de 8 sondas en la dirección predicha, p=%.4f\n",
            sum(ok2), n_ok, bt$p.value))
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_tronco.txt"))
