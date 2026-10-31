# 01c — Figures and tables of the article
#
# Reads the artifacts written by 01, 01b and 02 and draws the figures (TIFF, 1000 dpi) and the tables of
# the article without repeating any analysis.
# Paper section: Figs. 1-3 and Supplementary Fig. S1; Tables 1-3 and Supplementary Tables S1 and S11
# Input:  outputs/descriptivos_01b.RDS, T3_efectos_cohorte1.csv, T5_DMR_cohorte1.csv,
#         T6_interseccion_DMP_DMR.csv, S3_sensibilidad_replicas_tecnicas.csv, T0_comparabilidad_grupos.csv,
#         T4_ZG16B_correlacion_familiar.csv, tronco_ZG16B.RDS, cohorte_discordante.RDS
# Output: outputs/HernandezFig1-4.tif (Fig. 1 = HernandezFig1, Fig. 2 = HernandezFig3, Fig. 3 =
#         HernandezFig4, Supplementary Fig. S1 = HernandezFig2), F2_distribuciones.csv,
#         TABLA1_cohortes.csv, TABLA2_DMP.csv (Supplementary Table S1), TABLA3_DMR.csv (Table 2),
#         TABLA4_locus_tres_cohortes.csv (Table 3, complete statistics = Supplementary Table S11),
#         sessionInfo_01c.txt

# ---- 1. Packages, paths and locus -------------------------------------------------------------------
suppressMessages({
  library(data.table)
  library(minfi)                                             # getAnnotation
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(Gviz)                                              # genomic tracks
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)                 # hg19 gene models
  library(org.Hs.eg.db)                                      # gene symbol -> Entrez identifier
})
options(ucscChromosomeNames = FALSE)                         # Gviz without querying UCSC online
set.seed(42)
OUT <- "outputs"
stopifnot(dir.exists(OUT))                                   # 01b must have been run before
LOCUS <- data.table(
  cpg = c("cg01037823","cg16934444","cg03096107","cg06274007",   # TSS200
          "cg19128189","cg08829024","cg12475088","cg04543077"),  # TSS1500
  region   = rep(c("TSS200","TSS1500"), each = 4),
  dir_pred = rep(c(-1, 1), each = 4))
an <- as.data.table(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19),   # Illumina manifest
                    keep.rownames = "cpg")

# ---- 2. Artifacts of the analysis scripts -----------------------------------------------------------
for (f in c("descriptivos_01b.RDS", "T3_efectos_cohorte1.csv", "T5_DMR_cohorte1.csv",
            "T6_interseccion_DMP_DMR.csv", "S3_sensibilidad_replicas_tecnicas.csv",
            "T0_comparabilidad_grupos.csv", "tronco_ZG16B.RDS", "cohorte_discordante.RDS",
            "T4_ZG16B_correlacion_familiar.csv"))
  stopifnot(file.exists(file.path(OUT, f)))                  # every artifact must exist
de  <- readRDS(file.path(OUT, "descriptivos_01b.RDS"))       # genome-wide p values
T3  <- fread(file.path(OUT, "T3_efectos_cohorte1.csv"))      # the 22 significant probes
T5  <- fread(file.path(OUT, "T5_DMR_cohorte1.csv"))          # the three regions
T6  <- fread(file.path(OUT, "T6_interseccion_DMP_DMR.csv"))  # the intersection
S3  <- fread(file.path(OUT, "S3_sensibilidad_replicas_tecnicas.csv"))   # the locus, eight probes
T0  <- fread(file.path(OUT, "T0_comparabilidad_grupos.csv")) # group comparability
tr  <- readRDS(file.path(OUT, "tronco_ZG16B.RDS"))           # cohorts 1 and 2
dc  <- readRDS(file.path(OUT, "cohorte_discordante.RDS"))    # cohort 3
ttM <- data.table(cpg = de$ttM_cpg, p = de$ttM_p)            # the full EWAS
ttM[, fdr := p.adjust(p, "BH")]
lambdaM <- median(qchisq(ttM$p, 1, lower.tail = FALSE), na.rm = TRUE) / qchisq(0.5, 1)   # genomic inflation factor
cat(sprintf("Leído: %d sondas analizadas | %d a FDR<0.05 | λ = %.3f\n",
            nrow(ttM), sum(ttM$fdr < 0.05), lambdaM))        # 736065, 22, 1.033

# ---- 3. Figure 1: calibration (QQ) and genome-wide results (Manhattan) ------------------------------
tiff(file.path(OUT, "HernandezFig1.tif"), width = 6, height = 3.4, units = "in",
     res = 1000, compression = "lzw")                        # line art at 1000 dpi
par(mfrow = c(1, 2), mar = c(4, 4.2, 2.4, .8), family = "sans",
    cex.axis = .85, tcl = -0.3, mgp = c(2.4, .6, 0))
o <- -log10(sort(ttM$p)); e <- -log10(ppoints(nrow(ttM)))    # observed and expected
k <- unique(round(seq(1, length(o), length.out = 20000)))    # thinning of the QQ plot
plot(e[k], o[k], pch = 20, cex = .35, xlab = expression(Expected~-log[10](p)),
     ylab = expression(Observed~-log[10](p)), main = "(a) Quantile-quantile plot",
     cex.main = .95, font.main = 1)
abline(0, 1, col = "grey55")
legend("topleft", legend = bquote(lambda == .(sprintf("%.3f", lambdaM))),
       bty = "n", cex = .95)
mh <- merge(ttM, an[, .(cpg, chr, pos)], by = "cpg")[chr %in% paste0("chr", 1:22)]
mh[, cn := as.integer(sub("chr", "", chr))]
setorder(mh, cn, pos); mh[, idx := .I]                       # cumulative position
# No thinning: any sampling conditioned on p creates a density step at the cut-off that reads as a
# feature of the data; at this resolution all 736,065 probes fit
s2 <- mh
plot(s2$idx, -log10(s2$p), pch = 20, cex = .3,
     col = ifelse(s2$cn %% 2 == 0, "grey70", "grey35"),      # alternating chromosome shades
     xlab = "Chromosome", ylab = expression(-log[10](p)), xaxt = "n",
     main = "(b) Genome-wide results", cex.main = .95, font.main = 1)
ejes <- mh[, .(m = mean(idx)), by = cn]                      # center of each chromosome
axis(1, at = ejes[cn %% 2 == 1]$m, labels = ejes[cn %% 2 == 1]$cn,   # odd chromosomes above,
     tick = FALSE, cex.axis = .55, line = -.9)
axis(1, at = ejes[cn %% 2 == 0]$m, labels = ejes[cn %% 2 == 0]$cn,   # even below, so all 22 fit
     tick = FALSE, cex.axis = .55, line = -.35)
abline(h = -log10(max(ttM$p[ttM$fdr < 0.05])), col = "grey45")   # FDR threshold
zgp <- mh[cpg %in% LOCUS$cpg]                                # ZG16B probes highlighted
points(zgp$idx, -log10(zgp$p), pch = 1, cex = 1.2)
text(zgp[which.min(p)]$idx, -log10(min(zgp$p)), "ZG16B", pos = 3, cex = .7, offset = .4)   # label above the point
invisible(dev.off())
cat("Figura 1 escrita\n")

# ---- 4. Figure 2: distribution of the tested and the significant probes -----------------------------
# Gene region and CpG-island context of the analyzed background and of the significant set
cl <- list()
for (conj in c("fondo", "sig")) {
  a <- an[cpg %in% (if (conj == "fondo") ttM$cpg else T3$cpg)]
  grp <- a$UCSC_RefGene_Group
  g <- ifelse(gsub(";", "", grp) == "", "Intergenic",        # no gene -> intergenic
       ifelse(grepl("TSS|5'UTR|1stExon", grp), "Promoter",    # promoter takes precedence,
       ifelse(grepl("(^|;)Body(;|$)", grp), "Gene body",       # then gene body,
              "Other genic")))                                 # then 3'UTR and exon boundaries
  cl[[conj]] <- list(gen = factor(g, c("Promoter", "Gene body", "Other genic", "Intergenic")),
                     isla = factor(ifelse(a$Relation_to_Island == "", "OpenSea", a$Relation_to_Island),
                                   c("Island", "N_Shore", "S_Shore", "N_Shelf", "S_Shelf", "OpenSea")))
}
tiff(file.path(OUT, "HernandezFig2.tif"), width = 6, height = 3.4, units = "in",
     res = 1000, compression = "lzw")
par(mfrow = c(1, 2), mar = c(4, 5.4, 2.6, 1.2), oma = c(0, 0, 1.4, 0),   # left margin for category labels
    family = "sans", cex.axis = .76, tcl = -0.3, mgp = c(2.4, .6, 0))
porcentajes <- list()
for (panel in c("gen", "isla")) {                            # one horizontal bar panel each
  tit <- if (panel == "gen") "(a) Position relative to gene" else "(b) Relation to CpG island"
  m <- rbind(Found  = 100 * prop.table(table(cl$sig[[panel]])),   # significant probes, %
             Tested = 100 * prop.table(table(cl$fondo[[panel]])))  # background, %
  m <- m[, rev(seq_len(ncol(m))), drop = FALSE]              # top to bottom
  b <- barplot(m, beside = TRUE, horiz = TRUE, las = 1,      # horizontal: labels fit without rotation
               col = c("grey35", "grey75"), border = NA,
               xlim = c(0, max(m) * 1.30), xlab = "Probes (%)",
               main = tit, cex.main = .85, font.main = 1, cex.names = .75)
  text(m, b, sprintf("%.1f", m), pos = 4, cex = .62, xpd = NA, offset = .25)
  porcentajes[[panel]] <- m
}
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)   # canvas over the figure
plot(0, type = "n", axes = FALSE, xlab = "", ylab = "")      # for a single legend at the top
legend("top", c("Tested", "Found"), fill = c("grey75", "grey35"),
       border = NA, bty = "n", cex = .8, horiz = TRUE)
invisible(dev.off())
cat("Figura 2 escrita\n"); print(round(porcentajes$gen, 2)); print(round(porcentajes$isla, 2))
largas <- list()
for (panel in c("gen", "isla")) {                            # matrix to long format, both panels stacked
  largas[[panel]] <- melt(as.data.table(porcentajes[[panel]], keep.rownames = "set"),
                          variable.name = "category", value.name = "percent")[
                     , panel := if (panel == "gen") "Position relative to gene" else "Relation to CpG island"]
}
fwrite(rbind(largas$gen, largas$isla), file.path(OUT, "F2_distribuciones.csv"))

# ---- 5. Figure 3: genomic context of the locus (Gviz on hg19 TxDb) ----------------------------------
# Tracks: coordinates, gene model, the twelve probes by region, DMR, coefficients
eg  <- AnnotationDbi::select(org.Hs.eg.db, keys = "ZG16B",
                             keytype = "SYMBOL", columns = "ENTREZID")$ENTREZID[1]
gr_gen <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene, filter = list(gene_id = eg))
cat(sprintf("ZG16B en hg19: %s:%d-%d, hebra %s\n", as.character(seqnames(gr_gen)),
            start(gr_gen), end(gr_gen), as.character(strand(gr_gen))))
zg12 <- an[grepl("ZG16B", UCSC_RefGene_Name)][order(pos)]    # the twelve probes annotated to the gene
zg12[, grupo := ""]
for (i in seq_len(nrow(zg12)))                               # region, once per probe
  set(zg12, i, "grupo", paste(unique(strsplit(zg12$UCSC_RefGene_Group[i], ";")[[1]]), collapse = "/"))
zg12 <- merge(zg12[, .(cpg, pos, grupo)],
              S3[, .(cpg, M = M_principal)], by = "cpg", all.x = TRUE)
setorder(zg12, pos)
dmr <- T5[grepl("ZG16B", overlapping.genes)]                 # the region of the locus
ini <- min(zg12$pos) - 400; fin <- max(zg12$pos) + 400       # window to draw
pistas_sondas <- list()
for (nom in c("TSS1500", "TSS200", "Body / 3'UTR")) {        # one probe track per region
  gs <- if (nom == "Body / 3'UTR") c("Body", "3'UTR") else nom
  z <- zg12[grupo %in% gs]
  pistas_sondas[[nom]] <- AnnotationTrack(start = z$pos, width = 1, chromosome = "chr16", genome = "hg19",
                  name = nom, stacking = "dense", fill = "grey20", col = "grey20",
                  cex.title = .6, rotation.title = 0,        # horizontal title: rotated titles overlapped
                  background.title = "white", col.title = "black", fontcolor.title = "black")
}
zg_ana <- zg12[!is.na(M)]                                    # only the eight analyzed probes: drawing the
                                                             # filtered one as zero would make it look null
pistas <- list(
  GenomeAxisTrack(littleTicks = FALSE, cex = .6),            # 1. coordinates
  GeneRegionTrack(TxDb.Hsapiens.UCSC.hg19.knownGene, chromosome = "chr16",   # 2. gene model
                  start = ini, end = fin, name = "ZG16B", fill = "grey35",
                  col = "grey35", showId = FALSE, cex.title = .6, rotation.title = 0,
                  shape = "arrow",                           # exons as arrows: the plus strand is visible
                  background.title = "white", col.title = "black", fontcolor.title = "black"),
  pistas_sondas[["TSS1500"]],                                # 3. the four TSS1500 probes
  pistas_sondas[["TSS200"]],                                 # 4. the five TSS200 probes
  pistas_sondas[["Body / 3'UTR"]],                           # 5. the remaining three
  AnnotationTrack(start = dmr$start, end = dmr$end, chromosome = "chr16",    # 6. the DMR
                  genome = "hg19", name = "DMR", fill = "grey55", col = "grey55",
                  stacking = "dense", cex.title = .6, rotation.title = 0,
                  background.title = "white", col.title = "black", fontcolor.title = "black"),
  DataTrack(start = zg_ana$pos, width = 1, chromosome = "chr16",             # 7. the coefficients
            genome = "hg19", data = zg_ana$M, type = "h", lwd = 1.5,         # thin bars: two probes 18 bp apart
            col = "grey15", baseline = 0, col.baseline = "grey60",
            name = "Coefficient (M-value)", cex.title = .6, cex.axis = .55,
            showAxis = TRUE, col.axis = "black",            # otherwise Gviz uses the (white) title color
            yTicksAt = round(seq(min(zg_ana$M), max(zg_ana$M),
                                 length.out = 3), 2),        # numeric axis, so magnitude can be read
            rotation.title = 90,                             # tall track: the rotated title fits
            background.title = "white", col.title = "black", fontcolor.title = "black"))
tiff(file.path(OUT, "HernandezFig3.tif"), width = 6, height = 4.2, units = "in",
     res = 1000, compression = "lzw")
plotTracks(pistas, from = ini, to = fin, chromosome = "chr16",
           sizes = c(0.9, 0.8, 0.45, 0.45, 0.45, 0.4, 2.2),  # relative height of each track
           title.width = 1.6,
           fontfamily = "sans", fontfamily.title = "sans", margin = 6, innerMargin = 4)
invisible(dev.off())
cat("Figura 3 escrita\n")

# ---- 6. Figure 4: the locus in the three cohorts ----------------------------------------------------
paneles <- list(
  list(t = "(a) Discovery: buccal epithelium, continuous violence (n = 116)",
       d = merge(LOCUS, S3[, .(cpg, b = M_principal)], by = "cpg")),
  list(t = "(b) Replication: saliva, child maltreatment (n = 123)",
       d = merge(LOCUS, as.data.table(tr$cohorte2$locus)[, .(cpg, b = beta)], by = "cpg")),
  list(t = "(c) Third cohort: buccal epithelium, child maltreatment (n = 85)",
       d = merge(LOCUS, as.data.table(dc$locus)[, .(cpg, b = beta)], by = "cpg")))
pos8 <- an[cpg %in% LOCUS$cpg, .(cpg, pos)]                  # hg19 positions
tiff(file.path(OUT, "HernandezFig4.tif"), width = 6, height = 7.5, units = "in",
     res = 1000, compression = "lzw")
par(mfrow = c(3, 1), mar = c(4.2, 4.6, 2.6, 1), family = "sans",
    cex.axis = .95, tcl = -0.35, mgp = c(2.6, .7, 0))
for (p_ in paneles) {                                        # one panel per cohort
  dd <- merge(p_$d, pos8, by = "cpg")[order(pos)]
  dd[, kb := pos / 1000]
  n_ok <- sum(sign(dd$b) == dd$dir_pred)                     # concordance shown in the panel
  yl <- max(abs(dd$b)) * 1.5
  plot(dd$kb, dd$b, type = "n", ylim = c(-yl, yl),
       xlim = range(dd$pos)/1000 + c(-.05, .05),
       xlab = "Position on chromosome 16 (kb, hg19)",
       ylab = "Coefficient (M-value)", main = p_$t, cex.main = 1, font.main = 1)
  abline(h = 0, col = "grey60")
  abline(v = 2879.771, col = "grey75")                       # TSS1500/TSS200 boundary
  segments(dd$kb, 0, dd$kb, dd$b, col = "grey50")
  points(dd$kb, dd$b, pch = ifelse(dd$region == "TSS200", 19, 17), cex = 1.5)   # circles TSS200, triangles TSS1500
  legend("bottomleft", c("TSS1500", "TSS200"), pch = c(17, 19), bty = "n",
         cex = .9, horiz = TRUE)
  mtext(sprintf("%d of 8 probes in the predicted direction", n_ok),
        side = 3, line = -1.4, adj = 0.03, cex = .72)
}
invisible(dev.off())
cat("Figura 4 escrita\n")

# ---- 7. Tables of the article -----------------------------------------------------------------------
# Table 1: characteristics of the three cohorts
T1 <- data.table(
  Cohort    = c("Discovery (GSE226085)", "Replication (GSE239521)", "Third (GSE239520)"),
  Tissue    = c("Buccal epithelium", "Saliva", "Buccal epithelium"),
  Ancestry  = c("Syrian", "Japanese", "Japanese"),
  n         = c(116, 123, 85),
  Exposure  = c("War-related violence (0-19, continuous)",
                "Child maltreatment (binary)", "Child maltreatment (binary)"),
  Groups    = c("42 extended families", "61 exposed / 62 comparison", "36 exposed / 49 comparison"),
  `Sex (female/male)` = c("79/37", "47/76", "38/47"),
  `Age, years (M, SD)` = c("17.94 (14.80)", "13.80 (2.78)", "4.99 (2.33)"),
  `Median epithelial fraction` = c("0.815", "0.270", "0.862"))
fwrite(T1, file.path(OUT, "TABLA1_cohortes.csv")); print(T1)
# Supplementary Table S1: the differentially methylated probes
T2 <- T3[order(p), .(`Probe ID` = cpg, `Gene symbol` = ifelse(gene == "", "-", gene),
                     `Gene region` = ifelse(grupo == "", "-", grupo),
                     `Relation to island` = isla, Chr = sub("chr", "", chr), `Position hg19` = pos,
                     Dir = ifelse(M_est < 0, "Hypo", "Hyper"),
                     Coef = round(M_est, 3), `Δβ` = round(dBeta, 4),
                     p = signif(p, 3), `Adj. p` = signif(fdr, 3))]
fwrite(T2, file.path(OUT, "TABLA2_DMP.csv")); print(T2)
# Table 2: the regions
T3t <- T5[order(-no.cpgs), .(`Coord hg19` = sprintf("chr%s:%d-%d", sub("chr", "", seqnames),
                                                    start, end),
                             Width = width, `No. CpGs` = no.cpgs,
                             `Main gene assoc.` = overlapping.genes,
                             `Mean Δβ` = round(meandiff, 4),
                             `Min FDR` = signif(min_smoothed_fdr, 3))]
fwrite(T3t, file.path(OUT, "TABLA3_DMR.csv")); print(T3t)
# Table 3: the locus in the three cohorts, with both discovery FDRs side by side (primary model and
# within-family correlation model)
T4f <- fread(file.path(OUT, "T4_ZG16B_correlacion_familiar.csv"))
T4 <- merge(LOCUS, S3[, .(cpg, `Discovery coef` = M_principal, `Discovery p` = p_principal)],
            by = "cpg")
T4 <- merge(T4, T4f[, .(cpg, `FDR unadjusted` = fdr_primario,
                        `FDR family-adjusted` = fdr_familia)], by = "cpg")
T4 <- merge(T4, as.data.table(tr$cohorte2$locus)[, .(cpg, `Replication coef` = beta,
                                                     `Replication p` = p)], by = "cpg")
T4 <- merge(T4, as.data.table(dc$locus)[, .(cpg, `Third cohort coef` = beta,
                                            `Third cohort p` = p)], by = "cpg")
T4 <- merge(T4, an[, .(cpg, pos)], by = "cpg")[order(pos)]
T4 <- T4[, .(`Probe ID` = cpg, Region = region, `Position` = pos,
             `Sign` = ifelse(dir_pred < 0, "-", "+"),
             `Disc. coef` = round(`Discovery coef`, 4), `Disc. p` = signif(`Discovery p`, 3),
             `Disc. FDR` = sprintf("%.3f", `FDR unadjusted`),
             `Disc. FDR fam.` = sprintf("%.3f", `FDR family-adjusted`),
             `Repl. coef` = round(`Replication coef`, 4),
             `Repl. p` = signif(`Replication p`, 3),
             `Third coef` = round(`Third cohort coef`, 4),
             `Third p` = signif(`Third cohort p`, 3))]
fwrite(T4, file.path(OUT, "TABLA4_locus_tres_cohortes.csv")); print(T4)

# ---- 8. Session information and check of the files -----------------------------------------------
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_01c.txt"))
cat("\n\n=== FICHEROS PARA EL MANUSCRITO ===\n")
for (f in c("HernandezFig1.tif", "HernandezFig2.tif", "HernandezFig3.tif", "HernandezFig4.tif",
            "TABLA1_cohortes.csv", "TABLA2_DMP.csv", "TABLA3_DMR.csv",
            "TABLA4_locus_tres_cohortes.csv"))
  cat(sprintf("  %-32s %s\n", f, ifelse(file.exists(file.path(OUT, f)), "OK", "FALTA")))
