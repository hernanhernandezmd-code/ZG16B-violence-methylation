# 05 — Figure 4: the direction of the ZG16B association and age, in the three cohorts
#
# Draws the secondary analysis of age from the tables written by 03_exploratorio_edad.R, without
# refitting any model. Panel (a): for each of the eight promoter probes and each cohort, the
# exposure-by-age interaction expressed as a signed -log10(p), positive when the interaction moves the
# exposure effect toward the predicted sign (the discovery sign in the discovery cohort). Panel (b): the
# number of the eight probes that carry the predicted sign in each age group (discovery: under 10, 10 to
# 19 and 20 years or older; evaluation cohorts: halves split at the median age). The counts are checked
# against the values reported in the article (stopifnot).
# Paper section: Results (Secondary analysis of age); Discussion; Figure 4
# Input:  outputs/E1_interaccion_locus.csv, outputs/E1_interaccion_resumen.csv,
#         outputs/E1_edad_resumen.csv (written by 03), outputs/FigS2_una_por_familia.csv (probe positions,
#         written by 04)
# Output: outputs/HernandezFig_edad.tif (600 dpi), outputs/Fig4_edad.csv, outputs/sessionInfo_05.txt

# ---- 1. Packages and paths --------------------------------------------------------------------------
suppressMessages(library(data.table))
OUT <- "outputs"
COHORTES <- c("descubrimiento", "réplica (adolescentes, saliva)", "tercera (niños, bucal)")
ETIQUETAS <- c("Discovery (buccal, 1.5-64 y)", "Replication (saliva, 8.8-20 y)", "Third (buccal, 0.3-9.9 y)")

# ---- 2. Interaction by probe (panel a) --------------------------------------------------------------
il  <- fread(file.path(OUT, "E1_interaccion_locus.csv"))
ir  <- fread(file.path(OUT, "E1_interaccion_resumen.csv"))
pos <- unique(fread(file.path(OUT, "FigS2_una_por_familia.csv"))[, .(cpg, pos)])
il  <- merge(il, pos, by = "cpg")
il[, z := -log10(p_int) * ifelse(hacia_pred, 1, -1)]          # signed -log10(p): positive = toward the prediction
il[, orden := frank(pos, ties.method = "dense"), by = cohorte] # probes left to right by genomic position

# Control: the counts toward the prediction are those reported in Results (7, 8 and 8 of 8)
k <- il[, .(k = sum(hacia_pred)), by = cohorte][match(COHORTES, cohorte)]
print(k); stopifnot(identical(k$k, c(7L, 8L, 8L)))
stopifnot(identical(ir[match(COHORTES, cohorte), interaccion_hacia_pred], c("7 de 8", "8 de 8", "8 de 8")))

# ---- 3. Concordance by age group (panel b) ----------------------------------------------------------
er <- fread(file.path(OUT, "E1_edad_resumen.csv"))[grupo != "todos"]
er[, k := as.integer(sub(" de 8", "", concordancia))]
er[, etiqueta := c("< 10 y", "10-19 y", ">= 20 y", "< 14.1 y", ">= 14.1 y", "< 5.3 y", ">= 5.3 y")]
print(er[, .(cohorte, grupo, n, k)])
stopifnot(identical(er$k, c(8L, 8L, 8L, 0L, 8L, 1L, 0L)))    # values reported in Results
fwrite(rbind(il[, .(panel = "a", cohorte, cpg, region, pos, int, p_int, hacia_pred, z = round(z, 3))],
             er[, .(panel = "b", cohorte, cpg = grupo, region = NA, pos = NA, int = NA, p_int = NA,
                    hacia_pred = NA, z = k)], fill = TRUE),
       file.path(OUT, "Fig4_edad.csv"))

# ---- 4. Figure ---------------------------------------------------------------------------------------
# Same conventions as Fig. 3: circles TSS200, triangles TSS1500; filled symbols when p < .05
er[, etiqueta := c("<10", "10-19", "\u226520", "<14.1", "\u226514.1", "<5.3", "\u22655.3")]
for (dispositivo in "tiff") {                                # TIFF for the journal
  if (dispositivo == "tiff") tiff(file.path(OUT, "HernandezFig_edad.tif"), width = 7.5, height = 3.8, units = "in", res = 600, compression = "lzw")
  layout(matrix(1:2, 1), widths = c(1.2, 1))
  par(mar = c(4.8, 4.2, 2, 0.8), family = "sans", cex.axis = .75, cex.lab = .85, tcl = -0.3, mgp = c(2.3, .55, 0))

  # (a) signed -log10(p) of the interaction, one group of eight probes per cohort
  x  <- il[, (match(cohorte, COHORTES) - 1) * 9.5 + orden]
  yl <- c(-1.6, max(2.2, max(il$z) * 1.15))
  plot(x, il$z, type = "n", ylim = yl, xaxt = "n", xlab = "", xlim = c(0.3, 27.7),
       ylab = expression("Interaction with age, signed " * -log[10](italic(p))))
  abline(h = 0, col = "grey60"); abline(h = c(-1, 1) * -log10(.05), lty = 2, col = "grey75")
  segments(x, 0, x, il$z, col = "grey55")
  points(x, il$z, pch = ifelse(il$region == "TSS200", ifelse(il$p_int < .05, 19, 1), ifelse(il$p_int < .05, 17, 2)), cex = 1)
  centros <- (0:2) * 9.5 + 4.5
  mtext(c("Discovery", "Replication", "Third"), side = 1, at = centros, line = .6, cex = .62, font = 2)
  mtext(c("buccal, 1.5-64 y", "saliva, 8.8-20 y", "buccal, 0.3-9.9 y"), side = 1, at = centros, line = 1.5, cex = .58)
  text(centros, yl[2] * .97, paste0(k$k, " of 8"), cex = .7)
  legend("bottomright", c("TSS1500", "TSS200", "p < .05 (filled)"), pch = c(2, 1, 19), bg = "white", box.col = NA, cex = .6, inset = .02)
  mtext("(a)", side = 3, adj = -0.12, line = .5, font = 2, cex = .85)

  # (b) number of probes with the predicted sign by age group (dashed line: chance, 4 of 8)
  par(mar = c(4.8, 3.6, 2, 0.6))
  xb <- barplot(er$k, ylim = c(0, 9.6), space = c(0, 0, 0, .8, 0, .8, 0), col = "grey80", border = "grey35",
                ylab = "Probes with the predicted sign (of 8)", names.arg = rep("", nrow(er)),
                yaxt = "n")
  axis(2, at = c(0, 2, 4, 6, 8), las = 1)
  mtext(er$etiqueta, side = 1, at = xb, line = .4, cex = .5)   # one label per bar, none dropped
  abline(h = 4, lty = 2, col = "grey55")
  text(xb, er$k + .4, er$k, cex = .7)
  mtext(c("Discovery", "Replication", "Third"), side = 1, at = c(mean(xb[1:3]), mean(xb[4:5]), mean(xb[6:7])),
        line = 1.4, cex = .62, font = 2)
  mtext("Age (years)", side = 1, line = 2.5, cex = .6)
  mtext("(b)", side = 3, adj = -0.14, line = .5, font = 2, cex = .85)
  invisible(dev.off())
}
cat("\nFigura 4 escrita\n")
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_05.txt"))
