# =============================================================================
#  Figure 3, step 1: differential expression of RPH3A across diseases
#    Fig 3a  RPH3A on per-cohort UMAPs of RORB IT neurons
#    Fig 3b  RPH3A per donor (box plots) with edgeR results
#    Fig S1  the same UMAPs coloured by the original RORB IT subtype labels of each study
#    Fig S2  box plots and UMAPs for the six shared network genes
# =============================================================================
suppressPackageStartupMessages({ library(Matrix); library(data.table); library(ggplot2); library(patchwork); library(edgeR) })
set.seed(1234)
ARGS    <- commandArgs(trailingOnly = TRUE); DO_UMAP <- !("box" %in% ARGS)
PROJECT <- Sys.getenv("RPH3A_PROJECT", unset = if (length(ARGS) && dir.exists(ARGS[1])) ARGS[1] else getwd())
FIG3    <- file.path(PROJECT, "results", "Fig3"); CACHE <- file.path(FIG3, "cache")
OUT     <- file.path(FIG3, "Results"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(CACHE)) stop("Cache not found: ", CACHE, " (run step 0 first, or pass the project folder)")
GENE_MAIN  <- "RPH3A"
GENES_SUPP <- c("SOD1", "CSNK1A1", "PAFAH1B1", "MAPRE1", "CKAP5", "TAOK1")   # shared network genes (Fig S2)
GENES <- unique(c(GENE_MAIN, GENES_SUPP)); FLAG <- 10
N_HVG <- 2000; N_PCS <- 30; MAX_DRAW <- 3000
DS <- c(AD = "AD", sFTLD = "sFTLD", C9FTLD = "C9FTLD", sALS = "sALS", C9ALS = "C9ALS", PD_temporal = "PD temporal", PD_frontal = "PD frontal")
RL <- c("AD", "sFTLD", "C9FTLD", "C9ALS", "PD temporal")
dis_col <- c(AD = "#7B3294", sFTLD = "#D7301F", C9FTLD = "#FC8D59", sALS = "#4575B4", C9ALS = "#91BFDB", `PD temporal` = "#1B7837", `PD frontal` = "#7FBF7B")
dis_title <- c(AD = "AD\nMathys 2023", sFTLD = "sFTLD\nPineda 2024", C9FTLD = "C9FTLD\nPineda 2024", sALS = "sALS\nPineda 2024",
               C9ALS = "C9ALS\nPineda 2024", `PD temporal` = "PD temporal\nASAP 2025", `PD frontal` = "PD frontal\nASAP 2025")
say <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
save2 <- function(p, name, w, h) { ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white"); say("Saved ", name) }

# ── edgeR per dataset ─────────────────────────────────────────────────────────
run_edger <- function(M, info, min_nuclei = 1) {
  k <- info$n_nuclei >= min_nuclei; M <- M[, k]; info <- info[k]
  grp <- factor(ifelse(info$group == "Control", "Control", "Disease"), levels = c("Control", "Disease"))
  team <- if ("team" %in% names(info) && uniqueN(info$team) > 1) factor(info$team) else NULL
  design <- if (is.null(team)) model.matrix(~ grp) else model.matrix(~ team + grp)
  y0 <- DGEList(M); keep <- filterByExpr(y0, group = grp)
  y <- (if (exists("normLibSizes")) normLibSizes else calcNormFactors)(y0[keep, , keep.lib.sizes = FALSE]); y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE); tt <- topTags(glmQLFTest(fit, coef = "grpDisease"), n = Inf, sort.by = "none")$table
  yf <- y0; yf$samples$norm.factors <- 1; yf$samples$lib.size <- y$samples$lib.size * y$samples$norm.factors   # TMM-effective library sizes
  lc <- cpm(yf, log = TRUE, prior.count = 1)
  if (!is.null(team)) lc <- limma::removeBatchEffect(lc, batch = team, design = model.matrix(~ grp))
  list(de = data.table(gene = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM, F = tt$F, P = tt$PValue, FDR = tt$FDR),
       logcpm = lc, info = info[, grp := grp], n_donors = ncol(M), team_adjusted = !is.null(team))
}
RES <- list(); DE <- list(); VAL <- list()
for (d in names(DS)) {
  f <- file.path(CACHE, paste0(d, "_pseudobulk.rds")); if (!file.exists(f)) { message("Missing cache: ", f, " (run step 0)"); next }
  x <- readRDS(f); rownames(x$counts) <- toupper(rownames(x$counts)); x$counts <- x$counts[!duplicated(rownames(x$counts)), ]
  say(DS[[d]], ": edgeR on ", ncol(x$counts), " donors")
  r <- run_edger(x$counts, copy(x$donors)); r10 <- run_edger(x$counts, copy(x$donors), FLAG)
  RES[[d]] <- r
  DE[[d]] <- merge(r$de[, dataset := DS[[d]]], r10$de[, .(gene, logFC_min10 = logFC, FDR_min10 = FDR)], by = "gene", all.x = TRUE)
  g <- intersect(GENES, rownames(r$logcpm))
  VAL[[d]] <- data.table(dataset = DS[[d]], donor = rep(colnames(r$logcpm), each = length(g)), gene = g,
                         log2cpm = as.vector(r$logcpm[g, , drop = FALSE]))[r$info[, .(donor, grp, n_nuclei)], on = "donor"]
}
DE <- rbindlist(DE); VAL <- rbindlist(VAL)
fwrite(DE[order(dataset, P)], file.path(OUT, "Fig3_edgeR_all_genes.csv.gz"))
fwrite(DE[gene %in% GENES][order(gene, dataset)], file.path(OUT, "Fig3_edgeR_RPH3A_and_screen_hits.csv"))   # read by the forest plot and Table S
fwrite(VAL, file.path(OUT, "Fig3_donor_values.csv"))
print(DE[gene == GENE_MAIN, .(dataset, logFC = round(logFC, 2), FDR = signif(FDR, 2), logFC_min10 = round(logFC_min10, 2), FDR_min10 = signif(FDR_min10, 2))])
print(VAL[gene == GENE_MAIN, .(donors = .N, under_10 = sum(n_nuclei < FLAG)), by = .(dataset, grp)])

# ── Box plots ─────────────────────────────────────────────────────────────────
gf <- file.path(CACHE, "AD_grouping.txt")
AD_NOTE <- if (file.exists(gf) && readLines(gf)[1] == "pathology") "AD = pathologic AD vs no pathologic AD (public ROSMAP diagnosis)" else "AD = Mathys group 4 vs group 1"
star <- function(q) fcase(is.na(q), "", q < 0.001, "***", q < 0.01, "**", q < 0.05, "*", default = "n.s.")
LAB <- DE[gene %in% GENES, .(dataset, gene, label = sprintf("log2FC %+.2f\nFDR %s %s", logFC, formatC(FDR, format = "g", digits = 2), star(FDR)))]
box_plot <- function(genes) {
  d <- VAL[gene %in% genes][, `:=`(gene = factor(gene, levels = genes), dataset = factor(dataset, levels = DS),
                                   fillk = ifelse(grp == "Control", "Control", as.character(dataset)), low = n_nuclei < FLAG)]
  # robust axis range: 0.5th-99.5th percentile of donors with >= 10 nuclei; points outside are drawn at the edge
  d[, `:=`(lo = quantile(log2cpm[!low], 0.005), hi = quantile(log2cpm[!low], 0.995)), by = gene]
  d[, `:=`(y_pt = pmin(pmax(log2cpm, lo), hi), pt = factor(fcase(log2cpm < lo | log2cpm > hi, "out", low, "low", default = "ok"), levels = c("ok", "low", "out")))]
  l <- LAB[gene %in% genes][, `:=`(gene = factor(gene, levels = genes), dataset = factor(dataset, levels = DS))]
  ggplot(d, aes(grp, log2cpm)) +
    geom_boxplot(aes(fill = fillk), outliers = FALSE, width = 0.62, alpha = 0.75, linewidth = 0.35) +
    geom_point(aes(y = y_pt, size = n_nuclei, shape = pt), position = position_jitter(width = 0.14, height = 0, seed = 1), alpha = 0.65, stroke = 0.4) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 2.4, fill = "black", colour = "white") +
    geom_text(data = l, aes(x = 1.5, y = Inf, label = label), vjust = 1.15, size = 2.5, lineheight = 0.9, inherit.aes = FALSE) +
    facet_grid(gene ~ dataset, scales = "free_y", labeller = labeller(dataset = dis_title)) +
    scale_fill_manual(values = c(Control = "grey80", dis_col), guide = "none") +
    scale_shape_manual(values = c(ok = 16, low = 1, out = 6), drop = FALSE, name = NULL,
                       labels = c(ok = ">= 10 nuclei", low = "< 10 nuclei", out = "off-axis (drawn at edge)")) +
    scale_size_area(max_size = 2.6, breaks = c(10, 100, 1000), name = "Nuclei per donor") +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3))) +
    labs(x = NULL, y = "Pseudobulk log2 CPM (TMM), one point per donor",
         caption = paste("RORB IT neurons (all RORB IT subtypes per study). All donors with >= 1 nucleus; diamond = mean of donors.",
                         "\nLabel = edgeR QL on these donors (BH across genes within dataset).", AD_NOTE, "; FTLD/ALS vs shared PN controls;",
                         "PD frontal values shown with the team effect removed (team in the model).")) +
    theme_bw(base_size = 10) +
    theme(strip.text.y = element_text(face = "bold.italic", angle = 0), strip.text.x = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey96", colour = NA), panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
          legend.position = "bottom", plot.caption = element_text(hjust = 0, colour = "grey30", size = 7.5))
}
nd <- uniqueN(VAL$dataset)
save2(box_plot(GENE_MAIN), "Fig3b_RPH3A_boxplots", 1.75 * nd + 1.2, 4.4)
for (i in split(GENES_SUPP, ceiling(seq_along(GENES_SUPP) / 4)))
  save2(box_plot(i), paste0("FigS2_boxplots_", paste(i, collapse = "_")), 1.75 * nd + 1.2, 2 * length(i) + 1.6)
if (!DO_UMAP) quit(save = "no")

# ── UMAPs (one embedding per cohort, cached) ──────────────────────────────────
suppressPackageStartupMessages({ library(SingleCellExperiment); library(scuttle); library(scran); library(scater) })
rd <- function(d) { f <- file.path(CACHE, paste0(d, "_cells.rds")); if (file.exists(f)) readRDS(f) else NULL }
cohorts <- list(AD = list(rd("AD")), `FTLD/ALS` = NULL, `PD temporal` = list(rd("PD_temporal")), `PD frontal` = list(rd("PD_frontal")))
sf <- rd("sFTLD"); pin <- list(sf)
for (d in c("C9FTLD", "sALS", "C9ALS")) { x <- rd(d); if (is.null(x)) next; k <- x$group != "Control"
  pin[[d]] <- list(counts = x$counts[, k], donor = x$donor[k], group = x$group[k], subtype = x$subtype[k]) }
cohorts[["FTLD/ALS"]] <- pin
merge_parts <- function(parts) { parts <- Filter(Negate(is.null), parts); g <- Reduce(intersect, lapply(parts, function(p) rownames(p$counts)))
  list(counts = do.call(cbind, lapply(parts, function(p) p$counts[g, ])), donor = unlist(lapply(parts, `[[`, "donor")),
       group = unlist(lapply(parts, `[[`, "group")), subtype = unlist(lapply(parts, `[[`, "subtype")),
       team = unlist(lapply(parts, function(p) if (is.null(p$team)) rep(NA, length(p$donor)) else p$team))) }
embed <- function(x, tag) {
  f <- file.path(CACHE, paste0("UMAP_", gsub("[^A-Za-z]", "", tag), ".rds")); if (file.exists(f)) return(readRDS(f))
  rownames(x$counts) <- toupper(rownames(x$counts)); x$counts <- x$counts[!duplicated(rownames(x$counts)), ]
  sce <- SingleCellExperiment(list(counts = x$counts)); sce <- sce[rowSums(counts(sce) > 0) >= 10, ]
  sce$donor <- x$donor; if (!all(is.na(x$team))) sce$team <- x$team
  sce <- logNormCounts(sce); say("embedding ", tag, ": ", ncol(sce), " nuclei")
  sce <- runPCA(sce, subset_row = getTopHVGs(modelGeneVar(sce), n = N_HVG), ncomponents = N_PCS); dr <- "PCA"
  if (requireNamespace("harmony", quietly = TRUE)) { vars <- intersect(c("donor", "team"), names(colData(sce)))
    h <- tryCatch(harmony::RunHarmony(reducedDim(sce, "PCA"), meta_data = as.data.frame(colData(sce)[, vars, drop = FALSE]), vars_use = vars, verbose = FALSE),
                  error = function(e) NULL); if (!is.null(h)) { reducedDim(sce, "HARMONY") <- h; dr <- "HARMONY" } }
  sce <- runUMAP(sce, dimred = dr, n_neighbors = 30); um <- reducedDim(sce, "UMAP"); g <- intersect(GENES, rownames(sce))
  out <- data.table(UMAP1 = um[, 1], UMAP2 = um[, 2], donor = x$donor, group = x$group, subtype = x$subtype, as.matrix(t(logcounts(sce)[g, , drop = FALSE])))
  saveRDS(out, f); out
}
U <- lapply(names(cohorts), function(cn) embed(merge_parts(cohorts[[cn]]), cn)); names(U) <- names(cohorts)
lv_all <- c("Control", "AD", "sFTLD", "C9FTLD", "sALS", "C9ALS", "PD")
umap_plot <- function(gene, cn) {
  d <- copy(U[[cn]]); if (!gene %in% names(d)) return(NULL); d[, expr := get(gene)]
  lv <- intersect(lv_all, unique(d$group)); d[, group := factor(group, levels = lv)]
  st <- d[, .(n = .N, nd = uniqueN(donor), pct = 100 * mean(expr > 0)), by = group][order(group)]
  st[, lab := sprintf("%s\n%d donors, %.0f%% +", group, nd, pct)]; d[, facet := factor(st$lab[match(group, st$group)], levels = st$lab)]
  n <- min(min(st$n), MAX_DRAW); d <- d[, .SD[sample(.N, n)], by = facet][order(expr)]
  cap <- quantile(d$expr[d$expr > 0], 0.99, na.rm = TRUE)
  ggplot(d, aes(UMAP1, UMAP2, colour = pmin(expr, cap))) + geom_point(size = 0.35, stroke = 0) + facet_wrap(~ facet, nrow = 1) + coord_equal() +
    scale_colour_gradient(low = "grey88", high = "#08306B", name = paste0(gene, "\nlog-norm")) +
    labs(title = cn, subtitle = sprintf("%s nuclei drawn per group", format(n, big.mark = ","))) + theme_void(base_size = 9) +
    theme(strip.text = element_text(size = 7.5), plot.title = element_text(face = "bold"), plot.subtitle = element_text(colour = "grey40", size = 7))
}
subtype_plot <- function(cn) { d <- U[[cn]][sample(.N, min(.N, 20000))]
  ggplot(d, aes(UMAP1, UMAP2, colour = subtype)) + geom_point(size = 0.3, stroke = 0) + coord_equal() + labs(title = cn, colour = NULL) +
    theme_void(base_size = 8) + theme(plot.title = element_text(face = "bold"), legend.text = element_text(size = 6)) + guides(colour = guide_legend(override.aes = list(size = 2))) }
save2(wrap_plots(lapply(names(U), subtype_plot), nrow = 1), "FigS1_UMAP_RORB_IT_subtypes", 18, 3.6)   # Fig S1
for (g in GENES) {
  ps <- Filter(Negate(is.null), lapply(names(U), function(cn) umap_plot(g, cn)))
  w <- sapply(names(U), function(cn) uniqueN(U[[cn]]$group))
  ur <- wrap_plots(ps, nrow = 1, widths = w[seq_along(ps)])
  tag <- if (g == GENE_MAIN) "Fig3a_" else "FigS2_"
  save2(ur, paste0(tag, g, "_UMAPs"), 20, 3.6)
  panel <- ur / box_plot(g) + plot_layout(heights = c(1, 1.2)) +
    plot_annotation(title = sprintf("%s in RORB IT neurons across diseases", g), theme = theme(plot.title = element_text(face = "bold", size = 13)))
  save2(panel, paste0(tag, g, "_panel"), 20, 8.8)
}
say("Done: ", OUT)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_Fig3_step1.txt"))
