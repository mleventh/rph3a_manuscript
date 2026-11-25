# =============================================================================
#  Figure 4 panels 
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
ARGS    <- commandArgs(trailingOnly = TRUE)
PROJECT <- Sys.getenv("RPH3A_PROJECT", unset = if (length(ARGS) && dir.exists(ARGS[1])) ARGS[1] else getwd())
RES     <- file.path(PROJECT, "results")
OUT     <- file.path(RES, "Fig4"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RL <- c("AD", "sFTLD", "C9FTLD", "C9ALS", "PD temporal"); CT <- c("sALS", "PD frontal"); DIS <- c(RL, CT)
FOCUS <- c("NRXN", "PTPR", "NRG"); MIN_CHANGE <- 0.02
C_LOSS <- "#B4462F"; C_CTRL <- "#6E7F99"; C_GREY <- "grey78"
grp_col <- c(`RPH3A lost` = C_LOSS, `RPH3A not lost` = C_CTRL)
say <- function(...) message(format(Sys.time(), "%H:%M:%S"), "  ", ...)
th <- theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold", size = 11), plot.title.position = "plot",
        plot.subtitle = element_text(colour = "grey35", size = 8.5), axis.text = element_text(colour = "black"),
        strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "none")
save2 <- function(p, name, w, h) {
  tryCatch(ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h, device = cairo_pdf),
           error = function(e) ggsave(file.path(OUT, paste0(name, ".pdf")), p, width = w, height = h))
  ggsave(file.path(OUT, paste0(name, ".png")), p, width = w, height = h, dpi = 300, bg = "white"); say("Saved ", name) }
# y positions with a gap between the RPH3A-loss diseases and the two contrasts
ypos <- setNames(c(rev(seq_along(RL)) + length(CT) + 0.6, rev(seq_along(CT))), DIS)
y_axis <- scale_y_continuous(breaks = ypos, labels = names(ypos), expand = expansion(add = c(0.6, 1.1)))
grp_of <- function(d) factor(ifelse(d %in% RL, "RPH3A lost", "RPH3A not lost"), levels = names(grp_col))

lr <- fread(file.path(RES, "RPH3A_outgoing_LR_pairs.csv"))[disease %in% DIS]

# ── 4a. Total outgoing signaling from the primary RPH3A-loss neurons ─────────
A <- lr[, .(share_ctrl = sum(share_ctrl), share_dis = sum(share_dis), change = sum(change)), by = disease]
A[, `:=`(rel = 100 * change / share_ctrl, y = ypos[disease], grp = grp_of(disease))]
fwrite(A, file.path(OUT, "Fig4a_values.csv"))
lim <- max(abs(A$rel)) * 1.5
pa <- ggplot(A, aes(rel, y, fill = grp)) + geom_col(orientation = "y", width = 0.72) + geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = sprintf("%+.1f%%", rel), hjust = ifelse(rel < 0, 1.15, -0.15)), size = 3) +
  annotate("text", x = -lim, y = max(ypos) + 0.85, label = "RPH3A lost", hjust = 0, colour = C_LOSS, fontface = "bold", size = 3.2) +
  annotate("text", x = -lim, y = max(ypos[CT]) + 0.85, label = "RPH3A not lost", hjust = 0, colour = C_CTRL, fontface = "bold", size = 3.2) +
  scale_fill_manual(values = grp_col) + y_axis + scale_x_continuous(limits = c(-lim, lim)) +
  labs(title = "a  RPH3A-loss neurons send less", x = "Change in outgoing signaling share\n(% of control share)", y = NULL) + th
save2(pa, "Fig4a_total_outgoing", 3.8, 3.2)

# ── 4b. Pathway landscape ─────────────────────────────────────────────────────
pw <- lr[, .(share_ctrl = sum(share_ctrl), change = sum(change)), by = .(pathway = pathway_name, disease)]
pw[, rel := fifelse(share_ctrl > 0, 100 * change / share_ctrl, NA_real_)]
B <- pw[, .(n_lost = sum(disease %in% RL & change <= -MIN_CHANGE), rel_RL = mean(rel[disease %in% RL], na.rm = TRUE),
            rel_CT = mean(rel[disease %in% CT], na.rm = TRUE), base = mean(share_ctrl[disease %in% RL])), by = pathway][base >= 0.05]
B[, focus := pathway %in% FOCUS]
set.seed(1); B[, yj := n_lost + runif(.N, -0.18, 0.18)]
fwrite(B[order(rel_RL)], file.path(OUT, "Fig4b_pathway_landscape.csv"))
xl <- quantile(abs(B$rel_RL), 0.99, na.rm = TRUE) * 1.15
lab_b <- B[focus == TRUE][order(rel_RL)][, `:=`(lx = -xl * 0.55, ly = 5.35 - 0.38 * (seq_len(.N) - 1))]
pb <- ggplot(B, aes(rel_RL, yj)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 4.6, ymax = 5.45, fill = C_LOSS, alpha = 0.07) +
  geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_point(data = B[focus == FALSE], aes(size = base), colour = C_GREY) +
  geom_segment(data = lab_b, aes(x = lx, y = ly, xend = rel_RL, yend = yj), colour = C_LOSS, linewidth = 0.3) +
  geom_point(data = B[focus == TRUE], aes(size = base), fill = C_LOSS, colour = "black", shape = 21, stroke = 0.5) +
  geom_text(data = lab_b, aes(x = lx, y = ly, label = pathway), hjust = 1.1, colour = C_LOSS, fontface = "bold", size = 3.2) +
  scale_size_area(max_size = 6) + scale_y_continuous(breaks = 0:5, limits = c(-0.5, 5.6)) + scale_x_continuous(limits = c(-xl, xl)) +
  annotate("text", x = xl, y = -0.45, label = sprintf("%d pathways; dot size = baseline share", nrow(B)), hjust = 1, size = 2.4, colour = "grey45") +
  labs(title = "b  Only NRXN, PTPR and NRG fall in all five", x = "Mean change across the 5 RPH3A-loss diseases (% of control share)",
       y = "RPH3A-loss diseases with a loss\n(share down >= 0.02 points)") + th
save2(pb, "Fig4b_pathway_landscape", 4.6, 3.4)

# ── 4c. The three pathways, disease by disease ────────────────────────────────
C3 <- pw[pathway %in% FOCUS][, `:=`(y = ypos[disease], grp = grp_of(disease))]
bl <- B[pathway %in% FOCUS, .(pathway, base)]
C3 <- merge(C3, bl, by = "pathway")[, panel := factor(sprintf("%s  (control share %.1f%%)", pathway, base),
                                                       levels = sprintf("%s  (control share %.1f%%)", FOCUS, bl[match(FOCUS, pathway), base]))]
pc <- ggplot(C3, aes(rel, y, colour = grp)) + geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_segment(aes(x = 0, xend = rel, yend = y), linewidth = 1) + geom_point(size = 2.4) +
  facet_wrap(~ panel, nrow = 1, scales = "free_x") + scale_colour_manual(values = grp_col) + y_axis +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(title = "c  NRXN, PTPR and NRG, disease by disease", x = "Change (% of control share)", y = NULL) + th
save2(pc, "Fig4c_focus_pathways", 7.4, 2.9)
fwrite(C3[order(pathway, -y), .(pathway, disease, group = grp, share_ctrl, change_points = change, change_pct_of_control = rel)],
       file.path(OUT, "Fig4c_values.csv"))

# ── 4d. Ligand -> receptor pairs carrying the loss ────────────────────────────
D <- lr[pathway_name %in% FOCUS, .(change = sum(change)), by = .(interaction_name, pathway_name, disease)]
top <- D[disease %in% RL, .(m = mean(change)), by = interaction_name][order(m)][1:10, interaction_name]
D <- D[interaction_name %in% top]
# CellChat names one NRXN3 isoform by its UniProt ID (Q9HDB5); confirm in CellChatDB.human$interaction before publishing
D[, lab := sprintf("%s  (%s)", gsub("_", " → ", sub("^Q9HDB5", "NRXN3B", interaction_name)), pathway_name)]
D[, lab := factor(lab, levels = rev(unique(lab[match(top, interaction_name)])))]
Dm <- D[, .(RL = mean(change[disease %in% RL]), CT = mean(change[disease %in% CT])), by = lab]
fwrite(dcast(D, interaction_name + pathway_name ~ disease, value.var = "change"), file.path(OUT, "Fig4d_values.csv"))
pd <- ggplot() + geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_point(data = D[disease %in% RL], aes(change, lab), colour = C_LOSS, alpha = 0.35, size = 1, position = position_jitter(height = 0.12, seed = 1)) +
  geom_point(data = Dm, aes(RL, lab, shape = "RPH3A lost (mean of 5)"), colour = C_LOSS, size = 2.8) +
  geom_point(data = Dm, aes(CT, lab, shape = "sALS + PD frontal (mean)"), colour = C_CTRL, size = 2.8, stroke = 1.1) +
  scale_shape_manual(values = c(`RPH3A lost (mean of 5)` = 16, `sALS + PD frontal (mean)` = 1), name = NULL) +
  labs(title = "d  Carried by neurexin → LRRTM4 / neuroligin-1", x = "Change in share of total signaling (percentage points)", y = NULL) +
  th + theme(legend.position = c(0.02, 0.06), legend.justification = c(0, 0), legend.text = element_text(size = 7.5))
save2(pd, "Fig4d_LR_pairs", 5, 3.6)

# ── 4e. Receivers (step 5 output) ─────────────────────────────────────────────
f6 <- file.path(RES, "RPH3A_R6_receivers.csv")
if (file.exists(f6)) {
  r6 <- fread(f6)[disease %in% DIS, .(change = sum(change)), by = .(receiver, pathway, disease)]
  E <- r6[, .(change = mean(change)), by = .(receiver, pathway, grp = ifelse(disease %in% RL, "RPH3A lost", "not lost"))]
  ord <- E[grp == "RPH3A lost", .(s = sum(change)), by = receiver][order(-s), receiver]
  keep <- E[, .(m = max(abs(change))), by = receiver][m >= 0.005, receiver]
  E <- E[receiver %in% keep][, `:=`(receiver = factor(receiver, levels = ord), pathway = factor(pathway, levels = FOCUS),
                                    grp = factor(grp, levels = c("RPH3A lost", "not lost")),
                                    lab = ifelse(abs(change) >= 0.01, sprintf("%.2f", change), ""))]
  levels(E$grp) <- c("RPH3A lost", "Not lost\n(sALS, PD frontal)")
  fwrite(E, file.path(OUT, "Fig4e_values.csv"))
  v <- max(abs(E$change))
  pe <- ggplot(E, aes(grp, receiver, fill = change)) + geom_tile(colour = "white", linewidth = 0.6) + geom_text(aes(label = lab), size = 2.4) +
    facet_grid(. ~ pathway) + scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-v, v), name = "Mean change\n(% points)") +
    labs(title = "e  Excitatory → excitatory signaling is what is lost", subtitle = "Signals sent by RPH3A-loss neurons, by receiving class",
         x = NULL, y = NULL) + th + theme(legend.position = "right", axis.text.x = element_text(size = 7.5), axis.line = element_blank(),
                                         axis.ticks = element_blank())
  save2(pe, "Fig4e_receivers", 5.6, 0.3 * length(keep) + 1.9)
} else message("RPH3A_R6_receivers.csv not found: 4e skipped (run the LR-and-receivers step first)")

# ── S3. Every sender class (step 3 output) ───────────────────────────────────
# Two versions: percentage points (same units as R1) and % of control share (relative).
# Grey "absent" = class not present in that dataset's CellChat run (AD: no microglia or vascular
# cells were loaded; "Other Exc" exists only in the AD label mapping).
fs <- file.path(RES, "RPH3A_sender_share.csv")
if (file.exists(fs)) {
  classes <- c("RORB IT (target)", "Deep-layer Exc", "L2/3 IT", "Inh SST", "Other Exc", "Inh PV", "Inh LAMP5", "Inh VIP/CGE",
               "Astrocytes", "Microglia", "Oligodendrocytes", "OPC", "Vascular")
  S <- fread(fs)[disease %in% DIS]
  S <- merge(CJ(class = classes, disease = DIS), S, by = c("class", "disease"), all.x = TRUE)
  S[, absent := is.na(share_ctrl) | (share_ctrl == 0 & share_dis %in% c(0, NA))]
  S[, rel := fifelse(!absent & share_ctrl > 0, 100 * change / share_ctrl, NA_real_)]
  S[, `:=`(class = factor(class, levels = rev(classes)), disease = factor(disease, levels = DIS))]
  s4 <- function(var, unit, fmt, name) {
    S[, val := get(var)]; v <- quantile(abs(S$val), 0.97, na.rm = TRUE)
    ggplot(S, aes(disease, class)) + geom_tile(aes(fill = pmax(pmin(val, v), -v)), colour = "white") +
      geom_tile(data = S[absent == TRUE], fill = "grey92", colour = "white") +
      geom_text(aes(label = ifelse(absent, "absent", sprintf(fmt, val))), size = 2.3, colour = ifelse(S$absent, "grey55", "black")) +
      geom_vline(xintercept = length(RL) + 0.5, linewidth = 0.6) +
      geom_hline(yintercept = c(length(classes) - 1.5, length(classes) - 3.5), linetype = "dotted") +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-v, v), name = unit, na.value = "grey92") +
      labs(title = "S3  Outgoing signaling by sender class",
           subtitle = "Above the dotted lines: primary, then secondary RPH3A-loss classes",
           x = NULL, y = NULL) + th + theme(legend.position = "right", axis.text.x = element_text(angle = 35, hjust = 1), axis.line = element_blank(),
                                           axis.ticks = element_blank())
  }
  save2(s4("change", "Change\n(% points)", "%+.2f"), "FigS3_all_senders_points", 6, 4.6)
  save2(s4("rel", "Change\n(% of control)", "%+.0f"), "FigS3_all_senders_relative", 6, 4.6)
  fwrite(S[, .(class, disease, share_ctrl, share_dis, change, rel, absent)], file.path(OUT, "FigS3_values.csv"))
} else message("RPH3A_sender_share.csv not found: S3 skipped (run the outgoing-signaling step first)")

# ── Assembled Figure 4 (a-d, plus e if available) ─────────────────────────────
# row 1: a, b, d   row 2: c (three pathways) and e (receivers) if available
row2 <- if (exists("pe")) (pc | pe) + plot_layout(widths = c(1.5, 1)) else (pc | plot_spacer()) + plot_layout(widths = c(1.5, 1))
fig <- (pa | pb | pd) / row2 + plot_layout(heights = c(1, 1.15))
save2(fig, "Fig4_assembled", 16, 8.5)
say("Done: ", OUT)
writeLines(capture.output(sessionInfo()), file.path(OUT, "sessionInfo_Fig4.txt"))
