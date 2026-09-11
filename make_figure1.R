# Figure 1 - Study flow diagram (base R)
FIG_DIR <- "/Users/kuan-yuchu/Desktop/LOAD/figures"
dir.create(FIG_DIR, showWarnings = FALSE)

tiff(file.path(FIG_DIR, "Figure1_study_flow.tiff"),
     width = 8, height = 7, units = "in", res = 300,
     compression = "lzw")
op <- par(mar = c(0.5, 0.5, 0.5, 0.5))
plot(NA, xlim = c(0, 100), ylim = c(0, 100), axes = FALSE,
     xlab = "", ylab = "", asp = NA)

# helper: rounded-style rectangle with label
node <- function(x, y, w, h, text, fill = "#e3f2fd", border = "#1976d2",
                 text_col = "black", lwd = 1.5, cex = 0.95) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2,
       col = fill, border = border, lwd = lwd)
  text(x, y, text, cex = cex, col = text_col, font = 1)
}

# helper: arrow between nodes
arrow_d <- function(x1, y1, x2, y2, lbl = NULL) {
  arrows(x1, y1, x2, y2, length = 0.12, angle = 25, lwd = 1.5, col = "#333")
  if (!is.null(lbl)) text((x1 + x2)/2 + 18, (y1 + y2)/2,
                          lbl, cex = 0.75, col = "#777")
}

# Main column boxes (centred at x=40)
node(40, 92, 60, 9,
     "ADNI participants in the UPENN plasma file\n(GFAP_Q + pT217_F available)  N = 1,593",
     fill = "#e3f2fd")
node(40, 76, 60, 7,
     "Within 120-month follow-up window\nN = 1,470")
node(40, 60, 60, 7,
     "≥3 distinct plasma timepoints (GFAP and p-tau217)\nN = 111")
node(40, 44, 60, 7,
     "Baseline diagnosis available (CN / MCI / AD)\nN = 93")
node(40, 26, 60, 10,
     "Analysable cohort\nAge ≥65 at baseline\nN = 77",
     fill = "#c8e6c9", border = "#388e3c", lwd = 2)
node(40, 8, 60, 9,
     "CN: 42 (54.5%)   MCI: 31 (40.3%)   AD: 4 (5.2%)",
     fill = "#fff9c4", border = "#f57f17")

# arrows (vertical)
for (pair in list(c(87,80), c(72,72-3.5), c(56,48), c(40,32), c(20,12))) {
  arrows(40, pair[1], 40, pair[2], length = 0.12, angle = 25, lwd = 1.5, col = "#333")
}

# Exclusion side-boxes on the right
node(82, 76, 30, 9,
     "Excluded: follow-up\n>120 months\nor only 1–2 timepoints\n(n = 1,482)",
     fill = "#ffebee", border = "#999", cex = 0.78)
arrows(70, 76, 67, 76, length = 0.10, angle = 25, lwd = 1.2,
       col = "#999", lty = 2)

node(82, 35, 30, 6.5,
     "Excluded: age <65\n(n = 16)",
     fill = "#ffebee", border = "#999", cex = 0.78)
arrows(70, 35, 67, 35, length = 0.10, angle = 25, lwd = 1.2,
       col = "#999", lty = 2)

par(op); dev.off()
message("Saved: ", file.path(FIG_DIR, "Figure1_study_flow.tiff"))
