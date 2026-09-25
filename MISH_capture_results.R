# ============================================================================
#  MISH RESULTS CAPTURE
#  Run this AFTER MISH_PLS_SEM_analysis.R has completed in your session.
#  It collects every output needed to write Methodology / Analysis / Discussion
#  and writes it all to ./MISH_results_capture.txt
#
#  HOW TO USE
#  1. Run your main analysis script first so all objects exist in memory.
#  2. source("MISH_capture_results.R")
#  3. Send me the resulting MISH_results_capture.txt file.
# ============================================================================

out_file <- "MISH_results_capture.txt"
sink(out_file, split = FALSE)
on.exit(sink(), add = TRUE)

# Helper: print a labelled section header
section <- function(title) {
  cat("\n\n", strrep("=", 78), "\n", sep = "")
  cat("== ", title, "\n", sep = "")
  cat(strrep("=", 78), "\n\n", sep = "")
}

# Helper: safe print — won't crash the capture if an object doesn't exist
safe <- function(expr, label = NULL) {
  if (!is.null(label)) cat("\n--- ", label, " ---\n", sep = "")
  tryCatch(print(expr),
           error = function(e) cat("[not available: ", conditionMessage(e), "]\n", sep = ""))
}

# ---------------------------------------------------------------------------
# 1) METADATA
# ---------------------------------------------------------------------------
section("METADATA & ENVIRONMENT")
cat("Date generated   :", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
cat("R version        :", R.version.string, "\n")
cat("seminr version   :", as.character(packageVersion("seminr")), "\n")
cat("NCA version      :", as.character(packageVersion("NCA")), "\n")
cat("Working directory:", getwd(), "\n")

# ---------------------------------------------------------------------------
# 2) SAMPLE & DATA STRUCTURE
# ---------------------------------------------------------------------------
section("SAMPLE & DATA STRUCTURE")
cat("Sample size (n)    :", nrow(mish), "\n")
cat("Indicator count    :", ncol(mish), "\n")
cat("Missing cells      :", sum(is.na(mish)), "\n\n")

cat("Indicator blocks (after purification):\n")
for (nm in names(items)) cat(sprintf("  %-9s (%d items): %s\n",
                                     nm, length(items[[nm]]),
                                     paste(items[[nm]], collapse = ", ")))

cat("\nNote on indicator purification:\n")
cat("  SERVQUAL7: fully empty in the export (dropped before analysis).\n")
cat("  SERVQUAL3: 95.7% identical responses ('5') — zero-variance bootstrap\n")
cat("             instability + no covariance information. Dropped.\n")

# Show the raw context columns (role, modules, usage mode) if still in the env
if (exists("raw")) {
  cat("\nContext variable distributions (raw export):\n")
  if ("Označite kako koristite MISH:" %in% names(raw)) {
    cat("\nUsage mode:\n"); print(table(raw[["Označite kako koristite MISH:"]]))
  }
  if ("Radno mjesto" %in% names(raw)) {
    cat("\nRoles (top 15):\n")
    print(head(sort(table(raw[["Radno mjesto"]]), decreasing = TRUE), 15))
  }
}

# ---------------------------------------------------------------------------
# 3) DESCRIPTIVE STATISTICS
# ---------------------------------------------------------------------------
section("DESCRIPTIVE STATISTICS — ITEM LEVEL")
safe(round(desc_tbl, 3))

section("DESCRIPTIVE STATISTICS — CONSTRUCT LEVEL")
safe(con_summary, "Construct mean & SD (1–5 scale)")

cat("\nCronbach's alpha per construct:\n")
safe(round(alpha_tbl, 3))

cat("\nResponse distribution (counts across all items):\n")
safe(table(unlist(mish[, all_items])))

cat("\nPercentage ceiling per construct (% of '5' responses):\n")
for (nm in names(items)) {
  v <- as.matrix(mish[, items[[nm]]])
  cat(sprintf("  %-9s : %.1f%% '5'  |  mean = %.2f  |  sd = %.2f\n",
              nm, 100 * mean(v == 5, na.rm = TRUE),
              mean(v, na.rm = TRUE), sd(v, na.rm = TRUE)))
}

# ---------------------------------------------------------------------------
# COMMON METHOD BIAS
# ---------------------------------------------------------------------------
section("COMMON METHOD BIAS (CMB) CHECKS")

cat("Harman's single-factor test:\n")
cat("  Variance explained by first unrotated factor:",
    round(harman_first_factor_pct, 2), "%\n\n")

cat("Full-collinearity VIF values:\n")
print(round(full_vif, 3))

cat("\nMaximum full-collinearity VIF:",
    round(max_full_vif, 3), "\n")

# ---------------------------------------------------------------------------
# 4) MEASUREMENT MODEL ASSESSMENT (Model C)
# ---------------------------------------------------------------------------
section("MEASUREMENT MODEL — MODEL C (PRIMARY)")

cat("Structural specification (Model C):\n")
cat("  INFOQ,  SYSQUAL, SERVQUAL  ->  PU\n")
cat("  SERVQUAL                   ->  PEOU\n")
cat("  PU, PEOU                   ->  SAT\n")

safe(round(s$reliability, 3), "Reliability & convergent validity (alpha, rhoC, AVE, rhoA)")
safe(round(s$validity$htmt, 3), "Discriminant validity — HTMT matrix")
safe(round(s$loadings, 3),  "Standardised indicator loadings")
safe(round(s$weights, 3),   "Outer weights")
safe(round(s$validity$fl_criteria, 3), "Fornell–Larcker criterion (sqrt(AVE) on diagonal)")

# ---------------------------------------------------------------------------
# 5) STRUCTURAL MODEL — MODEL C
# ---------------------------------------------------------------------------
#Testira odnose između varijabla
section("STRUCTURAL MODEL — MODEL C (PRIMARY)")
safe(round(pm, 3),                     "Path coefficients matrix")
safe(round(s$total_effects, 3),        "Total effects")
safe(round(s$total_indirect_effects, 3), "Total indirect effects")
safe(round(s$fSquare, 3),              "f² effect sizes")
safe(s$vif_antecedents,                "VIF of antecedents (collinearity)")

cat("\nR² per endogenous construct:\n")
safe(round(r2_df, 3))

safe(boot_tbl, "Bootstrap path inference (5,000 resamples; 95% percentile CI)")

# ---------------------------------------------------------------------------
# 6) MEDIATION ANALYSIS — MODEL C
# ---------------------------------------------------------------------------
#Medijacija znači da jedna varijabla ne utječe izravno na drugu, nego preko posrednika (medijatora).
section("MEDIATION — MODEL C (indirect effects via PU/PEOU)")
safe(ind_df, "Indirect effects (product-of-coefficients bootstrap)")

# ---------------------------------------------------------------------------
# 7) MODEL B (FULL MEDIATION COMPARATOR)
# ---------------------------------------------------------------------------
#testira model u kojem INFOQ, SYSQUAL i SERVQUAL ne idu direktno na SAT, 
#nego djeluju preko medijatora PU i PEOU.

section("MODEL B — COMPARATOR (full mediation)")
cat("Structural specification (Model B):\n")
cat("  INFOQ, SYSQUAL, SERVQUAL  ->  PU and PEOU\n")
cat("  PU, PEOU                  ->  SAT\n\n")

safe(round(sB$reliability, 3), "Model B reliability")
safe(round(pmB, 3),            "Model B path coefficients + R²")
safe(boot_tbl_B,               "Model B bootstrap path inference")

# ---------------------------------------------------------------------------
# 8) MODEL COMPARISON
# ---------------------------------------------------------------------------
#uspoređuje koliko su Model C i Model B uspješni u objašnjavanju podataka i koliko su procijenjeni 
#efekti stabilni između specifikacija.
#znači da jedan model objašnjava 10 postotnih bodova više varijance od drugog za taj konstrukt.

section("MODEL COMPARISON — Model C vs Model B")
safe(as.data.frame(r2_cmp), "R² comparison with Δ")

cat("\nPaths into SAT — coefficient stability across specifications:\n")
safe(paths_cmp)

# ---------------------------------------------------------------------------
# 9) IPMA — IMPORTANCE-PERFORMANCE MAP ANALYSIS
# ---------------------------------------------------------------------------
#Koje konstrukte vrijedi unaprijediti da bi se najviše povećalo zadovoljstvo (SAT)?
#Ne gleda samo koliko je neki konstrukt važan, nego i koliko je trenutno dobro ocijenjen.

section("IPMA — target = SAT (Model C)")
safe(transform(as.data.frame(ipma),
               #To je ukupni efekt (total effect) konstrukta na ciljnu varijablu (SAT).
               #“Ako poboljšamo ovaj konstrukt, zadovoljstvo će se više povećati.”
               importance  = round(importance,  3),
               #To je prosječna razina konstrukta, transformirana na ljestvicu 0–100.
               performance = round(performance, 3)),
     "Importance (total effect) × Performance (0–100)")

cat("\nReference lines:\n")
#Ti prosjeci služe za crtanje okomitih i vodoravnih referentnih linija na IPMA grafu, 
#čime se dobivaju četiri kvadranta.
cat("  Mean importance  :", round(mean(ipma$importance,  na.rm = TRUE), 3), "\n")
cat("  Mean performance :", round(mean(ipma$performance, na.rm = TRUE), 3), "\n")

# ---------------------------------------------------------------------------
# 10) NCA — NECESSARY CONDITION ANALYSIS
# ---------------------------------------------------------------------------
#Koji uvjeti moraju postojati da bi     nužni uvjeti  zadovoljstvo bilo visoko?
#    prikazuje NCA effect size i p-vrijednosti,
section("NCA — necessary conditions for high SAT")
safe(nca_res$summaries, "NCA effect sizes (CE-FDH, CR-FDH) with permutation p-values")

cat("\nNCA bottleneck table (CE-FDH):\n")
#Kolika minimalna razina nekog uvjeta je potrebna za određenu razinu SAT-a.
tryCatch({
  #tryCatch() služi da se skripta ne sruši ako NCA bottleneck output ne radi.
  bt <- NCA::nca_output(nca_res, summaries = FALSE,
                         bottlenecks = TRUE, plots = FALSE)
  print(bt)
}, error = function(e) cat("[bottleneck output not available: ", conditionMessage(e), "]\n"))

# ---------------------------------------------------------------------------
# 11) PLOT INVENTORY
# ---------------------------------------------------------------------------
section("PLOT INVENTORY")
plot_dir <- "plots"
if (dir.exists(plot_dir)) {
  cat("Figures saved to ./", plot_dir, "/:\n", sep = "")
  for (f in sort(list.files(plot_dir, pattern = "\\.png$")))
    cat("  ", f, "\n", sep = "")
} else {
  cat("[plots directory not found]\n")
}

# ---------------------------------------------------------------------------
# 12) SESSION INFO
# ---------------------------------------------------------------------------
section("SESSION INFO")
print(sessionInfo())

# ---------------------------------------------------------------------------
# Close sink and notify
# ---------------------------------------------------------------------------
sink()
on.exit()  # clear the on.exit since sink is already closed
cat("\n=================================================================\n")
cat("Results captured to: ", normalizePath(out_file), "\n", sep = "")
cat("File size: ", round(file.info(out_file)$size / 1024, 1), " KB\n", sep = "")
cat("Send this file along with the ./plots/ folder.\n")
cat("=================================================================\n")
