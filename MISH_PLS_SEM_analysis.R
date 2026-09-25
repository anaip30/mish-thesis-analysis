# ============================================================================
#  MISH USER-EVALUATION STUDY  —  PLS-SEM + IPMA + NCA + Descriptives
#  Model family: DeLone & McLean IS-Success  ×  Technology Acceptance Model
#  Engine: seminr 2.5.0 (bootstrap via manual percentile loop) | ggplot2
# ----------------------------------------------------------------------------
#  STRUCTURAL MODELS IN THIS SCRIPT
#  Model C (PRIMARY — trimmed): quality → PU; SERVQUAL → PEOU; PU/PEOU → SAT
#    Retains only paths that showed bootstrap signal in the full model.
#    Reduces multicollinearity; stabilises SAT path coefficients.
#  Model B (COMPARATOR — full mediation): all quality → PU + PEOU; PU/PEOU → SAT
#    Quality reaches SAT only via beliefs. Used for R²/path comparison.
#
#  HOW TO RUN: source() top-to-bottom. First run installs missing packages.
#  All figures saved to ./plots/ as PNG.
# ============================================================================


# ---------------------------------------------------------------------------
# 0)  PACKAGES
# ---------------------------------------------------------------------------
required <- c("seminr", "NCA", "ggplot2", "dplyr", "tidyr", "tibble",
              "forcats", "scales", "ggrepel", "psych")
to_install <- setdiff(required, rownames(installed.packages()))
if (length(to_install)) install.packages(to_install, repos = "https://cloud.r-project.org")
invisible(lapply(required, library, character.only = TRUE))

set.seed(123)
dir.create("plots", showWarnings = FALSE)

# Shared ggplot theme
theme_mish <- theme_minimal(base_size = 12) +
  theme(panel.grid.minor  = element_blank(),
        plot.title        = element_text(face = "bold"),
        plot.subtitle     = element_text(colour = "grey35"),
        legend.position   = "bottom")

# ============================================================
# PROVJERA PLSpredict FUNKCIJA U seminr PAKETU
# ============================================================

packageVersion("seminr")

ls("package:seminr", pattern = "predict")
args(predict_pls)
help.search("PLSpredict", package = "seminr")

# ---------------------------------------------------------------------------
# HELPERS  (defined once, used throughout)
# ---------------------------------------------------------------------------

# Round only numeric columns of a matrix or data frame.
round_num <- function(x, digits = 3) {
  if (is.matrix(x)) return(round(x, digits))
  x <- as.data.frame(x)
  x[vapply(x, is.numeric, logical(1))] <-
    lapply(x[vapply(x, is.numeric, logical(1))], round, digits)
  x
}

# Safe R² lookup: returns NA for exogenous constructs; auto-detects row name.
r2_row_of <- function(m) {
  rn <- grep("^R\\^?2$", rownames(m), value = TRUE)
  if (length(rn)) rn[1] else NA_character_
}
r2_of <- function(m, cn) {
  rr <- r2_row_of(m)
  if (!is.na(rr) && cn %in% colnames(m)) as.numeric(m[rr, cn]) else NA_real_
}


# ---------------------------------------------------------------------------
# 1)  LOAD & PREPARE DATA
# ---------------------------------------------------------------------------
setwd("~/Desktop/MISH")
data_path <- file.choose()
raw <- read.csv2(data_path, check.names = FALSE, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
#data_path <- "MISH-2026-06-10.csv"   # edit if needed

raw  <- read.csv2(data_path, check.names = FALSE,
                  stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# Indicator blocks.
# SERVQUAL7 fully empty; SERVQUAL3 dropped (95.7% = "5" → zero-variance bootstrap).
items <- list(
  PU       = paste0("PU",       1:5),
  PEOU     = paste0("PEOU",     1:6),
  INFOQ    = paste0("INFOQ",    1:5),
  SYSQUAL  = paste0("SYSQUAL",  1:5),
  SERVQUAL = paste0("SERVQUAL", c(1,2,4,5,6)),
  SAT      = paste0("SAT",      1:5)
)
all_items <- unlist(items, use.names = FALSE)

mish        <- raw[, all_items]
mish[]      <- lapply(mish, function(x) as.numeric(as.character(x)))

cat("\n=== DATA CHECK ===\n")
cat("Respondents:", nrow(mish), "| Indicators:", ncol(mish),
    "| Missing:", sum(is.na(mish)), "\n")


# ---------------------------------------------------------------------------
# 2)  DESCRIPTIVE STATISTICS
# ---------------------------------------------------------------------------

## 2a) Per-item descriptives (skew/kurtosis flag the ceiling effect)
desc_tbl <- psych::describe(mish)[, c("n","mean","sd","median","min","max","skew","kurtosis")]
#skew i kurtosis pomažu vidjeti postoji li ceiling effect, tj. previše odgovora 5.
cat("\n=== ITEM-LEVEL DESCRIPTIVES ===\n"); print(round(desc_tbl, 2))

## 2b) Cronbach's alpha per construct
cat("\n=== CRONBACH'S ALPHA PER CONSTRUCT ===\n")
#Govori koliko su itemi unutar istog konstrukta međusobno konzistentni.
alpha_tbl <- sapply(items, function(cols)
  suppressWarnings(psych::alpha(mish[, cols], warnings = FALSE)$total$raw_alpha))
print(round(alpha_tbl, 3))

## 2c) Likert diverging plot
#Ovaj dio pretvara podatke iz širokog u dugi format:pivot_longer
long <- mish |>
  tibble::rownames_to_column("id") |>
  tidyr::pivot_longer(-id, names_to = "item", values_to = "resp") |>
  dplyr::filter(!is.na(resp)) |>
  dplyr::mutate(construct = forcats::fct_inorder(sub("[0-9]+$", "", item)),
                resp = factor(resp, levels = 1:5))

prop <- long |>
  dplyr::count(construct, item, resp, .drop = FALSE) |>
  dplyr::group_by(item) |>
  dplyr::mutate(p = n / sum(n)) |>
  dplyr::ungroup() |>
  dplyr::mutate(item = forcats::fct_inorder(item))

p_likert <- ggplot(prop, aes(x = p, y = forcats::fct_rev(item), fill = resp)) +
  geom_col(width = 0.8) +
  facet_grid(construct ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_continuous(labels = scales::percent_format(), expand = expansion()) +
  scale_fill_brewer(palette = "RdYlGn", name = "Response (1=low, 5=high)") +
  labs(title = "Response distribution per item (Likert)",
       subtitle = "Strong concentration on '5' = ceiling effect across all constructs",
       x = "Share of respondents", y = NULL) +
  theme_mish +
  theme(strip.placement = "outside",
        strip.text.y.left = element_text(angle = 0, face = "bold"))

## 2d) Construct mean ± SD
#Računa prosjek i standardnu devijaciju za svaki konstrukt.
con_summary <- long |>
  dplyr::mutate(resp = as.numeric(as.character(resp))) |>
  dplyr::group_by(construct) |>
  dplyr::summarise(mean = mean(resp), sd = sd(resp), .groups = "drop")

p_means <- ggplot(con_summary, aes(x = mean, y = reorder(construct, mean))) +
  geom_vline(xintercept = c(3, 4), linetype = 3, colour = "grey70") +
  geom_errorbarh(aes(xmin = mean - sd, xmax = mean + sd), height = 0.18, colour = "grey50") +
  geom_point(size = 4, colour = "#2c7fb8") +
  geom_text(aes(label = sprintf("%.2f", mean)), vjust = -1.1, size = 3.4) +
  scale_x_continuous(limits = c(1, 5), breaks = 1:5) +
  labs(title = "Construct means (± 1 SD)",
       subtitle = "All constructs sit near the top of the 1–5 scale",
       x = "Mean score", y = NULL) +
  theme_mish

ggsave("plots/01_likert.png", p_likert, width = 9, height = 8, dpi = 150)
ggsave("plots/02_means.png",  p_means,  width = 7, height = 4, dpi = 150)


# ---------------------------------------------------------------------------
# 3) COMMON METHOD BIAS (CMB) CHECKS
# ---------------------------------------------------------------------------

# 3a) Harman's single-factor test
# Unrotated exploratory factor analysis of all retained indicators.
# Reports the percentage of total variance explained by the first factor.

harman <- psych::fa(
  mish[, all_items],
  nfactors = 1,
  rotate = "none",
  fm = "pa"
)

harman_first_factor_pct <-
  as.numeric(harman$Vaccounted["Proportion Var", 1]) * 100

cat("\n=== HARMAN'S SINGLE-FACTOR TEST ===\n")
cat("Variance explained by first unrotated factor:",
    round(harman_first_factor_pct, 2), "%\n")


# 3b) Full-collinearity VIF
# Construct scores are calculated as means of their retained indicators.
# Each construct is regressed on all remaining constructs.

cmb_scores <- as.data.frame(
  lapply(items, function(cols) {
    rowMeans(mish[, cols, drop = FALSE], na.rm = TRUE)
  })
)

full_collinearity_vif <- function(target, data) {
  
  predictors <- setdiff(names(data), target)
  
  fit <- lm(
    reformulate(predictors, response = target),
    data = data
  )
  
  r2 <- summary(fit)$r.squared
  
  1 / (1 - r2)
}

full_vif <- sapply(
  names(cmb_scores),
  full_collinearity_vif,
  data = cmb_scores
)

cat("\n=== FULL COLLINEARITY VIF ===\n")
print(round(full_vif, 3))

max_full_vif <- max(full_vif, na.rm = TRUE)

cat("Maximum full-collinearity VIF:",
    round(max_full_vif, 3), "\n")

# ---------------------------------------------------------------------------
# 3)  PLS-SEM — MODEL C (PRIMARY, TRIMMED)
#
#  WHY THIS SPECIFICATION:
#  The original full model (all quality → PEOU + PU + SAT, PEOU → PU) was
#  estimated first. Bootstrap CIs revealed:
#   - Only SERVQUAL → PEOU was significant among quality → PEOU paths.
#   - INFOQ → PU and SYSQUAL → PU were significant; SERVQUAL → PU and
#     PEOU → PU were not.
#   - All direct quality → SAT and beliefs → SAT paths had CIs spanning
#     zero, with SEs up to 2.5 (severe multicollinearity from 5 correlated
#     predictors on n = 69).
#  Model C retains only the paths with signal and removes the redundant ones,
#  reducing the predictor count for SAT from 5 to 2 (PU, PEOU) and letting
#  their coefficients stabilise.
# ---------------------------------------------------------------------------

## 3a) Measurement model — all reflective (Mode A)
#Definira koji indikatori mjere koji konstrukt.

measurement_model <- constructs(
  reflective("PU",       items$PU),
  reflective("PEOU",     items$PEOU),
  reflective("INFOQ",    items$INFOQ),
  reflective("SYSQUAL",  items$SYSQUAL),
  reflective("SERVQUAL", items$SERVQUAL),
  reflective("SAT",      items$SAT)
)

## 3b) Structural model C — trimmed, empirically supported paths only
#Definira veze među konstruktima.
structural_model <- relationships(
  paths(from = c("INFOQ", "SYSQUAL", "SERVQUAL"), to = "PU"),    # quality → PU
  paths(from = "SERVQUAL",                         to = "PEOU"),  # only sig. quality → PEOU
  paths(from = c("PU", "PEOU"),                    to = "SAT")    # beliefs → SAT
)

## 3c) Estimate
#Procjenjuje PLS-SEM model.
pls_model <- estimate_pls(
  data              = mish,
  measurement_model = measurement_model,
  structural_model  = structural_model,
  inner_weights     = path_weighting
)
s  <- summary(pls_model)
#    pm = path coefficients + R² tablica
pm <- s$paths

## 3d) Self-contained percentile bootstrap
#  bootstrap_model() in seminr 2.5.0 returns a single fit, not resamples.
#  We resample manually: draw n rows with replacement, re-estimate, store
#  path_coef. Zero-variance resamples are caught and skipped silently.
run_boot <- function(model, data, R = 5000, seed = 123) {
  set.seed(seed)
  n      <- nrow(data)
  pc     <- model$path_coef
  paths  <- which(!is.na(pc) & pc != 0, arr.ind = TRUE)
  pnames <- paste0(rownames(pc)[paths[,1]], " -> ", colnames(pc)[paths[,2]])
  draws  <- matrix(NA_real_, nrow = R, ncol = nrow(paths),
                   dimnames = list(NULL, pnames))
  ok <- 0L
  for (i in seq_len(R * 2L)) {
    if (ok >= R) break
    idx <- sample.int(n, n, replace = TRUE)
    # Per-iteration sink to nullfile silences seminr's stdout chatter
    # ("All N observations are valid / Generating the seminr model").
    # Opening/closing inside the loop ensures no dangling sink on error.
    nf <- file(nullfile(), open = "wt")
    sink(nf)
    fit <- tryCatch(
      estimate_pls(data[idx, ],
                   measurement_model = model$measurement_model,
                   structural_model  = model$structural_model,
                   inner_weights     = path_weighting),
      error = function(e) NULL)
    sink()
    close(nf)
    if (is.null(fit) || !inherits(fit, "seminr_model")) next
    ok <- ok + 1L
    draws[ok, ] <- fit$path_coef[paths]
    if (ok %% 100 == 0 || ok == R)
      message(sprintf("Bootstrap: %4d / %d", ok, R))
  }
  draws[seq_len(ok), , drop = FALSE]
}

boot_summary <- function(model, draws) {
  pc    <- model$path_coef
  paths <- which(!is.na(pc) & pc != 0, arr.ind = TRUE)
  orig  <- pc[paths]
  lo95  <- apply(draws, 2, quantile, 0.025, na.rm = TRUE)
  hi95  <- apply(draws, 2, quantile, 0.975, na.rm = TRUE)
  se    <- apply(draws, 2, sd,       na.rm = TRUE)
  data.frame(path     = colnames(draws),
             original = round(orig,       3),
             boot_se  = round(se,         3),
             t_stat   = round(orig / se,  2),
             ci_lo    = round(lo95,       3),
             ci_hi    = round(hi95,       3),
             sig      = !(lo95 <= 0 & hi95 >= 0),
             row.names = NULL)
}

cat("\nRunning bootstrap for Model C (~2-3 min)...\n")
boot_draws <- run_boot(pls_model, mish, R = 5000, seed = 123)
boot_tbl   <- boot_summary(pls_model, boot_draws)

"ci_lo > 0 i ci_hi > 0  → značajno pozitivno
ci_lo < 0 i ci_hi < 0  → značajno negativno
ci_lo < 0 i ci_hi > 0  → nije značajno"

## 3e) Print assessment battery
cat("\n=== RELIABILITY & CONVERGENT VALIDITY (alpha, rhoC, AVE, rhoA) ===\n")
print(round(s$reliability, 3))
cat("\n=== DISCRIMINANT VALIDITY: HTMT (want < 0.85 / 0.90) ===\n")
print(round(s$validity$htmt, 3))
cat("\n=== INDICATOR LOADINGS (want > 0.708) ===\n")
print(round(s$loadings, 3))
cat("\n=== PATH COEFFICIENTS + R^2 ===\n")
print(round(pm, 3))
cat("\n=== f^2 EFFECT SIZES (0.02/0.15/0.35 = S/M/L) ===\n")
print(round(s$fSquare, 3))
cat("\n=== VIF of antecedents (want < 3, hard cap 5) ===\n")
print(s$vif_antecedents)
cat("\n=== BOOTSTRAPPED PATHS (95% percentile CI) ===\n")
print(boot_tbl)

## OPTIONAL Q² predict
#prikazuje prediktivnu sposobnost modela.
#try({
 # pred <- predict_pls(model = pls_model, technique = predict_DA, noFolds = 10, reps = 10)
  #cat("\n=== PLSpredict (PLS vs LM RMSE) ===\n"); print(summary(pred))
#}, silent = TRUE)

# ---------------------------------------------------------------------------

# 3f) PLSpredict — PREDIKTIVNA SPOSOBNOST MODELA C

# ---------------------------------------------------------------------------

set.seed(123)
pls_predict_C <- predict_pls(
  model     = pls_model,
  technique = predict_DA,
  noFolds   = 10,
  reps      = 10,
  cores     = 1
)
cat("\n=== PLSpredict — MODEL C ===\n")
summary(pls_predict_C)


# ---------------------------------------------------------------------------
# 4)  ggplot VISUALS — PLS-SEM RESULTS
# ---------------------------------------------------------------------------

## Helpers that use boot_tbl
bp       <- boot_tbl
bp$edge  <- bp$path

#Provjerava je li putanja značajna
is_sig <- function(from, to) {
  row <- bp[grepl(from, bp$edge) & grepl(to, bp$edge), , drop = FALSE]
  if (!nrow(row)) return(NA)
  !(row$ci_lo[1] <= 0 & row$ci_hi[1] >= 0)
}

#Izvlači koeficijent putanje iz matrice pm.
get_coef <- function(from, to) {
  if (from %in% rownames(pm) && to %in% colnames(pm)) as.numeric(pm[from, to]) else NA
}

## 4a) Structural path diagram -----------------------------------------------
# Node layout: quality left, beliefs centre, SAT right.
nodes <- tibble::tribble(
  ~construct,   ~x,  ~y,   ~role,
  "INFOQ",       0,  3.2, "exo",
  "SYSQUAL",     0,  2.2, "exo",
  "SERVQUAL",    0,  1.2, "exo",
  "PEOU",        2,  2.8, "endo",
  "PU",          2,  1.6, "endo",
  "SAT",         4,  2.2, "endo"
) |>
  dplyr::mutate(
    r2    = round(vapply(construct, function(cn) r2_of(pm, cn), numeric(1)), 2),
    label = ifelse(is.na(r2), construct,
                   paste0(construct, "\nR² = ", sprintf("%.2f", r2)))
  )

# Edges — mirrors structural_model above
edge_def <- tibble::tibble(
  from = c("INFOQ","SYSQUAL","SERVQUAL","SERVQUAL","PU","PEOU"),
  to   = c("PU",   "PU",     "PU",      "PEOU",    "SAT","SAT")
)

edges <- edge_def |>
  dplyr::rowwise() |>
  dplyr::mutate(coef = get_coef(from, to),
                sig  = is_sig(from, to)) |>
  dplyr::ungroup() |>
  dplyr::left_join(dplyr::select(nodes, construct, x1 = x, y1 = y),
                   by = c("from" = "construct")) |>
  dplyr::left_join(dplyr::select(nodes, construct, x2 = x, y2 = y),
                   by = c("to"   = "construct")) |>
  dplyr::mutate(midx = (x1 + x2) / 2,
                midy = (y1 + y2) / 2,
                sign = ifelse(coef >= 0, "positive", "negative"))

p_path <- ggplot() +
  geom_curve(data = edges,
             aes(x = x1, y = y1, xend = x2, yend = y2,
                 colour = sign, linewidth = abs(coef), linetype = sig),
             curvature = 0.12, alpha = 0.75,
             arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
  geom_label(data = edges,
             aes(x = midx, y = midy, label = sprintf("%.2f", coef)),
             size = 3, label.size = 0, fill = "white", colour = "grey20") +
  geom_label(data = nodes,
             aes(x = x, y = y, label = label, fill = role),
             size = 3.6, fontface = "bold",
             label.r = unit(0.4, "lines"), label.padding = unit(0.5, "lines")) +
  scale_fill_manual(values = c(exo = "#deebf7", endo = "#fee6ce"), guide = "none") +
  scale_colour_manual(values = c(positive = "#1b7837", negative = "#b2182b"),
                      name = "Path sign") +
  scale_linetype_manual(values   = c(`TRUE` = "solid", `FALSE` = "22"),
                        labels   = c(`TRUE` = "sig.", `FALSE` = "n.s."),
                        name     = "Significance",
                        na.value = "solid") +
  scale_linewidth(range = c(0.4, 2.2), guide = "none") +
  coord_cartesian(xlim = c(-0.6, 4.8), ylim = c(0.6, 3.8)) +
  labs(title    = "PLS-SEM structural model (Model C — trimmed)",
       subtitle = "Path coefficients on arrows; R² inside endogenous constructs") +
  theme_void(base_size = 12) +
  theme(plot.title    = element_text(face = "bold"),
        plot.subtitle = element_text(colour = "grey35"),
        legend.position = "bottom")

## 4b) Path forest plot -------------------------------------------------------
#Crta koeficijente putanja s 95% bootstrap intervalima pouzdanosti.
p_forest <- ggplot(boot_tbl, aes(x = original, y = reorder(path, original),
                                 colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.7) +
  scale_colour_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "significant", `FALSE` = "n.s."),
                      name = NULL) +
  labs(title    = "Structural paths — Model C (bootstrap 95% CIs)",
       subtitle = "CIs crossing 0 are not significant",
       x = "Path coefficient", y = NULL) +
  theme_mish

## 4c) R² bar chart -----------------------------------------------------------
#Crta koliko model objašnjava varijance za endogene konstrukte:    PU,    PEOU,    SAT
r2_df <- tibble::tibble(
  construct = colnames(pm),
  r2 = vapply(colnames(pm), function(cn) r2_of(pm, cn), numeric(1))
) |> dplyr::filter(!is.na(r2))

p_r2 <- ggplot(r2_df, aes(x = r2, y = reorder(construct, r2))) +
  geom_col(fill = "#3182bd", width = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r2)), hjust = -0.2, size = 3.6) +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Explained variance (R²) — Model C",
       x = expression(R^2), y = NULL) +
  theme_mish

## 4d) HTMT heatmap -----------------------------------------------------------
#Crta matricu HTMT za diskriminantnu valjanost.

htmt      <- s$validity$htmt
htmt_long <- as.data.frame(as.table(as.matrix(htmt))) |>
  setNames(c("c1","c2","value")) |>
  dplyr::filter(!is.na(value), value != 0)

p_htmt <- ggplot(htmt_long, aes(c1, c2, fill = value)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient2(low = "#2166ac", mid = "#f7f7f7", high = "#b2182b",
                       midpoint = 0.85, limits = c(0,1), name = "HTMT") +
  labs(title    = "Discriminant validity — HTMT matrix",
       subtitle = "Values above ~0.85–0.90 flag concerns",
       x = NULL, y = NULL) +
  theme_mish + theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("plots/03_path_diagram.png", p_path,   width = 9, height = 6, dpi = 150)
ggsave("plots/04_path_forest.png",  p_forest, width = 8, height = 5, dpi = 150)
ggsave("plots/05_r2.png",           p_r2,     width = 7, height = 4, dpi = 150)
ggsave("plots/06_htmt.png",         p_htmt,   width = 7, height = 6, dpi = 150)




# ---------------------------------------------------------------------------
# 5)  IPMA — Importance–Performance Map Analysis (target = SAT)
#
#  IMPORTANCE  = total effect of each predictor on SAT (direct + indirect).
#  PERFORMANCE = construct item-mean rescaled to 0–100.
#  Both computed from Model C.
# ---------------------------------------------------------------------------
#Koji konstrukti su važni za zadovoljstvo i koliko su trenutno dobro ocijenjeni?
target     <- "SAT"
predictors <- setdiff(names(items), target)
te         <- s$total_effects

importance <- sapply(predictors, function(p) {
  if (p %in% rownames(te) && target %in% colnames(te))
    as.numeric(te[p, target]) else NA_real_
})

performance <- sapply(predictors, function(p) {
  m <- mean(as.matrix(mish[, items[[p]]]), na.rm = TRUE)
  (m - 1) / 4 * 100
})

ipma <- tibble::tibble(construct  = predictors,
                       importance = importance,
                       performance = performance)
cat("\n=== IPMA TABLE (target = SAT) ===\n")
print(transform(as.data.frame(ipma),
                importance  = round(importance,  3),
                performance = round(performance, 3)))

p_ipma <- ggplot(ipma, aes(importance, performance)) +
  geom_hline(yintercept = mean(ipma$performance, na.rm = TRUE),
             linetype = 2, colour = "grey60") +
  geom_vline(xintercept = mean(ipma$importance,  na.rm = TRUE),
             linetype = 2, colour = "grey60") +
  geom_point(size = 4, colour = "#d95f02") +
  ggrepel::geom_text_repel(aes(label = construct), size = 4, seed = 1) +
  labs(title    = "IPMA — drivers of Satisfaction (SAT) — Model C",
       subtitle = "Bottom-right = high importance but lower performance = priority to improve",
       x = "Importance (total effect on SAT)",
       y = "Performance (0–100)") +
  theme_mish
ggsave("plots/07_ipma.png", p_ipma, width = 8, height = 6, dpi = 150)


# ---------------------------------------------------------------------------
# 6)  NCA — Necessary Condition Analysis
#  Tests whether each antecedent is a NECESSARY (bottleneck) condition for
#  high SAT. Complements PLS-SEM's average-effects logic.
# ---------------------------------------------------------------------------
#Koji uvjet mora postojati da bi SAT bio visok?
#Za svaki konstrukt računa prosjek njegovih itema.
cscore  <- as.data.frame(
  sapply(items, function(cols) rowMeans(mish[, cols], na.rm = TRUE)))
#Testira jesu li ti konstrukti nužni uvjeti za visoki SAT.
nca_x   <- c("PU", "PEOU", "INFOQ", "SYSQUAL", "SERVQUAL")

nca_res <- NCA::nca_analysis(cscore, x = nca_x, y = "SAT",
                             ceilings  = c("ce_fdh", "cr_fdh"),
                             test.rep  = 1000)
cat("\n=== NCA SUMMARY (d ~0.1 small / 0.3 medium / 0.5 large) ===\n")
print(nca_res$summaries)

# CE-FDH ceiling helper
#To je pomoćna funkcija koja ručno računa CE-FDH ceiling line.
#pronalazi granicu iznad koje nema opažanja — tzv. impossible zone.

ce_fdh <- function(x, y) {
  o   <- order(x); xs <- x[o]; ys <- y[o]; cmy <- cummax(ys)
  peer <- ys == cmy & c(TRUE, diff(cmy) > 0)
  list(line  = data.frame(x = xs, y = cmy),
       peers = data.frame(x = xs[peer], y = ys[peer]))
}
#Crta NCA graf za svaki prediktor prema SAT.
nca_plot <- function(xvar) {
  d  <- data.frame(x = cscore[[xvar]], y = cscore[["SAT"]])
  ce <- ce_fdh(d$x, d$y)
  cr <- lm(y ~ x, data = ce$peers)
  ggplot(d, aes(x, y)) +
    geom_ribbon(data = ce$line,
                aes(x = x, ymin = y, ymax = max(d$y)),
                inherit.aes = FALSE, fill = "#fdae6b", alpha = 0.25) +
    geom_jitter(width = 0.04, height = 0.04, alpha = 0.6, colour = "grey30") +
    geom_step(data = ce$line, aes(x, y), colour = "#e6550d", linewidth = 0.9) +
    geom_abline(intercept = coef(cr)[1], slope = coef(cr)[2],
                colour = "#3182bd", linetype = 2) +
    labs(title    = paste0("NCA ceiling: ", xvar, " \u2192 SAT"),
         subtitle = "Orange step = CE-FDH; blue dashed = CR-FDH; shaded = impossible zone",
         x = xvar, y = "SAT") +
    theme_mish
}
" 
* sive točke = ispitanici,
* narančasta stepenasta linija = CE-FDH ceiling,
* plava isprekidana linija = CR-FDH ceiling,
* zasjenjeno područje = impossible zone(Za određenu nisku vrijednost X-a nije moguće postići visoki SAT.)."

nca_plots <- lapply(nca_x, nca_plot); names(nca_plots) <- nca_x
for (xv in nca_x)
  ggsave(sprintf("plots/08_nca_%s.png", xv), nca_plots[[xv]],
         width = 6, height = 5, dpi = 150)

cat("\n=== NCA BOTTLENECK TABLE (CE-FDH) ===\n")
#Za SAT od 80%, koliko minimalno treba PU, PEOU, INFOQ itd.?
try(print(NCA::nca_output(nca_res, summaries = FALSE,
                          bottlenecks = TRUE, plots = FALSE)), silent = TRUE)


# ---------------------------------------------------------------------------
# 7)  MODEL B — COMPARATOR (full mediation, all quality → PU + PEOU)- standardizirani koeficijenti putanja,
#  All quality dimensions feed both beliefs; beliefs → SAT only.
#  Compare its R² against Model C to quantify what the richer specification
#  of quality antecedents adds (or doesn't add).
# ---------------------------------------------------------------------------
"radi se o modelu potpune medijacije (full mediation),što znači da kvaliteta ne utječe izravno na 
zadovoljstvo (SAT), nego isključivo preko percipirane korisnosti (PU) 
i percipirane jednostavnosti korištenja (PEOU)."

"Važno je primijetiti da ne postoje izravne veze:

* INFOQ → SAT
* SYSQUAL → SAT
* SERVQUAL → SAT

Zbog toga se govori o potpunoj medijaciji, jer kvaliteta utječe na zadovoljstvo samo preko PU i PEOU."

structural_model_B <- relationships(
  paths(from = c("INFOQ","SYSQUAL","SERVQUAL"), to = c("PU","PEOU")),
  paths(from = c("PU","PEOU"),                  to = "SAT")
)

pls_B <- estimate_pls(data              = mish,
                      measurement_model = measurement_model,
                      structural_model  = structural_model_B,
                      inner_weights     = path_weighting)
sB   <- summary(pls_B)
pmB  <- sB$paths

"Bootstrap omogućuje procjenu:

* standardnih pogrešaka,
* t-vrijednosti,
* p-vrijednosti,
* intervala pouzdanosti.

Na temelju tih rezultata utvrđuje se jesu li pojedini odnosi među konstruktima statistički značajni."

cat("\nRunning bootstrap for Model B (~2-3 min)...\n")
boot_draws_B <- run_boot(pls_B, mish, R = 5000, seed = 123)
boot_tbl_B   <- boot_summary(pls_B, boot_draws_B)

cat("\n=== MODEL B PATHS + R^2 ===\n");    print(round(pmB, 3))
cat("\n=== MODEL B BOOTSTRAPPED PATHS ===\n"); print(boot_tbl_B)


# ---------------------------------------------------------------------------
# 8)  MEDIATION ANALYSIS (Model C)
#  For each quality → PU → SAT chain: compute indirect effect as product of
#  bootstrap draws (path1_b × path2_b) and derive percentile CI.
#  If indirect CI excludes 0 and direct quality → SAT path doesn't exist (or
#  is n.s.), the effect is fully mediated through PU.
# ---------------------------------------------------------------------------
cat("\n=== INDIRECT EFFECTS via PU (bootstrap product-of-coefficients) ===\n")

#Provjerava prenose li se učinci kvalitete na zadovoljstvo korisnika (SAT) neizravno,
#preko PU i djelomično preko PEOU.
#zOvaj kod provjerava preko kojih mehanizama kvaliteta utječe na zadovoljstvo.

indirect_boot <- function(draws, path1, path2) {
  # draws columns named "A -> B"; multiply the two legs
  c1 <- draws[, path1, drop = TRUE]
  c2 <- draws[, path2, drop = TRUE]
  c1 * c2
}

quality_dims <- c("INFOQ","SYSQUAL","SERVQUAL")
ind_results  <- lapply(quality_dims, function(q) {
  p1 <- paste0(q,   " -> PU")
  p2 <- "PU -> SAT"
  if (!p1 %in% colnames(boot_draws) || !p2 %in% colnames(boot_draws))
    return(NULL)
  ind <- indirect_boot(boot_draws, p1, p2)
  data.frame(path    = paste0(q, " -> PU -> SAT"),
             indirect = round(mean(ind, na.rm = TRUE), 3),
             ci_lo    = round(quantile(ind, 0.025, na.rm = TRUE), 3),
             ci_hi    = round(quantile(ind, 0.975, na.rm = TRUE), 3),
             sig      = !(quantile(ind, 0.025) <= 0 & quantile(ind, 0.975) >= 0),
             row.names = NULL)
})
ind_df <- do.call(rbind, Filter(Negate(is.null), ind_results))
print(ind_df)
#spituje i utjecaj kvalitete usluge na zadovoljstvo preko jednostavnosti korištenja.
# Also include SERVQUAL → PEOU → SAT chain
p1s <- "SERVQUAL -> PEOU"; p2s <- "PEOU -> SAT"
if (p1s %in% colnames(boot_draws) && p2s %in% colnames(boot_draws)) {
  ind_s  <- indirect_boot(boot_draws, p1s, p2s)
  ind_df <- rbind(ind_df,
                  data.frame(path     = "SERVQUAL -> PEOU -> SAT",
                             indirect = round(mean(ind_s, na.rm = TRUE), 3),
                             ci_lo    = round(quantile(ind_s, 0.025, na.rm = TRUE), 3),
                             ci_hi    = round(quantile(ind_s, 0.975, na.rm = TRUE), 3),
                             sig      = !(quantile(ind_s, 0.025) <= 0 & quantile(ind_s, 0.975) >= 0),
                             row.names = NULL))
}

p_mediation <- ggplot(ind_df, aes(x = indirect, y = reorder(path, indirect),
                                  colour = sig)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.7) +
  scale_colour_manual(values = c(`TRUE` = "#1b7837", `FALSE` = "grey60"),
                      labels = c(`TRUE` = "significant", `FALSE` = "n.s."),
                      name = NULL) +
  labs(title    = "Indirect (mediated) effects on SAT — Model C",
       subtitle = "CI excludes 0 = significant indirect path",
       x = "Indirect effect (95% CI)", y = NULL) +
  theme_mish
ggsave("plots/09_mediation.png", p_mediation, width = 8, height = 5, dpi = 150)

# ---------------------------------------------------------------------------
# 8b) VAF — dodatni medijacijski model s izravnim putanjama prema SAT
# ---------------------------------------------------------------------------

structural_model_VAF <- relationships(
  paths(from = c("INFOQ", "SYSQUAL", "SERVQUAL"), to = "PU"),
  paths(from = "SERVQUAL", to = "PEOU"),
  paths(from = c("PU", "PEOU"), to = "SAT"),
  
  # Direct effects potrebni za VAF
  paths(from = c("INFOQ", "SYSQUAL", "SERVQUAL"), to = "SAT")
)

pls_VAF <- estimate_pls(
  data              = mish,
  measurement_model = measurement_model,
  structural_model  = structural_model_VAF,
  inner_weights     = path_weighting
)

summary(pls_VAF)$paths

# ============================================================
# VAF — Variance Accounted For
# ============================================================

pmVAF <- summary(pls_VAF)$paths

# Izravni učinci na SAT
direct_INFOQ    <- as.numeric(pmVAF["INFOQ", "SAT"])
direct_SYSQUAL  <- as.numeric(pmVAF["SYSQUAL", "SAT"])
direct_SERVQUAL <- as.numeric(pmVAF["SERVQUAL", "SAT"])

# Indirektni učinci
ind_INFOQ <- as.numeric(pmVAF["INFOQ", "PU"]) *
  as.numeric(pmVAF["PU", "SAT"])

ind_SYSQUAL <- as.numeric(pmVAF["SYSQUAL", "PU"]) *
  as.numeric(pmVAF["PU", "SAT"])

# SERVQUAL ima dva indirektna puta prema SAT
ind_SERVQUAL_PU <- as.numeric(pmVAF["SERVQUAL", "PU"]) *
  as.numeric(pmVAF["PU", "SAT"])

ind_SERVQUAL_PEOU <- as.numeric(pmVAF["SERVQUAL", "PEOU"]) *
  as.numeric(pmVAF["PEOU", "SAT"])

ind_SERVQUAL_total <- ind_SERVQUAL_PU + ind_SERVQUAL_PEOU

# Ukupni učinci
total_INFOQ    <- direct_INFOQ + ind_INFOQ
total_SYSQUAL  <- direct_SYSQUAL + ind_SYSQUAL
total_SERVQUAL <- direct_SERVQUAL + ind_SERVQUAL_total

# VAF
VAF_INFOQ    <- ind_INFOQ / total_INFOQ * 100
VAF_SYSQUAL  <- ind_SYSQUAL / total_SYSQUAL * 100
VAF_SERVQUAL <- ind_SERVQUAL_total / total_SERVQUAL * 100

vaf_tbl <- data.frame(
  odnos = c("INFOQ -> SAT",
            "SYSQUAL -> SAT",
            "SERVQUAL -> SAT"),
  direct = c(direct_INFOQ,
             direct_SYSQUAL,
             direct_SERVQUAL),
  indirect = c(ind_INFOQ,
               ind_SYSQUAL,
               ind_SERVQUAL_total),
  total = c(total_INFOQ,
            total_SYSQUAL,
            total_SERVQUAL),
  VAF_percent = c(VAF_INFOQ,
                  VAF_SYSQUAL,
                  VAF_SERVQUAL)
)

cat("\n=== VAF MEDIATION ANALYSIS ===\n")
print(round_num(vaf_tbl, 3))
# ---------------------------------------------------------------------------
# 9)  MODEL C vs MODEL B — SIDE-BY-SIDE COMPARISON
# ---------------------------------------------------------------------------
"Cilj je provjeriti objašnjava li jednostavniji, “trimmed” Model C jednako 
dobro podatke kao složeniji Model B."

pmA  <- s$paths
endo <- intersect(colnames(pmA), colnames(pmB))

r2_cmp <- tibble::tibble(
  construct            = endo,
  `Model C (trimmed)`  = vapply(endo, function(cn) r2_of(pmA, cn), numeric(1)),
  `Model B (full med.)`= vapply(endo, function(cn) r2_of(pmB, cn), numeric(1))
)
"Ako je Δ R² blizu nule, znači da Model C objašnjava gotovo jednako varijance kao Model B.
Ako je Δ R² pozitivan, Model C ima veći R².
Ako je Δ R² negativan, Model B ima veći R²."

"Gubimo li nešto ako koristimo jednostavniji Model C umjesto složenijeg Modela B?"

r2_cmp$`Δ R² (C − B)` <- round(
  r2_cmp$`Model C (trimmed)` - r2_cmp$`Model B (full med.)`, 3)
cat("\n=== R² COMPARISON: Model C vs Model B ===\n")
print(as.data.frame(r2_cmp))

r2_long <- tidyr::pivot_longer(r2_cmp[,1:3], -construct,
                               names_to = "model", values_to = "r2") |>
  dplyr::filter(!is.na(r2))

p_r2_cmp <- ggplot(r2_long, aes(x = construct, y = r2, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.62) +
  geom_text(aes(label = sprintf("%.2f", r2)),
            position = position_dodge(width = 0.7), vjust = -0.4, size = 3.2) +
  scale_fill_manual(values = c("Model C (trimmed)"   = "#3182bd",
                               "Model B (full med.)" = "#fd8d3c"), name = NULL) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.08))) +
  labs(title    = "R² comparison: Model C (trimmed) vs Model B (full mediation)",
       subtitle = "Similar R² = trimmed model loses nothing; better stability",
       x = NULL, y = expression(R^2)) +
  theme_mish

# Paths into SAT — compare bootstrap CIs across models
paths_cmp <- dplyr::bind_rows(
  boot_tbl   |> dplyr::mutate(model = "Model C (trimmed)")   |>
    dplyr::filter(grepl("SAT$", path)),
  boot_tbl_B |> dplyr::mutate(model = "Model B (full med.)") |>
    dplyr::filter(grepl("SAT$", path))
)

p_paths_cmp <- ggplot(paths_cmp,
                      aes(x = original, y = path, colour = model)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi),
                  position = position_dodge(width = 0.55), linewidth = 0.7) +
  scale_colour_manual(values = c("Model C (trimmed)"   = "#3182bd",
                                 "Model B (full med.)" = "#fd8d3c"), name = NULL) +
  labs(title    = "Paths into SAT: Model C vs Model B",
       subtitle = "Narrower CIs in Model C = less multicollinearity",
       x = "Path coefficient (95% CI)", y = NULL) +
  theme_mish

ggsave("plots/10_r2_compare.png",    p_r2_cmp,    width = 8, height = 5, dpi = 150)
ggsave("plots/11_paths_compare.png", p_paths_cmp, width = 8, height = 6, dpi = 150)

# ---------------------------------------------------------------------------
# f² EFFECT SIZE PLOT — MODEL C
# ---------------------------------------------------------------------------

f2_plot_df <- data.frame(
  Path = c(
    "INFOQ → PU",
    "SYSQUAL → PU",
    "SERVQUAL → PU",
    "SERVQUAL → PEOU",
    "PU → SAT",
    "PEOU → SAT"
  ),
  f2 = c(
    0.406,
    0.665,
    0.235,
    1.075,
    1.456,
    0.232
  )
)

# Klasifikacija veličine učinka
f2_plot_df$Effect <- cut(
  f2_plot_df$f2,
  breaks = c(-Inf, 0.02, 0.15, 0.35, Inf),
  labels = c(
    "Zanemariv",
    "Mali",
    "Srednji",
    "Veliki"
  )
)

# Poredaj putanje prema veličini f²
f2_plot_df$Path <- reorder(
  f2_plot_df$Path,
  f2_plot_df$f2
)

# Graf
p_f2 <- ggplot(
  f2_plot_df,
  aes(x = f2, y = Path)
) +
  geom_col(width = 0.65) +
  
  # Referentne granice za f²
  geom_vline(
    xintercept = c(0.02, 0.15, 0.35),
    linetype = "dashed",
    linewidth = 0.5
  ) +
  
  # Vrijednost f² uz svaki stupac
  geom_text(
    aes(label = sprintf("%.3f", f2)),
    hjust = -0.15,
    size = 3.8
  ) +
  
  scale_x_continuous(
    limits = c(0, 1.65),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    title = "Veličine učinaka (f²) — Model C",
    subtitle = "Referentne vrijednosti: 0,02 = mali; 0,15 = srednji; 0,35 = veliki učinak",
    x = expression(f^2),
    y = NULL
  ) +
  
  theme_mish

# Spremi graf
ggsave(
  "plots/12_f2_effect_sizes.png",
  p_f2,
  width = 9,
  height = 5.5,
  dpi = 300
)
# ---------------------------------------------------------------------------
# 10)  DISPLAY ALL PLOTS
# ---------------------------------------------------------------------------
print(p_likert)
print(p_means)
print(p_path)
print(p_forest)
print(p_r2)
print(p_htmt)
print(p_ipma)
print(nca_plots[["INFOQ"]])
print(p_mediation)
print(p_r2_cmp)
print(p_paths_cmp)
print(p_f2) 

cat("\nDone. All figures saved to ./plots/\n",
    "Model C (trimmed, primary):  structural_model  in section 3b\n",
    "Model B (full med, compare): structural_model_B in section 7\n")
# ============================================================================