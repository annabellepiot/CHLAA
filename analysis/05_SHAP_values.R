# =========================================================================
# Shapley decomposition of intervention contributions to averted burden
# =========================================================================
#
# For each health zone and each posterior draw, builds the fully-on fitted
# parameter set (all six levers at their historically observed windows,
# effects, capacities and vaccine schedule), simulates all 2^6 = 64 on/off
# combinations under ONE shared seed (common random numbers), and computes
# the averted burden V(S) = burden(none) - burden(S) for every combination.
# The Shapley value for each lever is assembled from that lookup table
# (no re-simulation per lever/subset pair, so there is no hidden
# approximation). Contributions sum EXACTLY to the total averted burden
# (the efficiency axiom, asserted below via verdict()); leave-one-out and
# add-one-in are reported as the two bounds Shapley averages over, plus a
# redundancy/synergy diagnostic (sum(add-one-in) vs total vs
# sum(leave-one-out)).
#
# The empty set is "no response"; the full set is the real deployed
# programme (inherited unmodified from the fit) - so this partitions the
# burden the actual programme averted, at its ACTUAL observed timing (not
# a re-optimised schedule). A lever not deployed in a given zone (e.g.
# vaccination where no campaign ran) correctly receives ~0.
#
# Both outcomes (cases, deaths) are derived from the SAME 64 x N_DRAWS
# simulations per zone: chlaa_simulate() already returns cum_symptoms and
# cum_deaths together in one call, so reading off both columns costs
# nothing extra and guarantees identical CRN trajectories underlie both
# decompositions.
#
# Horizon (HORIZON_EXTRA = 182, weekly grid) and dt (0.25) match
# 02_02_scenario_analysis_all_HZs.R so this analysis is on the same
# footing as the scenario / intervention-contribution figures. Negative
# Shapley contributions are legitimate (not clipped) - cholera control
# levers are frequently substitutes, and a lever can look harmful in
# combination while helping alone, or vice versa.
#
# Designed for PBS array job submission (one HZ per job, compute + cache
# only) or sequential interactive use (also rebuilds the aggregate
# figures at the end) - mirrors the pattern of 02_02_scenario_analysis_all_HZs.R.
#
# =========================================================================

source("00_setup.R") #/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA/analysis/00_setup.R

suppressMessages({
    library(ggplot2)
    library(dplyr)
})

dir.create(RDS_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TAB_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

## ---- Fixed method settings -----------------------------------------
N_DRAWS <- 60 # posterior draws -> parameter-uncertainty CrIs
N_PART <- 30 # particles per simulation (averaged within a draw)
BURNIN <- 0.25
HORIZON_EXTRA <- 182 # days beyond last observed week (matches 02_02)
DT <- 0.25 # matches 02_02_scenario_analysis_all_HZs.R
N_THREADS <- 4 # matches the ncpus=4 PBS resource request
SEED <- 1
OUTCOMES <- c("cases", "deaths")

INTS <- c("CTC", "ORC", "CATI", "Hygiene", "Chlorination", "Vax1")
k <- length(INTS)

intervention_labels <- c(
    CTC = "CTC", ORC = "ORC", CATI = "CATI", Hygiene = "Hygiene",
    Chlorination = "Chlorination", Vax1 = "Vaccination"
)

# Display-friendly health zone names (e.g. "ngiri_ngiri" -> "Ngiri Ngiri").
# Same helper as in 02_02_scenario_analysis_all_HZs.R.
hz_label <- function(x) {
    x <- gsub("_", " ", x)
    gsub("(?<=^|\\s)([a-z])", "\\U\\1", x, perl = TRUE)
}

# For each lever, the parameter overrides that switch it OFF. Everything not
# listed keeps its HISTORICALLY OBSERVED value carried in the fitted pars, so
# the FULL set = the real deployed programme and the EMPTY set = no response.
off_overrides <- list(
    CTC          = list(ctc_start = 0,   ctc_end = 0,   ctc_capacity = 0),
    ORC          = list(orc_start = 0,   orc_end = 0,   orc_capacity = 0),
    CATI         = list(cati_start = 0,  cati_end = 0,  cati_effect = 0),
    Hygiene      = list(hyg_start = 0,   hyg_end = 0,   hyg_effect = 0),
    Chlorination = list(chlor_start = 0, chlor_end = 0, chlor_effect = 0),
    Vax1         = list(vax1_start = 0,  vax1_end = 0,  vax1_total_doses = 0,
                        vax1_schedule_time = c(0L, 1L),
                        vax1_schedule_doses = c(0, 0),
                        n_vax1_schedule = 2L)
)

## ---- Helpers -------------------------------------------------------

# Build the parameter list for a subset S (a logical vector over INTS):
# start from the fully-on fitted pars for this draw, then switch OFF every
# lever NOT in S.
make_subset_pars <- function(full_pars, S) {
    p <- full_pars
    for (j in which(!S)) p <- utils::modifyList(p, off_overrides[[INTS[j]]])
    chlaa_parameters_validate(p)
    p
}

# Burden at end of horizon for a config, averaged over CRN particles, for
# BOTH outcomes at once (one chlaa_simulate() call).
burden_multi <- function(pars, tvec, seed) {
    s <- chlaa_simulate(
        pars = pars, time = tvec, n_particles = N_PART, dt = DT,
        seed = seed, n_threads = N_THREADS, deterministic = FALSE
    ) # <-- CRN: same seed for all 64 subsets in a draw
    end <- s[s$time == max(s$time), ]
    c(cases = mean(end$cum_symptoms), deaths = mean(end$cum_deaths))
}

# Shapley weight for a subset of size s (excluding lever i), k levers total.
shap_w <- function(s, k) factorial(s) * factorial(k - s - 1) / factorial(k)

# Quantiles reported: 2.5/25/50/75/97.5%, mirroring the box+whisker
# convention used throughout this project's forecast figures.
q5 <- function(x) stats::quantile(x, c(.025, .25, .5, .75, .975), na.rm = TRUE)

## ---- Enumerate all 2^k subsets once (as bitmasks) ------------------
masks <- 0:(2^k - 1)
as_S <- function(m) as.logical(bitwAnd(m, 2^(0:(k - 1))) > 0) # length-k logical
full_mask <- 2^k - 1

## ---- Main per-zone computation (cached) -----------------------------

compute_shapley_hz <- function(hz_name, n_draws = N_DRAWS, n_part = N_PART, burnin = BURNIN,
                                seed = SEED, horizon_extra = HORIZON_EXTRA, dt = DT,
                                force = FALSE, verbose = TRUE) {
    out_path <- file.path(RDS_DIR, sprintf("%s_shapley.rds", hz_name))
    if (file.exists(out_path) && !force) {
        if (verbose) cat("Using cached Shapley results for", hz_name, "\n")
        return(readRDS(out_path))
    }

    if (verbose) cat("\n--- Shapley decomposition:", hz_name, "---\n")

    fo <- load_fit_obj(hz_name)
    fit <- fo$fit
    tvec <- seq(7, max(fo$observed$time) + horizon_extra, by = 7)

    # Per-chain burn-in then combine (NOT chlaa_fit_draws() +
    # chlaa_fit_select_iterations(), which stacks all chains BEFORE burn-in
    # and so only discards the first `burnin` fraction of the concatenated
    # stack - i.e. effectively burns in just the first chain and keeps 100%
    # of every later one). This is the same internal helper
    # chlaa_forecast_scenarios_from_fit() itself uses.
    dr <- chlaa:::.chlaa_fit_selected_draws_matrix(fit, burnin = burnin, thin = 1)
    set.seed(seed)
    draw_idx <- sample.int(nrow(dr), n_draws, replace = n_draws > nrow(dr))

    phi_arr <- loo_arr <- aoi_arr <- array(
        NA_real_,
        dim = c(n_draws, k, length(OUTCOMES)),
        dimnames = list(NULL, INTS, OUTCOMES)
    )
    tot_mat <- matrix(NA_real_, n_draws, length(OUTCOMES), dimnames = list(NULL, OUTCOMES))

    for (d in seq_len(n_draws)) {
        seed_d <- 10000 + d # CRN seed shared by all 64 configs this draw
        full_pars <- chlaa:::.chlaa_update_pars_from_theta(dr[draw_idx[d], ], fo$pars_warm, fit)

        # 1) burden for every subset (both outcomes at once)
        b <- matrix(NA_real_, length(masks), length(OUTCOMES), dimnames = list(as.character(masks), OUTCOMES))
        for (m in masks) {
            b[as.character(m), ] <- burden_multi(make_subset_pars(full_pars, as_S(m)), tvec, seed_d)
        }
        # V(S) = burden(empty) - burden(S) (averted), per outcome column
        V <- b
        for (o in OUTCOMES) V[, o] <- b["0", o] - b[, o]
        tot_mat[d, ] <- V[as.character(full_mask), ]

        # 2) assemble Shapley + the two bounds from the V-lookup (no extra sims)
        for (i in seq_len(k)) {
            bit_i <- 2^(i - 1)
            others <- setdiff(masks, masks[bitwAnd(masks, bit_i) > 0]) # subsets WITHOUT i
            for (o in OUTCOMES) {
                phi <- 0
                for (m in others) {
                    s <- sum(as_S(m))
                    phi <- phi + shap_w(s, k) * (V[as.character(m + bit_i), o] - V[as.character(m), o])
                }
                phi_arr[d, i, o] <- phi
                aoi_arr[d, i, o] <- V[as.character(bit_i), o] # add-one-in (marginal first)
                loo_arr[d, i, o] <- V[as.character(full_mask), o] - V[as.character(full_mask - bit_i), o] # leave-one-out (marginal last)
            }
        }
    }

    ## ---- Summarise per outcome ---------------------------------------
    summarise_outcome <- function(o) {
        phi_mat <- phi_arr[, , o]
        loo_mat <- loo_arr[, , o]
        aoi_mat <- aoi_arr[, , o]
        tot_vec <- tot_mat[, o]

        # Per-draw shares: sum_i phi_mat[d, i] == tot_vec[d] EXACTLY for every
        # draw d (the efficiency axiom holds per draw), so dividing row-wise
        # before summarising gives shares that sum to exactly 100% for every
        # draw, and a properly propagated CI on the share itself - unlike
        # normalising by the sum of independently-computed per-lever medians
        # (which are not additive).
        share_mat <- 100 * sweep(phi_mat, 1, tot_vec, "/")

        phi_q <- t(apply(phi_mat, 2, q5))
        share_q <- t(apply(share_mat, 2, q5))

        res <- data.frame(
            intervention      = INTS,
            shapley_q0p025    = phi_q[, 1], shapley_q0p25 = phi_q[, 2], shapley_q0p5 = phi_q[, 3],
            shapley_q0p75     = phi_q[, 4], shapley_q0p975 = phi_q[, 5],
            share_pct_q0p025  = share_q[, 1], share_pct_q0p25 = share_q[, 2], share_pct_q0p5 = share_q[, 3],
            share_pct_q0p75   = share_q[, 4], share_pct_q0p975 = share_q[, 5],
            loo_med           = apply(loo_mat, 2, stats::median, na.rm = TRUE),
            aoi_med           = apply(aoi_mat, 2, stats::median, na.rm = TRUE)
        )

        sum_phi <- mean(rowSums(phi_mat))
        total <- mean(tot_vec)
        eff_pass <- abs(sum_phi - total) < 1e-6 * max(1, abs(total))
        verdict(eff_pass, sprintf(
            "[%s/%s] Shapley values sum to the total averted burden (%.1f vs %.1f).",
            hz_name, o, sum_phi, total
        ))

        sum_aoi <- mean(rowSums(aoi_mat))
        sum_loo <- mean(rowSums(loo_mat))
        interaction_note <- if (sum_aoi > total && total > sum_loo) {
            "substitutes (redundant levers - strong alone, less so together)"
        } else if (sum_loo > total && total > sum_aoi) {
            "complements (synergistic levers)"
        } else {
            "mixed substitute/complement pattern"
        }
        cat(sprintf(
            "[%s/%s] sum(add-one-in)=%.1f  >=  total=%.1f  >=  sum(leave-one-out)=%.1f  ->  %s\n",
            hz_name, o, sum_aoi, total, sum_loo, interaction_note
        ))

        list(
            table = res, phi_mat = phi_mat, loo_mat = loo_mat, aoi_mat = aoi_mat,
            share_mat = share_mat, tot_vec = tot_vec,
            efficiency_pass = eff_pass, efficiency_gap = sum_phi - total,
            sum_aoi = sum_aoi, sum_loo = sum_loo, total = total,
            interaction_note = interaction_note
        )
    }

    outcome_results <- setNames(lapply(OUTCOMES, summarise_outcome), OUTCOMES)

    for (o in OUTCOMES) {
        write.csv(
            outcome_results[[o]]$table,
            file.path(TAB_DIR, sprintf("shapley_%s_%s.csv", hz_name, o)),
            row.names = FALSE
        )
    }

    result <- list(
        hz_name = hz_name, outcomes = outcome_results,
        n_draws = n_draws, n_part = n_part, burnin = burnin, seed = seed,
        dt = dt, horizon_extra = horizon_extra, timestamp = Sys.time()
    )
    saveRDS(result, out_path)
    if (verbose) cat("Saved:", out_path, "\n")
    result
}

## ---- Main Execution --------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
    # Array job mode: run single HZ specified by argument
    hz_to_run <- args[1]

    cat("\n", rep("=", 60), "\n", sep = "")
    cat("Running Shapley decomposition for:", hz_to_run, "\n")
    cat(rep("=", 60), "\n\n", sep = "")

    result <- tryCatch(
        {
            compute_shapley_hz(hz_to_run, verbose = TRUE)
        },
        error = function(e) {
            cat("\nERROR in Shapley decomposition for", hz_to_run, ":\n")
            cat(conditionMessage(e), "\n")
            saveRDS(list(
                hz_name = hz_to_run,
                error = conditionMessage(e),
                timestamp = Sys.time()
            ), file.path(RDS_DIR, sprintf("%s_shapley_FAILED.rds", hz_to_run)))
            NULL
        }
    )

    if (!is.null(result)) {
        cat("\nSUCCESS: Completed Shapley decomposition for", hz_to_run, "\n")
    } else {
        cat("\nFAILURE: Could not complete Shapley decomposition for", hz_to_run, "\n")
        quit(status = 1)
    }
} else {
    # Interactive mode: run all HZs sequentially (cache-backed, so already
    # computed zones return instantly)
    cat("No HZ specified. Running in sequential mode for all HZs.\n")
    cat("For production use, submit as array job with HZ name as argument.\n\n")

    all_fits <- list.files(RDS_DIR, pattern = "_fit\\.rds$", full.names = FALSE)
    all_hzs <- sort(gsub("_fit\\.rds$", "", all_fits[!grepl("comparative", all_fits)]))

    failed_fits <- list.files(RDS_DIR, pattern = "_shapley_FAILED\\.rds$", full.names = FALSE)
    failed_hzs <- gsub("_shapley_FAILED\\.rds$", "", failed_fits)
    all_hzs <- setdiff(all_hzs, failed_hzs)

    if (length(all_hzs) == 0) stop("No completed fit artifacts found in: ", RDS_DIR)

    cat("Health zones with completed fits:", paste(all_hzs, collapse = ", "), "\n\n")

    results <- list()
    for (hz_name in all_hzs) {
        results[[hz_name]] <- tryCatch(
            {
                compute_shapley_hz(hz_name, verbose = TRUE)
            },
            error = function(e) {
                cat("\nERROR in Shapley decomposition for", hz_name, ":\n")
                cat(conditionMessage(e), "\n\n")
                NULL
            }
        )
    }

    n_success <- sum(sapply(results, function(x) !is.null(x)))
    n_failed <- length(results) - n_success

    cat("\n", rep("=", 60), "\n", sep = "")
    cat("SHAPLEY ANALYSIS SUMMARY\n")
    cat(rep("=", 60), "\n", sep = "")
    cat("Total HZs:", length(results), "\n")
    cat("Successful:", n_success, "\n")
    cat("Failed:", n_failed, "\n")
    if (n_failed > 0) {
        failed_now <- names(results)[sapply(results, is.null)]
        cat("Failed HZs:", paste(failed_now, collapse = ", "), "\n")
    }
}

## =========================================================================
## Aggregate Shapley figures across all health zones
## =========================================================================
#
# Runs every time this script is invoked (array task or interactive), same
# as the composite/contribution figures in 02_02_scenario_analysis_all_HZs.R
# - cheap once all zones are cached, and the last array task to finish
# leaves complete figures behind.

cat("\n", rep("=", 60), "\n", sep = "")
cat("Building Shapley relative-contribution figures\n")
cat(rep("=", 60), "\n\n", sep = "")

shapley_rds_files <- list.files(RDS_DIR, pattern = "_shapley\\.rds$", full.names = TRUE)
shapley_rds_files <- shapley_rds_files[!grepl("FAILED", shapley_rds_files)]
shapley_objs <- lapply(shapley_rds_files, readRDS)
names(shapley_objs) <- vapply(shapley_objs, `[[`, character(1), "hz_name")

if (length(shapley_objs) == 0) {
    cat("No Shapley results available - skipping figures.\n")
} else {
    shap_dat <- lapply(names(shapley_objs), function(hz) {
        obj <- shapley_objs[[hz]]
        lapply(OUTCOMES, function(o) {
            tbl <- obj$outcomes[[o]]$table
            tbl$hz <- hz
            tbl$variable <- o
            tbl$interaction_note <- obj$outcomes[[o]]$interaction_note
            tbl
        }) %>% bind_rows()
    }) %>% bind_rows()

    shap_dat <- shap_dat %>%
        mutate(
            intervention = factor(intervention_labels[intervention], levels = unname(intervention_labels[INTS])),
            variable = factor(variable, levels = OUTCOMES),
            num_label = paste0(round(share_pct_q0p5), "% (", round(share_pct_q0p025), " to ", round(share_pct_q0p975), ")"),
            # Anchored at the 25/75% box edge (not the 95% whisker end), same
            # reasoning as 02_02's intervention-contribution figure: a few
            # small-population HZs have wide enough 95% intervals to overflow
            # any sensible shared x-axis.
            label_x = ifelse(share_pct_q0p5 >= 0, share_pct_q0p75, share_pct_q0p25),
            label_hjust = ifelse(share_pct_q0p5 >= 0, -0.08, 1.08)
        )

    write.csv(shap_dat, file.path(TAB_DIR, "shapley_all_hz.csv"), row.names = FALSE)

    # HZ ordering: reuse the same "no-intervention excess" magnitude ranking
    # as the composite scenario / intervention-contribution figures, so all
    # of this project's figures are directly comparable panel-for-panel.
    hz_rank_all <- tryCatch(
        {
            scenario_rds_files4 <- list.files(RDS_DIR, pattern = "_scenarios\\.rds$", full.names = TRUE)
            scenario_rds_files4 <- scenario_rds_files4[!grepl("FAILED", scenario_rds_files4)]
            scenario_objs4 <- lapply(scenario_rds_files4, readRDS)
            names(scenario_objs4) <- vapply(scenario_objs4, `[[`, character(1), "hz_name")

            denom_dat4 <- lapply(names(scenario_objs4), function(hz) {
                fc <- scenario_objs4[[hz]]$scenario_forecasts
                snap_time <- fc$time[which.min(abs(unique(fc$time) - 365))]
                fc %>%
                    filter(
                        type == "difference", scenario == "no_interventions",
                        variable == "cum_symptoms", time == snap_time
                    ) %>%
                    transmute(hz = hz, denom = q0p5)
            }) %>% bind_rows()

            r <- denom_dat4 %>% arrange(denom) %>% pull(hz)
            c(intersect(r, names(shapley_objs)), setdiff(names(shapley_objs), r))
        },
        error = function(e) sort(names(shapley_objs))
    )

    shap_dat <- shap_dat %>% mutate(hz = factor(hz, levels = hz_rank_all))

    caption_txt_shapley <- paste(
        "Shapley decomposition: exact partition of the burden averted by the full historical response",
        "(all 6 levers at fitted timing/effects) vs. no response, using common random numbers\n",
        "across all 64 lever on/off combinations per posterior draw. Shares sum to exactly 100% for every draw.\n",
        "KNOWN ISSUE: Vaccination currently shows ~0% for every zone due to a bug in how vax1_start/vax1_end/vax1_total_doses were computed at the fitting stage\n",
        "(fixed in 01_02_fitting_all_HZs.R but not yet reflected in the saved fits) -",
        "re-run the fit for the 7 affected zones before trusting this lever."
    )

    ## ---- Forest-plot figure (uncertainty, one facet per intervention) ----

    build_shapley_forest_plot <- function(var_name, plot_title, xlim_clip = NULL) {
        dat <- shap_dat %>% filter(variable == var_name)

        n_clipped <- 0
        if (!is.null(xlim_clip)) {
            n_clipped <- sum(dat$share_pct_q0p025 < xlim_clip[1] | dat$share_pct_q0p975 > xlim_clip[2], na.rm = TRUE)
            x_pad <- data.frame(
                x = xlim_clip,
                hz = factor(hz_rank_all[1], levels = hz_rank_all),
                intervention = factor(levels(dat$intervention)[1], levels = levels(dat$intervention))
            )
            inner_margin <- diff(xlim_clip) * 0.02
            dat <- dat %>%
                mutate(
                    label_x_clamped = ifelse(
                        share_pct_q0p5 >= 0,
                        pmin(label_x, xlim_clip[2] - inner_margin),
                        pmax(label_x, xlim_clip[1] + inner_margin)
                    ),
                    was_clamped = label_x_clamped != label_x,
                    label_hjust = ifelse(was_clamped, -label_hjust + 1, label_hjust),
                    label_x = label_x_clamped
                )
        } else {
            x_range <- range(c(dat$share_pct_q0p025, dat$share_pct_q0p975), na.rm = TRUE)
            x_pad_amt <- diff(x_range) * 0.4
            x_pad <- data.frame(
                x = c(x_range[1] - x_pad_amt, x_range[2] + x_pad_amt),
                hz = factor(hz_rank_all[1], levels = hz_rank_all),
                intervention = factor(levels(dat$intervention)[1], levels = levels(dat$intervention))
            )
        }

        caption_txt <- caption_txt_shapley
        if (n_clipped > 0) {
            caption_txt <- paste0(
                caption_txt, sprintf(
                    " %d interval(s) extend beyond the axis range shown (x-axis clipped to %g to %g%% for readability).",
                    n_clipped, xlim_clip[1], xlim_clip[2]
                )
            )
        }

        p <- ggplot(dat, aes(y = hz)) +
            geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
            geom_blank(data = x_pad, aes(x = x, y = hz)) +
            geom_crossbar(
                aes(x = share_pct_q0p5, xmin = share_pct_q0p25, xmax = share_pct_q0p75, fill = intervention),
                orientation = "y", width = 0.65, alpha = 0.8, colour = "grey30",
                linewidth = 0.3, middle.linewidth = 0.6
            ) +
            geom_errorbar(
                aes(xmin = share_pct_q0p025, xmax = share_pct_q0p975),
                orientation = "y", width = 0.35, linewidth = 0.4, colour = "grey30"
            ) +
            geom_label(
                aes(x = label_x, label = num_label, hjust = label_hjust),
                size = 2.3, family = "Helvetica", colour = "black",
                fill = "white", linewidth = 0, label.padding = unit(0.12, "lines")
            ) +
            facet_wrap(~intervention, ncol = 2) +
            scale_y_discrete(labels = hz_label) +
            scale_fill_brewer(palette = "Set2", guide = "none") +
            labs(
                x = "Shapley share of total averted burden (%)",
                y = NULL,
                title = plot_title,
                subtitle = "Exact partition of averted burden across levers (efficiency axiom); shares sum to 100% for every posterior draw",
                caption = caption_txt
            )

        if (!is.null(xlim_clip)) p <- p + coord_cartesian(xlim = xlim_clip)

        p +
            theme_minimal(base_family = "Helvetica", base_size = 12) +
            theme(
                panel.grid       = element_blank(),
                panel.background = element_rect(fill = "white", colour = "grey70"),
                panel.border     = element_rect(fill = NA, colour = "grey70", linewidth = 0.5),
                plot.background  = element_rect(fill = "white", colour = NA),
                strip.text       = element_text(face = "bold", size = 11),
                axis.ticks.x     = element_line(colour = "grey40"),
                axis.ticks.y     = element_blank(),
                plot.title       = element_text(face = "bold", size = 15),
                plot.subtitle    = element_text(size = 9.5, colour = "grey40"),
                plot.caption     = element_text(size = 7.5, colour = "grey40", hjust = 0),
                panel.spacing    = unit(1, "lines")
            )
    }

    ## ---- Stacked-bar figure (medians only, visualises the 100% partition) ----

    build_shapley_stacked_plot <- function(var_name, plot_title) {
        dat <- shap_dat %>% filter(variable == var_name)

        stack_caption <- paste0(
            caption_txt_shapley, "\n",
            "Segments show per-draw-median shares only (medians are not strictly additive across levers, so per-zone\n",
            "totals may deviate slightly from 100%; see shapley_all_hz.csv for exact per-draw sums, and the companion forest-plot figure for uncertainty)."
        )

        ggplot(dat, aes(x = hz, y = share_pct_q0p5, fill = intervention)) +
            geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
            geom_hline(yintercept = 100, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
            coord_flip() +
            scale_x_discrete(labels = hz_label) +
            scale_fill_brewer(palette = "Set2", name = "Intervention") +
            labs(
                x = NULL, y = "Median Shapley share of total averted burden (%)",
                title = plot_title,
                subtitle = "Shapley shares sum to exactly 100% of the total averted burden for every posterior draw",
                caption = stack_caption
            ) +
            theme_minimal(base_family = "Helvetica", base_size = 12) +
            theme(
                panel.grid       = element_blank(),
                panel.background = element_rect(fill = "white", colour = "grey70"),
                panel.border     = element_rect(fill = NA, colour = "grey70", linewidth = 0.5),
                plot.background  = element_rect(fill = "white", colour = NA),
                axis.ticks.x     = element_line(colour = "grey40"),
                axis.ticks.y     = element_blank(),
                plot.title       = element_text(face = "bold", size = 15),
                plot.subtitle    = element_text(size = 9.5, colour = "grey40"),
                plot.caption     = element_text(size = 7.5, colour = "grey40", hjust = 0),
                legend.position  = "bottom"
            )
    }

    p_forest_cases <- build_shapley_forest_plot("cases", "Shapley contribution of individual interventions - Cases")
    ggsave(file.path(FIG_DIR, "shapley_forest_cases.png"), p_forest_cases, width = 12, height = 10, dpi = 300)
    ggsave(file.path(FIG_DIR, "shapley_forest_cases.pdf"), p_forest_cases, width = 12, height = 10)

    p_forest_deaths <- build_shapley_forest_plot("deaths", "Shapley contribution of individual interventions - Deaths")
    ggsave(file.path(FIG_DIR, "shapley_forest_deaths.png"), p_forest_deaths, width = 12, height = 10, dpi = 300)
    ggsave(file.path(FIG_DIR, "shapley_forest_deaths.pdf"), p_forest_deaths, width = 12, height = 10)

    p_stacked_cases <- build_shapley_stacked_plot("cases", "Relative contribution (Shapley) by health zone - Cases")
    ggsave(file.path(FIG_DIR, "shapley_stacked_cases.png"), p_stacked_cases, width = 12, height = 9, dpi = 300)
    ggsave(file.path(FIG_DIR, "shapley_stacked_cases.pdf"), p_stacked_cases, width = 12, height = 9)

    p_stacked_deaths <- build_shapley_stacked_plot("deaths", "Relative contribution (Shapley) by health zone - Deaths")
    ggsave(file.path(FIG_DIR, "shapley_stacked_deaths.png"), p_stacked_deaths, width = 12, height = 9, dpi = 300)
    ggsave(file.path(FIG_DIR, "shapley_stacked_deaths.pdf"), p_stacked_deaths, width = 12, height = 9)

    cat("\nShapley figures saved to:\n")
    cat("  ", file.path(FIG_DIR, "shapley_forest_cases.png"), "\n")
    cat("  ", file.path(FIG_DIR, "shapley_forest_deaths.png"), "\n")
    cat("  ", file.path(FIG_DIR, "shapley_stacked_cases.png"), "\n")
    cat("  ", file.path(FIG_DIR, "shapley_stacked_deaths.png"), "\n")
}
