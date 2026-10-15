# Does an intervention change the trajectory of weekly cholera cases?
# Hierarchical GAM (partial pooling across HZs; Pedersen et al. 2019, model "GS"),
# one model per intervention group x lag:
#
#   cases ~ s(time_index)                              global baseline epidemic curve
#         + s(time_index, hz, bs = "fs")               HZ deviations from the baseline
#         + s(time_index, by = intervention)           global intervention effect f(t)
#         + s(time_index, hz, by = intervention, "fs") HZ deviations from f(t)
#   family = nb(), method = "REML"
#
# intervention = group intensity c_g(t) in [0, 1] from intervention_intensity.R (run that first).
# f(t) is on the log scale: the log rate ratio of cases at full intensity (1) vs none (0)
# in week t. f(t) below 0 = fewer cases when the intervention is running.
#
# Inputs / assumptions:
# - Missing case weeks are left out (not set to 0): mgcv does not need consecutive rows.
# - Intervention intensity is 0 before the first activity and after the last report.
# - Lag 0: intervention in week t acting on cases in week t. Lag 1: intervention in week t-1.
# - The HZ-deviation term for the intervention is only included when at least
#   min_hz_term4 HZs have at least min_wk_term4 active weeks; groups with fewer than
#   min_hz_fit active HZs are not fitted.
# - Within one HZ, its baseline deviation and its intervention deviation are both
#   functions of time only, so they are weakly identified. The global effect is the
#   result to interpret; HZ curves are shrunk deviations, not separate HZ effects.
#
# Outputs:
#   analysis/data/intervention_effect_gam.csv           one row per group x lag x smooth term
#   figures/gam_effect_<group>.png                      global f(t), lag 0 and lag 1 panels
#   figures/gam_effect_<group>_hz.png                   global + HZ deviation per HZ (term-4 groups)
#   output/intervention_effect_gam/gam_check_<group>_lag<l>.txt / .png

library(dplyr)
library(tidyr)
library(ggplot2)
library(mgcv)

data_dir <- "Z:/acp25_DIDE/CHLAA/analysis/data"
fig_dir <- "Z:/acp25_DIDE/CHLAA/figures"
chk_dir <- "Z:/acp25_DIDE/CHLAA/output/intervention_effect_gam"
dir.create(chk_dir, showWarnings = FALSE, recursive = TRUE)
week1_date <- as.Date("2024-12-30") # IDSR date of 2025 week 1

lags <- c(0, 1)
lag_labels <- c("0" = "Lag 0: intervention in same week", "1" = "Lag 1: intervention in previous week")
k_fs <- 5 # basis size per HZ in the factor-smooth terms
min_hz_fit <- 3 # fit a group only if at least this many HZs are active
min_hz_term4 <- 5 # include HZ deviations of the effect if at least this many HZs...
min_wk_term4 <- 8 # ...have at least this many active weeks

group_colours <- c(
    ctc = "#00BFC4", cati = "#7B68EE", hyg = "#FF69B4",
    chlor = "#32CD32", lat = "#f86d6d"
)

term_labels <- c(
    "s(time_index)" = "baseline (global)",
    "s(time_index,hz)" = "baseline (HZ deviations)",
    "s(time_index):intervention" = "intervention (global)",
    "s(time_index,hz):intervention" = "intervention (HZ deviations)"
)

# ---- data ----
cases <- read.csv(file.path(data_dir, "IDSR_AA_indicators.csv"), stringsAsFactors = FALSE) %>%
    filter(!is.na(cases)) %>%
    select(hz, week, cases)
hz_levels <- sort(unique(cases$hz))

intensity <- read.csv(file.path(data_dir, "intervention_intensity.csv"), stringsAsFactors = FALSE) %>%
    pivot_longer(c(intensity, intensity_lag1), names_to = "lag", values_to = "intervention") %>%
    mutate(lag = ifelse(lag == "intensity", 0, 1)) %>%
    complete(hz = hz_levels, week = 1:52, group, lag, fill = list(intervention = 0, active = FALSE))

active_wk <- intensity %>%
    filter(lag == 0) %>%
    group_by(group, hz) %>%
    summarise(active_weeks = sum(active), .groups = "drop")

group_plan <- active_wk %>%
    group_by(group) %>%
    summarise(
        hz_active = sum(active_weeks > 0),
        hz_term4 = sum(active_weeks >= min_wk_term4),
        .groups = "drop"
    ) %>%
    mutate(model = case_when(
        hz_active < min_hz_fit ~ "not fitted",
        hz_term4 >= min_hz_term4 ~ "full (terms 1-4)",
        TRUE ~ "global effect only (terms 1-3)"
    ))
cat("Model per group (active = any indicator reported > 0):\n")
print(as.data.frame(group_plan), row.names = FALSE)

# ---- helpers ----

# fitted smooth (sum of the named terms) with 95% CI on newdata, via the lpmatrix
term_curve <- function(fit, nd, terms) {
    X <- predict(fit, nd, type = "lpmatrix")
    idx <- unlist(lapply(fit$smooth, function(s) {
        if (s$label %in% terms) s$first.para:s$last.para
    }))
    Xs <- X[, idx, drop = FALSE]
    est <- drop(Xs %*% coef(fit)[idx])
    se <- sqrt(rowSums((Xs %*% fit$Vp[idx, idx]) * Xs))
    data.frame(nd, est = est, lower = est - 1.96 * se, upper = est + 1.96 * se)
}

fit_one <- function(g, l, full) {
    d <- cases %>%
        inner_join(intensity %>% filter(group == g, lag == l) %>% select(hz, week, intervention),
            by = c("hz", "week")
        ) %>%
        transmute(hz = factor(hz, levels = hz_levels), time_index = week, cases, intervention)
    f <- cases ~ s(time_index, bs = "tp") + s(time_index, hz, bs = "fs", k = k_fs) +
        s(time_index, by = intervention, bs = "tp")
    if (full) f <- update(f, . ~ . + s(time_index, hz, by = intervention, bs = "fs", k = k_fs))
    fit <- withCallingHandlers(
        gam(f, family = nb(), data = d, method = "REML"),
        warning = function(w) {
            if (grepl("repeated 1-d smooths", conditionMessage(w))) invokeRestart("muffleWarning")
        }
    )
    list(fit = fit, data = d)
}

# ---- fit all groups x lags ----
fits <- list()
rows <- list()
for (g in group_plan$group[group_plan$model != "not fitted"]) {
    full <- group_plan$model[group_plan$group == g] == "full (terms 1-4)"
    for (l in lags) {
        key <- paste0(g, "_lag", l)
        cat("\nFitting", g, "lag", l, if (full) "(full)" else "(global effect only)", "\n")
        res <- tryCatch(fit_one(g, l, full), error = function(e) {
            cat("  failed:", conditionMessage(e), "\n")
            NULL
        })
        if (is.null(res)) next
        fits[[key]] <- res
        fit <- res$fit

        # gam.check: text + residual plots
        png(file.path(chk_dir, paste0("gam_check_", key, ".png")), width = 1600, height = 1600, res = 200)
        chk <- capture.output(gam.check(fit))
        dev.off()
        writeLines(c(capture.output(summary(fit)), "", chk), file.path(chk_dir, paste0("gam_check_", key, ".txt")))

        # missing case weeks while the intervention was running (any HZ)
        miss <- intensity %>%
            filter(group == g, lag == l, intervention > 0) %>%
            anti_join(cases, by = c("hz", "week")) %>%
            nrow()
        # concurvity of each term with the terms of the other component (baseline vs
        # intervention): near 1 = intervention and baseline curves cannot be separated.
        # Global vs HZ-deviation terms of the same component overlap by construction.
        # Fails when HZs without the intervention leave all-zero columns in term 4.
        labs <- sapply(fit$smooth, `[[`, "label")
        is_int <- grepl("intervention", labs)
        conc <- tryCatch({
            cw <- concurvity(fit, full = FALSE)$worst[labs, labs]
            setNames(sapply(seq_along(labs), function(i) max(cw[i, is_int != is_int[i]])), labs)
        }, error = function(e) setNames(rep(NA_real_, length(labs)), labs))

        st <- summary(fit)$s.table
        rows[[key]] <- data.frame(
            group = g, lag = l, model = group_plan$model[group_plan$group == g],
            term = term_labels[rownames(st)], edf = round(st[, "edf"], 2),
            chi_sq = round(st[, "Chi.sq"], 1), p_value = signif(st[, "p-value"], 3),
            concurvity_vs_other = round(conc[rownames(st)], 2),
            dev_expl = round(summary(fit)$dev.expl, 3),
            theta = round(fit$family$getTheta(TRUE), 2),
            n_rows = nrow(res$data),
            hz_term4 = if (full) group_plan$hz_term4[group_plan$group == g] else 0,
            missing_case_weeks_during_intervention = miss,
            row.names = NULL
        )
    }
}

results <- bind_rows(rows)
write.csv(results, file.path(data_dir, "intervention_effect_gam.csv"), row.names = FALSE)
cat("\nSmooth terms (edf near 1 = effect close to linear in time; concurvity_vs_other near 1 =",
    "intervention and baseline curves hard to separate):\n")
print(results %>% select(group, lag, term, edf, p_value, concurvity_vs_other, dev_expl, n_rows,
    missing_case_weeks_during_intervention), row.names = FALSE)

# ---- figures ----
caption_txt <- paste(
    "f(t): log rate ratio of cases at full intervention intensity (1) vs none (0); below 0 = fewer cases.",
    "Band = 95% CI. Shown only for weeks when the intervention was running.",
    "\nIntensity is 0 before first activity and after the last report; missing case weeks are left out of the fit.",
    "Negative binomial hierarchical GAM, REML."
)

for (g in unique(results$group)) {
    keys <- paste0(g, "_lag", lags)
    keys <- keys[keys %in% names(fits)]

    # global effect: weeks where at least one HZ had intervention > 0 (in the data used for the fit)
    glob <- bind_rows(lapply(keys, function(key) {
        res <- fits[[key]]
        on <- res$data %>%
            filter(intervention > 0) %>%
            count(time_index, name = "n_hz")
        nd <- data.frame(time_index = on$time_index, hz = factor(hz_levels[1], levels = hz_levels), intervention = 1)
        term_curve(res$fit, nd, "s(time_index):intervention") %>%
            mutate(n_hz = on$n_hz, lag = sub(".*_lag", "", key))
    })) %>%
        mutate(date = week1_date + 7 * (time_index - 1), lag = lag_labels[lag])

    p <- ggplot(glob, aes(date, est)) +
        geom_hline(yintercept = 0, colour = "grey40") +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = group_colours[g], alpha = 0.3) +
        geom_line(colour = group_colours[g], linewidth = 0.9) +
        geom_text(aes(y = -Inf, label = n_hz), vjust = -0.5, size = 2.2, colour = "grey30") +
        facet_wrap(~lag, ncol = 1) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        labs(
            title = paste0(g, " \u2013 global intervention effect on weekly cases"),
            subtitle = paste0("Model: ", unique(results$model[results$group == g]),
                ". Numbers along the bottom: HZs with the intervention running that week."),
            x = NULL, y = "f(t), log rate ratio", caption = caption_txt
        ) +
        theme_bw(base_size = 9) +
        theme(panel.grid.minor = element_blank())
    ggsave(file.path(fig_dir, paste0("gam_effect_", g, ".png")), p, width = 9, height = 6.5, dpi = 200)

    # HZ curves (global + HZ deviation) for groups with the HZ-deviation term
    if (!any(results$term[results$group == g] == "intervention (HZ deviations)")) next
    hzc <- bind_rows(lapply(keys, function(key) {
        res <- fits[[key]]
        on <- res$data %>% filter(intervention > 0)
        if (nrow(on) == 0) return(NULL)
        nd <- data.frame(time_index = on$time_index, hz = on$hz, intervention = 1)
        term_curve(res$fit, nd, c("s(time_index):intervention", "s(time_index,hz):intervention")) %>%
            mutate(lag = sub(".*_lag", "", key))
    })) %>%
        mutate(date = week1_date + 7 * (time_index - 1), lag = lag_labels[lag])

    # weeks with the intervention running but no case report
    miss <- intensity %>%
        filter(group == g, intervention > 0, lag %in% lags) %>%
        anti_join(cases, by = c("hz", "week")) %>%
        mutate(date = week1_date + 7 * (week - 1), lag = lag_labels[as.character(lag)],
            hz = factor(hz, levels = hz_levels))

    p_hz <- ggplot(hzc, aes(date, est)) +
        geom_hline(yintercept = 0, colour = "grey40") +
        geom_ribbon(aes(ymin = lower, ymax = upper), fill = group_colours[g], alpha = 0.3) +
        geom_line(colour = group_colours[g], linewidth = 0.8) +
        geom_point(data = miss, aes(date, y = -Inf), shape = 4, colour = "red", size = 1.6,
            inherit.aes = FALSE) +
        facet_grid(hz ~ lag) +
        scale_x_date(date_breaks = "1 month", date_labels = "%b") +
        labs(
            title = paste0(g, " \u2013 intervention effect per HZ (global + shrunk HZ deviation)"),
            subtitle = paste("HZ deviations are weakly identified (confounded with the HZ baseline curve):",
                "read as a check on the global effect, not as separate HZ effects. Red x: week with no case report."),
            x = NULL, y = "log rate ratio", caption = caption_txt
        ) +
        theme_bw(base_size = 8) +
        theme(panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0))
    ggsave(file.path(fig_dir, paste0("gam_effect_", g, "_hz.png")), p_hz,
        width = 10, height = 1.2 + 1.1 * n_distinct(hzc$hz), dpi = 200, limitsize = FALSE
    )
}
