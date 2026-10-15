# Weekly intervention intensity c_g(t) in [0, 1] per HZ and intervention group,
# built from the AA tracking indicators. Used as the `intervention` covariate in
# intervention_effect_gam.R.
#
# Outputs:
#   analysis/data/intervention_intensity.csv  weekly profiles (hz, week, date, group,
#                                              intensity, intensity_lag1, active)
#   figures/intensity_<hz>.png                 scaled indicator exposures and group intensity
#   printed table: backlog spikes
#
# Method:
# 1. Indicators: model groups ctc, cati, hyg, chlor, lat (orc has only case-driven
#    counts, i.e. demand, so it is left out). Keep unit types weekly_amount, level,
#    level_probable, people_reached; drop rate (unknown reference amounts),
#    rolling_4wk, case_driven_count, no_data; drop duplicates id46 (= id44) and
#    id18 (subset of id17).
# 2. Backlog spikes (> spike_mult x median non-zero value of that HZ x indicator)
#    are spread evenly over the preceding run of missing weeks (if any).
# 3. Exposure: weekly_amount (deliveries) -> decaying stock
#    S_t = S_{t-1} exp(-1/tau) + x_t (missing weeks inside the reporting period
#    count as no delivery); level / people_reached -> as reported. Each series is
#    scaled to 0-1 by its maximum within the HZ (1 = that HZ's peak activity).
# 4. Group intensity: mean of the scaled indicator exposures each week, linearly
#    interpolated over weeks with no indicator reported, rescaled to max 1 within
#    the HZ. 0 before the first activity and after the last reported week.
# 5. intensity_lag1 = intensity of the previous week (intervention at t-1 acting on
#    cases at t). active = any indicator of the group reported > 0 that week.

library(dplyr)
library(tidyr)
library(ggplot2)

data_dir <- "Z:/acp25_DIDE/CHLAA/analysis/data"
fig_dir <- "Z:/acp25_DIDE/CHLAA/figures"
week1_date <- as.Date("2024-12-30") # IDSR date of 2025 week 1

tau <- 4 # weeks a delivery keeps acting (decay time constant)
spike_mult <- 10

# ---- groups (as in interventions_correlations.R), intensity groups only ----
group_ids <- list(
    ctc   = c(2, 4, 13, 14, 24, 25, 26, 37, 42, 43, 45),
    cati  = c(11, 16, 32, 33),
    hyg   = c(9, 36, 38),
    chlor = c(10, 39),
    lat   = 21
)
group_colours <- c(
    ctc = "#00BFC4", cati = "#7B68EE", hyg = "#FF69B4",
    chlor = "#32CD32", lat = "#f86d6d"
)
group_levels <- names(group_ids)

keep_types <- c("weekly_amount", "level", "level_probable", "people_reached")
drop_dups <- c("id46", "id18")

codebook <- read.csv(file.path(data_dir, "AA_indicator_codebook.csv"),
    stringsAsFactors = FALSE, encoding = "UTF-8"
)

ind_key <- bind_rows(lapply(names(group_ids), function(g) {
    data.frame(id = paste0("id", group_ids[[g]]), group = g, stringsAsFactors = FALSE)
})) %>%
    left_join(codebook[, c("id", "unit_type")], by = "id")

cat("Indicators excluded by unit type:\n")
print(ind_key %>% filter(!unit_type %in% keep_types), row.names = FALSE)
ind_key <- ind_key %>% filter(unit_type %in% keep_types, !id %in% drop_dups)

# ---- weekly indicator data, completed to weeks 1-52 ----
dat <- read.csv(file.path(data_dir, "IDSR_AA_indicators.csv"), stringsAsFactors = FALSE) %>%
    select(-date) %>%
    complete(hz, week = 1:52)

long <- dat %>%
    select(hz, week, all_of(ind_key$id)) %>%
    pivot_longer(-c(hz, week), names_to = "id", values_to = "value") %>%
    left_join(ind_key, by = "id") %>%
    group_by(hz, id) %>%
    filter(any(!is.na(value))) %>%
    arrange(week, .by_group = TRUE) %>%
    ungroup()

# ---- helpers ----

# spread values > mult x median non-zero over the preceding run of missing weeks
spread_spikes <- function(v, mult) {
    flag <- rep(FALSE, length(v))
    nz <- v[!is.na(v) & v > 0]
    if (length(nz) < 3) {
        return(list(v = v, flag = flag))
    }
    thr <- mult * median(nz)
    for (i in which(!is.na(v) & v > thr)) {
        flag[i] <- TRUE
        j <- i - 1
        while (j >= 1 && is.na(v[j])) j <- j - 1
        if (j >= 1 && j < i - 1) v[(j + 1):i] <- v[i] / (i - j) # only gaps inside the reporting period
    }
    list(v = v, flag = flag)
}

# deliveries -> decaying stock over the reporting period; other types as reported
make_exposure <- function(v, type, tau) {
    if (type != "weekly_amount") {
        return(v)
    }
    rep_wk <- which(!is.na(v))
    s <- rep(NA_real_, length(v))
    stock <- 0
    for (t in min(rep_wk):max(rep_wk)) {
        stock <- stock * exp(-1 / tau) + ifelse(is.na(v[t]), 0, v[t])
        s[t] <- stock
    }
    s
}

# group intensity for one HZ x group on weeks 1:52
group_intensity <- function(d) {
    m <- d %>%
        group_by(week) %>%
        summarise(
            e = if (all(is.na(exposure))) NA_real_ else mean(exposure, na.rm = TRUE),
            active = any(value_adj > 0, na.rm = TRUE), .groups = "drop"
        ) %>%
        complete(week = 1:52, fill = list(active = FALSE))
    out <- data.frame(week = 1:52, intensity = 0, active = m$active)
    act <- m$week[!is.na(m$e) & m$e > 0]
    if (length(act) == 0) {
        return(out)
    }
    wk_range <- min(act):max(m$week[!is.na(m$e)])
    obs <- !is.na(m$e)
    mu <- approx(m$week[obs], m$e[obs], xout = wk_range, rule = 2)$y
    out$intensity[wk_range] <- mu / max(mu)
    out
}

# ---- exposure per HZ x indicator ----
expo <- long %>%
    group_by(hz, id) %>%
    group_modify(function(d, key) {
        sp <- spread_spikes(d$value, spike_mult)
        d$value_adj <- sp$v
        d$spike <- sp$flag
        d$exposure <- make_exposure(sp$v, d$unit_type[1], tau)
        mx <- max(d$exposure, na.rm = TRUE)
        d$exposure <- if (mx > 0) d$exposure / mx else d$exposure
        d
    }) %>%
    ungroup()

cat("\nBacklog spikes (>", spike_mult, "x median non-zero):\n")
print(as.data.frame(expo %>% filter(spike) %>%
    transmute(hz, id, group, week,
        reported = value,
        spread = ifelse(value_adj < value, "spread over preceding gap", "kept (no preceding gap)")
    )), row.names = FALSE)

# ---- group intensity, lag 1 ----
weekly <- expo %>%
    group_by(hz, group) %>%
    group_modify(function(d, key) group_intensity(d)) %>%
    arrange(week, .by_group = TRUE) %>%
    mutate(intensity_lag1 = lag(intensity, default = 0)) %>%
    ungroup() %>%
    mutate(date = week1_date + 7 * (week - 1)) %>%
    select(hz, week, date, group, intensity, intensity_lag1, active)

write.csv(weekly, file.path(data_dir, "intervention_intensity.csv"), row.names = FALSE)
cat("\nWrote", nrow(weekly), "weekly rows to intervention_intensity.csv\n")

# ---- figures ----
windows <- read.csv(file.path(data_dir, "hz_parameters.csv"), stringsAsFactors = FALSE) %>%
    filter(grepl("^(ctc|cati|hyg|chlor|lat)_(start|end)$", parameter)) %>%
    separate(parameter, into = c("group", "bound"), sep = "_") %>%
    mutate(value = as.Date(na_if(value, ""))) %>%
    pivot_wider(names_from = bound, values_from = value) %>%
    filter(!is.na(start), !is.na(end))

x_scale <- scale_x_date(
    limits = range(weekly$date) + c(-4, 4), date_breaks = "1 month", date_labels = "%b",
    expand = c(0, 0), oob = scales::oob_squish
)

for (h in sort(unique(weekly$hz))) {
    w <- weekly %>%
        filter(hz == h) %>%
        group_by(group) %>%
        filter(any(intensity > 0)) %>%
        ungroup()
    if (nrow(w) == 0) next
    pts <- expo %>%
        filter(hz == h, group %in% w$group, !is.na(exposure)) %>%
        mutate(date = week1_date + 7 * (week - 1))
    win <- windows %>% filter(hz == h, group %in% w$group)
    lv <- intersect(group_levels, unique(w$group))
    w$group <- factor(w$group, levels = lv)
    pts$group <- factor(pts$group, levels = lv)
    win$group <- factor(win$group, levels = lv)

    p <- ggplot(w, aes(date)) +
        geom_rect(
            data = win, aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = group),
            inherit.aes = FALSE, alpha = 0.12
        ) +
        geom_step(aes(y = intensity, colour = group), linewidth = 0.9, direction = "mid") +
        geom_point(data = pts, aes(y = exposure, shape = id), size = 1.4, colour = "grey20", alpha = 0.7) +
        facet_wrap(~group, ncol = 1) +
        scale_colour_manual(values = group_colours, guide = "none") +
        scale_fill_manual(values = group_colours, guide = "none") +
        scale_shape_manual(values = c(16, 17, 15, 3, 4, 8, 1, 2, 0)) +
        scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1)) +
        x_scale +
        labs(
            title = paste0(h, " \u2013 intervention intensity c(t)"),
            x = NULL, y = "Relative intensity (0-1)", shape = "Indicator",
            caption = paste0(
                "Points: indicator exposure scaled to its HZ max (deliveries as decaying stock, tau = ", tau,
                " wk). Line: group mean, rescaled to max 1.\n",
                "Shaded: model window (hz_parameters.csv). 0 before first activity and after last report."
            )
        ) +
        theme_bw(base_size = 9) +
        theme(legend.position = "bottom", panel.grid.minor = element_blank())

    ggsave(file.path(fig_dir, paste0("intensity_", h, ".png")), p,
        width = 11, height = 1.2 + 1.9 * length(lv), dpi = 200
    )
}
