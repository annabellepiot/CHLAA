#=====================================================================
#04_descriptive_map_of_healthzones.R
#---------------------------------------------------------------------
#Map of the 2025 cumulative attack rate (incidence) for the 12 CHLAA
#study health zones.
#
#   Attack rate = (sum of 2025 cases) / population * 100,000
#   Study window: 2025 epi-year (IDSR `year == 2025`)
#   Population  : per-zone constant `population` from IDSR_dataset.csv
#
#Design: one DRC province locator (province borders only, the 4 study
#provinces highlighted + labelled) with circular "porthole" insets -
#one per study province - connected by arrows. Each inset shows that
#province's health zones clipped to a circle, study zones coloured on a
#shared viridis scale, surrounding zones light grey for context. The
#highlighted health-zone names are printed above each inset.
#
#The whole figure is a single ggplot on a synthetic canvas: geometries
#are affine-transformed (translate + uniform scale) into the DRC box or
#each inset circle, so arrows can be drawn across the entire canvas.
#
#Requires `sf` in the Renv conda env:
#   conda install -n Renv -c conda-forge r-sf
#=====================================================================

suppressMessages({
    library(sf)
    library(dplyr)
    library(tidyr)
    library(readr)
    library(ggplot2)
    library(scales)
    library(stringr)
})
sf_use_s2(FALSE) #planar ops are fine at this scale and avoid s2 edge cases

##---- Paths ---------------------------------------------------------
base_dir <- "/rds/general/user/acp25/home/MIMIC/Clean_data/Proj_2/CHLAA"
data_dir <- file.path(base_dir, "analysis", "data")
fig_dir <- file.path(base_dir, "figures")
out_dir <- file.path(base_dir, "analysis", "output")
hz_shp <- file.path(data_dir, "Healthzone_shapefiles", "RDC_Zones de santé.shp")
prov_shp <- file.path(data_dir, "cod_admin_boundaries.shp", "cod_admin1.shp")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

##---- Zone lookup: hz slug -> shapefile `Nom` + province ------------
hz_lookup <- tibble::tribble(
    ~hz,           ~nom,           ~province,
    "alunguli",    "Alunguli",     "Maniema",
    "bumbu",       "Bumbu",        "Kinshasa",
    "goma",        "Goma",         "Nord-Kivu",
    "kindu",       "Kindu",        "Maniema",
    "kingabwa",    "Kingabwa",     "Kinshasa",
    "kirotshe",    "Kirotshe",     "Nord-Kivu",
    "kokolo",      "Kokolo",       "Kinshasa",
    "limete",      "Limete",       "Kinshasa",
    "maluku_i",    "Maluku 1",     "Kinshasa",
    "ngiri_ngiri", "Ngir Ngiri",   "Kinshasa",
    "nyiragongo",  "Nyirangongo",  "Nord-Kivu",
    "opala",       "Opala",        "Tshopo"
)
study_provinces <- c("Kinshasa", "Nord-Kivu", "Maniema", "Tshopo")

##---- Attack-rate table from IDSR -----------------------------------
idsr <- read_csv(file.path(data_dir, "IDSR_dataset.csv"), show_col_types = FALSE)

ar_tbl <- idsr %>%
    filter(hz %in% hz_lookup$hz, year == 2025) %>%
    group_by(hz) %>%
    summarise(
        cases_2025 = sum(cases, na.rm = TRUE),
        population = dplyr::first(population), .groups = "drop"
    ) %>%
    mutate(ar_100k = 1e5 * cases_2025 / population) %>%
    right_join(hz_lookup, by = "hz") %>%
    arrange(desc(ar_100k))

cat("\n2025 cumulative attack rate (per 100,000):\n")
print(as.data.frame(ar_tbl[, c("hz", "province", "cases_2025", "population", "ar_100k")]),
    digits = 5, row.names = FALSE
)
write_csv(ar_tbl, file.path(out_dir, "healthzone_attack_rate_2025.csv"))

##---- Geometry (all WGS84) ------------------------------------------
zones <- st_read(hz_shp, quiet = TRUE) %>%
    st_transform(4326) %>%
    st_make_valid()
prov <- st_read(prov_shp, quiet = TRUE) %>%
    st_transform(4326) %>%
    st_make_valid()

#Attach attack rate to every health zone (NA for non-study zones)
zones <- zones %>%
    left_join(ar_tbl %>% select(nom, hz, ar_100k), by = c("Nom" = "nom")) %>%
    mutate(is_study = !is.na(ar_100k))

#Integrity check
stopifnot(sum(zones$is_study) == nrow(hz_lookup))
cat("\nMatched", sum(zones$is_study), "of", nrow(hz_lookup), "study zones to geometry.\n\n")

#=====================================================================
#Canvas compositing helpers
#=====================================================================
#Affine transform of an sf object's geometry: translate source centre
#(cx,cy) to target (tx,ty) and scale uniformly by s. Returns an sf with
#CRS dropped (canvas units).
#CANVAS_CRS is a projected CRS attached to canvas-space geometries only
#so coord_sf treats them as planar (equal aspect) and never reprojects.
CANVAS_CRS <- 3857
aff <- function(x, cx, cy, s, tx, ty) {
    g <- st_geometry(x)
    g2 <- st_sfc(lapply(g, function(p) (p - c(cx, cy)) * s + c(tx, ty)), crs = CANVAS_CRS)
    st_set_geometry(x, g2)
}

#Circle as an sf polygon (geographic), for clipping
circ_poly <- function(cx, cy, r, n = 240) {
    a <- seq(0, 2 * pi, length.out = n)
    ring <- cbind(cx + r * cos(a), cy + r * sin(a))
    ring[n, ] <- ring[1, ] #force exact closure
    st_sfc(st_polygon(list(ring)), crs = 4326)
}
#Circle outline as a data frame (canvas units), for drawing borders
circ_df <- function(cx, cy, r, grp, n = 180) {
    a <- seq(0, 2 * pi, length.out = n)
    data.frame(x = cx + r * cos(a), y = cy + r * sin(a), grp = grp)
}
FONT <- "Helvetica"
#Display form of a health-zone slug: drop underscores, Title Case
pretty_hz <- function(v) {
    v <- gsub("_", " ", v)
    gsub("(^|\\s)([a-z])", "\\1\\U\\2", v, perl = TRUE)
}
wrap_names <- function(v, width = 22) {
    paste(strwrap(paste(pretty_hz(sort(v)), collapse = ", "), width = width), collapse = "\n")
}

##---- DRC locator transform (province polygons -> top-left box) -----
#Canvas is 0..100 in both axes; coord_sf keeps aspect 1:1 (round circles)
DRC_BOX <- c(left = 0, right = 46, bottom = 42, top = 96)
bd <- st_bbox(prov)
sd <- min(
    (DRC_BOX["right"] - DRC_BOX["left"]) / (bd["xmax"] - bd["xmin"]),
    (DRC_BOX["top"] - DRC_BOX["bottom"]) / (bd["ymax"] - bd["ymin"])
)
dcx <- mean(bd[c("xmin", "xmax")])
dcy <- mean(bd[c("ymin", "ymax")])
dtx <- mean(DRC_BOX[c("left", "right")])
dty <- mean(DRC_BOX[c("bottom", "top")])
to_drc <- function(x) aff(x, dcx, dcy, sd, dtx, dty)

prov$is_study <- prov$adm1_name %in% study_provinces
prov_t <- to_drc(prov)

#Province label / arrow anchor points (point-on-surface, on canvas)
prov_pts <- prov %>% filter(is_study)
pp_geo <- suppressWarnings(st_point_on_surface(st_geometry(prov_pts)))
pp_can <- st_coordinates(st_geometry(to_drc(st_set_geometry(prov_pts, pp_geo))))
prov_anchor <- tibble(province = prov_pts$adm1_name, px = pp_can[, 1], py = pp_can[, 2])

#Manual label offsets so the three central provinces don't collide
prov_anchor <- prov_anchor %>%
    left_join(tibble::tribble(
        ~province, ~ndx, ~ndy,
        "Kinshasa", -1, 3,
        "Tshopo", -3, 4.5,
        "Maniema", -7, -2,
        "Nord-Kivu", -1, -5
    ), by = "province")

##---- Inset layout (edit these to refine placement) -----------------
#Arranged to echo geography and pack tightly around the locator:
#Kinshasa (far west) bottom-left; Maniema (S) bottom-centre; Nord-Kivu
#(E) mid-right; Tshopo (N) top-right. `curv` bows each connector arrow.
insets <- tibble::tribble(
    ~province, ~cx, ~cy, ~r, ~curv,
    "Kinshasa", 18, 22, 15, 0.20,
    "Maniema", 50, 21, 15, 0.18,
    "Nord-Kivu", 72, 45, 16, -0.18,
    "Tshopo", 64, 81, 15, -0.20
)
MARGIN <- 1.6 #circle radius = half study-zone extent * MARGIN

#Second-level "urban zoom" panel: a circular inset LEFT of the Kinshasa
#locator inset that enlarges just the 5 clustered urban zones (too small
#to read in the province inset). Its own circle, grey context map, shared
#colour scale; fed by an arrow from a small "lens" circle drawn around the
#cluster inside the Kinshasa province inset.
UZ_cx <- -15 #new-panel centre (canvas units): left of, and above, the Kinshasa
UZ_cy <- 40 #inset, so the connector arrow curves up-left (compact, saves width)
UZ_r <- 13 #new-panel circle radius
MARGIN_URBAN <- 2.0 #circle radius = max anchor distance * MARGIN_URBAN

#Zone labels sit at each territory INSIDE its circle. In the Kinshasa
#province inset only Maluku 1 is labelled; the 5 urban zones are labelled
#individually on the urban-zoom panel instead.
KIN_URBAN <- c("bumbu", "kingabwa", "kokolo", "limete", "ngiri_ngiri")
#Fine-tune individual zone-label anchors (canvas units), keyed by id. The
#`uz_*` ids nudge the 5 urban labels apart on the zoom panel.
#The 5 urban labels are pushed off their (clustered) anchors into the
#surrounding space and connected by thin leader lines. Kokolo/Ngiri Ngiri/
#Bumbu fan out left; Limete/Kingabwa up-right.
zone_nudge <- tibble::tribble(
    ~id, ~dx, ~dy,
    "uz_kokolo", -2.3, 2.8,
    "uz_ngiri_ngiri", 4.5, -6.3,
    "uz_bumbu", -2.9, -3.3,
    "uz_limete", 0.9, 6.4,
    "uz_kingabwa", 2.85, 3.9,
    "kirotshe", 0, 2.5, #nudge the Nord-Kivu label up slightly
    "goma", -2, -5 #push Goma's label off its tiny territory (leader added)
)

##---- Build each inset (clip province context to a circle) ----------
inset_zones <- list()
disc_df <- list()
border_df <- list()
label_df <- list()
arrow_df <- list()

for (i in seq_len(nrow(insets))) {
    pr <- insets$province[i]
    cx <- insets$cx[i]
    cy <- insets$cy[i]
    r <- insets$r[i]
    sz <- zones %>% filter(is_study, PROVINCE == pr)
    bb <- st_bbox(sz)
    fcx <- mean(bb[c("xmin", "xmax")])
    fcy <- mean(bb[c("ymin", "ymax")])
    fr <- 0.5 * max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"]) * MARGIN

    circ <- circ_poly(fcx, fcy, fr)
    clip <- suppressWarnings(st_intersection(zones, circ))
    clip <- clip[!st_is_empty(st_geometry(clip)), ]
    clip <- suppressWarnings(st_collection_extract(clip, "POLYGON"))
    clip_t <- aff(clip, fcx, fcy, r / fr, cx, cy) #scale circle radius fr -> r
    inset_zones[[i]] <- clip_t

    disc_df[[i]] <- circ_df(cx, cy, r, grp = pr) #white backing disc
    border_df[[i]] <- circ_df(cx, cy, r, grp = pr) #dark outline

    #Per-zone labels, anchored at each territory (point-on-surface of the
    #clipped, transformed polygon). In the Kinshasa inset only Maluku 1 is
    #labelled - the 5 urban zones are labelled on the urban-zoom panel, and
    #a "lens" circle round their cluster feeds the arrow to that panel.
    study_clip <- clip_t %>% filter(is_study)
    if (pr == "Kinshasa") {
        groups <- list(list(id = "maluku_i", hz = "maluku_i", text = "Maluku 1"))
        #Lens circle (Kinshasa-inset canvas coords) enclosing the 5 urban zones
        urb <- study_clip %>% filter(hz %in% KIN_URBAN)
        ub <- st_bbox(urb)
        lens_cx <- mean(ub[c("xmin", "xmax")])
        lens_cy <- mean(ub[c("ymin", "ymax")])
        lens_r <- 0.5 * max(ub["xmax"] - ub["xmin"], ub["ymax"] - ub["ymin"]) * 1.35
    } else {
        groups <- lapply(sort(unique(study_clip$hz)), function(h)
            list(id = h, hz = h, text = pretty_hz(h)))
    }
    label_df[[i]] <- bind_rows(lapply(groups, function(g) {
        gg <- study_clip %>% filter(hz %in% g$hz)
        if (nrow(gg) == 0) {
            return(NULL)
        }
        pt <- st_coordinates(st_point_on_surface(st_union(st_geometry(gg))))
        tibble(id = g$id, x = pt[1], y = pt[2], label = g$text)
    }))

    #arrow: province anchor -> nearest point on the inset circle
    pa <- prov_anchor %>% filter(province == pr)
    dx <- cx - pa$px
    dy <- cy - pa$py
    d <- sqrt(dx^2 + dy^2)
    arrow_df[[i]] <- tibble(
        x = pa$px, y = pa$py,
        xend = cx - dx / d * r, yend = cy - dy / d * r,
        curv = insets$curv[i]
    )
}

##---- Urban-zoom panel (second-level zoom on the 5 Kinshasa zones) ---
#Same build as an inset, but focused on just KIN_URBAN so the tiny city
#zones fill the circle. Surrounding zones stay grey; the 5 are coloured
#on the shared scale and labelled individually.
urban <- zones %>% filter(is_study, PROVINCE == "Kinshasa", hz %in% KIN_URBAN)
#Centre on the 5 zones' label anchors, not the union bbox: a thin southern
#spur inflates the bbox and would push the cluster up into a corner. This
#keeps the five zones centred and filling the circle.
uanch <- suppressWarnings(st_coordinates(st_point_on_surface(st_geometry(urban))))
ufcx <- mean(uanch[, 1])
ufcy <- mean(uanch[, 2])
ufr <- max(sqrt((uanch[, 1] - ufcx)^2 + (uanch[, 2] - ufcy)^2)) * MARGIN_URBAN
ucirc <- circ_poly(ufcx, ufcy, ufr)
uclip <- suppressWarnings(st_intersection(zones, ucirc))
uclip <- uclip[!st_is_empty(st_geometry(uclip)), ]
uclip <- suppressWarnings(st_collection_extract(uclip, "POLYGON"))
uclip_t <- aff(uclip, ufcx, ufcy, UZ_r / ufr, UZ_cx, UZ_cy)

j <- length(inset_zones) + 1L
inset_zones[[j]] <- uclip_t
disc_df[[j]] <- circ_df(UZ_cx, UZ_cy, UZ_r, grp = "urban_zoom") #white backing
border_df[[j]] <- circ_df(UZ_cx, UZ_cy, UZ_r, grp = "urban_zoom") #dark outline

#Individual labels for the 5 urban zones (point-on-surface of each polygon)
uz_study <- uclip_t %>% filter(is_study, hz %in% KIN_URBAN)
label_df[[j]] <- bind_rows(lapply(KIN_URBAN, function(h) {
    gg <- uz_study %>% filter(hz == h)
    if (nrow(gg) == 0) {
        return(NULL)
    }
    pt <- st_coordinates(st_point_on_surface(st_union(st_geometry(gg))))
    tibble(id = paste0("uz_", h), x = pt[1], y = pt[2], label = pretty_hz(h))
}))

#Lens circle round the cluster in the Kinshasa inset (drawn as a border).
border_df[[j + 1L]] <- circ_df(lens_cx, lens_cy, lens_r, grp = "urban_lens")
#Connector arrow from the lens edge to the urban-zoom panel edge. This one
#starts INSIDE the Kinshasa disc, so it is drawn as its own late layer (over
#the zones) rather than via arrow_df, which is drawn behind the discs.
adx <- UZ_cx - lens_cx
ady <- UZ_cy - lens_cy
ad <- sqrt(adx^2 + ady^2)
urban_arrow <- tibble(
    x = lens_cx + adx / ad * lens_r, y = lens_cy + ady / ad * lens_r,
    xend = UZ_cx - adx / ad * UZ_r, yend = UZ_cy - ady / ad * UZ_r,
    curv = 0.30
)

inset_zones <- do.call(rbind, inset_zones)
disc_df <- do.call(rbind, disc_df)
border_df <- do.call(rbind, border_df)
label_df <- do.call(rbind, label_df) %>%
    left_join(zone_nudge, by = "id") %>%
    mutate(
        ax = x, ay = y, #keep anchor for leader lines
        x = x + coalesce(dx, 0), y = y + coalesce(dy, 0)
    )
#Leader lines from anchor -> nudged label: the 5 urban zones, plus Goma (a
#tiny zone whose label would otherwise cover its whole territory).
LEADER_IDS <- c(paste0("uz_", KIN_URBAN), "goma")
leader_df <- label_df %>% filter(id %in% LEADER_IDS)
urban_lab <- label_df %>% filter(grepl("^uz_", id))
arrow_df <- do.call(rbind, arrow_df)

#Curved connector arrows: one geom_curve per arrow so curvature can vary
arrow_layers <- lapply(seq_len(nrow(arrow_df)), function(i) {
    geom_curve(
        data = arrow_df[i, ], aes(x = x, y = y, xend = xend, yend = yend),
        curvature = arrow_df$curv[i], ncp = 12, colour = "grey40",
        linewidth = 0.55, arrow = arrow(length = unit(0.22, "cm"), type = "closed")
    )
})

#Lens -> urban-zoom arrow, same style but drawn as a late layer (its tail
#sits inside the Kinshasa disc, so it must render over the zones).
urban_arrow_layer <- geom_curve(
    data = urban_arrow, aes(x = x, y = y, xend = xend, yend = yend),
    curvature = urban_arrow$curv, ncp = 12, colour = "grey40",
    linewidth = 0.55, arrow = arrow(length = unit(0.22, "cm"), type = "closed")
)

##---- Shared colour scale (ColorBrewer YlGnBu, middle 80%) ----------
#Sequential yellow -> green -> blue (low -> high). Trim the lightest and
#darkest 10% so the extremes aren't near-white / near-black.
ylgnbu <- scales::brewer_pal(palette = "YlGnBu", direction = 1)(9)
ylgnbu_mid <- colorRampPalette(ylgnbu)(100)[11:90]
fill_scale <- scale_fill_gradientn(
    colours = ylgnbu_mid, name = "Cumulative attack rate\n(per 100,000)",
    labels = comma, limits = c(0, max(zones$ar_100k, na.rm = TRUE)),
    guide = guide_colourbar(
        barwidth = 15, barheight = 1.0, title.position = "top",
        title.hjust = 0.5
    )
)

##---- Compose -------------------------------------------------------
CTX_GREY <- "grey85" #shared context grey (DRC locator + inside circles)

p <- ggplot() +
    #DRC locator (same context grey as the circles; darker borders so the
    #province outlines stay visible against that grey)
    geom_sf(
        data = prov_t %>% filter(!is_study), fill = CTX_GREY,
        colour = "grey45", linewidth = 0.3
    ) +
    geom_sf(
        data = prov_t %>% filter(is_study), fill = "#d7e0e6",
        colour = "grey25", linewidth = 0.5
    ) +
    geom_label(
        data = prov_anchor, aes(px + ndx, py + ndy, label = province),
        size = 5.6, fontface = "bold", colour = "grey15", family = FONT,
        fill = "white", linewidth = 0.25, label.padding = unit(0.16, "lines")
    ) +
    #curved connector arrows
    arrow_layers +
    #inset porthole discs (white backing so gaps/water read clean)
    geom_polygon(data = disc_df, aes(x, y, group = grp), fill = "white") +
    geom_sf(
        data = inset_zones %>% filter(!is_study), fill = CTX_GREY,
        colour = "grey65", linewidth = 0.4
    ) +
    geom_sf(
        data = inset_zones %>% filter(is_study), aes(fill = ar_100k),
        colour = "grey25", linewidth = 0.55
    ) +
    geom_polygon(
        data = border_df, aes(x, y, group = grp), fill = NA,
        colour = "grey30", linewidth = 0.6
    ) +
    #lens -> urban-zoom arrow, over the zones so it reads from the lens
    urban_arrow_layer +
    #leader lines connecting the 5 urban labels to their (clustered) zones
    geom_segment(
        data = leader_df, aes(x = ax, y = ay, xend = x, yend = y),
        colour = "grey40", linewidth = 0.6
    ) +
    #province / larger-zone labels
    geom_label(
        data = label_df %>% filter(!grepl("^uz_", id)), aes(x, y, label = label),
        size = 4.9, fontface = "bold", colour = "grey15", family = FONT,
        fill = alpha("white", 0.82), linewidth = 0.2, lineheight = 0.9,
        label.padding = unit(0.12, "lines")
    ) +
    #5 urban-zone labels (slightly smaller; they sit close together)
    geom_label(
        data = urban_lab, aes(x, y, label = label),
        size = 4.9, fontface = "bold", colour = "grey15", family = FONT,
        fill = alpha("white", 0.9), linewidth = 0.2, lineheight = 0.9,
        label.padding = unit(0.1, "lines")
    ) +
    fill_scale +
    coord_sf(
        xlim = c(-30, 90), ylim = c(2, 96), expand = FALSE,
        crs = CANVAS_CRS, datum = NA
    ) +
    theme_void(base_size = 13, base_family = FONT) +
    theme(
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        legend.position = "inside",
        legend.position.inside = c(0.575, 0.44),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
        plot.margin = margin(8, 8, 8, 8)
    )

out_png <- file.path(fig_dir, "healthzone_attack_rate_map_2025.png")
ggsave(out_png, p, width = 15, height = 12, dpi = 600, bg = "white")
cat("Saved map to:", out_png, "\n")
cat("Saved table to:", file.path(out_dir, "healthzone_attack_rate_2025.csv"), "\n")
