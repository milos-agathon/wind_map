###########################
# Mapping wind data with R
# Milos Popovic 2022/08/28
###########################
# libraries we need
libs <- c(
    "tidyverse", "sf", "giscoR",
    "lubridate", "classInt",
    "rWind", "metR", "oce"
)

# install missing libraries
installed_libs <- libs %in% rownames(installed.packages())
if (any(installed_libs == F)) {
    install.packages(libs[!installed_libs])
}

# load libraries
invisible(lapply(libs, library, character.only = T))

sysfonts::font_add_google("Montserrat", "Montserrat")
showtext::showtext_auto()

# 1. GET WIND DATA
#------------------

get_wind_data <- function(time_range, mean_wind_data, eur_wind_df) {
    time_range <- seq(ymd_hms(paste(2022, 8, 27, 00, 00, 00, sep = "-")),
        ymd_hms(paste(2022, 8, 28, 00, 00, 00, sep = "-")),
        by = "1 hours"
    )

    mean_wind_data <- rWind::wind.dl_2(time_range, -28.5, 58.5, 34.0, 73.5) %>%
        rWind::wind.mean()

    eur_wind_df <- as.data.frame(mean_wind_data)
    return(eur_wind_df)
}

eur_wind_df <- get_wind_data()

# 2. TRANSFORM WIND DATA
#----------------------

get_wind_points <- function() {
    # convert wind data into points
    eur_wind_pts <- eur_wind_df %>%
        sf::st_as_sf(coords = c("lon", "lat")) %>%
        sf::st_set_crs(4326)
    return(eur_wind_pts)
}

eur_wind_pts <- get_wind_points()

# 3. GRIDDED WIND DATA
#----------------------

# divide wind points into grids
get_wind_grid <- function() {
    eur_wind_grid <- eur_wind_pts %>%
        sf::st_make_grid(n = c(80, 100)) %>%
        sf::st_sf() %>%
        dplyr::mutate(id = row_number())
    return(eur_wind_grid)
}

eur_wind_grid <- get_wind_grid()

get_wind_grid_aggregated <- function() {
    eur_wind_grid_agg <-
        sf::st_join(eur_wind_pts, eur_wind_grid,
            join = sf::st_within
        ) %>%
        sf::st_drop_geometry() %>%
        dplyr::group_by(id) %>%
        dplyr::summarise(
            n = n(), u = mean(ugrd10m),
            v = mean(vgrd10m), speed = mean(speed)
        ) %>%
        dplyr::inner_join(eur_wind_grid, by = "id") %>%
        dplyr::select(n, u, v, speed, geometry) %>%
        sf::st_as_sf() %>%
        na.omit()
    return(eur_wind_grid_agg)
}

eur_wind_grid_agg <- get_wind_grid_aggregated()

# 4. INTERPOLATE WIND DATA
#-------------------------

## obtain the centroid coordinates from the grid as table
get_wind_coords <- function() {
    coords <- eur_wind_grid_agg %>%
        st_centroid() %>%
        st_coordinates() %>%
        as_tibble() %>%
        rename(lon = X, lat = Y)
    return(coords)
}

coords <- get_wind_coords()

## remove the geometry from the simple feature of wind data
get_wind_df <- function() {
    eur_df <- coords %>%
        bind_cols(sf::st_drop_geometry(eur_wind_grid_agg))
    return(eur_df)
}

eur_df <- get_wind_df()


## interpolate the U component
get_u_interpolation <- function() {
    wu <- oce::interpBarnes(
        x = eur_df$lon,
        y = eur_df$lat,
        z = eur_df$u
    )
    return(wu)
}

wu <- get_u_interpolation()

## get the width and length for tranforming wide into long table
dimension <- data.frame(lon = wu$xg, wu$zg) %>% dim()

## make a U component data table from interpolated matrix
get_u_table <- function() {
    udf <- data.frame(
        lon = wu$xg,
        wu$zg
    ) %>%
        gather(key = "lata", value = "u", 2:dimension[2]) %>%
        mutate(lat = rep(wu$yg, each = dimension[1])) %>%
        select(lon, lat, u) %>%
        as_tibble()

    return(udf)
}

udf <- get_u_table()

## interpolate the V component
get_u_interpolation <- function() {
    wv <- oce::interpBarnes(
        x = eur_df$lon,
        y = eur_df$lat,
        z = eur_df$v
    )
    return(wv)
}

wv <- get_u_interpolation()

## make the V component data table from interpolated matrix
get_v_table <- function() {
    vdf <- data.frame(lon = wv$xg, wv$zg) %>%
        gather(key = "lata", value = "v", 2:dimension[2]) %>%
        mutate(lat = rep(wv$yg, each = dimension[1])) %>%
        select(lon, lat, v) %>%
        as_tibble()
    return(vdf)
}

vdf <- get_v_table()

## merge the V and U component tables and compute velocity
get_final_table <- function() {
    df <- udf %>%
        bind_cols(vdf %>% select(v)) %>%
        mutate(vel = sqrt(u^2 + v^2))
    return(df)
}

df <- get_final_table()

# 5. MAP WIND DATA
#-------------------------

get_europe_sf <- function() {
    eur_sf <- giscoR::gisco_get_countries(
        year = "2016", epsg = "4326",
        resolution = "10", region = c("Europe", "Asia")
    )

    return(eur_sf)
}

# bounding box
crsLONGLAT <- "+proj=longlat +datum=WGS84 +no_defs"

get_bounding_box <- function(bbox, bb) {
    bbox <- st_sfc(
        st_polygon(list(cbind(
            c(-25, 48.5, 48.5, -25, -25),
            c(32.000, 32.000, 69.5, 69.5, 32.000)
        ))),
        crs = crsLONGLAT
    )

    bb <- sf::st_bbox(bbox)

    return(bb)
}

# colors
cols <- c(
    "#feebe2", "#d84594", "#bc2b8a", "#7a0177"
)

newcol <- colorRampPalette(cols)
ncols <- 6
cols2 <- newcol(ncols)

# breaks
vmin <- min(df$vel, na.rm = T)
vmax <- max(df$vel, na.rm = T)

brk <- classInt::classIntervals(df$vel,
    n = 6,
    style = "fisher"
)$brks %>%
    head(-1) %>%
    tail(-1) %>%
    append(vmax)

breaks <- c(vmin, brk)

make_wind_map <- function(eur_sf, bb) {
    eur_sf <- get_europe_sf()
    bb <- get_bounding_box()

    p <- df %>%
        ggplot() +
        metR::geom_streamline(
            data = df,
            aes(
                x = lon, y = lat, dx = u, dy = v,
                color = sqrt(..dx..^2 + ..dy..^2)
            ),
            L = 2, res = 2, n = 60,
            arrow = NULL, lineend = "round",
            alpha = .85
        ) +
        geom_sf(
            data = eur_sf,
            fill = NA,
            color = "#07CFF7",
            size = .25,
            alpha = .99
        ) +
        coord_sf(
            crs = crsLONGLAT,
            xlim = c(bb["xmin"], bb["xmax"]),
            ylim = c(bb["ymin"], bb["ymax"])
        ) +
        scale_color_gradientn(
            name = "Average speed (m/s)",
            colours = cols2,
            breaks = breaks,
            labels = round(breaks, 1),
            limits = c(vmin, vmax)
        ) +
        guides(
            fill = "none",
            color = guide_legend(
                override.aes = list(size = 3, alpha = 1, shape = 15),
                direction = "horizontal",
                keyheight = unit(2.5, units = "mm"),
                keywidth = unit(15, units = "mm"),
                title.position = "top",
                title.hjust = .5,
                label.hjust = .5,
                nrow = 1,
                byrow = T,
                reverse = F,
                label.position = "bottom"
            )
        ) +
        theme_bw() +
        theme(
            text = element_text(family = "Montserrat"),
            axis.line = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            legend.position = c(.5, .05),
            legend.text = element_text(size = 60, color = "white"),
            legend.title = element_text(size = 80, color = "white"),
            legend.key = element_blank(),
            plot.margin = unit(c(t = 1, r = -2, b = -1, l = -2), "lines"),
            panel.grid.major = element_line(color = "#070C33", size = 0.2),
            panel.grid.minor = element_blank(),
            plot.title = element_text(
                face = "bold", size = 120,
                color = "white", hjust = .5
            ),
            plot.caption = element_text(
                size = 50, color = "grey80",
                hjust = .25, vjust = -10
            ),
            plot.subtitle = element_blank(),
            plot.background = element_rect(fill = "#070C33", color = NA),
            panel.background = element_rect(fill = "#070C33", color = NA),
            legend.background = element_rect(fill = "#070C33", color = NA),
            legend.spacing.y = unit(.5, "pt"),
            panel.border = element_blank()
        ) +
        labs(
            x = "",
            y = "",
            title = "Average wind speed in Europe (27 August 2022)",
            subtitle = "",
            caption = ""
        )
    return(p)
}

p <- make_wind_map()

ggsave(
    filename = "eur_wind_27august2022.png",
    width = 8.5, height = 7, dpi = 600, device = "png", p
)
