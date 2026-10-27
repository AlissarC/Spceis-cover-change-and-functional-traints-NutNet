
install.packages(c("rnaturalearth", "rnaturalearthdata", "sf"))
install.packages("officer") ## to export table in word
install.packages("flextable")## to export table in word
install.packages("plotbiomes") ## to export table in word

library(ggplot2)
library(rnaturalearth) ## for world map
library(rnaturalearthdata)
library(sf)
library(officer)
library(flextable)
library(dplyr)

###################### Load data ######################################################################
summary = read.csv("./output/final_df.csv")
names(summary)

summary_clean <- summary %>% filter(!is.na(Nmass_NutNet))
nrow(summary_clean) # 37967

length(unique(summary_clean$site_code)) #94

# Number of years for each site
years_per_site <- summary_clean %>%
  distinct(site_name, site_code, n_years) %>%
  arrange(site_code)
##################### extract Table S2 for the paper ######################
sites_summary <- summary_clean %>%
  group_by(site_code) %>%
  summarise(
    Latitude = first(latitude),
    Longitude = first(longitude),
    Continent = first(continent),
    Tg = sprintf("%.2f", first(tmp_GS)),
    MI =sprintf("%.2f", first(aridity)),
    MAP = round(first(pre_annual)),
    MAT = round(first(tmp_annual)),
    .groups = "drop"
  )

sites_summary
names(sites_summary)
write.csv(sites_summary, "./Output/sites_summary.csv")

####################### function to transform lat and long in degrees and seconds, N, S , E , W ###########################
convert_to_dms <- function(decimal, is_lat = TRUE) {
  dir <- if (is_lat) {
    ifelse(decimal >= 0, "N", "S")
  } else {
    ifelse(decimal >= 0, "E", "W")
  }
  
  abs_dec <- abs(decimal)
  degrees <- floor(abs_dec)
  minutes_full <- (abs_dec - degrees) * 60
  minutes <- floor(minutes_full)
  seconds <- round((minutes_full - minutes) * 60, 2)
  
  paste0(degrees, "°", minutes, "'", seconds, "\" ", dir)
}

names(sites_summary)
sites_summary_new <- sites_summary %>%
  mutate(
    Latitude_dms = sapply(Latitude, convert_to_dms, is_lat = TRUE),
    Longitude_dms = sapply(Longitude, convert_to_dms, is_lat = FALSE)
  )
sites_summary_new

sites_summary_final <- sites_summary_new %>%
  dplyr::select(-Latitude, -Longitude)


sites_summary_year <- sites_summary_final %>% 
  left_join(years_per_site, by = c("site_code"))
sites_summary_year
sites_summary_year = as.data.frame(sites_summary_year)

names(sites_summary_year)

length(unique(sites_summary_year$site_code)) # 94

write.csv(sites_summary_year, "./Tables/sites_summary_year.csv")

####################### Export for Supplementary information Table ########################
site_summary_clean <- sites_summary_year %>%
  rename(`Number of years` = n_years,
         `Latitude` = Latitude_dms,
         `Longitude` = Longitude_dms)
names(site_summary_clean)

site_summary_clean2 <- site_summary_clean %>%
  dplyr::select(site_code, Continent, Latitude, Longitude, Tg, MAP,MI,`Number of years`) 
names(site_summary_clean2)

# Convert to flextable
ft <- flextable(site_summary_clean2) %>%
  autofit() %>%
  width(width = rep(0.85, ncol(site_summary_clean2)))  # 8 columns × 0.8 inches = 6.4 inches

# Export to Word
doc <- read_docx() %>%
  body_add_par("Site Summary Table", style = "heading 1") %>%
  body_add_flextable(value = ft)

print(doc, target = "./Tables/TableS1.docx")


############################## World map #####################
world <- ne_countries(scale = "medium", returnclass = "sf")

sites_sf <- sites_summary |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

sites_sf <- sites_summary |>
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326)

world_map <- ggplot() +
  
  # Land
  geom_sf(
    data = world,
    fill = "#f2f2f2",        # soft grey land
    color = "#d0d0d0",       # subtle borders
    linewidth = 0.2
  ) +
  
  # Sites with white halo (very important visually)
  geom_sf(
    data = sites_sf,
    shape = 21,
    fill = "#2c7fb8",        # professional blue
    color = "white",         # halo effect
    size = 3,
    stroke = 0.6,
    alpha = 0.95
  ) +
  
  # Better projection
  coord_sf(
    crs = "+proj=robin",
    expand = FALSE
  ) +
  
  theme_void() +
  theme(
    plot.background  = element_rect(fill = "#dbe9f4", color = NA),
    panel.background = element_rect(fill = "#dbe9f4", color = NA)
  )

print(world_map)

ggsave( "./Fig/FigS1.png", plot = world_map,width = 10,height = 10,units = "in",dpi = 800)

######################################## Fig S2 ######################################################

theme_FigS2 <- theme_classic(base_size = 12) +
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8))

p1 <- dominance_Nmass_MI_NutNet + theme_FigS2 + labs(tag = "(a)")
p2 <- dominance_logNmass_Tg_NutNet + theme_FigS2 + labs(tag = "(b)")
p3 <- dominance_chi_MI_NutNet + theme_FigS2 + labs(tag = "(c)")
p4 <- dominance_chi_Tg_NutNet + theme_FigS2 + labs(tag = "(d)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 4) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 4) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 3.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4)) +
  plot_layout(heights = c(0.12, 1, 1))


# Final figure 

FigS2 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.6, 1.4))
FigS2

# save figure
ggsave("./Fig/FigS2_final.pdf",FigS2,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/FigS2_final.tiff",FigS2,width = 10,height = 8.5,dpi=600, units = "in",bg = "white")

######################################## Fig S3 ######################################################

theme_FigS3 <- theme_classic(base_size = 12) +
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8))

p1 <- lambda_SLA_MI_NutNet + theme_FigS3 + labs(tag = "(a)")
p2 <- lambda_logSLA_Tg_NutNet + theme_FigS3 + labs(tag = "(b)")
p3 <- dominance_SLA_MI_NutNet + theme_FigS3 + labs(tag = "(b)")
p4 <- dominance_logSLA_Tg_NutNet + theme_FigS3 + labs(tag = "(d)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 4) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 4) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 3.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4)) +
  plot_layout(heights = c(0.12, 1, 1))


# Final figure 

FigS3 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.6, 1.4))
FigS3

# save figure
ggsave("./Fig/FigS3_final.pdf",FigS3,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/FigS3_final.tiff",FigS3,width = 10,height = 8.5,dpi=600, units = "in",bg = "white")

##########################################################################################
######################################## Fig S4 ######################################################
theme_FigS4 <- theme_classic(base_size = 12) +
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8))

p1 <- dominance_LDMC_MI_NutNet + theme_FigS4 + labs(tag = "(a)")
p2 <- dominance_logLDMC_Tg_NutNet + theme_FigS4 + labs(tag = "(b)")
p3 <- dominance_SRL_MI_NutNet + theme_FigS4 + labs(tag = "(c)")
p4 <- dominance_logSRL_Tg_NutNet + theme_FigS4 + labs(tag = "(d)")
p5 <- dominance_Hmax_MI_NutNet + theme_FigS4 + labs(tag = "(e)")
p6 <- dominance_logHmax_Tg_NutNet + theme_FigS4 + labs(tag = "(f)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 4) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 4) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 3.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4) /
                  (p5 | p6)) +
  plot_layout(heights = c(0.12, 1, 1, 1))

# Final figure 

FigS4 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.6, 1.4))
FigS4


# save figure
ggsave("./Fig/FigS4_final.pdf",FigS4,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/FigS4_final.tiff",FigS4,width = 10,height = 8.5,dpi=600, units = "in",bg = "white")

##################### Figure S5 Heat-map Nmass + chi ###############################################
theme_heat_panel <- theme_minimal(base_size = 9) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 8, angle = 25, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 8, hjust = 0.5),
    legend.text  = element_text(size = 7),
    legend.key.height = unit(7, "mm"),
    legend.key.width  = unit(3.5, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.margin = margin(3, 3, 3, 3)
  )
# Scales
Nmass_limits <- c(-0.06, 0.15)
chi_limits   <- c(-0.7, 0.18)

Nmass_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = Nmass_limits,
  oob = squish,
  name = expression(N[mass]),
  breaks = seq(-0.05, 0.15, by = 0.03),
  labels = number_format(accuracy = 0.05),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")
  )
)

chi_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = chi_limits,
  oob = squish,
  name = expression(chi),
  breaks = seq(-0.7, 0.18, by = 0.1),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")
  )
)

# Panels 
reldom_Nmass_N0_nolabel <- reldom_Nmass_N0 +
  labs(title = NULL, subtitle = NULL) +
  Nmass_fill_scale +
  theme_heat_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

reldom_Nmass_N_nolabel <- reldom_Nmass_N +
  labs(title = NULL, subtitle = NULL) +
  Nmass_fill_scale +
  theme_heat_panel +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )

reldom_chi_N0_nolabel <- reldom_chi_N0 +
  labs(title = NULL, subtitle = NULL) +
  chi_fill_scale +
  theme_heat_panel

reldom_chi_N_nolabel <- reldom_chi_N +
  labs(title = NULL, subtitle = NULL) +
  chi_fill_scale +
  theme_heat_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

row_title_Nmass <- ggplot() + theme_void() +
  labs(title = expression(N[mass])) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0)
  )

row_title_chi <- ggplot() + theme_void() +
  labs(title = expression(chi)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0)
  )


###################################### assemble all pannels #####################
## legend
get_legend_grob <- function(p) {
  cowplot::get_legend(p + theme(legend.position = "right"))
}

make_trait_row <- function(trait_label_plot, p_N0, p_N) {
  lg <- get_legend_grob(p_N0)
  p_N0_noleg <- p_N0 + theme(legend.position = "none")
  p_N_noleg  <- p_N  + theme(legend.position = "none")
  (trait_label_plot | p_N0_noleg | p_N_noleg | wrap_elements(lg)) +
    plot_layout(widths = c(0.09, 1, 1, 0.16))
}


# column headers row
col_headers <- (
  plot_spacer() |
    (ggplot() + theme_void() + labs(title = "Ambient N") +
       theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 14))) |
    (ggplot() + theme_void() + labs(title = "Added N") +
       theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 14))) |
    plot_spacer()
) + plot_layout(widths = c(0.09, 1, 1, 0.16))


############### Fig S5 ###########################
FigS5 <-
  col_headers /
  make_trait_row(row_title_Nmass, reldom_Nmass_N0_nolabel, reldom_Nmass_N_nolabel) /
  make_trait_row(row_title_chi,   reldom_chi_N0_nolabel,   reldom_chi_N_nolabel) +
  plot_layout(heights = c(0.10, 1, 1))
FigS5

ggsave("./Fig/FigS5_final.pdf",FigS5,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/FigS5_final.tiff",FigS5,width = 12,height = 7.5, units = "in",bg = "white")
ggsave("./Fig/FigS5_word.png",FigS5,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")



############################ Figure S6 ############################################
##################### SLA + LDMC ###############################################

# Scales
SLA_limits <-  c(-0.04, 0.06)
LDMC_limits   <-  c(-0.21, 0.13)

SLA_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = SLA_limits,
  oob = squish,
  name = expression(SLA),
  breaks = seq(-0.04, 0.06, by = 0.03),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")
  )
)

LDMC_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = LDMC_limits,
  oob = squish,
  name = expression(LDMC),
  breaks = seq(-0.21, 0.13, by = 0.05),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")
  )
)

# Panels 
reldom_SLA_N0_nolabel <- reldom_SLA_N0 +
  labs(title = NULL, subtitle = NULL) +
  SLA_fill_scale +
  theme_heat_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
reldom_SLA_N0_nolabel

reldom_SLA_N_nolabel <- reldom_SLA_N +
  labs(title = NULL, subtitle = NULL) +
  SLA_fill_scale +
  theme_heat_panel +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )
reldom_SLA_N_nolabel

reldom_LDMC_N0_nolabel <- reldom_LDMC_N0 +
  labs(title = NULL, subtitle = NULL) +
  LDMC_fill_scale +
  theme_heat_panel+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

reldom_LDMC_N_nolabel <- reldom_LDMC_N +
  labs(title = NULL, subtitle = NULL) +
  LDMC_fill_scale +
  theme_heat_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

row_title_SLA <- ggplot() + theme_void() +
  labs(title = expression(SLA)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0)
  )

row_title_LDMC <- ggplot() + theme_void() +
  labs(title = expression(LDMC)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0)
  )

#####################  SRL + H ###############################################

# Scales
SRL_limits <-  c(-0.02, 0.04)
H_limits   <-  c(0.03, 0.14)

SRL_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = SRL_limits,
  oob = squish,
  name = expression(SRL),
  breaks = seq(-0.02, 0.04, by = 0.01),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")
  )
)

H_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = H_limits,
  oob = squish,
  name = expression(H),
  breaks = seq(0.03, 0.14, by = 0.05),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")
  )
)

# Panels 
reldom_SRL_N0_nolabel <- reldom_SRL_N0 +
  labs(title = NULL, subtitle = NULL) +
  SRL_fill_scale +
  theme_heat_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
reldom_SRL_N0_nolabel

reldom_SRL_N_nolabel <- reldom_SRL_N +
  labs(title = NULL, subtitle = NULL) +
  SRL_fill_scale +
  theme_heat_panel +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank()
  )
reldom_SRL_N_nolabel

reldom_H_N0_nolabel <- reldom_Hmax_N0 +
  labs(title = NULL, subtitle = NULL) +
  H_fill_scale +
  theme_heat_panel

reldom_H_N_nolabel <- reldom_Hmax_N +
  labs(title = NULL, subtitle = NULL) +
  H_fill_scale +
  theme_heat_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

row_title_SRL <- ggplot() + theme_void() +
  labs(title = expression(SRL)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0)
  )

row_title_H <- ggplot() + theme_void() +
  labs(title = expression(H[max])) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0)
  )

###################################### assemble all pannels #####################
## legend
get_legend_grob <- function(p) {
  cowplot::get_legend(p + theme(legend.position = "right"))
}

make_trait_row <- function(trait_label_plot, p_N0, p_N) {
  lg <- get_legend_grob(p_N0)
  p_N0_noleg <- p_N0 + theme(legend.position = "none")
  p_N_noleg  <- p_N  + theme(legend.position = "none")
  (trait_label_plot | p_N0_noleg | p_N_noleg | wrap_elements(lg)) +
    plot_layout(widths = c(0.09, 1, 1, 0.16))
}


# column headers row
col_headers <- (
  plot_spacer() |
    (ggplot() + theme_void() + labs(title = "Ambient N") +
       theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 14))) |
    (ggplot() + theme_void() + labs(title = "Added N") +
       theme(plot.title = element_text(face = "italic", hjust = 0.5, size = 14))) |
    plot_spacer()
) + plot_layout(widths = c(0.09, 1, 1, 0.16))

FigS6 <-
  col_headers /
  make_trait_row(row_title_LDMC,  reldom_LDMC_N0_nolabel,  reldom_LDMC_N_nolabel) /
  make_trait_row(row_title_SRL,   reldom_SRL_N0_nolabel,   reldom_SRL_N_nolabel) /
  make_trait_row(row_title_H,     reldom_H_N0_nolabel,     reldom_H_N_nolabel) +
  plot_layout(heights = c(0.10, 1, 1, 1))

FigS6

ggsave("./Fig/FigS6_final.pdf",FigS6,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/FigS6_final.tiff",FigS6,width = 12,height = 7.5, units = "in",bg = "white")
ggsave("./Fig/FigS6_word.png",FigS6,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")


##########################################################################################
######################################## Fig S7 ######################################################

theme_figS7 <- theme_classic(base_size = 12) +
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8))

p1 <- lambda_MI_Nmass_rawdata + theme_figS7 + labs(tag = "(a)")
p2 <- lambda_Nmass_rawdata_Tg + theme_figS7 + labs(tag = "(b)")
p3 <- lambda_MI_chi_rawdata + theme_figS7 + labs(tag = "(c)")
p4 <- lambda_chi_rawdata_Tg + theme_figS7 + labs(tag = "(d)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 5) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 5) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 4.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4)) +
  plot_layout(heights = c(0.12, 1, 1))


# Final figure 

FigS7 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.3, 1.7))
FigS7

ggsave("./Fig/FigS7.png",plot = FigS7,width = 12, height = 12, dpi = 300, bg = "white")
ggsave("./Fig/FigS7.pdf",FigS7,width  = 180/25.4, height = 220/25.4, device = cairo_pdf)
ggsave("./Fig/FigS7_word.png",FigS7,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")

##################################### Fig S8 ##################################################
p1 <- dominance_MI_Nmass_rawdata + theme_figS7 + labs(tag = "(a)")
p2 <- reldom_Nmass_rawdata_Tg + theme_figS7 + labs(tag = "(b)")
p3 <- dominance_MI_chi_rawdata + theme_figS7 + labs(tag = "(c)")
p4 <- reldom_chi_rawdata_Tg + theme_figS7 + labs(tag = "(d)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 5) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 5) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 4.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4)) +
  plot_layout(heights = c(0.12, 1, 1))


# Final figure 

FigS8 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.3, 1.7))
FigS8

ggsave("./Fig/FigS8.png",plot = FigS8,width = 12, height = 12, dpi = 300, bg = "white")
ggsave("./Fig/FigS8.pdf",FigS8,width  = 180/25.4, height = 220/25.4, device = cairo_pdf)
ggsave("./Fig/FigS8_word.png",FigS8,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")

######################################## Fig S9 ######################################################

theme_FigS3 <- theme_classic(base_size = 12) +
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8))

p1 <- lambda_MI_SLA_rawdata + theme_FigS3 + labs(tag = "(a)")
p2 <- lambda_SLA_rawdata_Tg + theme_FigS3 + labs(tag = "(b)")
p3 <- dominance_MI_SLA_rawdata + theme_FigS3 + labs(tag = "(c)")
p4 <- reldom_SLA_rawdata_Tg + theme_FigS3 + labs(tag = "(d)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 4) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 4) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 3.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4)) +
  plot_layout(heights = c(0.12, 1, 1))


# Final figure 

FigS9 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.6, 1.4))
FigS9

# save figure
ggsave("./Fig/FigS9.png",plot = FigS9,width = 12, height = 12, dpi = 300, bg = "white")
ggsave("./Fig/FigS9.pdf",FigS9,width  = 180/25.4, height = 220/25.4, device = cairo_pdf)
ggsave("./Fig/FigS9_word.png",FigS9,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")

############################### Fig S10 #################################
p1 <- lambda_MI_LDMC_rawdata + theme_FigS3 + labs(tag = "(a)")
p2 <- lambda_LDMC_rawdata_Tg + theme_FigS3 + labs(tag = "(b)")
p3 <- lambda_MI_SRL_rawdata + theme_FigS3 + labs(tag = "(c)")
p4 <- lambda_SRL_rawdata_Tg + theme_FigS3 + labs(tag = "(d)")
p5 <- lambda_MI_Hmax_rawdata + theme_FigS3 + labs(tag = "(e)")
p6 <- lambda_Hmax_rawdata_Tg + theme_FigS3 + labs(tag = "(f)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 4) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 4) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 3.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4) /
                  (p5 | p6)) +
  plot_layout(heights = c(0.12, 1, 1, 1))

# Final figure 

FigS10 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.6, 1.4))
FigS10

# save figure
ggsave("./Fig/FigS10.png",plot = FigS10,width = 12, height = 12, dpi = 300, bg = "white")
ggsave("./Fig/FigS10.pdf",FigS10,width  = 180/25.4, height = 220/25.4, device = cairo_pdf)
ggsave("./Fig/FigS10_word.png",FigS10,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")

########################## Fig S11 #############################################################

p1 <- dominance_MI_LDMC_rawdata + theme_FigS3 + labs(tag = "(a)")
p2 <- reldom_LDMC_rawdata_Tg + theme_FigS3 + labs(tag = "(b)")
p3 <- dominance_MI_SRL_rawdata + theme_FigS3 + labs(tag = "(c)")
p4 <- reldom_SRL_rawdata_Tg + theme_FigS3 + labs(tag = "(d)")
p5 <- dominance_MI_Hmax_rawdata + theme_FigS3 + labs(tag = "(e)")
p6 <- reldom_Hmax_rawdata_Tg + theme_FigS3 + labs(tag = "(f)")

# combined MI / Tg color legend
legend_df <- data.frame(
  y = c(0.43, 0.32, 0.21),
  colour = c("#7F0000", "#228B22", "#08306B"),
  label = c(
    "MI = 0.2 / Tg = 6 °C",
    "MI = 0.7 / Tg = 12 °C",
    "MI = 1.8 / Tg = 20 °C"))

# One legend
legend_panel <- ggplot() +
  # Soil nitrogen legend
  annotate("text", x = 0.05, y = 0.88,label = "Soil N",
           hjust = 0, fontface = "bold", size = 4.5) +
  annotate( "segment",x = 0.05, xend = 0.27, y = 0.79, yend = 0.79,
            linetype = "dashed", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.79,label = "Ambient N", hjust = 0, size = 4) +
  annotate("segment", x = 0.05, xend = 0.27, y = 0.70, yend = 0.70, linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 0.33, y = 0.70, label = "Added N", hjust = 0, size = 4) +
  # Combined climate legend
  annotate("text", x = 0.05, y = 0.54,label = "Climate",hjust = 0, fontface = "bold", size = 4.5) +
  geom_segment(data = legend_df,aes(  x = 0.05, xend = 0.27, y = y, yend = y,
                                      colour = colour),linewidth = 2.5,show.legend = FALSE) +
  geom_text(data = legend_df,aes(x = 0.33, y = y, label = label),hjust = 0,size = 3.5) +
  scale_colour_identity() + coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off") +
  theme_void() +
  theme(plot.margin = margin(20, 15, 8, 5))

#  column titles

title_MI <- ggplot() +
  annotate( "text", x = 0.5, y = 0.5, label = "MI", fontface = "bold", size = 5) +
  theme_void()

title_Tg <- ggplot() +annotate( "text", x = 0.5, y = 0.5, label = "Tg", fontface = "bold", size = 5) +
  theme_void()

# titles and panels
main_panels <- ((title_MI | title_Tg) /
                  (p1 | p2) /
                  (p3 | p4) /
                  (p5 | p6)) +
  plot_layout(heights = c(0.12, 1, 1, 1))

# Final figure 

FigS11 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.6, 1.4))
FigS11

# save figure
ggsave("./Fig/FigS11.png",plot = FigS11,width = 12, height = 12, dpi = 300, bg = "white")
ggsave("./Fig/FigS11.pdf",FigS11,width  = 180/25.4, height = 220/25.4, device = cairo_pdf)
ggsave("./Fig/FigS11_word.png",FigS11,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")
