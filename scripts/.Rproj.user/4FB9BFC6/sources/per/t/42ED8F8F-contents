install.packages(c("magick", "ggplot2", "patchwork"))

library(magick)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)
library(grid)

##################### Conceptual figure 1 #################################

col_arid  <- "#7F0000"  # Low MI, low Tg
col_humid <- "#08306B"  # High MI, high Tg

make_trait_panel <- function(
    xlab,
    tag,
    slopes,
    intercept_arid = 0.45,
    intercept_humid = 0.45) {
  
  df <- data.frame(
    x = rep(c(0, 1), 4),
    y = c(
      intercept_arid  + slopes["arid_ambient"]  * c(0, 1),
      intercept_humid + slopes["humid_ambient"] * c(0, 1),
      intercept_arid  + slopes["arid_added"]    * c(0, 1),
      intercept_humid + slopes["humid_added"]   * c(0, 1)),
    clim = rep(
      c("Low MI, Low Tg", "High MI, High Tg"),
      times = 2,
      each = 2),
    Ntrt = rep(
      c("Ambient N", "Added N"),
      each = 4))
  
  ggplot(df,aes(x = x,y = y,colour = clim,linetype = Ntrt)) +
  geom_line(linewidth = 1.5) +
  scale_colour_manual(name = NULL,values = c("Low MI, Low Tg" = col_arid,"High MI, High Tg" = col_humid),
      breaks = c("Low MI, Low Tg", "High MI, High Tg"),
      guide = guide_legend(order = 1,
      override.aes = list(linewidth = 1.8, linetype = "solid"))) +
    
    scale_linetype_manual(name = NULL,values = c("Ambient N" = "dashed","Added N" = "solid"),
      breaks = c("Ambient N", "Added N"),
      guide = guide_legend(order = 2,
    override.aes = list(linewidth = 1.8, colour = "grey25"))) +
    labs( x = xlab, y = expression(lambda),tag = tag) +
    coord_cartesian(xlim = c(0, 1),ylim = c(0, 1.15), clip = "off") +
    theme_classic(base_size = 18) +
    theme(axis.ticks = element_blank(),axis.text  = element_blank(),
      axis.line = element_line( linewidth = 1.1, colour = "black"),
      axis.title = element_text(size = 20,colour = "black"),
      axis.title.x = element_text(margin = margin(t = 8)),
      axis.title.y = element_text(margin = margin(r = 8)),
      plot.tag = element_text(face = "bold",size = 20),
      plot.margin = margin(12, 12, 12, 12))
}

p_Nmass <- make_trait_panel(
  xlab = expression(N[mass]),
  tag = "(a)",
  slopes = c(
    arid_ambient  = 0.50,
    humid_ambient = -0.25,
    arid_added    = 0.25,
    humid_added   = -0.10))

p_chi <- make_trait_panel(
  xlab = expression(italic(chi)),
  tag = "(b)",
  slopes = c(
    arid_ambient  = -0.30,
    humid_ambient = 0.50,
    arid_added    = -0.10,
    humid_added   = 0.25))

p_LDMC <- make_trait_panel(
  xlab = expression(LDMC),
  tag = "(c)",
  slopes = c(
    arid_ambient  = 0.60,
    humid_ambient = -0.20,
    arid_added    = 0.35,
    humid_added   = -0.10))


p_SRL <- make_trait_panel(
  xlab = expression(SRL),
  tag = "(d)",
  slopes = c(
    arid_ambient  = 0.50,
    humid_ambient = -0.25,
    arid_added    = 0.25,
    humid_added   = -0.10))

p_Hmax <- make_trait_panel(
  xlab = expression(H[max]),
  tag = "(e)",
  slopes = c(
    arid_ambient  = -0.10,
    humid_ambient = 0.45,
    arid_added    = -0.05,
    humid_added   = 0.65))

# Combine panels + legend 
Fig1 <- ( p_Nmass + p_chi + p_LDMC +
    p_SRL   + p_Hmax + guide_area()) +
  plot_layout(ncol = 3,guides = "collect") &
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    legend.key.width  = unit(2, "cm"),
    legend.key.height = unit(0.7, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    legend.margin = margin(5, 5, 5, 5),
    plot.tag = element_text(
      face = "bold",
      size = 20),
    plot.margin = margin(15, 15, 15, 15))

Fig1

ggsave("./Fig/Fig1_final.pdf",Fig1,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/Fig1_final.tiff",Fig1,width = 12.5,height = 7,units = "in",bg = "white")

######################################## Fig 2 ######################################################
theme_fig2 <- theme_classic(base_size = 12) +
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8))

p1 <- lambda_Nmass_MI_NutNet + theme_fig2 + labs(tag = "(a)")
p2 <- lambda_logNmass_Tg_NutNet + theme_fig2 + labs(tag = "(b)")
p3 <- lambda_chi_MI_NutNet + theme_fig2 + labs(tag = "(c)")
p4 <- lambda_chi_Tg_NutNet + theme_fig2 + labs(tag = "(d)")

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

Fig2 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.3, 1.7))
Fig2

ggsave("./Fig/Fig2_final.tiff", plot = Fig2,width = 12,height = 7.5,units = "in",dpi = 600,compression = "lzw",bg = "white")
ggsave("./Fig/Fig2_final.pdf", plot = Fig2,width = 10,height = 8.5,units = "in",dpi = 600,bg = "white")# save figure
ggsave("./Fig/Fig2_final.png", plot = Fig2,width = 10,height = 8.5,units = "in",dpi = 600,bg = "white")


##########################################################################################
######################################## Fig 3 ######################################################
theme_fig2 <- theme_classic(base_size = 12) +
  theme(
    plot.tag.position = c(0.02, 0.98),
    plot.tag = element_text(face = "bold", size = 14),
    legend.position = "none",
    axis.title.y = element_text(margin = margin(r = 6)),
    axis.title.x = element_text(margin = margin(t = 6)),
    axis.text = element_text(size = 11),
    plot.margin = margin(8, 10, 8, 8))

p1 <- lambda_LDMC_MI_NutNet + theme_fig2 + labs(tag = "(a)")
p2 <- lambda_logLDMC_Tg_NutNet + theme_fig2 + labs(tag = "(b)")
p3 <- lambda_SRL_MI_NutNet + theme_fig2 + labs(tag = "(c)")
p4 <- lambda_logSRL_Tg_NutNet + theme_fig2 + labs(tag = "(d)")
p5 <- lambda_Hmax_MI_NutNet + theme_fig2 + labs(tag = "(e)")
p6 <- lambda_logHmax_Tg_NutNet + theme_fig2 + labs(tag = "(f)")

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

Fig3 <- (main_panels | legend_panel) +
  plot_layout(widths = c(3.6, 1.4))
Fig3


# save figure
ggsave("./Fig/Fig3_final.pdf",Fig3,width = 12,height = 7.5,units = "in",bg = "white")
ggsave("./Fig/Fig3_final.tiff",Fig3,width = 10,height = 8.5,dpi=600, units = "in",bg = "white")

##################### Figure 4 Heat-map Nmass + chi ###############################################
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
    plot.margin = margin(3, 3, 3, 3))

# Scales
Nmass_limits <- c(-0.1, 0.4)
chi_limits   <- c(-1.5, 0.5)

Nmass_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = Nmass_limits,
  oob = squish,
  name = expression(N[mass]),
  breaks = seq(-0.1, 0.4, by = 0.07),
  labels = number_format(accuracy = 0.1),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")))

chi_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = chi_limits,
  oob = squish,
  name = expression(chi),
  breaks = seq(-1.5, 0.5, by = 0.5),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")))

# Panels 
HM_Nmass_N0_nolabel <- lambda_Nmass_N0 +
  labs(title = NULL, subtitle = NULL) +
  Nmass_fill_scale +
  theme_heat_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

HM_Nmass_N_nolabel <- lambda_Nmass_N +
  labs(title = NULL, subtitle = NULL) +
  Nmass_fill_scale +
  theme_heat_panel +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank())

HM_chi_N0_nolabel <- lambda_chi_N0 +
  labs(title = NULL, subtitle = NULL) +
  chi_fill_scale +
  theme_heat_panel

HM_chi_N_nolabel <- lambda_chi_N +
  labs(title = NULL, subtitle = NULL) +
  chi_fill_scale +
  theme_heat_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

row_title_Nmass <- ggplot() + theme_void() +
  labs(title = expression(N[mass])) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0))

row_title_chi <- ggplot() + theme_void() +
  labs(title = expression(chi)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0))

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
    plot_spacer()) + plot_layout(widths = c(0.09, 1, 1, 0.16))

############### Fig 4 ###########################
Fig4 <-
  col_headers /
  make_trait_row(row_title_Nmass, HM_Nmass_N0_nolabel, HM_Nmass_N_nolabel) /
  make_trait_row(row_title_chi,   HM_chi_N0_nolabel,   HM_chi_N_nolabel) +
  plot_layout(heights = c(0.10, 1, 1))
Fig4

ggsave("./Fig/Fig4_final.pdf",Fig4,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/Fig4_final.tiff",Fig4,width = 12,height = 7.5, units = "in",bg = "white")
ggsave("./Fig/Fig4_word.png",Fig4,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")

############################ Figure 5 ############################################
##################### SLA + LDMC ###############################################

# Scales
SLA_limits <-  c(-0.04, 0.06)
LDMC_limits   <-  c(-0.4, 0.3)

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
    barwidth  = unit(3.5, "mm")))

LDMC_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = LDMC_limits,
  oob = squish,
  name = expression(LDMC),
  breaks = seq(-0.4, 0.3, by = 0.1),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")))

# Panels 
HM_SLA_N0_nolabel <- lambda_SLA_N0 +
  labs(title = NULL, subtitle = NULL) +
  SLA_fill_scale +
  theme_heat_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
HM_SLA_N0_nolabel

HM_SLA_N_nolabel <- lambda_SLA_N +
  labs(title = NULL, subtitle = NULL) +
 SLA_fill_scale +
  theme_heat_panel +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank())
HM_SLA_N_nolabel

HM_LDMC_N0_nolabel <- lambda_LDMC_N0 +
  labs(title = NULL, subtitle = NULL) +
  LDMC_fill_scale +
  theme_heat_panel+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

HM_LDMC_N_nolabel <- lambda_LDMC_N +
  labs(title = NULL, subtitle = NULL) +
  LDMC_fill_scale +
  theme_heat_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())+
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

row_title_SLA <- ggplot() + theme_void() +
  labs(title = expression(SLA)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0))

row_title_LDMC <- ggplot() + theme_void() +
  labs(title = expression(LDMC)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0))

#####################  SRL + H ###############################################
# Scales
SRL_limits <-  c(-0.1, 0.15)
H_limits   <-  c(0, 0.4)

SRL_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = SRL_limits,
  oob = squish,
  name = expression(SRL),
  breaks = seq(-0.1, 0.15, by = 0.05),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")))

H_fill_scale <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  limits = H_limits,
  oob = squish,
  name = expression(H[max]),
  breaks = seq(0, 0.4, by = 0.05),
  labels = number_format(accuracy = 0.01),
  guide = guide_colorbar(
    title.position = "top",
    title.hjust = 0.5,
    ticks = TRUE,
    frame.colour = "grey60",
    barheight = unit(25, "mm"),
    barwidth  = unit(3.5, "mm")))

# Panels 
HM_SRL_N0_nolabel <- lambda_SRL_N0 +
  labs(title = NULL, subtitle = NULL) +
  SRL_fill_scale +
  theme_heat_panel +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
HM_SRL_N0_nolabel

HM_SRL_N_nolabel <- lambda_SRL_N +
  labs(title = NULL, subtitle = NULL) +
  SRL_fill_scale +
  theme_heat_panel +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    axis.text.y = element_blank(), axis.ticks.y = element_blank())
HM_SRL_N_nolabel

HM_H_N0_nolabel <- lambda_Hmax_N0 +
  labs(title = NULL, subtitle = NULL) +
  H_fill_scale +
  theme_heat_panel

HM_H_N_nolabel <- lambda_Hmax_N +
  labs(title = NULL, subtitle = NULL) +
  H_fill_scale +
  theme_heat_panel +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())

row_title_SRL <- ggplot() + theme_void() +
  labs(title = expression(SRL)) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0))

row_title_H <- ggplot() + theme_void() +
  labs(title = expression(H[max])) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0, margin = margin(b = 1)),
    plot.margin = margin(0, 0, 0, 0))

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

Fig5 <-
  col_headers /
  make_trait_row(row_title_LDMC,  HM_LDMC_N0_nolabel,  HM_LDMC_N_nolabel) /
  make_trait_row(row_title_SRL,   HM_SRL_N0_nolabel,   HM_SRL_N_nolabel) /
  make_trait_row(row_title_H,     HM_H_N0_nolabel,     HM_H_N_nolabel) +
  plot_layout(heights = c(0.10, 1, 1, 1))

Fig5
ggsave("./Fig/Fig5_final.pdf",Fig5,width = 12.5,height = 7,units = "in",bg = "white")
ggsave("./Fig/Fig5_final.tiff",Fig5,width = 12,height = 7.5, units = "in",bg = "white")
ggsave("./Fig/Fig5_word.png",Fig5,width  = 200/25.4,height = 220/25.4,dpi= 600,bg= "white")

