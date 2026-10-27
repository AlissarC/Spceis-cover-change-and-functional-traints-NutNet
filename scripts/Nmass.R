#########################################
# necessary packages 
#########################################
R.version.string
library(tidyverse)  # Data manipulation, string processing, and plots
library(ggplot2)    # Figures, heatmaps, and graphical customization
library(scales)     # Color transparency and numerical formatting
########################################################
# Statistics
########################################################
library(lme4)       # Linear mixed-effects models
library(car)        # ANOVA and model diagnostics
library(emmeans)    # Estimated marginal means and slopes
library(multcomp)   # Multiple comparisons and CLD
library(MuMIn)      # Marginal and conditional R²
########################################################
# 3.Tables Export and formatting
########################################################

library(officer)    # Create and edit Word documents
library(flextable)  # Format and export statistical tables

###################################### Nmass ##########################################################
stat_Nmass = read.csv("./output/final_df.csv")
names(stat_Nmass)
table(stat_Nmass$growth_form)
length(unique(stat_Nmass$species)) #2059

n_distinct(stat_Nmass$site_code) # 94
length(unique(stat_Nmass$species[!is.na(stat_Nmass$Nmass_NutNet)])) # 1514

stat_Nmass %>%
  distinct(site_code, block) %>%
  count(site_code, name = "n_blocks") # some sites have only one block ! ex: pape.de, ufrec.us, sevis.us

stat_Nmass %>%
  distinct(site_code, trt) %>%
  count(site_code, name = "n_trt") # some sites have only one block ! ex: pape.de, ufrec.us, sevis.us

################### Natural log transforming of traits ##########################
plot(stat_Nmass$aridity, stat_Nmass$tmp_GS)
stat_Nmass$logNmass = log(stat_Nmass$Nmass_NutNet)
hist(stat_Nmass$logNmass)

max(stat_Nmass$tmp_GS) # 27.625
min(stat_Nmass$tmp_GS) # 3.402083

max(stat_Nmass$aridity) # 2.43
min(stat_Nmass$aridity) # 0.0644

hist(stat_Nmass$aridity)

######################################### clean ###############################################################
names(stat_Nmass)
nrow(stat_Nmass) # 44419

stat_Nmass_clean <- stat_Nmass %>% 
  filter(!is.na(logNmass), !is.na(mean_log_lambda_cover_all))
nrow(stat_Nmass_clean) #  37967

################################################################### LMERs for Nmass ############################################
######################################################## Nmass NutNet analysis: log lambda corrected by K ##############################################
names(stat_Nmass_clean)

stat_Nmass_clean %>%
  dplyr::distinct(species, growth_form) %>%
  dplyr::count(growth_form)
#forb 1119
#graminoid  395

lambda_Nmass_NutNet <- lmer(mean_log_lambda_cover_all~logNmass + logNmass:(aridity*Ntrt+ tmp_GS*Ntrt) + aridity + tmp_GS + 
                            logNmass:(Ntrt*Ptrt*Ktrt) + n_years +
                            (1|site_code) + (1|site_code:block),
                          data=stat_Nmass_clean)
Anova(lambda_Nmass_NutNet) # effects: Nmass , Ntrt*Nmass, Ptrt*Nmass, aridity*Nmass*Ntrt
vif(lambda_Nmass_NutNet)

confint(lambda_Nmass_NutNet, method = "Wald")
summary(lambda_Nmass_NutNet)
r.squaredGLMM(lambda_Nmass_NutNet)              
#R2m       R2c
#0.005402666 0.1005177
#emm_options(pbkrtest.limit = 51099)
AIC(lambda_Nmass_NutNet) #   132001.1

res <- residuals(lambda_Nmass_NutNet)
qqnorm(res)
qqline(res, col = "red")

hist(res, breaks = 20, main = "Histogram of Residuals") ## good

plot(fitted(lambda_Nmass_NutNet), res, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity :OK

densityPlot(res)

ks.test(res, "pnorm", mean = mean(res), sd = sd(res))

######################### Export Anova table ######################################

# Define the p-value formatting function
fmt_p <- function(p) {
  ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", ifelse(p < 0.01, "<0.01", ifelse(p < 0.05, "<0.05",sprintf("%.3f", p)))))
}

lambda_Nmass_table = as.data.frame(Anova(lambda_Nmass_NutNet))
lambda_Nmass_model_NutNet <- data.frame(
  Var   = c(
    'Nmass', 'MI', 'Tg', 'number of years','Nmass x MI', 
    'Nmass x Soil N', 'Nmass x Tg','Nmass x Soil P',
    'Nmass x Soil K', 'Nmass x MI x Soil N', 'Nmass x Tg x Soil N', 'Nmass x Soil N x Soil P',
    'Nmass x Soil N x Soil K','Nmass x Soil P x Soil K', 'Nmass x Soil N x Soil P x Soil K '),
  chisq = round(lambda_Nmass_table[1:15, 1], 2),
  df    = lambda_Nmass_table[1:15, 2],
  P     = fmt_p(lambda_Nmass_table[1:15, 3])   
)

# Create flextable
ft <- flextable(lambda_Nmass_model_NutNet) %>%
  autofit()

# Export to Word
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "./Tables/lambda_Nmass_model_NutNet.docx")

################################ Extract aridity slopes and intercepts ################################
################ get slopes and of they are significantly different from zero##### 
#  MI levels from the observed distribution (quantiles): 10%, 50%, 95%
aridity_levels_quantile <- quantile(
  stat_Nmass_clean$aridity,
  probs = c(0.1, 0.5, 0.96),
  na.rm = TRUE)

aridity_levels_num <- as.numeric(aridity_levels_quantile)
aridity_levels <- round(aridity_levels_num, 1) 

lambda_slopes_MI_Nmass_NutNet  <- summary(
  emtrends(lambda_Nmass_NutNet , ~ aridity * Ntrt, 
           var = "logNmass", 
           at = list(aridity = aridity_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
lambda_slopes_MI_Nmass_NutNet 

########### comp between slopes ################
lambda_slope_comp_MI_Nmass_NutNet  <- emtrends(lambda_Nmass_NutNet , ~ aridity * Ntrt, 
                                               var = "logNmass",
                                               at = list(aridity = aridity_levels, Ntrt = c(0, 1)))
lambda_slope_comp_MI_Nmass_NutNet 

lambda_comp_slopes_MI_Nmass_NutNet  <- cld(lambda_slope_comp_MI_Nmass_NutNet , Letters = letters, adjust = "tukey")
lambda_comp_slopes_MI_Nmass_NutNet 


slopes_lambda_Ntrt_Nmass_NutNet  <- summary(
  emtrends(lambda_Nmass_NutNet , ~ Ntrt, 
           var = "logNmass", infer = c(TRUE, TRUE)))
slopes_lambda_Ntrt_Nmass_NutNet 


slopes_lambda_Ptrt_Nmass_NutNet  <- summary(
  emtrends(lambda_Nmass_NutNet , ~ Ptrt, 
           var = "logNmass", infer = c(TRUE, TRUE)))
slopes_lambda_Ptrt_Nmass_NutNet 

#################### export tables ############################################
##############################################################################

################# Slopes ##############################
lambda_slopes_MI_Nmass_NutNet_df <- as.data.frame(lambda_slopes_MI_Nmass_NutNet) %>%
  
  transmute(
    aridity,
    Ntrt,
    slope_est = logNmass.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
lambda_slopes_MI_Nmass_NutNet_df

lambda_comp_slopes_MI_Nmass_NutNet_groups <- as.data.frame(lambda_comp_slopes_MI_Nmass_NutNet) %>%
  transmute(aridity, Ntrt, slope_grp = str_squish(.group))
lambda_comp_slopes_MI_Nmass_NutNet_groups
# Merge slope stats and letters
lambda_slopes_Nmass_NutNet_all <- left_join(lambda_slopes_MI_Nmass_NutNet_df, lambda_comp_slopes_MI_Nmass_NutNet_groups, by = c("aridity","Ntrt"))
##################################################################################################################

##### formating ###########################
#Keep helper columns for composing
Table_lambda_Nmass_NutNet_ft <- lambda_slopes_Nmass_NutNet_all %>%
  transmute(
    MI = aridity,
    N  = Ntrt,
    slope_est = as.numeric(slope_est),
    slope_SE  = as.numeric(slope_SE),
    slope_p   = as.numeric(slope_p),
    slope_grp = as.character(slope_grp)
    
  ) %>%
  mutate(
    `Slopes ± SE`     = "",
    `Slopes P-value`        = fmt_p(slope_p),
     )
Table_lambda_Nmass_NutNet_ft


## check format and columns
Table_lambda_Nmass_NutNet_out <- Table_lambda_Nmass_NutNet_ft %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)

  ) %>%
  
relocate(MI, N, `Slopes ± SE`, `Slopes P-value`)
Table_lambda_Nmass_NutNet_out

ft <- flextable(Table_lambda_Nmass_NutNet_out)

for (i in seq_len(nrow(Table_lambda_Nmass_NutNet_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_lambda_Nmass_NutNet_out$slope_est[i])),
      as_chunk(Table_lambda_Nmass_NutNet_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_lambda_Nmass_NutNet_out$slope_SE[i]))))}

helper_cols <- c("slope_est","slope_SE","slope_p","slope_grp")

if ("delete_columns" %in% getNamespaceExports("flextable")) {
  ft <- delete_columns(ft, j = helper_cols)
} else {
  keep <- c("MI","N","Slopes ± SE","Slopes P-value")
  ft$col_keys <- keep
}

ft <- bold(ft, part = "header")
ft <- fontsize(ft, size = 12, part = "all")
ft <- align(ft, align = "center", part = "all")
ft <- theme_box(ft)
ft <- autofit(ft)

# Export
doc <- read_docx()
doc <- body_add_par(doc, "Nmass", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_lambda_Nmass_NutNet.docx")

###############################################################################################################################
############################### Figure lambda cover_Nmass_MI #############################################
Nmass_grid <- seq(min(stat_Nmass_clean$logNmass, na.rm = TRUE),
  max(stat_Nmass_clean$logNmass, na.rm = TRUE),
  length.out = 300)

emm_df_Nmass <- as.data.frame(
  emmeans(
    lambda_Nmass_NutNet,
    ~ logNmass | aridity * Ntrt,
    at = list(
      logNmass = Nmass_grid,
      aridity  = aridity_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_Nmass

pred_df_Nmass <- emm_df_Nmass %>%
  mutate(
    lambda_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    MI = factor(
      aridity,
      levels = aridity_levels,
      labels = paste0("MI = ", aridity_levels)))
names(pred_df_Nmass)

# Center where polygons meet: at the median
x_med <- median(stat_Nmass_clean$logNmass, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_Nmass <- pred_df_Nmass %>%
  group_by(Ntrt, MI) %>%
  arrange(logNmass, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logNmass),
    x_max   = last(logNmass),
    fit_min = first(lambda_fit),
    se_min  = first(SE),
    fit_max = last(lambda_fit),
    se_max  = last(SE),
    fit_med = approx(x = logNmass, y = lambda_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_Nmass <- pmap_dfr(poly_df_Nmass, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                        fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})

poly_right_Nmass <- pmap_dfr(poly_df_Nmass, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                         fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})

poly_plot_Nmass <- bind_rows(poly_left_Nmass, poly_right_Nmass) %>%
  group_by(Ntrt, MI, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_Nmass
dplyr::count(poly_plot_Nmass, id)
# MI colors by level (low/mid/high) 
MI_cols <- c("#7F0000", "#228B22","#08306B")

names(MI_cols) <- levels(pred_df_Nmass$MI)

# Plot: polygons + 6 lines (3 MI x 2 Ntrt)
lambda_Nmass_MI_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_Nmass,
    aes(x = x, y = y, group = id, fill = MI),
    alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_Nmass,
            aes(x = logNmass, y = lambda_fit,
                color = MI,
                linetype = Ntrt,
                linewidth = Ntrt,
                group = interaction(MI, Ntrt))) +
  
  scale_linewidth_manual(
    values = c("0" = 1, "1" = 1),
    guide = "none"
  ) +
  
  scale_linetype_manual(
    name   = "Soil N",
    values = c("0" = "dashed", "1" = "solid"),
    labels = c("0" = "Ambient N", "1" = "Added N")) +
  
  scale_color_manual(name = "MI", values = MI_cols) +
  scale_fill_manual (name = "MI", values = MI_cols) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["cover"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

lambda_Nmass_MI_NutNet

############################### Extract Tg  slopes ################################

################ get slopes and of they are significantly different from zero##### 
tmp_levels_quantile <- quantile(
  stat_Nmass_clean$tmp_GS,
  probs = c(0.1, 0.5, 0.96),
  na.rm = TRUE)
tmp_levels_num <- as.numeric(tmp_levels_quantile)
tmp_levels <- round(tmp_levels_num, 0) 


lambda_slopes_tmp_logNmass_NutNet  <- summary(
  emtrends(lambda_Nmass_NutNet , ~ tmp_GS * Ntrt, 
           var = "logNmass", 
           at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
lambda_slopes_tmp_logNmass_NutNet 

########### comp between slopes ################
lambda_slope_comp_tmp_logNmass_NutNet  <- emtrends(lambda_Nmass_NutNet , ~ tmp_GS * Ntrt, 
                                              var = "logNmass",
                                              at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1)))
lambda_slope_comp_tmp_logNmass_NutNet 

lambda_comp_slopes_tmp_logNmass_NutNet  <- cld(lambda_slope_comp_tmp_logNmass_NutNet , Letters = letters, adjust = "tukey")
lambda_comp_slopes_tmp_logNmass_NutNet 


#################### export tables ############################################
##############################################################################

################# Slopes ##############################
lambda_slopes_tmp_logNmass_NutNet_df <- as.data.frame(lambda_slopes_tmp_logNmass_NutNet) %>%
  
  transmute(
    tmp_GS,
    Ntrt,
    slope_est = logNmass.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
lambda_slopes_tmp_logNmass_NutNet_df

lambda_comp_slopes_tmp_logNmass_NutNet_groups <- as.data.frame(lambda_comp_slopes_tmp_logNmass_NutNet) %>%
  transmute(tmp_GS , Ntrt, slope_grp = str_squish(.group))
lambda_comp_slopes_tmp_logNmass_NutNet_groups
# Merge slope stats and letters
lambda_slopes_logNmass_NutNet_all <- left_join(lambda_slopes_tmp_logNmass_NutNet_df, lambda_comp_slopes_tmp_logNmass_NutNet_groups, by = c("tmp_GS","Ntrt"))
##################################################################################################################
### formating ###########################
#Keep helper columns for composing
Table_lambda_logNmass_NutNet_ft_tmp <- lambda_slopes_logNmass_NutNet_all %>%
  transmute(
    Tg = tmp_GS,
    N  = Ntrt,
    slope_est = as.numeric(slope_est),
    slope_SE  = as.numeric(slope_SE),
    slope_p   = as.numeric(slope_p),
    slope_grp = as.character(slope_grp)
  ) %>%
  mutate(
    `Slopes ± SE`     = "",
    `Slopes P-value`        = fmt_p(slope_p))
Table_lambda_logNmass_NutNet_ft_tmp

## check format and columns
Table_lambda_logNmass_NutNet_ft_tmp_out <- Table_lambda_logNmass_NutNet_ft_tmp %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)) %>%
  
  relocate(Tg, N, `Slopes ± SE`, `Slopes P-value`)
Table_lambda_logNmass_NutNet_ft_tmp_out

# build flextable
ft <- flextable(Table_lambda_logNmass_NutNet_ft_tmp_out)

for (i in seq_len(nrow(Table_lambda_logNmass_NutNet_ft_tmp_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_lambda_logNmass_NutNet_ft_tmp_out$slope_est[i])),
      as_chunk(Table_lambda_logNmass_NutNet_ft_tmp_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_lambda_logNmass_NutNet_ft_tmp_out$slope_SE[i]))
    )
  )
}

helper_cols <- c("slope_est","slope_SE","slope_p","slope_grp")

if ("delete_columns" %in% getNamespaceExports("flextable")) {
  ft <- delete_columns(ft, j = helper_cols)
} else {
  keep <- c("Tg","N","Slopes ± SE","Slopes P-value")
  ft$col_keys <- keep
}

# Minimal styling:
ft <- bold(ft, part = "header")
ft <- fontsize(ft, size = 12, part = "all")
ft <- align(ft, align = "center", part = "all")
ft <- theme_box(ft)
ft <- autofit(ft)

# Export
doc <- read_docx()
doc <- body_add_par(doc, "logNmass", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_lambda_tmp_logNmass_NutNet.docx")

##################################################################################################
logNmass_grid <- seq(min(stat_Nmass_clean$logNmass, na.rm = TRUE),
                max(stat_Nmass_clean$logNmass, na.rm = TRUE),
                length.out = 300)
logNmass_grid

emm_df_logNmass_tg <- as.data.frame(
  emmeans(
    lambda_Nmass_NutNet,
    ~ logNmass | tmp_GS * Ntrt,
    at = list(
      logNmass = logNmass_grid,
      tmp_GS  = tmp_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_logNmass_tg  

# prediction dataframe : fit + Tg factor (3 levels)
pred_df_logNmass_tg   <- emm_df_logNmass_tg   %>%
  mutate(
    lambda_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    Tg = factor(
      tmp_GS,
      levels = tmp_levels,
      labels = paste0("Tg = ", tmp_levels)))
names(pred_df_logNmass_tg)

# Center where polygons meet: at the median
x_med <- median(stat_Nmass_clean$logNmass, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_logNmass_tg <- pred_df_logNmass_tg %>%
  group_by(Ntrt, Tg) %>%
  arrange(logNmass, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logNmass),
    x_max   = last(logNmass),
    fit_min = first(lambda_fit),
    se_min  = first(SE),
    fit_max = last(lambda_fit),
    se_max  = last(SE),
    fit_med = approx(x = logNmass, y = lambda_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_logNmass_tg <- pmap_dfr(poly_df_logNmass_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                      fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})
poly_left_logNmass_tg

poly_right_logNmass_tg <- pmap_dfr(poly_df_logNmass_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                       fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})
poly_right_logNmass_tg


poly_plot_logNmass_tg <- bind_rows(poly_left_logNmass_tg, poly_right_logNmass_tg) %>%
  group_by(Ntrt, Tg, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_logNmass_tg


# Tg colors by level (low/mid/high) 
#Tg_cols <- c("#2A9D8F", "#E9C46A", "#F4A261")
Tg_cols <- c("#7F0000", "#228B22","#08306B")

names(Tg_cols) <- levels(pred_df_logNmass_tg$tmp_GS)

# Plot: polygons + 6 lines (3 Tg x 2 Ntrt)
lambda_logNmass_Tg_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_logNmass_tg,
               aes(x = x, y = y, group = id, fill = Tg),
               alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_logNmass_tg,
            aes(x = logNmass, y = lambda_fit,
                color = Tg, linetype = Ntrt,
                group = interaction(Tg, Ntrt)),
            linewidth = 1) +
  
  scale_linetype_manual(
    name   = "Soil N",
    values = c("0" = "dashed", "1" = "solid"),
    labels = c("0" = "Ambient N", "1" = "Added N")) +
  
  scale_color_manual(name = "Tg", values = Tg_cols) +
  scale_fill_manual (name = "Tg", values = Tg_cols) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["cover"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

lambda_logNmass_Tg_NutNet

###################################################################################

######################################################## Nmass NutNet analysis: log lambda relative dominance corrected by K ##############################################
names(stat_Nmass_clean)

dominance_Nmass_NutNet <- lmer(mean_log_lambda_Ci_all~logNmass + logNmass:(aridity*Ntrt+ tmp_GS*Ntrt) + aridity + tmp_GS + 
                              logNmass:(Ntrt*Ptrt*Ktrt) + n_years +
                              (1|site_code) + (1|site_code:block),
                            data=stat_Nmass_clean)
Anova(dominance_Nmass_NutNet)

summary(dominance_Nmass_NutNet)
r.squaredGLMM(dominance_Nmass_NutNet)              
#R2m       R2c
#0.004957285 0.1264718
AIC(dominance_Nmass_NutNet) # 72688.7

res <- residuals(dominance_Nmass_NutNet)
qqnorm(res)
qqline(res, col = "red")

hist(res, breaks = 20, main = "Histogram of Residuals") ## good

plot(fitted(dominance_Nmass_NutNet), res, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity :OK

densityPlot(res)

ks.test(res, "pnorm", mean = mean(res), sd = sd(res))

######################### Export Anova table #######################################

# Define the p-value formatting function
fmt_p <- function(p) {
  ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", ifelse(p < 0.01, "<0.01", ifelse(p < 0.05, "<0.05",sprintf("%.3f", p)))))
}
dominance_Nmass_NutNet_table = as.data.frame(Anova(dominance_Nmass_NutNet))

dominance_Nmass_model_NutNet <- data.frame(
  Var   = c(
    'Nmass', 'MI', 'Tg', 'Number of years', 'Nmass x MI', 
    'Nmass x Soil N', 'Nmass x Tg','Nmass x Soil P',
    'Nmass x Soil K', 'Nmass x MI x Soil N', 'Nmass x Tg x Soil N', 'Nmass x Soil N x Soil P',
    'Nmass x Soil N x Soil K','Nmass x Soil P x Soil K', 'Nmass x Soil N x Soil P x Soil K '),
  chisq = round(dominance_Nmass_NutNet_table[1:15, 1], 2),
  df    = dominance_Nmass_NutNet_table[1:15, 2],
  P     = fmt_p(dominance_Nmass_NutNet_table[1:15, 3])   
)

# Create flextable
ft <- flextable(dominance_Nmass_model_NutNet) %>%
  autofit()

# Export to Word
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "./Tables/dominance_Nmass_model_NutNet.docx")

################################ Extract aridity slopes and intercepts ################################

################ get slopes and of they are significantly different from zero##### 
dominance_slopes_MI_Nmass_NutNet  <- summary(
  emtrends(dominance_Nmass_NutNet , ~ aridity * Ntrt, 
           var = "logNmass", 
           at = list(aridity = aridity_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
dominance_slopes_MI_Nmass_NutNet 

########### comp between slopes ################
dominance_slope_comp_MI_Nmass_NutNet  <- emtrends(dominance_Nmass_NutNet , ~ aridity * Ntrt, 
                                               var = "logNmass",
                                               at = list(aridity = aridity_levels, Ntrt = c(0, 1)))
dominance_slope_comp_MI_Nmass_NutNet 

dominance_comp_slopes_MI_Nmass_NutNet  <- cld(dominance_slope_comp_MI_Nmass_NutNet , Letters = letters, adjust = "tukey")
dominance_comp_slopes_MI_Nmass_NutNet 


slopes_dominance_Ntrt_Nmass_NutNet  <- summary(
  emtrends(dominance_Nmass_NutNet , ~ Ntrt, 
           var = "logNmass", infer = c(TRUE, TRUE)))
slopes_dominance_Ntrt_Nmass_NutNet 


slopes_dominance_Ptrt_Nmass_NutNet  <- summary(
  emtrends(dominance_Nmass_NutNet , ~ Ptrt, 
           var = "logNmass", infer = c(TRUE, TRUE)))
slopes_dominance_Ptrt_Nmass_NutNet 

#################### export tables ############################################
##############################################################################
dominance_slopes_MI_Nmass_NutNet_df <- as.data.frame(dominance_slopes_MI_Nmass_NutNet) %>%
  
  transmute(
    aridity,
    Ntrt,
    slope_est = logNmass.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
dominance_slopes_MI_Nmass_NutNet_df

dominance_comp_slopes_MI_Nmass_NutNet_groups <- as.data.frame(dominance_comp_slopes_MI_Nmass_NutNet) %>%
  transmute(aridity, Ntrt, slope_grp = str_squish(.group))
dominance_comp_slopes_MI_Nmass_NutNet_groups
# Merge slope stats and letters
dominance_slopes_Nmass_NutNet_all <- left_join(dominance_slopes_MI_Nmass_NutNet_df, dominance_comp_slopes_MI_Nmass_NutNet_groups, by = c("aridity","Ntrt"))
##################################################################################################################
### formating ###########################
#Keep helper columns for composing
Table_dominance_Nmass_NutNet_ft <- dominance_slopes_Nmass_NutNet_all %>%
  transmute(
    MI = aridity,
    N  = Ntrt,
    slope_est = as.numeric(slope_est),
    slope_SE  = as.numeric(slope_SE),
    slope_p   = as.numeric(slope_p),
    slope_grp = as.character(slope_grp)
  ) %>%
  mutate(
    `Slopes ± SE`     = "",
    `Slopes P-value`        = fmt_p(slope_p)
  )
Table_dominance_Nmass_NutNet_ft

fmt_p <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

## check format and columns
Table_dominance_Nmass_NutNet_out <- Table_dominance_Nmass_NutNet_ft %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)
   
  ) %>%
  
relocate(MI, N, `Slopes ± SE`, `Slopes P-value`)
Table_dominance_Nmass_NutNet_out

# build flextable
ft <- flextable(Table_dominance_Nmass_NutNet_out)

for (i in seq_len(nrow(Table_dominance_Nmass_NutNet_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_dominance_Nmass_NutNet_out$slope_est[i])),
      as_chunk(Table_dominance_Nmass_NutNet_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_dominance_Nmass_NutNet_out$slope_SE[i]))
    )
  )

}

#Drop helper columns from the flextable
helper_cols <- c("slope_est","slope_SE","slope_p","slope_grp")

if ("delete_columns" %in% getNamespaceExports("flextable")) {
  ft <- delete_columns(ft, j = helper_cols)
} else {
  keep <- c("MI","N","Slopes ± SE","Slopes P-value")
  ft$col_keys <- keep
}

# Minimal styling:
ft <- bold(ft, part = "header")
ft <- fontsize(ft, size = 12, part = "all")
ft <- align(ft, align = "center", part = "all")
ft <- theme_box(ft)
ft <- autofit(ft)

# Export
doc <- read_docx()
doc <- body_add_par(doc, "Nmass", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_dominance_Nmass_NutNet.docx")
##################################################################################################

###############################################################################################################################
############################### Figure lambda relative dominance_Nmass_MI #############################################
logNmass_grid <- seq(min(stat_Nmass_clean$logNmass, na.rm = TRUE),
                 max(stat_Nmass_clean$logNmass, na.rm = TRUE),
                 length.out = 300)

emm_df_Nmass_dominance <- as.data.frame(
  emmeans(
    dominance_Nmass_NutNet,
    ~ logNmass | aridity * Ntrt,
    at = list(
      logNmass = logNmass_grid,
      aridity  = aridity_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_Nmass_dominance

# prediction dataframe : fit + MI factor (3 levels)
pred_df_Nmass_dominance <- emm_df_Nmass_dominance %>%
  mutate(
    lambda_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    MI = factor(
      aridity,
      levels = aridity_levels,
      labels = paste0("MI = ", aridity_levels)))
names(pred_df_Nmass_dominance)

# Center where polygons meet: at the median
x_med <- median(stat_Nmass_clean$logNmass, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_Nmass_dominance <- pred_df_Nmass_dominance %>%
  group_by(Ntrt, MI) %>%
  arrange(logNmass, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logNmass),
    x_max   = last(logNmass),
    fit_min = first(lambda_fit),
    se_min  = first(SE),
    fit_max = last(lambda_fit),
    se_max  = last(SE),
    fit_med = approx(x = logNmass, y = lambda_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_Nmass_dominance <- pmap_dfr(poly_df_Nmass_dominance, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                    fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})

poly_right_Nmass_dominance <- pmap_dfr(poly_df_Nmass_dominance, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                     fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})

poly_plot_Nmass_dominance <- bind_rows(poly_left_Nmass_dominance, poly_right_Nmass_dominance) %>%
  group_by(Ntrt, MI, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_Nmass_dominance

# MI colors by level (low/mid/high) 
MI_cols <- c("#7F0000", "#228B22","#08306B")

names(MI_cols) <- levels(pred_df_Nmass$MI)

# Plot: polygons + 6 lines (3 MI x 2 Ntrt)
dominance_Nmass_MI_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_Nmass_dominance,
               aes(x = x, y = y, group = id, fill = MI),
               alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_Nmass_dominance,
            aes(x = logNmass, y = lambda_fit,
                color = MI, linetype = Ntrt,
                group = interaction(MI, Ntrt)),
            linewidth = 1) +
  
  scale_linetype_manual(
    name   = "Soil N",
    values = c("0" = "dashed", "1" = "solid"),
    labels = c("0" = "Ambient N", "1" = "Added N")) +
  
  scale_color_manual(name = "MI", values = MI_cols) +
  scale_fill_manual (name = "MI", values = MI_cols) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["reldom"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

dominance_Nmass_MI_NutNet


############################### Extract Tg  slopes ################################

################ get slopes and of they are significantly different from zero##### 
tmp_levels_quantile <- quantile(
  stat_Nmass_clean$tmp_GS,
  probs = c(0.1, 0.5, 0.96),
  na.rm = TRUE)
tmp_levels_num <- as.numeric(tmp_levels_quantile)
tmp_levels <- round(tmp_levels_num, 0) 


dominance_slopes_tmp_logNmass_NutNet  <- summary(
  emtrends(dominance_Nmass_NutNet , ~ tmp_GS * Ntrt, 
           var = "logNmass", 
           at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
dominance_slopes_tmp_logNmass_NutNet 

########### comp between slopes ################
dominance_slope_comp_tmp_logNmass_NutNet  <- emtrends(dominance_Nmass_NutNet , ~ tmp_GS * Ntrt, 
                                                   var = "logNmass",
                                                   at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1)))
dominance_slope_comp_tmp_logNmass_NutNet 

dominance_comp_slopes_tmp_logNmass_NutNet  <- cld(dominance_slope_comp_tmp_logNmass_NutNet , Letters = letters, adjust = "tukey")
dominance_comp_slopes_tmp_logNmass_NutNet 


#################### export tables ############################################
##############################################################################

################# Slopes ##############################
dominance_slopes_tmp_logNmass_NutNet_df <- as.data.frame(dominance_slopes_tmp_logNmass_NutNet) %>%
  
  transmute(
    tmp_GS,
    Ntrt,
    slope_est = logNmass.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
dominance_slopes_tmp_logNmass_NutNet_df

dominance_comp_slopes_tmp_logNmass_NutNet_groups <- as.data.frame(dominance_comp_slopes_tmp_logNmass_NutNet) %>%
  transmute(tmp_GS , Ntrt, slope_grp = str_squish(.group))
dominance_comp_slopes_tmp_logNmass_NutNet_groups
# Merge slope stats and letters
dominance_slopes_logNmass_NutNet_all <- left_join(dominance_slopes_tmp_logNmass_NutNet_df, dominance_comp_slopes_tmp_logNmass_NutNet_groups, by = c("tmp_GS","Ntrt"))
##################################################################################################################
### formating ###########################
#Keep helper columns for composing
Table_dominance_logNmass_NutNet_ft_tmp <- dominance_slopes_logNmass_NutNet_all %>%
  transmute(
    Tg = tmp_GS,
    N  = Ntrt,
    slope_est = as.numeric(slope_est),
    slope_SE  = as.numeric(slope_SE),
    slope_p   = as.numeric(slope_p),
    slope_grp = as.character(slope_grp)
  ) %>%
  mutate(
    `Slopes ± SE`     = "",
    `Slopes P-value`        = fmt_p(slope_p))
Table_dominance_logNmass_NutNet_ft_tmp

## check format and columns
Table_dominance_logNmass_NutNet_ft_tmp_out <- Table_dominance_logNmass_NutNet_ft_tmp %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)) %>%
  
  relocate(Tg, N, `Slopes ± SE`, `Slopes P-value`)
Table_dominance_logNmass_NutNet_ft_tmp_out

# build flextable
ft <- flextable(Table_dominance_logNmass_NutNet_ft_tmp_out)

for (i in seq_len(nrow(Table_dominance_logNmass_NutNet_ft_tmp_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_dominance_logNmass_NutNet_ft_tmp_out$slope_est[i])),
      as_chunk(Table_dominance_logNmass_NutNet_ft_tmp_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_dominance_logNmass_NutNet_ft_tmp_out$slope_SE[i]))
    )
  )
}

helper_cols <- c("slope_est","slope_SE","slope_p","slope_grp")

if ("delete_columns" %in% getNamespaceExports("flextable")) {
  ft <- delete_columns(ft, j = helper_cols)
} else {
  keep <- c("Tg","N","Slopes ± SE","Slopes P-value")
  ft$col_keys <- keep
}

# Minimal styling:
ft <- bold(ft, part = "header")
ft <- fontsize(ft, size = 12, part = "all")
ft <- align(ft, align = "center", part = "all")
ft <- theme_box(ft)
ft <- autofit(ft)

# Export
doc <- read_docx()
doc <- body_add_par(doc, "logNmass", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_dominance_tmp_logNmass_NutNet.docx")

##################################################################################################
logNmass_grid <- seq(min(stat_Nmass_clean$logNmass, na.rm = TRUE),
                     max(stat_Nmass_clean$logNmass, na.rm = TRUE),
                     length.out = 300)
logNmass_grid

emm_df_logNmass_tg <- as.data.frame(
  emmeans(
    dominance_Nmass_NutNet,
    ~ logNmass | tmp_GS * Ntrt,
    at = list(
      logNmass = logNmass_grid,
      tmp_GS  = tmp_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_logNmass_tg  

# prediction dataframe : fit + Tg factor (3 levels)
pred_df_logNmass_tg   <- emm_df_logNmass_tg   %>%
  mutate(
    dominance_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    Tg = factor(
      tmp_GS,
      levels = tmp_levels,
      labels = paste0("Tg = ", tmp_levels)))
names(pred_df_logNmass_tg)

# Center where polygons meet: at the median
x_med <- median(stat_Nmass_clean$logNmass, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_logNmass_tg <- pred_df_logNmass_tg %>%
  group_by(Ntrt, Tg) %>%
  arrange(logNmass, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logNmass),
    x_max   = last(logNmass),
    fit_min = first(dominance_fit),
    se_min  = first(SE),
    fit_max = last(dominance_fit),
    se_max  = last(SE),
    fit_med = approx(x = logNmass, y = dominance_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_logNmass_tg <- pmap_dfr(poly_df_logNmass_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                                fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})
poly_left_logNmass_tg

poly_right_logNmass_tg <- pmap_dfr(poly_df_logNmass_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                                 fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})
poly_right_logNmass_tg


poly_plot_logNmass_tg <- bind_rows(poly_left_logNmass_tg, poly_right_logNmass_tg) %>%
  group_by(Ntrt, Tg, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_logNmass_tg


# Tg colors by level (low/mid/high) 
Tg_cols <- c("#7F0000", "#228B22","#08306B")

names(Tg_cols) <- levels(pred_df_logNmass_tg$tmp_GS)

# Plot: polygons + 6 lines (3 Tg x 2 Ntrt)
dominance_logNmass_Tg_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_logNmass_tg,
               aes(x = x, y = y, group = id, fill = Tg),
               alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_logNmass_tg,
            aes(x = logNmass, y = dominance_fit,
                color = Tg, linetype = Ntrt,
                group = interaction(Tg, Ntrt)),
            linewidth = 1) +
  
  scale_linetype_manual(
    name   = "Soil N",
    values = c("0" = "dashed", "1" = "solid"),
    labels = c("0" = "Ambient N", "1" = "Added N")) +
  
  scale_color_manual(name = "Tg", values = Tg_cols) +
  scale_fill_manual (name = "Tg", values = Tg_cols) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["reldom"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

dominance_logNmass_Tg_NutNet

############################################# Heat maps ###############################################
## common code
## make tmp bins
make_tmp_bin <- function(x){
  cut(
    x,
    breaks = c(-Inf, 5, 10, 15, 20, 25),
    include.lowest = TRUE,
    right = FALSE,
    labels = c("Cold (<5)", "Cool (5–10)", "Mild (10–15)", "Warm (15–20)", "Hot (20–25)")
  )
}

## make tmp bins
make_aridity_class <- function(a){
  cut(
    a,
    breaks = c(-Inf, 0.20, 0.50, 0.65, Inf),
    include.lowest = TRUE,
    right = TRUE,
    labels = c("Arid (<0.2)", "Semi-arid (0.20–0.50)", "Dry sub-humid (0.50–0.65)",
               "Humid (>0.65)"))}

stat_Nmass_clean <- stat_Nmass_clean %>%
  mutate(
    tmp_bin = make_tmp_bin(tmp_GS),
    aridity_class = make_aridity_class(aridity)
  )

# representative values per bin/class (midpoints)
tmp_mid <- tibble(
  tmp_bin = factor(c("Cold (<5)","Cool (5–10)","Mild (10–15)","Warm (15–20)","Hot (20–25)"),
                   levels = c("Cold (<5)","Cool (5–10)","Mild (10–15)","Warm (15–20)","Hot (20–25)")),
  tmp_val = c(2.5, 7.5, 12.5, 17.5, 22.5)
)

arid_mid <- tibble(
  aridity_class = factor(c("Arid (<0.20)","Semi-arid (0.20–0.50)","Dry sub-humid (0.50–0.65)", "Humid (>0.65)"),
                         levels = c("Arid (<0.20)","Semi-arid (0.20–0.50)","Dry sub-humid (0.50–0.65)", "Humid (>0.65)")),
  arid_val = c(0.1, 0.35, 0.575, 0.9)
)

# Build grid including P x K to calculate the average 
grid_pk <- tidyr::crossing(
  tmp_mid,
  arid_mid,
  Ntrt = c(0, 1),
  Ptrt = c(0, 1),
  Ktrt = c(0, 1)
)

######################################### Lambda cover ###############################################
# Compute slope of logNmass from fixed effects for lambda cover
# Slope = d/d(logNmass) of linear predictor

b_lambda_cover <- lme4::fixef(lambda_Nmass_NutNet)
getb <- function(nm) if (nm %in% names(b_lambda_cover)) unname(b_lambda_cover[nm]) else 0

# helper:
grid_pk <- grid_pk %>%
  mutate(
    slope_logNmass_cover =
      getb("logNmass") +
      getb("logNmass:aridity") * arid_val +
      getb("logNmass:tmp_GS") * tmp_val +
      getb("logNmass:Ntrt") * Ntrt +
      getb("logNmass:Ptrt") * Ptrt +
      getb("logNmass:Ktrt") * Ktrt +
      getb("logNmass:aridity:Ntrt") * arid_val * Ntrt +
      getb("logNmass:Ntrt:tmp_GS") * tmp_val * Ntrt +
      getb("logNmass:Ntrt:Ptrt") * Ntrt * Ptrt +
      getb("logNmass:Ntrt:Ktrt") * Ntrt * Ktrt +
      getb("logNmass:Ptrt:Ktrt") * Ptrt * Ktrt +
      getb("logNmass:Ntrt:Ptrt:Ktrt") * Ntrt * Ptrt * Ktrt
  )

# Average across plots within Ntrt groups 
#Ntrt=0: average across Control, P, K, PK (Ptrt/Ktrt combos)
#Ntrt=1: average across N, NP, NK, NPK (Ptrt/Ktrt combos)

grid_avg <- grid_pk %>%
  group_by(tmp_bin, tmp_val, aridity_class, arid_val, Ntrt) %>%
  summarise(slope_logNmass_cover = mean(slope_logNmass_cover), .groups = "drop") %>%
  mutate(
    Ntrt_lab = ifelse(Ntrt == 0, "No N (Control, P, K, PK)", "+N (N, NP, NK, NPK)")
  )

plot_heat_lambda <- function(df_oneN, title_text,
                             shared_scale = NULL,
                             digits = 2) {
  
  p <- ggplot(df_oneN,
              aes(x = aridity_class,
                  y = tmp_bin,
                  fill = slope_logNmass_cover)) +
    
    # Tiles
    geom_tile(color = "white", linewidth = 0.5) +
    
    # Text labels — readable everywhere
    geom_label(
      aes(label = sprintf(paste0("%.", digits, "f"),
                          slope_logNmass_cover)),
      color = "black",
      fill = scales::alpha("white", 0.70),
      label.size = 0,
      label.padding = unit(0.08, "lines"),
      size = 3.0,
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    
    labs(
      x = "MI",
      y = "Tg",
      title = title_text
    ) +
    
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      
      axis.title = element_text(size = 9),
      axis.text.x = element_text(
        size = 8,
        angle = 25,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = 8),
      
      plot.title = element_text(
        face = "bold",
        size = 10,
        hjust = 0.5
      ),
      
      plot.margin = margin(3, 3, 3, 3),
      
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      
      # Better legend proportions
      legend.key.height = unit(18, "mm"),
      legend.key.width  = unit(5, "mm")
    )
  
  if (!is.null(shared_scale)) {
    p <- p + shared_scale
  } else {
    p <- p + scale_fill_viridis_c(
      option = "cividis",
      direction = -1,
      begin = 0.18,
      end = 0.98,
      limits = c(-0.15, 0.40),
      oob = scales::squish,
      name = expression("Slopes of " * N[mass]),
      breaks = seq(-0.1, 0.4, by = 0.1),
      labels = scales::number_format(accuracy = 0.01)
    )
  }
  
  p
}

# Shared legend scale (fixed)
scale_fill_shared <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  begin = 0.18,
  end = 0.98,
  limits = c(-0.15, 0.40),
  oob = scales::squish,
  name = expression("Slopes of " * N[mass]),
  breaks = seq(-0.1, 0.4, by = 0.1),
  labels = scales::number_format(accuracy = 0.01)
)

# Plot Ambient N
lambda_Nmass_N0 <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 0),
  expression(lambda[cover] ~ "\u2013 Ambient N"),
  shared_scale = scale_fill_shared,
  digits = 2
) +
  theme(legend.position = "right")

lambda_Nmass_N0

lambda_Nmass_N  <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 1),
  expression(lambda[cover] ~ "\u2013 Added N"),
  shared_scale = scale_fill_shared,
  digits = 2
) + theme(legend.position = "right")
lambda_Nmass_N

#########################################################################################

######################################### Relative dominance ###############################################

b_lambda_reldom <- lme4::fixef(dominance_Nmass_NutNet)
getb <- function(nm) if (nm %in% names(b_lambda_reldom)) unname(b_lambda_reldom[nm]) else 0

# helper:
grid_pk <- grid_pk %>%
  mutate(
    slope_logNmass_reldom =
      getb("logNmass") +
      getb("logNmass:aridity") * arid_val +
      getb("logNmass:tmp_GS") * tmp_val +
      getb("logNmass:Ntrt") * Ntrt +
      getb("logNmass:Ptrt") * Ptrt +
      getb("logNmass:Ktrt") * Ktrt +
      getb("logNmass:aridity:Ntrt") * arid_val * Ntrt +
      getb("logNmass:Ntrt:tmp_GS") * tmp_val * Ntrt +
      getb("logNmass:Ntrt:Ptrt") * Ntrt * Ptrt +
      getb("logNmass:Ntrt:Ktrt") * Ntrt * Ktrt +
      getb("logNmass:Ptrt:Ktrt") * Ptrt * Ktrt +
      getb("logNmass:Ntrt:Ptrt:Ktrt") * Ntrt * Ptrt * Ktrt
  )

# Average across plots within Ntrt groups 
#Ntrt=0: average across Control, P, K, PK (Ptrt/Ktrt combos)
#Ntrt=1: average across N, NP, NK, NPK (Ptrt/Ktrt combos)

grid_avg <- grid_pk %>%
  group_by(tmp_bin, tmp_val, aridity_class, arid_val, Ntrt) %>%
  summarise(slope_logNmass_reldom = mean(slope_logNmass_reldom), .groups = "drop") %>%
  mutate(
    Ntrt_lab = ifelse(Ntrt == 0, "No N (Control, P, K, PK)", "+N (N, NP, NK, NPK)")
  )

plot_heat_lambda <- function(df_oneN, title_text,
                             shared_scale = NULL,
                             digits = 2) {
  
  p <- ggplot(df_oneN,
              aes(x = aridity_class,
                  y = tmp_bin,
                  fill = slope_logNmass_reldom)) +
    
    # Tiles
    geom_tile(color = "white", linewidth = 0.5) +
    
    # Text labels — readable everywhere
    geom_label(
      aes(label = sprintf(paste0("%.", digits, "f"),
                          slope_logNmass_reldom)),
      color = "black",
      fill = scales::alpha("white", 0.70),
      label.size = 0,
      label.padding = unit(0.08, "lines"),
      size = 3.0,
      lineheight = 0.95,
      show.legend = FALSE
    ) +
    
    labs(
      x = "MI",
      y = "Tg",
      title = title_text
    ) +
    
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      
      axis.title = element_text(size = 9),
      axis.text.x = element_text(
        size = 8,
        angle = 25,
        hjust = 1,
        vjust = 1
      ),
      axis.text.y = element_text(size = 8),
      
      plot.title = element_text(
        face = "bold",
        size = 10,
        hjust = 0.5
      ),
      
      plot.margin = margin(3, 3, 3, 3),
      
      legend.title = element_text(size = 9),
      legend.text  = element_text(size = 8),
      
      # Better legend proportions
      legend.key.height = unit(18, "mm"),
      legend.key.width  = unit(5, "mm")
    )
  
  if (!is.null(shared_scale)) {
    p <- p + shared_scale
  } else {
    p <- p + scale_fill_viridis_c(
      option = "cividis",
      direction = -1,
      begin = 0.18,
      end = 0.98,
      limits = c(-0.15, 0.40),
      oob = scales::squish,
      name = expression("Slopes of " * N[mass]),
      breaks = seq(-0.1, 0.4, by = 0.1),
      labels = scales::number_format(accuracy = 0.01)
    )
  }
  
  p
}

# Shared legend scale (fixed)
scale_fill_shared <- scale_fill_viridis_c(
  option = "cividis",
  direction = -1,
  begin = 0.18,
  end = 0.98,
  limits = c(-0.15, 0.40),
  oob = scales::squish,
  name = expression("Slopes of " * N[mass]),
  breaks = seq(-0.1, 0.4, by = 0.1),
  labels = scales::number_format(accuracy = 0.01)
)

# Plot Ambient N
reldom_Nmass_N0 <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 0),
  expression(lambda[cover] ~ "\u2013 Ambient N"),
  shared_scale = scale_fill_shared,
  digits = 2
) +
  theme(legend.position = "right")

reldom_Nmass_N0

reldom_Nmass_N  <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 1),
  expression(lambda[cover] ~ "\u2013 Added N"),
  shared_scale = scale_fill_shared,
  digits = 2
) + theme(legend.position = "right")
reldom_Nmass_N

###############################################################################################

##################### Nmass-lambda cover raw data #################################################
names(stat_Nmass_clean)

MI_limits <- range(stat_Nmass_clean$aridity, na.rm = TRUE)

lambda_MI_Nmass_rawdata <- ggplot(
  stat_Nmass_clean,
  aes(x = logNmass, y = mean_log_lambda_cover_all)
) +
  geom_point(
    aes(color = aridity, shape = factor(Ntrt)),
    size = 2.4, alpha = 0.25
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    labels = c("0" = "Ambient N", "1" = "Added N"),
    name = NULL
  ) +
  scale_color_gradientn(
    colours = c("#7F0000", "#228B22", "#08306B"),
    limits  = MI_limits,
    name    = "MI"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["cover"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

lambda_MI_Nmass_rawdata


##################### Nmass-lambda reldom raw data #################################################
names(stat_Nmass_clean)

MI_limits <- range(stat_Nmass_clean$aridity, na.rm = TRUE)

dominance_MI_Nmass_rawdata <- ggplot(
  stat_Nmass_clean,
  aes(x = logNmass, y = mean_log_lambda_Ci_all)
) +
  geom_point(
    aes(color = aridity, shape = factor(Ntrt)),
    size = 2.4, alpha = 0.25
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    labels = c("0" = "Ambient N", "1" = "Added N"),
    name = NULL
  ) +
  scale_color_gradientn(
    colours = c("#7F0000", "#228B22", "#08306B"),
    limits  = MI_limits,
    name    = "MI"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["reldom"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

dominance_MI_Nmass_rawdata


##################### Nmass-lambda cover raw data Tg #################################################

Tg_limits <- range(stat_Nmass_clean$tmp_GS, na.rm = TRUE)

lambda_Nmass_rawdata_Tg <- ggplot(
  stat_Nmass_clean,
  aes(x = logNmass, y = mean_log_lambda_cover_all)
) +
  geom_point(
    aes(color = tmp_GS, shape = factor(Ntrt)),
    size = 2.4, alpha = 0.25
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    labels = c("0" = "Ambient N", "1" = "Added N"),
    name = NULL
  ) +
  scale_color_gradientn(
    colours = c("#7F0000", "#228B22", "#08306B"),
    limits  = Tg_limits,
    name    = "Tg"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["cover"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

lambda_Nmass_rawdata_Tg

##################### chi-lambda reldom Tg raw data #################################################

Tg_limits <- range(stat_Nmass_clean$tmp_GS, na.rm = TRUE)

reldom_Nmass_rawdata_Tg <- ggplot(
  stat_Nmass_clean,
  aes(x = logNmass, y = mean_log_lambda_Ci_all)
) +
  geom_point(
    aes(color = tmp_GS, shape = factor(Ntrt)),
    size = 2.4, alpha = 0.25
  ) +
  scale_shape_manual(
    values = c("0" = 16, "1" = 17),
    labels = c("0" = "Ambient N", "1" = "Added N"),
    name = NULL
  ) +
  scale_color_gradientn(
    colours = c("#7F0000", "#228B22", "#08306B"),
    limits  = Tg_limits,
    name    = "Tg"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = alpha("white", 0.75), colour = NA),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
    axis.title = element_text(colour = "black"),
    axis.text  = element_text(colour = "black"),
    axis.title.y = element_text(size = rel(2.2), colour = "black"),
    axis.title.x = element_text(size = rel(2.2), colour = "black"),
    axis.text.x  = element_text(size = rel(2.0), colour = "black"),
    axis.text.y  = element_text(size = rel(2.0), colour = "black")
  ) +
  ylab(expression(lambda["reldom"])) +
  xlab(expression(log~italic(N)[mass]~"(gN g"^-1*")"))

reldom_Nmass_rawdata_Tg

