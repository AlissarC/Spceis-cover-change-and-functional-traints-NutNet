#########################################
# necessary packages 
#########################################

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

###################################### SLA ##########################################################
stat_SLA_clean = stat_Nmass_clean

################### Natural log transforming of traits ##########################
plot(stat_SLA_clean$aridity, stat_SLA_clean$tmp_GS)
stat_SLA_clean$logSLA = log(stat_SLA_clean$SLA_NutNet)
hist(stat_SLA_clean$logSLA)
################################################################### LMERs for SLA ############################################
######################################################## SLA NutNet analysis: log lambda corrected by K ##############################################
names(stat_SLA_clean)

lambda_SLA_NutNet <- lmer(mean_log_lambda_cover_all~logSLA + logSLA:(aridity*Ntrt+ tmp_GS*Ntrt) + aridity + tmp_GS + 
                            logSLA:(Ntrt*Ptrt*Ktrt) + n_years +
                            (1|site_code) + (1|site_code:block),
                          data=stat_SLA_clean)
Anova(lambda_SLA_NutNet) 
vif(lambda_SLA_NutNet)

confint(lambda_SLA_NutNet, method = "Wald")
summary(lambda_SLA_NutNet)
r.squaredGLMM(lambda_SLA_NutNet)              
#R2m       R2c
#0.004857121 0.100169
#emm_options(pbkrtest.limit = 51099)
AIC(lambda_SLA_NutNet) #   132012.8

res <- residuals(lambda_SLA_NutNet)
qqnorm(res)
qqline(res, col = "red")

hist(res, breaks = 20, main = "Histogram of Residuals") ## good

plot(fitted(lambda_SLA_NutNet), res, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity :OK

densityPlot(res)

ks.test(res, "pnorm", mean = mean(res), sd = sd(res))

######################### Export Anova table ######################################

# Define the p-value formatting function
fmt_p <- function(p) {
  ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", ifelse(p < 0.01, "<0.01", ifelse(p < 0.05, "<0.05",sprintf("%.3f", p)))))
}

lambda_SLA_table = as.data.frame(Anova(lambda_SLA_NutNet))
lambda_SLA_model_NutNet <- data.frame(
  Var   = c(
    'SLA', 'MI', 'Tg', 'number of years','SLA x MI', 
    'SLA x Soil N', 'SLA x Tg','SLA x Soil P',
    'SLA x Soil K', 'SLA x MI x Soil N', 'SLA x Tg x Soil N', 'SLA x Soil N x Soil P',
    'SLA x Soil N x Soil K','SLA x Soil P x Soil K', 'SLA x Soil N x Soil P x Soil K '),
  chisq = round(lambda_SLA_table[1:15, 1], 2),
  df    = lambda_SLA_table[1:15, 2],
  P     = fmt_p(lambda_SLA_table[1:15, 3])   
)

# Create flextable
ft <- flextable(lambda_SLA_model_NutNet) %>%
  autofit()

# Export to Word
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "./Tables/lambda_SLA_model_NutNet.docx")

################################ Extract aridity slopes and intercepts ################################
################ get slopes and of they are significantly different from zero##### 
#  MI levels from the observed distribution (quantiles): 10%, 50%, 95%
aridity_levels_quantile <- quantile(
  stat_SLA_clean$aridity,
  probs = c(0.1, 0.5, 0.96),
  na.rm = TRUE)

aridity_levels_num <- as.numeric(aridity_levels_quantile)
aridity_levels <- round(aridity_levels_num, 1) 

lambda_slopes_MI_SLA_NutNet  <- summary(
  emtrends(lambda_SLA_NutNet , ~ aridity * Ntrt, 
           var = "logSLA", 
           at = list(aridity = aridity_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
lambda_slopes_MI_SLA_NutNet 

########### comp between slopes ################
lambda_slope_comp_MI_SLA_NutNet  <- emtrends(lambda_SLA_NutNet , ~ aridity * Ntrt, 
                                               var = "logSLA",
                                               at = list(aridity = aridity_levels, Ntrt = c(0, 1)))
lambda_slope_comp_MI_SLA_NutNet 

lambda_comp_slopes_MI_SLA_NutNet  <- cld(lambda_slope_comp_MI_SLA_NutNet , Letters = letters, adjust = "tukey")
lambda_comp_slopes_MI_SLA_NutNet 


slopes_lambda_Ntrt_SLA_NutNet  <- summary(
  emtrends(lambda_SLA_NutNet , ~ Ntrt, 
           var = "logSLA", infer = c(TRUE, TRUE)))
slopes_lambda_Ntrt_SLA_NutNet 


slopes_lambda_Ptrt_SLA_NutNet  <- summary(
  emtrends(lambda_SLA_NutNet , ~ Ptrt, 
           var = "logSLA", infer = c(TRUE, TRUE)))
slopes_lambda_Ptrt_SLA_NutNet 

#################### export tables ############################################
##############################################################################

################# Slopes ##############################
lambda_slopes_MI_SLA_NutNet_df <- as.data.frame(lambda_slopes_MI_SLA_NutNet) %>%
  
  transmute(
    aridity,
    Ntrt,
    slope_est = logSLA.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
lambda_slopes_MI_SLA_NutNet_df

lambda_comp_slopes_MI_SLA_NutNet_groups <- as.data.frame(lambda_comp_slopes_MI_SLA_NutNet) %>%
  transmute(aridity, Ntrt, slope_grp = str_squish(.group))
lambda_comp_slopes_MI_SLA_NutNet_groups
# Merge slope stats and letters
lambda_slopes_SLA_NutNet_all <- left_join(lambda_slopes_MI_SLA_NutNet_df, lambda_comp_slopes_MI_SLA_NutNet_groups, by = c("aridity","Ntrt"))
##################################################################################################################

##### formating ###########################
#Keep helper columns for composing
Table_lambda_SLA_NutNet_ft <- lambda_slopes_SLA_NutNet_all %>%
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
Table_lambda_SLA_NutNet_ft


## check format and columns
Table_lambda_SLA_NutNet_out <- Table_lambda_SLA_NutNet_ft %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)

  ) %>%
  
relocate(MI, N, `Slopes ± SE`, `Slopes P-value`)
Table_lambda_SLA_NutNet_out

ft <- flextable(Table_lambda_SLA_NutNet_out)

for (i in seq_len(nrow(Table_lambda_SLA_NutNet_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_lambda_SLA_NutNet_out$slope_est[i])),
      as_chunk(Table_lambda_SLA_NutNet_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_lambda_SLA_NutNet_out$slope_SE[i]))))}

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
doc <- body_add_par(doc, "SLA", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_lambda_SLA_NutNet.docx")

###############################################################################################################################
############################### Figure lambda cover_SLA_MI #############################################
SLA_grid <- seq(min(stat_SLA_clean$logSLA, na.rm = TRUE),
  max(stat_SLA_clean$logSLA, na.rm = TRUE),
  length.out = 300)

emm_df_SLA <- as.data.frame(
  emmeans(
    lambda_SLA_NutNet,
    ~ logSLA | aridity * Ntrt,
    at = list(
      logSLA = SLA_grid,
      aridity  = aridity_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_SLA

pred_df_SLA <- emm_df_SLA %>%
  mutate(
    lambda_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    MI = factor(
      aridity,
      levels = aridity_levels,
      labels = paste0("MI = ", aridity_levels)))
names(pred_df_SLA)

# Center where polygons meet: at the median
x_med <- median(stat_SLA_clean$logSLA, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_SLA <- pred_df_SLA %>%
  group_by(Ntrt, MI) %>%
  arrange(logSLA, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logSLA),
    x_max   = last(logSLA),
    fit_min = first(lambda_fit),
    se_min  = first(SE),
    fit_max = last(lambda_fit),
    se_max  = last(SE),
    fit_med = approx(x = logSLA, y = lambda_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_SLA <- pmap_dfr(poly_df_SLA, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                        fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})

poly_right_SLA <- pmap_dfr(poly_df_SLA, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                         fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})

poly_plot_SLA <- bind_rows(poly_left_SLA, poly_right_SLA) %>%
  group_by(Ntrt, MI, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_SLA
dplyr::count(poly_plot_SLA, id)
# MI colors by level (low/mid/high) 
MI_cols <- c("#7F0000", "#228B22","#08306B")

names(MI_cols) <- levels(pred_df_SLA$MI)

# Plot: polygons + 6 lines (3 MI x 2 Ntrt)
lambda_SLA_MI_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_SLA,
    aes(x = x, y = y, group = id, fill = MI),
    alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_SLA,
            aes(x = logSLA, y = lambda_fit,
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))

lambda_SLA_MI_NutNet

############################### Extract Tg  slopes ################################

################ get slopes and of they are significantly different from zero##### 
tmp_levels_quantile <- quantile(
  stat_SLA_clean$tmp_GS,
  probs = c(0.1, 0.5, 0.96),
  na.rm = TRUE)
tmp_levels_num <- as.numeric(tmp_levels_quantile)
tmp_levels <- round(tmp_levels_num, 0) 


lambda_slopes_tmp_logSLA_NutNet  <- summary(
  emtrends(lambda_SLA_NutNet , ~ tmp_GS * Ntrt, 
           var = "logSLA", 
           at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
lambda_slopes_tmp_logSLA_NutNet 

########### comp between slopes ################
lambda_slope_comp_tmp_logSLA_NutNet  <- emtrends(lambda_SLA_NutNet , ~ tmp_GS * Ntrt, 
                                              var = "logSLA",
                                              at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1)))
lambda_slope_comp_tmp_logSLA_NutNet 

lambda_comp_slopes_tmp_logSLA_NutNet  <- cld(lambda_slope_comp_tmp_logSLA_NutNet , Letters = letters, adjust = "tukey")
lambda_comp_slopes_tmp_logSLA_NutNet 


#################### export tables ############################################
##############################################################################

################# Slopes ##############################
lambda_slopes_tmp_logSLA_NutNet_df <- as.data.frame(lambda_slopes_tmp_logSLA_NutNet) %>%
  
  transmute(
    tmp_GS,
    Ntrt,
    slope_est = logSLA.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
lambda_slopes_tmp_logSLA_NutNet_df

lambda_comp_slopes_tmp_logSLA_NutNet_groups <- as.data.frame(lambda_comp_slopes_tmp_logSLA_NutNet) %>%
  transmute(tmp_GS , Ntrt, slope_grp = str_squish(.group))
lambda_comp_slopes_tmp_logSLA_NutNet_groups
# Merge slope stats and letters
lambda_slopes_logSLA_NutNet_all <- left_join(lambda_slopes_tmp_logSLA_NutNet_df, lambda_comp_slopes_tmp_logSLA_NutNet_groups, by = c("tmp_GS","Ntrt"))
##################################################################################################################
### formating ###########################
#Keep helper columns for composing
Table_lambda_logSLA_NutNet_ft_tmp <- lambda_slopes_logSLA_NutNet_all %>%
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
Table_lambda_logSLA_NutNet_ft_tmp

## check format and columns
Table_lambda_logSLA_NutNet_ft_tmp_out <- Table_lambda_logSLA_NutNet_ft_tmp %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)) %>%
  
  relocate(Tg, N, `Slopes ± SE`, `Slopes P-value`)
Table_lambda_logSLA_NutNet_ft_tmp_out

# build flextable
ft <- flextable(Table_lambda_logSLA_NutNet_ft_tmp_out)

for (i in seq_len(nrow(Table_lambda_logSLA_NutNet_ft_tmp_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_lambda_logSLA_NutNet_ft_tmp_out$slope_est[i])),
      as_chunk(Table_lambda_logSLA_NutNet_ft_tmp_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_lambda_logSLA_NutNet_ft_tmp_out$slope_SE[i]))
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
doc <- body_add_par(doc, "logSLA", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_lambda_tmp_logSLA_NutNet.docx")

##################################################################################################
logSLA_grid <- seq(min(stat_SLA_clean$logSLA, na.rm = TRUE),
                max(stat_SLA_clean$logSLA, na.rm = TRUE),
                length.out = 300)
logSLA_grid

emm_df_logSLA_tg <- as.data.frame(
  emmeans(
    lambda_SLA_NutNet,
    ~ logSLA | tmp_GS * Ntrt,
    at = list(
      logSLA = logSLA_grid,
      tmp_GS  = tmp_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_logSLA_tg  

# prediction dataframe : fit + Tg factor (3 levels)
pred_df_logSLA_tg   <- emm_df_logSLA_tg   %>%
  mutate(
    lambda_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    Tg = factor(
      tmp_GS,
      levels = tmp_levels,
      labels = paste0("Tg = ", tmp_levels)))
names(pred_df_logSLA_tg)

# Center where polygons meet: at the median
x_med <- median(stat_SLA_clean$logSLA, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_logSLA_tg <- pred_df_logSLA_tg %>%
  group_by(Ntrt, Tg) %>%
  arrange(logSLA, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logSLA),
    x_max   = last(logSLA),
    fit_min = first(lambda_fit),
    se_min  = first(SE),
    fit_max = last(lambda_fit),
    se_max  = last(SE),
    fit_med = approx(x = logSLA, y = lambda_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_logSLA_tg <- pmap_dfr(poly_df_logSLA_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                      fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})
poly_left_logSLA_tg

poly_right_logSLA_tg <- pmap_dfr(poly_df_logSLA_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                       fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})
poly_right_logSLA_tg


poly_plot_logSLA_tg <- bind_rows(poly_left_logSLA_tg, poly_right_logSLA_tg) %>%
  group_by(Ntrt, Tg, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_logSLA_tg


# Tg colors by level (low/mid/high) 
Tg_cols <- c("#7F0000", "#228B22","#08306B")


names(Tg_cols) <- levels(pred_df_logSLA_tg$tmp_GS)

# Plot: polygons + 6 lines (3 Tg x 2 Ntrt)
lambda_logSLA_Tg_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_logSLA_tg,
               aes(x = x, y = y, group = id, fill = Tg),
               alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_logSLA_tg,
            aes(x = logSLA, y = lambda_fit,
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))

lambda_logSLA_Tg_NutNet

###################################################################################

######################################################## SLA NutNet analysis: log lambda relative dominance corrected by K ##############################################
names(stat_SLA_clean)

dominance_SLA_NutNet <- lmer(mean_log_lambda_Ci_all~logSLA + logSLA:(aridity*Ntrt+ tmp_GS*Ntrt) + aridity + tmp_GS + 
                              logSLA:(Ntrt*Ptrt*Ktrt) + n_years +
                              (1|site_code) + (1|site_code:block),
                            data=stat_SLA_clean)
Anova(dominance_SLA_NutNet)

summary(dominance_SLA_NutNet)
r.squaredGLMM(dominance_SLA_NutNet)              
#R2m       R2c
#0.004411358 0.1262298
AIC(dominance_SLA_NutNet) # 72705.53

res <- residuals(dominance_SLA_NutNet)
qqnorm(res)
qqline(res, col = "red")

hist(res, breaks = 20, main = "Histogram of Residuals") ## good

plot(fitted(dominance_SLA_NutNet), res, xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs. Fitted Values")  # heteroscedasticity :OK

densityPlot(res)

ks.test(res, "pnorm", mean = mean(res), sd = sd(res))

######################### Export Anova table #######################################

# Define the p-value formatting function
fmt_p <- function(p) {
  ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", ifelse(p < 0.01, "<0.01", ifelse(p < 0.05, "<0.05",sprintf("%.3f", p)))))
}
dominance_SLA_NutNet_table = as.data.frame(Anova(dominance_SLA_NutNet))

dominance_SLA_model_NutNet <- data.frame(
  Var   = c(
    'SLA', 'MI', 'Tg', 'Number of years', 'SLA x MI', 
    'SLA x Soil N', 'SLA x Tg','SLA x Soil P',
    'SLA x Soil K', 'SLA x MI x Soil N', 'SLA x Tg x Soil N', 'SLA x Soil N x Soil P',
    'SLA x Soil N x Soil K','SLA x Soil P x Soil K', 'SLA x Soil N x Soil P x Soil K '),
  chisq = round(dominance_SLA_NutNet_table[1:15, 1], 2),
  df    = dominance_SLA_NutNet_table[1:15, 2],
  P     = fmt_p(dominance_SLA_NutNet_table[1:15, 3])   
)

# Create flextable
ft <- flextable(dominance_SLA_model_NutNet) %>%
  autofit()

# Export to Word
doc <- read_docx()
doc <- body_add_flextable(doc, value = ft)
print(doc, target = "./Tables/dominance_SLA_model_NutNet.docx")

################################ Extract aridity slopes and intercepts ################################

################ get slopes and of they are significantly different from zero##### 
dominance_slopes_MI_SLA_NutNet  <- summary(
  emtrends(dominance_SLA_NutNet , ~ aridity * Ntrt, 
           var = "logSLA", 
           at = list(aridity = aridity_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
dominance_slopes_MI_SLA_NutNet 

########### comp between slopes ################
dominance_slope_comp_MI_SLA_NutNet  <- emtrends(dominance_SLA_NutNet , ~ aridity * Ntrt, 
                                               var = "logSLA",
                                               at = list(aridity = aridity_levels, Ntrt = c(0, 1)))
dominance_slope_comp_MI_SLA_NutNet 

dominance_comp_slopes_MI_SLA_NutNet  <- cld(dominance_slope_comp_MI_SLA_NutNet , Letters = letters, adjust = "tukey")
dominance_comp_slopes_MI_SLA_NutNet 


slopes_dominance_Ntrt_SLA_NutNet  <- summary(
  emtrends(dominance_SLA_NutNet , ~ Ntrt, 
           var = "logSLA", infer = c(TRUE, TRUE)))
slopes_dominance_Ntrt_SLA_NutNet 


slopes_dominance_Ptrt_SLA_NutNet  <- summary(
  emtrends(dominance_SLA_NutNet , ~ Ptrt, 
           var = "logSLA", infer = c(TRUE, TRUE)))
slopes_dominance_Ptrt_SLA_NutNet 

#################### export tables ############################################
##############################################################################
dominance_slopes_MI_SLA_NutNet_df <- as.data.frame(dominance_slopes_MI_SLA_NutNet) %>%
  
  transmute(
    aridity,
    Ntrt,
    slope_est = logSLA.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
dominance_slopes_MI_SLA_NutNet_df

dominance_comp_slopes_MI_SLA_NutNet_groups <- as.data.frame(dominance_comp_slopes_MI_SLA_NutNet) %>%
  transmute(aridity, Ntrt, slope_grp = str_squish(.group))
dominance_comp_slopes_MI_SLA_NutNet_groups
# Merge slope stats and letters
dominance_slopes_SLA_NutNet_all <- left_join(dominance_slopes_MI_SLA_NutNet_df, dominance_comp_slopes_MI_SLA_NutNet_groups, by = c("aridity","Ntrt"))
##################################################################################################################
### formating ###########################
#Keep helper columns for composing
Table_dominance_SLA_NutNet_ft <- dominance_slopes_SLA_NutNet_all %>%
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
Table_dominance_SLA_NutNet_ft

fmt_p <- function(p) ifelse(is.na(p), NA, ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))

## check format and columns
Table_dominance_SLA_NutNet_out <- Table_dominance_SLA_NutNet_ft %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)
   
  ) %>%
  
relocate(MI, N, `Slopes ± SE`, `Slopes P-value`)
Table_dominance_SLA_NutNet_out

# build flextable
ft <- flextable(Table_dominance_SLA_NutNet_out)

for (i in seq_len(nrow(Table_dominance_SLA_NutNet_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_dominance_SLA_NutNet_out$slope_est[i])),
      as_chunk(Table_dominance_SLA_NutNet_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_dominance_SLA_NutNet_out$slope_SE[i]))
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
doc <- body_add_par(doc, "SLA", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_dominance_SLA_NutNet.docx")
##################################################################################################

###############################################################################################################################
############################### Figure lambda relative dominance_SLA_MI #############################################
logSLA_grid <- seq(min(stat_SLA_clean$logSLA, na.rm = TRUE),
                 max(stat_SLA_clean$logSLA, na.rm = TRUE),
                 length.out = 300)

emm_df_SLA_dominance <- as.data.frame(
  emmeans(
    dominance_SLA_NutNet,
    ~ logSLA | aridity * Ntrt,
    at = list(
      logSLA = logSLA_grid,
      aridity  = aridity_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_SLA_dominance

# prediction dataframe : fit + MI factor (3 levels)
pred_df_SLA_dominance <- emm_df_SLA_dominance %>%
  mutate(
    lambda_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    MI = factor(
      aridity,
      levels = aridity_levels,
      labels = paste0("MI = ", aridity_levels)))
names(pred_df_SLA_dominance)

# Center where polygons meet: at the median
x_med <- median(stat_SLA_clean$logSLA, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_SLA_dominance <- pred_df_SLA_dominance %>%
  group_by(Ntrt, MI) %>%
  arrange(logSLA, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logSLA),
    x_max   = last(logSLA),
    fit_min = first(lambda_fit),
    se_min  = first(SE),
    fit_max = last(lambda_fit),
    se_max  = last(SE),
    fit_med = approx(x = logSLA, y = lambda_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_SLA_dominance <- pmap_dfr(poly_df_SLA_dominance, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                    fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})

poly_right_SLA_dominance <- pmap_dfr(poly_df_SLA_dominance, function(Ntrt, MI, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                     fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    MI   = MI,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})

poly_plot_SLA_dominance <- bind_rows(poly_left_SLA_dominance, poly_right_SLA_dominance) %>%
  group_by(Ntrt, MI, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_SLA_dominance

# MI colors by level (low/mid/high) 
MI_cols <- c("#7F0000", "#228B22","#08306B")

names(MI_cols) <- levels(pred_df_SLA$MI)

# Plot: polygons + 6 lines (3 MI x 2 Ntrt)
dominance_SLA_MI_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_SLA_dominance,
               aes(x = x, y = y, group = id, fill = MI),
               alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_SLA_dominance,
            aes(x = logSLA, y = lambda_fit,
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))

dominance_SLA_MI_NutNet

############################### Extract Tg  slopes ################################

################ get slopes and of they are significantly different from zero##### 
tmp_levels_quantile <- quantile(
  stat_SLA_clean$tmp_GS,
  probs = c(0.1, 0.5, 0.96),
  na.rm = TRUE)
tmp_levels_num <- as.numeric(tmp_levels_quantile)
tmp_levels <- round(tmp_levels_num, 0) 


dominance_slopes_tmp_logSLA_NutNet  <- summary(
  emtrends(dominance_SLA_NutNet , ~ tmp_GS * Ntrt, 
           var = "logSLA", 
           at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1))), infer = c(TRUE, TRUE))
dominance_slopes_tmp_logSLA_NutNet 

########### comp between slopes ################
dominance_slope_comp_tmp_logSLA_NutNet  <- emtrends(dominance_SLA_NutNet , ~ tmp_GS * Ntrt, 
                                                   var = "logSLA",
                                                   at = list(tmp_GS = tmp_levels, Ntrt = c(0, 1)))
dominance_slope_comp_tmp_logSLA_NutNet 

dominance_comp_slopes_tmp_logSLA_NutNet  <- cld(dominance_slope_comp_tmp_logSLA_NutNet , Letters = letters, adjust = "tukey")
dominance_comp_slopes_tmp_logSLA_NutNet 


#################### export tables ############################################
##############################################################################

################# Slopes ##############################
dominance_slopes_tmp_logSLA_NutNet_df <- as.data.frame(dominance_slopes_tmp_logSLA_NutNet) %>%
  
  transmute(
    tmp_GS,
    Ntrt,
    slope_est = logSLA.trend,
    slope_SE  = SE,
    slope_df  = df,
    slope_t   = z.ratio,
    slope_p   = p.value
  )
dominance_slopes_tmp_logSLA_NutNet_df

dominance_comp_slopes_tmp_logSLA_NutNet_groups <- as.data.frame(dominance_comp_slopes_tmp_logSLA_NutNet) %>%
  transmute(tmp_GS , Ntrt, slope_grp = str_squish(.group))
dominance_comp_slopes_tmp_logSLA_NutNet_groups
# Merge slope stats and letters
dominance_slopes_logSLA_NutNet_all <- left_join(dominance_slopes_tmp_logSLA_NutNet_df, dominance_comp_slopes_tmp_logSLA_NutNet_groups, by = c("tmp_GS","Ntrt"))
##################################################################################################################
### formating ###########################
#Keep helper columns for composing
Table_dominance_logSLA_NutNet_ft_tmp <- dominance_slopes_logSLA_NutNet_all %>%
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
Table_dominance_logSLA_NutNet_ft_tmp

## check format and columns
Table_dominance_logSLA_NutNet_ft_tmp_out <- Table_dominance_logSLA_NutNet_ft_tmp %>%
  mutate(
    `Slopes ± SE`     = "",                          # to be filled
    `Slopes P-value`        = fmt_p(slope_p)) %>%
  
  relocate(Tg, N, `Slopes ± SE`, `Slopes P-value`)
Table_dominance_logSLA_NutNet_ft_tmp_out

# build flextable
ft <- flextable(Table_dominance_logSLA_NutNet_ft_tmp_out)

for (i in seq_len(nrow(Table_dominance_logSLA_NutNet_ft_tmp_out))) {
  # slope cell
  ft <- compose(
    ft, i = i, j = "Slopes ± SE",
    value = as_paragraph(
      as_chunk(sprintf("%.3f", Table_dominance_logSLA_NutNet_ft_tmp_out$slope_est[i])),
      as_chunk(Table_dominance_logSLA_NutNet_ft_tmp_out$slope_grp[i], props = fp_text(vertical.align = "subscript")),
      as_chunk(sprintf(" \u00B1 %.3f", Table_dominance_logSLA_NutNet_ft_tmp_out$slope_SE[i]))
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
doc <- body_add_par(doc, "logSLA", style = "Normal")  # just a raw indicator
doc <- body_add_flextable(doc, ft)
print(doc, target = "./Tables/slopes_dominance_tmp_logSLA_NutNet.docx")

##################################################################################################
logSLA_grid <- seq(min(stat_SLA_clean$logSLA, na.rm = TRUE),
                     max(stat_SLA_clean$logSLA, na.rm = TRUE),
                     length.out = 300)
logSLA_grid

emm_df_logSLA_tg <- as.data.frame(
  emmeans(
    dominance_SLA_NutNet,
    ~ logSLA | tmp_GS * Ntrt,
    at = list(
      logSLA = logSLA_grid,
      tmp_GS  = tmp_levels,
      Ntrt     = c(0, 1)), infer = c(TRUE, TRUE)))  # I will use SE for polygons
emm_df_logSLA_tg  

# prediction dataframe : fit + Tg factor (3 levels)
pred_df_logSLA_tg   <- emm_df_logSLA_tg   %>%
  mutate(
    dominance_fit = emmean,
    Ntrt = factor(Ntrt, levels = c(0, 1)),
    Tg = factor(
      tmp_GS,
      levels = tmp_levels,
      labels = paste0("Tg = ", tmp_levels)))
names(pred_df_logSLA_tg)

# Center where polygons meet: at the median
x_med <- median(stat_SLA_clean$logSLA, na.rm = TRUE)

# Use plus or minus 1 SE polygons 
k <- 1

poly_df_logSLA_tg <- pred_df_logSLA_tg %>%
  group_by(Ntrt, Tg) %>%
  arrange(logSLA, .by_group = TRUE) %>%
  summarise(
    x_min   = first(logSLA),
    x_max   = last(logSLA),
    fit_min = first(dominance_fit),
    se_min  = first(SE),
    fit_max = last(dominance_fit),
    se_max  = last(SE),
    fit_med = approx(x = logSLA, y = dominance_fit, xout = x_med, rule = 2)$y,
    .groups = "drop"
  ) %>%
  mutate(
    upr_min = fit_min + k * se_min,
    lwr_min = fit_min - k * se_min,
    upr_max = fit_max + k * se_max,
    lwr_max = fit_max - k * se_max
  )

poly_left_logSLA_tg <- pmap_dfr(poly_df_logSLA_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                                fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "left",
    x    = c(x_min, x_min, x_med),
    y    = c(upr_min, lwr_min, fit_med)
  )
})
poly_left_logSLA_tg

poly_right_logSLA_tg <- pmap_dfr(poly_df_logSLA_tg, function(Ntrt, Tg, x_min, x_max, fit_min, se_min, fit_max, se_max,
                                                                 fit_med, upr_min, lwr_min, upr_max, lwr_max) {
  tibble(
    Ntrt = Ntrt,
    Tg   = Tg,
    side = "right",
    x    = c(x_max, x_max, x_med),
    y    = c(upr_max, lwr_max, fit_med)
  )
})
poly_right_logSLA_tg


poly_plot_logSLA_tg <- bind_rows(poly_left_logSLA_tg, poly_right_logSLA_tg) %>%
  group_by(Ntrt, Tg, side) %>%
  mutate(id = dplyr::cur_group_id()) %>%
  ungroup()
poly_plot_logSLA_tg


# Tg colors by level (low/mid/high) 
Tg_cols <- c("#7F0000", "#228B22","#08306B")


names(Tg_cols) <- levels(pred_df_logSLA_tg$tmp_GS)

# Plot: polygons + 6 lines (3 Tg x 2 Ntrt)
dominance_logSLA_Tg_NutNet <- ggplot() +
  geom_polygon(data = poly_plot_logSLA_tg,
               aes(x = x, y = y, group = id, fill = Tg),
               alpha = 0.22, color = NA) +
  
  geom_line(data = pred_df_logSLA_tg,
            aes(x = logSLA, y = dominance_fit,
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))

dominance_logSLA_Tg_NutNet

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

stat_SLA_clean <- stat_SLA_clean %>%
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
# Compute slope of logSLA from fixed effects for lambda cover
# Slope = d/d(logSLA) of linear predictor

b_lambda_cover <- lme4::fixef(lambda_SLA_NutNet)
getb <- function(nm) if (nm %in% names(b_lambda_cover)) unname(b_lambda_cover[nm]) else 0

# helper:
grid_pk <- grid_pk %>%
  mutate(
    slope_logSLA_cover =
      getb("logSLA") +
      getb("logSLA:aridity") * arid_val +
      getb("logSLA:tmp_GS") * tmp_val +
      getb("logSLA:Ntrt") * Ntrt +
      getb("logSLA:Ptrt") * Ptrt +
      getb("logSLA:Ktrt") * Ktrt +
      getb("logSLA:aridity:Ntrt") * arid_val * Ntrt +
      getb("logSLA:Ntrt:tmp_GS") * tmp_val * Ntrt +
      getb("logSLA:Ntrt:Ptrt") * Ntrt * Ptrt +
      getb("logSLA:Ntrt:Ktrt") * Ntrt * Ktrt +
      getb("logSLA:Ptrt:Ktrt") * Ptrt * Ktrt +
      getb("logSLA:Ntrt:Ptrt:Ktrt") * Ntrt * Ptrt * Ktrt
  )

# Average across plots within Ntrt groups 
#Ntrt=0: average across Control, P, K, PK (Ptrt/Ktrt combos)
#Ntrt=1: average across N, NP, NK, NPK (Ptrt/Ktrt combos)

grid_avg <- grid_pk %>%
  group_by(tmp_bin, tmp_val, aridity_class, arid_val, Ntrt) %>%
  summarise(slope_logSLA_cover = mean(slope_logSLA_cover), .groups = "drop") %>%
  mutate(
    Ntrt_lab = ifelse(Ntrt == 0, "No N (Control, P, K, PK)", "+N (N, NP, NK, NPK)")
  )

plot_heat_lambda <- function(df_oneN, title_text,
                             shared_scale = NULL,
                             digits = 2) {
  
  p <- ggplot(df_oneN,
              aes(x = aridity_class,
                  y = tmp_bin,
                  fill = slope_logSLA_cover)) +
    
    # Tiles
    geom_tile(color = "white", linewidth = 0.5) +
    
    # Text labels — readable everywhere
    geom_label(
      aes(label = sprintf(paste0("%.", digits, "f"),
                          slope_logSLA_cover)),
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
      name = expression("Slopes of " * SLA),
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
  name = expression("Slopes of " * SLA),
  breaks = seq(-0.1, 0.4, by = 0.1),
  labels = scales::number_format(accuracy = 0.01)
)

# Plot Ambient N
lambda_SLA_N0 <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 0),
  expression(lambda[cover] ~ "\u2013 Ambient N"),
  shared_scale = scale_fill_shared,
  digits = 2
) +
  theme(legend.position = "right")

lambda_SLA_N0

lambda_SLA_N  <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 1),
  expression(lambda[cover] ~ "\u2013 Added N"),
  shared_scale = scale_fill_shared,
  digits = 2
) + theme(legend.position = "right")
lambda_SLA_N

###############################################################################################

######################################### Relative dominance ###############################################

b_lambda_reldom <- lme4::fixef(dominance_SLA_NutNet)
getb <- function(nm) if (nm %in% names(b_lambda_reldom)) unname(b_lambda_reldom[nm]) else 0

# helper:
grid_pk <- grid_pk %>%
  mutate(
    slope_logSLA_reldom =
      getb("logSLA") +
      getb("logSLA:aridity") * arid_val +
      getb("logSLA:tmp_GS") * tmp_val +
      getb("logSLA:Ntrt") * Ntrt +
      getb("logSLA:Ptrt") * Ptrt +
      getb("logSLA:Ktrt") * Ktrt +
      getb("logSLA:aridity:Ntrt") * arid_val * Ntrt +
      getb("logSLA:Ntrt:tmp_GS") * tmp_val * Ntrt +
      getb("logSLA:Ntrt:Ptrt") * Ntrt * Ptrt +
      getb("logSLA:Ntrt:Ktrt") * Ntrt * Ktrt +
      getb("logSLA:Ptrt:Ktrt") * Ptrt * Ktrt +
      getb("logSLA:Ntrt:Ptrt:Ktrt") * Ntrt * Ptrt * Ktrt
  )

# Average across plots within Ntrt groups 
#Ntrt=0: average across Control, P, K, PK (Ptrt/Ktrt combos)
#Ntrt=1: average across N, NP, NK, NPK (Ptrt/Ktrt combos)

grid_avg <- grid_pk %>%
  group_by(tmp_bin, tmp_val, aridity_class, arid_val, Ntrt) %>%
  summarise(slope_logSLA_reldom = mean(slope_logSLA_reldom), .groups = "drop") %>%
  mutate(
    Ntrt_lab = ifelse(Ntrt == 0, "No N (Control, P, K, PK)", "+N (N, NP, NK, NPK)")
  )

plot_heat_lambda <- function(df_oneN, title_text,
                             shared_scale = NULL,
                             digits = 2) {
  
  p <- ggplot(df_oneN,
              aes(x = aridity_class,
                  y = tmp_bin,
                  fill = slope_logSLA_reldom)) +
    
    # Tiles
    geom_tile(color = "white", linewidth = 0.5) +
    
    # Text labels — readable everywhere
    geom_label(
      aes(label = sprintf(paste0("%.", digits, "f"),
                          slope_logSLA_reldom)),
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
      name = expression("Slopes of " * SLA),
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
  name = expression("Slopes of " * SLA),
  breaks = seq(-0.1, 0.4, by = 0.1),
  labels = scales::number_format(accuracy = 0.01)
)

# Plot Ambient N
reldom_SLA_N0 <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 0),
  expression(lambda[cover] ~ "\u2013 Ambient N"),
  shared_scale = scale_fill_shared,
  digits = 2
) +
  theme(legend.position = "right")

reldom_SLA_N0

reldom_SLA_N  <- plot_heat_lambda(
  grid_avg %>% dplyr::filter(Ntrt == 1),
  expression(lambda[cover] ~ "\u2013 Added N"),
  shared_scale = scale_fill_shared,
  digits = 2
) + theme(legend.position = "right")
reldom_SLA_N

###############################################################################################

##################### SLA-lambda cover raw data #################################################
names(stat_SLA_clean)

MI_limits <- range(stat_SLA_clean$aridity, na.rm = TRUE)

lambda_MI_SLA_rawdata <- ggplot(
  stat_SLA_clean,
  aes(x = logSLA, y = mean_log_lambda_cover_all)
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))

lambda_MI_SLA_rawdata


##################### SLA-lambda reldom raw data #################################################
names(stat_SLA_clean)

MI_limits <- range(stat_SLA_clean$aridity, na.rm = TRUE)

dominance_MI_SLA_rawdata <- ggplot(
  stat_SLA_clean,
  aes(x = logSLA, y = mean_log_lambda_Ci_all)
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))

dominance_MI_SLA_rawdata


##################### SLA-lambda cover raw data Tg #################################################

Tg_limits <- range(stat_SLA_clean$tmp_GS, na.rm = TRUE)

lambda_SLA_rawdata_Tg <- ggplot(
  stat_SLA_clean,
  aes(x = logSLA, y = mean_log_lambda_cover_all)
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))
lambda_SLA_rawdata_Tg

##################### chi-lambda reldom Tg raw data #################################################

Tg_limits <- range(stat_SLA_clean$tmp_GS, na.rm = TRUE)

reldom_SLA_rawdata_Tg <- ggplot(
  stat_SLA_clean,
  aes(x = logSLA, y = mean_log_lambda_Ci_all)
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
  xlab(expression(log~italic(SLA)~"(mm"^2~"mg"^-1*")"))

reldom_SLA_rawdata_Tg

