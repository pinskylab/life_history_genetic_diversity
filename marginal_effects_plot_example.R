######### Abslat figures #######

abslat_model_pi <- lm(logpi ~ abslat_scale, data = mtdna_small_pi, 
                      na.action = "na.fail")

#### Predict ####
#marginal effects (effects package)
abslat_eff <- plot_model(abslat_model_pi, type = "pred",
                         terms = "abslat_scale [all]")

#pull out marginal effects dataframe
abslat_eff_data <- as.data.frame(abslat_eff$data)

#unscale lat
#use same scaled:center & scaled:scale from original data
abslat_scale <- scale(mtdna_small_pi$abslat) #bc had to convert to numeric to run model/calculate marginal effects
abslat_eff_data$abslat <- (abslat_eff_data$x * attr(abslat_scale, "scaled:scale")) + 
  attr(abslat_scale, "scaled:center")

#unlog pi
abslat_eff_data$unlog_pi <- 10^(abslat_eff_data$predicted)
abslat_eff_data$unlog_conf.low <- 10^(abslat_eff_data$conf.low)
abslat_eff_data$unlog_conf.high <- 10^(abslat_eff_data$conf.high)

#### Plot abslat ####
#for legend
colors <- c("10-degree binned means" = "#3F6DAA", "Regression" = "black")

#plot
mtdna_pi_abslat_plot_both <- ggplot() + 
  geom_line(data = abslat_eff_data, 
            aes(x = abslat, y = unlog_pi, color = "Regression"), size = 6) + 
  geom_ribbon(data = abslat_eff_data, 
              aes(x = abslat, ymin = unlog_conf.low, ymax = unlog_conf.high, color = "Regression"), alpha = 0.1) + 
  geom_point(data = abslat_pi_binned_means, 
             aes(x = X, y = pi_mean, color = "10-degree binned means", size = pi_count), shape = "square") + 
  geom_errorbar(data = abslat_pi_binned_means, 
                aes(x = X, ymin = mean_lowerSE, ymax = mean_upperSE, 
                    color = "10-degree binned means"), width = 1.75, size = 3) + 
  xlab("absolute latitude") + ylab("mtDNA pi") + labs(color = "Legend") + 
  scale_y_continuous(limits = c(0, 0.0075)) + 
  scale_color_manual(values = colors) +
  scale_size_continuous(breaks = c(5, 50, 100, 200, 500), 
                        range = c(10, 20))
mtdna_pi_abslat_plot_annotated_both <- mtdna_pi_abslat_plot_both + theme_bw() + coord_flip() + 
  theme(panel.border = element_rect(size = 1), axis.title = element_text(size = 110), 
        axis.ticks = element_line(color = "black", size = 1), 
        axis.text = element_text(size = 100, color = "black"), 
        axis.line = element_line(size = 2, color = "black"),
        legend.position = "top", legend.box = "vertical", 
        legend.text = element_text(size = 100), 
        legend.key.size = unit(3, "cm"),
        legend.title = element_blank())
mtdna_pi_abslat_plot_annotated_both