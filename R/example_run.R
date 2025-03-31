# # library(Budgeon)
# library(truncnorm)
# library(rlang)
# library(dplyr)
# library(ggplot2)
# library(gridExtra)
# df <- generate_synth_data(13, 3, 2)
# PET_pos_threshold <- 0.92
# value_name <- "Marker_level"
# time_name <- "Time"
# id_name <- "ID"
# 
# individual_coeffs <- get_individ_lm(df, value_name, time_name, id_name)
# 
# 
# individual_coeffs <- individual_coeffs[individual_coeffs$Slope >= 0,]
# 
# eval_at_midpoint = eval_lmm(df[df[[id_name]] %in% individual_coeffs[[id_name]],], id_name = id_name, time_name = time_name, individual_coeffs) 
# 
# df_tau <- get_integration_estimates(eval_at_midpoint, id_name = id_name, degree)
# 
# a <- ggplot(df, aes(x = Time, y = Marker_level, group = ID)) + geom_point(alpha = 0.8) + 
#   geom_line(alpha = 0.4) + theme_bw() + xlab("Time") + ylab("Simulated Marker Level") +
#   ggtitle("A. Simulated Longitudinal Data")
# 
# b <- ggplot(eval_at_midpoint, aes(x = eval_at_midpoint, y = Slope)) + geom_point() + 
#   theme_bw() + geom_smooth(method = "gam", se = FALSE, colour = "#782223") +
#   ggtitle("B. Results of Mixed Effects Models") + xlab("Estimated Mean by Participant") +
#   ylab("Estimated Slope by Participant")
# 
# c <- ggplot(df_tau, aes(x = mu_tau, y = tau)) + geom_point() + theme_bw() +
#   geom_smooth(method = "gam", se = FALSE, colour = "#782223") +
#   ggtitle("C. Integrate reciprocal of B") + xlab("Marker Level") + ylab("Disease Progression Time")
# 
# d <- ggplot(df_tau, aes(y = mu_tau, x = tau)) + geom_point() + theme_bw() +
#   geom_smooth(method = "gam", se = FALSE, colour = "#782223") + 
#   ggtitle("D. Inversion of C") + xlab("Disease Progression Time") + ylab("Marker Level")
# 
# grid.arrange(a, b, c, d, nrow = 2, ncol = 2)
# 
# hor_adj <- get_horizontal_adjustment(df_tau, PET_pos_threshold)
# 
# df_tau$Time_to_Positivity <- df_tau$tau - hor_adj
