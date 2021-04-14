library(ggplot2)
load("./data/sim_bayes_data.RData")
load("./data/sim_bayes_fit_results.RData")
load("./data/sim_bayes_fit_fix_results.RData")
N_warmup = 500
N_sample = 500
N_chain = 4

df <- data.frame(tausq = c(c(fit_results["tausq", , ]), 
                           c(fit_results_fixphi["0.2phi", "tausq", , ]),
                           c(fit_results_fixphi["0.5phi", "tausq", , ]),
                           c(fit_results_fixphi["phi", "tausq", , ]),
                           c(fit_results_fixphi["2phi", "tausq", , ]),
                           c(fit_results_fixphi["5phi", "tausq", , ])),
                 sigmasq = c(c(fit_results["sigmasq", , ]),
                             c(fit_results_fixphi["0.2phi", "sigmasq", , ]),
                             c(fit_results_fixphi["0.5phi", "sigmasq", , ]),
                             c(fit_results_fixphi["phi", "sigmasq", , ]),
                             c(fit_results_fixphi["2phi", "sigmasq", , ]),
                             c(fit_results_fixphi["5phi", "sigmasq", , ])),
                 kappa = c(c(fit_results["kappa", , ]),
                           c(fit_results_fixphi["0.2phi", "kappa", , ]),
                           c(fit_results_fixphi["0.5phi", "kappa", , ]),
                           c(fit_results_fixphi["phi", "kappa", , ]),
                           c(fit_results_fixphi["2phi", "kappa", , ]),
                           c(fit_results_fixphi["5phi", "kappa", , ])),
                 sample_size = rep(rep(c(1, 2, 3), 
                                       N_chain * N_sample), 6),
                 type = rep(c(3, 1, 2, 4, 5, 6), 
                            each = length(c(fit_results["kappa", , ]))))

df$sample_size <- factor(df$sample_size, levels = c(1, 2, 3), 
                  labels = c("400", "900", "1600"))

df$type <- factor(df$type, levels = c(1, 2, 3, 4, 5, 6), 
                  labels = c("0.2phi", "0.5phi", "unknown", 
                             "phi", "2phi", "5phi"))
                  #labels = c( "0.2phi", "MCMC", "phi", "5phi"))


p_sigmasq_compar <- ggplot(data = df, 
                       aes(y = sigmasq, x = type, 
                           fill = sample_size, color = sample_size)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("grey90", "grey70", "grey50")) + 
  scale_color_manual(values=c("black", "black", "black"))+ 
  geom_hline(yintercept = sigmasq, linetype = "dashed") + 
  scale_y_log10() + 
  scale_x_discrete(labels = c('0.2phi' = expression(0.2 * phi[0]),
                              '0.5phi' = expression(0.5 * phi[0]),
                              'unknown' = 'unknown',
                              'phi' = expression(phi[0]),
                              '2phi' = expression(2 * phi[0]),
                              '5phi' = expression(5 * phi[0]))) +
  ylab(expression(paste("Posterior of ", sigma^{2}, " (log scale)"))) + 
  xlab(expression(phi)) + labs(fill = "Sample Size", color = "Sample Size") +
  theme_bw(base_size = 12) + theme(legend.position='none')

width <- 8.0
height <- 3.0

setEPS()
postscript("./pics/simBayes_sigmasq.eps", width = width, height = height)
print(p_sigmasq_compar)
dev.off()

p_tausq_compar <- ggplot(data = df, 
                           aes(y = tausq, x = type, 
                               fill = sample_size, color = sample_size)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("grey90", "grey70", "grey50")) + 
  scale_color_manual(values=c("black", "black", "black"))+ 
  geom_hline(yintercept = tausq, linetype = "dashed") + 
  scale_y_log10() + 
  scale_x_discrete(labels = c('0.2phi' = expression(0.2 * phi[0]),
                              '0.5phi' = expression(0.5 * phi[0]),
                              'unknown' = 'unknown',
                              'phi' = expression(phi[0]),
                              '2phi' = expression(2 * phi[0]),
                              '5phi' = expression(5 * phi[0]))) +
  ylab(expression(paste("Posterior of ", tau^{2}, " (log scale)"))) + 
  xlab(expression(phi)) + labs(fill = "Sample Size", color = "Sample Size") +
  theme_bw(base_size = 12) + theme(legend.position='none')

setEPS()
postscript("./pics/simBayes_tausq.eps", width = width, height = height)
print(p_tausq_compar)
dev.off()


p_kappa_compar <- ggplot(data = df, 
                           aes(y = kappa, x = type, 
                               fill = sample_size, color = sample_size)) + 
  geom_boxplot() + 
  scale_fill_manual(values=c("grey90", "grey70", "grey50")) + 
  scale_color_manual(values=c("black", "black", "black"))+ 
  geom_hline(yintercept = sigmasq * phi, linetype = "dashed") + 
  ylab(expression(paste("Posterior of ", kappa, " (log scale)"))) + 
  xlab(expression(phi)) + labs(fill = "Sample Size", color = "Sample Size") +
  scale_y_log10() + 
  scale_x_discrete(labels = c('0.2phi' = expression(0.2 * phi[0]),
                              '0.5phi' = expression(0.5 * phi[0]),
                              'unknown' = 'unknown',
                              'phi' = expression(phi[0]),
                              '2phi' = expression(2 * phi[0]),
                              '5phi' = expression(5 * phi[0]))) +
  theme_bw(base_size = 12) + theme(legend.position='none') 

setEPS()
postscript("./pics/simBayes_kappa.eps", width = width, height = height)
print(p_kappa_compar)
dev.off()



