#' A function to draw diganostic plot
#'
#' @param sfpca_data: The prepared data from prepare_data() function (list)
#' @param sfpca_model: The optimal sfpca model
#' @import ggplot2
#' @import rstan
#' @import bayesplot

plot_diagnostic <- function(sfpca_data, sfpca_model){
  N <- sfpca_data$num_subjects
  loo_best <- sfpca_model$looic
  pkdf <- data.frame(pk = loo_best$diagnostics$pareto_k, id = 1:N)
  print(ggplot2::ggplot(pkdf, aes(x = id,y = pk)) +
          geom_point(shape = 3, color = "blue") +
          labs(x = "Observation left out", y = "Pareto shape k") +
          geom_hline(yintercept = 0.7, linetype = 2, color = "red", size = 0.2) +
          ggtitle("PSIS-LOO diagnostics") +
          theme_classic() +
          theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"),
                axis.text.x = element_text(size=10, face="bold"),
                axis.text.y = element_text(size=10, face="bold"),
                axis.title.x = element_text(size=12, face="bold"),
                axis.title.y = element_text(size=12, face="bold")) )

  sa <- sfpca_model$sa
  Nsamples <- sfpca_model$Nsamples
  Nchains <- sfpca_model$Nchains
  Ynew <- rstan::extract(sa, "Ynew", permuted = FALSE)
  V <- sfpca_data$visits.vector
  Ynew_transform <- matrix(rep(0, Nsamples / 2 * Nchains * sum(V)), ncol = sum(V))
  ind <- 0
  for (i in 1:(Nsamples / 2)){
    for (j in 1:Nchains){
      ind <- ind + 1
      Ynew_transform[ind, ] < Ynew[i, j, ]
    }
  }
  Ynew_mean <- colMeans(Ynew_transform)
  bayesplot::color_scheme_set("brightblue")
  k <- sfpca_model$pc
  d <- sfpca_model$knot
  print(bayesplot::ppc_dens_overlay(sfpca_data$data$response, Ynew_transform) +
          ggplot2::ggtitle(paste(k, 'pc_', d, 'knot')) )
}
