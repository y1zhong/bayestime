#' plot group mean of fpc scores
#'
#' @param output: The model output list from output_results() function
#' @param group_name One column name of interested group variables in data
#' @export
plot_fpc_group_mean <- function(output, group_name){
  data <- output$df
  K <- output$rotation$npcs

  fpcs <- c()
  for (k in 1:K){
    fpcs <- c(fpcs, paste('fpc', k, sep = ''))
    data_temp <- data[, c('ID', 'delivery', 'fpc1')]
    classes <- factor(data_temp$delivery)
    data_wide <- data_temp[data_temp[, 'delivery'] == classes[1], ]
  }
  ## effect of mean of two groups
  for (k in 1:K){
    #pdf(paste(paste('shannon_sfpca/FPCs_mean_PC_GroupMean_shannon_fecal', k, sep=''), 'pdf', sep='.'), width=5, height=5)
    scores_mu_g1 = mean(df[df$birth_mode_ms %in% 'CS', paste('fpc', k, sep='')])
    scores_mu_g2 = mean(df[df$birth_mode_ms %in% 'CSseed', paste('fpc', k, sep='')])
    scores_mu_g3 = mean(df[df$birth_mode_ms %in% 'Vag', paste('fpc', k, sep='')])
    plot(time_cont*max(dat$Time), Mu_functions*sigma_y + mu_y, type="n",
         lwd=2,col=1, xlab='Time', ylab='shannon', font.lab=2, cex.lab=1.2, ylim=c(1.7,3.5))

    lines(time_cont*max(dat$Time), (Mu_functions + FPC_mean[,k]*scores_mu_g1)*sigma_y + mu_y,type="l",lwd=3,lty=2,col=2) # red
    lines(time_cont*max(dat$Time), (Mu_functions + FPC_mean[,k]*scores_mu_g2)*sigma_y + mu_y,type="l",lwd=3,lty=2,col=3) # green
    lines(time_cont*max(dat$Time), (Mu_functions + FPC_mean[,k]*scores_mu_g3)*sigma_y + mu_y,type="l",lwd=3,lty=2,col=4) #
    title(main=paste(paste('PC', k, sep=' '), ' (', prop_var_avg[k], ' )', sep=''))
    #axis(1, font=2) # make x-axis ticks label bold
    legend('bottomright', c('CS', 'CSseed', 'Vaginal'), lty=c(2,2,2), lwd=c(3,3,3), col=c(2,3,4), bty='n', cex=0.5)
    #dev.off()
  }

}
