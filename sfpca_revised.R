# set the script directory as the working directly

# Three main functions
# 1. prepare_data
# 2. visu_data: time plots, by group also (to be added)
# 3. get spline basis

### warning: in basis_setup_sparse, K = # basis, Q= # PCs
### warning: need to check whether my basis results are the same as in basis_setup_sparse

### standarize.y and center.y cannot do separately; has to choose from one of them; thus revise code on 08/16/2019
### see application to MG's code

### prepare data for sfpca model
prepare_data = function(data, unique_subject_id, time_var, response, transform.y='standardize', scale.time=FALSE, group.var=NA){
	# data: target longitudinal data for analysis (must be a data frame)
	# unique_subject_id: the column name corresponding to unique subject id in the data (string)
	# time_var: the column name corresponding to the time variable in the data (string)
	# response: the column name of the intersted response variable
	# standardize.y: the option of whether or not to standardize response variable (True/False) with mean 0 and sd 1
	# scale.time: the option of whether or not to scale the time variable to be within [0, 1] (True/False)

  #print warnings
  data_check = data[, c(as.character(unique_subject_id), as.character(time_var))]
  if (sum(duplicated(data_check)) != 0) return("Each subject need to have unique measurement at each time point ")
  
	# create new ID
	data$ID = as.numeric(as.numeric(data[, unique_subject_id]))
	N = length(unique(data$ID)) # total number of unique subjects

	# convert group id to be numeric
	if (!is.na(group.var)){
	data[, group.var] = as.numeric(as.factor(data[, group.var]))
	}

	# create time 
	if (scale.time == TRUE){
		data$time = (data[, time_var] - min(data[, time_var])) / (max(data[, time_var]) - min(data[, time_var]))
	} else{
		data$time = data[, time_var]
	}

	T = length(unique(data$time)) # total number of sampling time points

	# transform response (code updated on 08/16/2019)
	if (transform.y == 'standardize'){
		data$response = (data[, response] - mean(data[, response], na.rm=T)) / sd(data[, response], na.rm=T)
	} else if (transform.y == 'center'){
		data$response = (data[, response] - mean(data[, response], na.rm=T)) 
	} else {
		data$response = data[, response]
	}	

	
	# re-order the data by ID and time
	data = data[order(data$ID, data$time), ]

	# create visits vector, response and time matrix
	ID.list = unique(data$ID)
	#visits.vector = matrix(rep(0, N*1), nrow=N)
	visits.vector = vector(mode = "numeric", length = N)
	response.list = NULL
	time.matrix = matrix(rep(0, N*T), nrow=N)
	# response.matrix = time.matrix = matrix(rep(0, N*T), nrow=N)
	# rownames(visits.vector) = rownames(response.matrix) = rownames(time.matrix) = ID.list
	# colnames(visits.vector) = 'n_visits'
	# rownames(response.matrix) = rownames(time.matrix) = ID.list
	# colnames(response.matrix) = colnames(time.matrix) = seq(1, T, 1)

	# visits index for each individual when stacking the data
	visits.stop = vector(mode = "numeric", length = N)

	# size index for each individual in covariance matrix
	cov.start = cov.stop = vector(mode = "numeric", length = N)

	# group id based on interested group
	id_group = vector(mode = "numeric", length = N)

	for(i in 1:N){ 
		# visits vector
		subject_i = data[data$ID==ID.list[i],]
		subject_i$n_visits = dim(subject_i)[1]	
		visits.vector[i] = unique(subject_i$n_visits)

		# visits index
		visits.stop[i] = sum(visits.vector)

		# covariance size index
		cov.stop[i] = sum(visits.vector^2)

	    # response matrix
	    # response.matrix[i, ] = c(subject_i$response, rep(0, T - unique(subject_i$n_visits)))
	    response.list = c(response.list, subject_i$response)

		# time matrix
		time.matrix[i, ] = c(subject_i$time, rep(0, T - unique(subject_i$n_visits)))

		# group id based on interested group
		if (!is.na(group.var)){
		id_group[i] = unique(subject_i[, group.var])
		}

		rm(subject_i)
	}	
	visits.start = c(1, visits.stop[-N] + 1)
	cov.start = c(1, cov.stop[-N] + 1)
	cov.size = sum(visits.vector^2)

	prepared_data = list(data=data, num_subjects=N, num_times=T, response.list=response.list, time.matrix=time.matrix,
		                 visits.vector=visits.vector, visits.start=visits.start, visits.stop=visits.stop,
		                 cov.start=cov.start, cov.stop=cov.stop, cov.size=cov.size, id_group=id_group)
	return(prepared_data)
}


### set up spline basis for sparse data
basis_setup_sparse = function(prepared_data, nknots, orth=TRUE, delta=1/10000){
	# prepared_data: longitudinal data after applying prepared_data() function
	# knots: user-defined number of knots
	# orth: default setting for orth should be TRUE (after discussed with Wes on 02/13/2019)

	time_var = prepared_data$data$time
	num_subjects = prepared_data$num_subjects
	num_times = prepared_data$num_times
	S = prepared_data$time.matrix
	V = prepared_data$visits.vector

	# continuous time interval
	time_unique = sort(unique(time_var))
	time_min = min(time_unique)
	time_max = max(time_unique)
	time_cont = seq(time_min, time_max / delta) * delta # chop the entire time interval into many small subintervals
	time_cont = round(time_cont / delta)*delta # to avoid rounding error?

	# specify placement of knots
	qs = 1/(nknots + 1)
	knots = quantile(time_unique, qs)
	if(nknots > 1){
		for(q in 2:nknots){
			knots = c(knots, q*quantile(time_unique,qs))
		}
	}

	knots = as.vector(knots)


	# obtain cubic spline basis
	library('splines')

	##### option 1: force all matrices to have the same size
	# ## 1. for densely sampled time points
	# phi_t_cont=list()
	# phi_t_cont = bs(time_cont, knots=knots, degree=3,intercept=TRUE) # cubic spline, degree=spline_degree

	# # Gram-Schmidt Orthonormalization
	# temp = phi_t_cont
	# orth = TRUE
	# K = nknots + 4 # num of spline basis 

	# for(k in 1:K){
	# 	if(orth==TRUE){
	# 		if(k > 1){
	# 			for(q in 1:(k-1)){
	# 				temp[,k]=temp[,k]-(sum(temp[,k]*temp[,k-q])/
	# 					sum(temp[,k-q]^2))*temp[,k-q];
	# 			}
	# 		}
	# 	}		
	#     temp[,k]=temp[,k]/sqrt(sum(temp[,k]*temp[,k]))
	# }

	# phi_t_cont=t(sqrt(1/delta)*temp)

	# ## 2. for sparsely sampled time points
	# phi_t=list()
	# for(i in 1:num_subjects){
	# 	phi_t[[i]] = array(0,dim=c(K, V[i])) # phi_t: K (number of basis function) * number of total visit for each subject

	# 	for(k in 1:K){
	# 		for(t in 1:V[i]){
	# 			phi_t[[i]][k, t] = phi_t_cont[k, abs(time_cont - S[i, t]) == min(abs(time_cont - S[i, t]))]
	# 		}
	# 	}

	# 	## fill up zeros to make matrices same size
	# 	miss_visits = num_times - V[i]
	# 	fill_zeros = matrix(rep(0, K*miss_visits), nrow=K)
	# 	phi_t[[i]] = cbind(phi_t[[i]], fill_zeros)

	# }

	#### option 2: stack subjects and visits
	## 1. for densely sampled time points
	phi_t_cont=list()
	phi_t_cont = bs(time_cont, knots=knots, degree=3,intercept=TRUE) # cubic spline, degree=spline_degree

	### the same as in setup_basis_sparse so far

	# Gram-Schmidt Orthonormalization
	temp = phi_t_cont
	K = nknots + 4 # num of spline basis 

	for(k in 1:K){
		if(orth==TRUE){
			if(k > 1){
				for(q in 1:(k-1)){
					temp[,k]=temp[,k]-(sum(temp[,k]*temp[,k-q])/
						sum(temp[,k-q]^2))*temp[,k-q];
				}
			}
		}		
	    temp[,k]=temp[,k]/sqrt(sum(temp[,k]*temp[,k]))
	}

	phi_t_cont=t(sqrt(1/delta)*temp)

	## 2. for sparsely sampled time points
	phi_t_stacked=NULL
	phi_t=list()
	for(i in 1:num_subjects){
		phi_t[[i]] = array(0,dim=c(K, V[i])) # phi_t: K (number of basis function) * number of total visit for each subject

		for(k in 1:K){
			for(t in 1:V[i]){
				phi_t[[i]][k, t] = phi_t_cont[k, abs(time_cont - S[i, t]) == min(abs(time_cont - S[i, t]))]
			}
		}

		# stack subjects and visits: number of visits as rows, and number of basis as columns
		phi_t_stacked = rbind(phi_t_stacked, t(phi_t[[i]]))
	}

	results_basis = list(knot_place=knots, time_cont=time_cont, orth_spline_basis_sparse=phi_t, 
						 orth_spline_basis_sparse_stacked=phi_t_stacked, orth_spline_basis_cont=phi_t_cont)
	return(results_basis)
}

###perform post hoc rotation
post_hoc_rotation <- function(prepared_data, model_index, npcs, nknots, sa, Nchains, Nsamples){
  Sigma = extract(sa,"Sigma",permuted=FALSE)
  W = extract(sa,"W",permuted=FALSE)
  sigma_eps = extract(sa,"sigma_eps",permuted=FALSE)
  theta_mu = extract(sa,"theta_mu",permuted=FALSE)
  alpha = extract(sa,"alpha",permuted=FALSE)
  Theta = extract(sa,"Theta",permuted=FALSE)
  
  ## Reshape parameters and reorient loadings with PCA rotation 
  N = prepared_data$num_subjects
  K = npcs
  Q = nknots + 4
  
  theta_mu_new = array(0, dim=c(Q, Nchains*Nsamples/2))
  alpha_old = alpha_new = array(0, dim=c(K, N, Nchains*Nsamples/2)) 
  Theta_old = Theta_new = array(0, dim=c(Q, K, Nchains*Nsamples/2))
  W_old = array(0, dim=c(Q, Q, Nchains*Nsamples/2)) 
  
  ind = 0
  prop_var = NULL
  for(i in 1:dim(W)[1]){
   for(j in 1:dim(W)[2]){
      print(ind)
      ind = ind + 1
      theta_mu_new[,ind] = array(theta_mu[i,j,])
      alpha_old[,,ind] = t(array(alpha[i,j,],dim=c(N, K)))
      Theta_old[,,ind] = array(Theta[i,j,],dim=c(Q, K))
      W_old[,,ind] = array(W[i,j,],dim=c(Q,Q)) 
      
      eigen_temp_sigma=eigen(W_old[,,ind])
      v_temp=eigen_temp_sigma$vectors
      d_temp=eigen_temp_sigma$values 
      prop_var = rbind(prop_var, d_temp/sum(d_temp)) # proportion of variance explained by each PC
      
      for(com in 1:length(d_temp)){
        if(!(d_temp[com]-Re(d_temp[com])==0)){
          d_temp[com]=-1*10^5
        }
      }
      pos_temp=array(0,dim=c(K,1))
      for(pos in 1:K){
        pos_temp[pos]=(1:length(d_temp))[max(d_temp)==d_temp]
        d_temp[pos_temp[pos]]=-1e+5
      }
      
      Theta_new[,,ind]=v_temp[,pos_temp]
      for(k in 1:K){
        Theta_new[, k, ind]=sign(Theta_new[1,k,ind]) * Theta_new[,k,ind]
      }
      
      alpha_new[,, ind] = t(Theta_new[,,ind]) %*% Theta_old[,,ind] %*% alpha_old[,,ind]
    }
  }
  prop_var_avg_origin = colMeans(prop_var)
  (prop_var_avg = paste(round(colMeans(prop_var)*100, 2), '%', sep=''))
  #rename Q
  li <- list(num_subjects = N, npcs = K, Q = Q, alpha_new = alpha_new, 
            theta_mu_new = theta_mu_new, Theta_new = Theta_new, prop_var_avg_origin = prop_var_avg_origin, 
            prop_var_avg = prop_var_avg)
  return(li)
}



output_results <- function(prepared_data, npcs, vars_select, results_list, results_basis){
  ALPHA_array = results_list$alpha_new
  MU_array = results_list$theta_mu_new
  THETA_array = results_list$Theta_new
  phi_t_cont = results_basis$orth_spline_basis_cont
  phi_t = results_basis$orth_spline_basis_sparse
  time_cont = results_basis$time_cont
  
  nloop=dim(ALPHA_array)[3]
  first=1
  last=nloop
  
  MU_mean = MU_array[, first] #mean function across sampling sessions
  ALPHA_mean = ALPHA_array[,,first] # mean factor scores
  THETA_mean = THETA_array[,,first] # mean factor loading
  
  for(iter in 2:nloop){
    MU_mean = MU_mean + MU_array[, iter]
    ALPHA_mean = ALPHA_mean + ALPHA_array[,,iter]
    THETA_mean = THETA_mean + THETA_array[,,iter]
  }
  
  MU_mean=cbind(MU_mean/(last-first+1))
  ALPHA_mean=cbind(ALPHA_mean/(last-first+1))
  THETA_mean=cbind(THETA_mean/(last-first+1))
  
  Mu_functions = t(bdiag(cbind(phi_t_cont)))%*%MU_mean
  FPC_mean=t(phi_t_cont)%*%THETA_mean
  
  if(npcs == 1){
    ### create data frame containing needed information ####
    df = prepared_data$data[, vars_select]
    Y_sparse = list()
    time_sparse = list()
    scores = data.frame(t(ALPHA_mean)) 
    df$fpc1 = 0 # principle component scores
    
    i = 0
    for (pid in unique(df$ID)){
      i = i + 1
      Y_sparse[[i]] = df$response[df$ID == pid]
      time_sparse[[i]] = df$time[df$ID == pid]
      df$fpc1[df$ID == pid] = scores[i]
    }
    df$fpc1 = as.numeric(df$fpc1) # data type issue 
    
    Fits_sparse=list()
    for(i in 1:N){
      Fits_sparse[[i]] = t(phi_t[[i]]) %*% MU_mean + t(phi_t[[i]]) %*% THETA_mean %*% ALPHA_mean[i]
    }
    
    df$Y_sparse = unlist(Y_sparse) 
    df$Fits_sparse = unlist(Fits_sparse)
    df$residuals = df$Y_sparse - df$Fits_sparse
  } else {
    ### create data frame containing needed information ####
    df = prepared_data$data[, vars_select]
    Y_sparse = list()
    time_sparse = list()
    scores = data.frame(t(ALPHA_mean)) 
    for(k in 1:npcs){
      names(scores)[k] = paste('fpc', k, sep = '')
      df[names(scores)[k]] = 0 # principle component scores  # it depends of PCs (better to choose number of PCs as input)
    }
    
    i = 0
    for (pid in unique(df$ID)){
      i = i + 1
      Y_sparse[[i]] = df$response[df$ID == pid]
      time_sparse[[i]] = df$time[df$ID == pid]
      for(k in 1:npcs){
        df[,names(scores)[k]][df$ID == pid] = scores[i, k]
      }
    }
    
    Fits_sparse=list()
    for(i in 1:N){
      Fits_sparse[[i]] = t(phi_t[[i]]) %*% MU_mean + t(phi_t[[i]]) %*% THETA_mean %*% ALPHA_mean[, i]
    }
    
    df$Y_sparse = unlist(Y_sparse) # check: sum(df$Y_sparse != df$response) == 0
    df$Fits_sparse = unlist(Fits_sparse)
    df$residuals = df$Y_sparse - df$Fits_sparse
  }
  
  return_list = list(df = df, Mu_functions = Mu_functions, time_sparse = time_sparse,
                       Y_sparse = Y_sparse, FPC_mean = FPC_mean)
  return(return_list)
}

plot_qqplot <- function(response){
  par(mfrow=c(2,2))
  qqPlot(response)
  qqPlot(log(response))
  qqPlot(sqrt(response))
  qqPlot(response)
}

plot_bygroup <- function(df, xaxis, yaxis, unique_id, groupby){
  ggplot(df, aes(x=xaxis, y=yaxis, group=unique_id, color=groupby)) + 
    geom_line(alpha=0.2) + 
    geom_smooth(se=TRUE, size=2, aes(group=groupby, fill=groupby), level=0.95) +
    xlab('xaxis') + ylab("yaxis") + ggtitle('Gut') + 
    theme(axis.title.x = element_text(face="bold"),
          axis.title.y = element_text(face="bold"),
          legend.position="top") 
}

plot_k_diagnostic <- function(pkdf, id, pk){
  ggplot(pkdf, aes(x=id,y=pk)) + geom_point(shape=3, color="blue") +
    labs(x="Observation left out", y="Pareto shape k") +
    geom_hline(yintercept = 0.7, linetype=2, color="red", size=0.2) +
    ggtitle("PSIS-LOO diagnostics") + theme_classic() + 
    theme(plot.title = element_text(hjust = 0.5, size=15, face="bold"),
          axis.text.x= element_text(size=10, face="bold"),
          axis.text.y= element_text(size=10, face="bold"),
          axis.title.x= element_text(size=12, face="bold"),
          axis.title.y= element_text(size=12, face="bold")) 
}


plot_mean_curve <- function(time_cont, data, Mu_functions, sigma_y, mu_y, time_sparse, Y_sparse){
  plot(time_cont*(max(data$Time) - min(data$Time)) + min(data$Time), Mu_functions*sigma_y + mu_y, type="l",ylim=c(-1, 4),
       xlab='Days of life', ylab='shannon diversity', lwd=5, col=4, font.lab=2, cex.lab=1.2)
  for(i in 1:N){
    lines(time_sparse[[i]]*(max(data$Time) - min(data$Time)) + min(data$Time),Y_sparse[[i]]*sigma_y + mu_y,type="l",lwd=.25)
  }
}


