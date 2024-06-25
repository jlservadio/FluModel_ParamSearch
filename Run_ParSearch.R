
rm(list = ls())

library(Rcpp)
library(tictoc)

setwd('~/work/Param_Search/RR_2L_I')

source('~/work/Param_Search/Clean_Res.R')
source('~/work/Param_Search/Incidence_checks.R')


dt = Sys.Date()

#
# Set up model configurations (population(s), R classes)
#


group_labels = c('E1', 'E2')
populations = c(500000, 500000)

num.r = 2


starts = c('S' = 0.6, 'I' = 10/1000000, 'II' = 10/1000000, 'R' = round(0.25 / num.r / 3, 3))
starts[1] = 1 - sum(starts[2:3]*3) - (starts[4]*num.r*3)


file.name = 'Code/ModelCode_RR_2L_I.txt'



#
# Flu characteristics
#



tic()
# write model file

source('Code/write_model.R') 
toc()






#
# Single Run
#

tic()
source(file.name) 
toc()



tic()
out = sirout_groups(beta_H1 = 0.22, beta_B = 0.21, beta_H3 = 0.24, 
	h_H1 = 0, h_B = 0, h_H3 = 0, 
	d_H1 = 0, d_B = 0, d_H3 = 0, 
	nu = 1/(5/2), rho = 1/(400/num.r), 
	sigma12 = 1 - 0.27, sigma13 = 1 - 0.3, sigma23 = 1 - 0.17, 
	eta12 = 0.7, eta21 = 0.5,
	N0 = 1000000, tf = 365*20, 
	days_import = 100)
toc()
names(out) = all.cols

mod.out = as.data.frame(out)	
summary(mod.out$N)




#
# Parameter search
#


keep.pars = NULL

n.its = 100000
source('~/work/Param_Search/set_nits.R')


tic()
for(i in 1:n.its) {
	
	# draw betas
	m.bs = round(runif(1, 0.1, 0.5), 4)
	betas = sort(round(rnorm(3, m.bs, sd = 0.01), 4))
	new.params = betas[c(2, 1, 3)]
	
	# draw duration
	#new.params = c(new.params, round(runif(1, 120, 365*6)))
	new.params = c(new.params, round(runif(1, 720, 365*6)))
	
	# draw cross immunity
	for (j in 1:3) { new.params = c(new.params, round(runif(1, 0.05, 0.95), 4)) }
	
	# draw mixing
	for (j in 1:2) { new.params = c(new.params, round(runif(1, 0.05, 1), 4)) }
	
	# draw importation freq
	#new.params = c(new.params, round(runif(1, 30, 365)))
	new.params = c(new.params, round(runif(1, 130, 365)))
	
	
	names(new.params) = c('beta_H1', 'beta_B', 'beta_H3', 'duration', 'sigma12', 'sigma13', 'sigma23', 
		'eta12', 'eta21', 'impfreq')
		
	
	out = sirout_groups(beta_H1 = new.params['beta_H1'], beta_B = new.params['beta_B'], beta_H3 = new.params['beta_H3'], 
		h_H1 = 0, h_B = 0, h_H3 = 0, 
		d_H1 = 0, d_B = 0, d_H3 = 0, 
		nu = 1/(5/2), rho = 1/(new.params['duration']/num.r), 
		sigma12 = 1 - new.params['sigma12'], sigma13 = 1 - new.params['sigma13'], sigma23 = 1 - new.params['sigma23'], 
		eta12 = new.params['eta12'], eta21 = new.params['eta21'],
		N0 = 1000000, tf = 365*20, 
		days_import = new.params['impfreq'])
	names(out) = all.cols
	
	mod.out = as.data.frame(out)	
	
	out.inc = get.inc(mod.out, 365*10)
	
	
	checks = check.par.traj(out.inc)
	
	if (i == 1 || sum(checks) > 6) {
		keep.pars = rbind(keep.pars, c(new.params, c(checks, sum(checks))))
		if (!is.null(dim(keep.pars))) { if (dim(keep.pars)[1] %% 10 == 0) { save(keep.pars, file = paste('Model_Results/Params_', dt, '_of_', n.its, '.Rdata', sep = '')) } } 
	}
	
	
}
toc()

if (!is.null(dim(keep.pars))) { save(keep.pars, file = paste('Model_Results/Params_', dt, '_of_', n.its, '.Rdata', sep = '')) }


dim(keep.pars)

cat('\nRR_2L_I \t', 
	 as.character(dt), '\t', dim(keep.pars)[1], ' of ', n.its, '\t\t', round((proc.time()-start.time)[3]/60/60, 2), ' hours', sep = '', 
	file = write.file, append = TRUE)
