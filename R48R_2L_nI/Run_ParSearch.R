

Sys.sleep(120*9)

rm(list = ls())

tryCatch(library(Rcpp), error = function(e) { install.packages('Rcpp', repos = 'http://cran.us.r-projct.org'); library(Rcpp) } )
tryCatch(library(parallel), error = function(e) { install.packages('parallel', repos = 'http://cran.us.r-projct.org'); library(Rcpp) } )

setwd('~/work/Param_Search/R48R_2L_nI')

source('~/work/Param_Search/Code/Clean_Res.R')
source('~/work/Param_Search/Code/Incidence_checks.R')

dt = Sys.Date()
date()


#
# Set up model configurations (population(s), R classes)
#


group_labels = c('E', 'e')
populations = c(500000, 500000)

num.r = 48


starts = c('S' = 0.6, 'I' = 10/1000000, 'II' = 10/1000000, 'R' = round(0.25 / num.r / 3, 3))
starts[1] = 1 - sum(starts[2:3]*3) - (starts[4]*num.r*3)


file.name = 'Code/ModelCode_R48R_2L_nI.txt'



#
# Flu characteristics
#


# write model file

source('Code/write_model.R') 




#
# Single Run
#

source(file.name) 


out = sirout_groups(beta_H1 = 0.22, beta_B = 0.21, beta_H3 = 0.24, 
	h_H1 = 0, h_B = 0, h_H3 = 0, 
	d_H1 = 0, d_B = 0, d_H3 = 0, 
	nu = 1/(5/2), rho = 1/(400/num.r), 
	sigma12 = 1 - 0.27, sigma13 = 1 - 0.3, sigma23 = 1 - 0.17, 
	eta12 = 1, eta21 = 1,
	N0 = 1000000, tf = 365*20)
names(out) = all.cols

mod.out = as.data.frame(out)	
summary(mod.out$N)




#
# Parameter search
#


keep.pars = NULL

n.par = 10
tot.time = 60*60*24.5
start.time = proc.time(); end.time = proc.time(); cur.time = (end.time - start.time)[3]

it = 0
while (cur.time < tot.time) {

	it = it+1

	it.pars = list()
	it.pars = mclapply(c(1:(n.par*500)), function(x) {
	
	# draw betas
	m.bs = round(runif(1, 0.1, 0.5), 4)
	betas = sort(round(rnorm(3, m.bs, sd = 0.01), 4))
	new.params = betas[c(2, 1, 3)]
	
	# draw duration
	new.params = c(new.params, round(runif(1, 120, 365*6)))
	
	# draw cross immunity
	for (j in 1:3) { new.params = c(new.params, round(runif(1, 0.05, 0.95), 4)) }
	
	# draw mixing
	for (j in 1:2) { new.params = c(new.params, round(runif(1, 0.05, 1), 4)) }

	# draw importation freq
#	new.params = c(new.params, round(runif(1, 130, 365)))
	
	names(new.params) = c('beta_H1', 'beta_B', 'beta_H3', 'duration', 'sigma12', 'sigma13', 'sigma23', 
		'eta12', 'eta21')
		
	
	out = sirout_groups(beta_H1 = new.params['beta_H1'], beta_B = new.params['beta_B'], beta_H3 = new.params['beta_H3'], 
		h_H1 = 0, h_B = 0, h_H3 = 0, 
		d_H1 = 0, d_B = 0, d_H3 = 0, 
		nu = 1/(5/2), rho = 1/(new.params['duration']/num.r), 
		sigma12 = 1 - new.params['sigma12'], sigma13 = 1 - new.params['sigma13'], sigma23 = 1 - new.params['sigma23'], 
		eta12 = new.params['eta12'], eta21 = new.params['eta21'],
		N0 = 1000000, tf = 365*20)
	names(out) = all.cols
	
	mod.out = as.data.frame(out)	
	
	out.inc = get.inc(mod.out, 365*10)
	
	
	checks = check.par.traj(out.inc)

	return(c(new.params, checks, sum(checks)))

	}, mc.cores = n.par)

	for (i in 1:length(it.pars)) {
		if (it.pars[[i]][length(it.pars[[i]])] > 6 || is.null(dim(keep.pars))) {
			keep.pars = rbind(keep.pars, it.pars[[i]])
		}
	}

	end.time = proc.time(); cur.time = (end.time - start.time)[3]
	
	if (!is.null(dim(keep.pars)) && nrow(keep.pars) > 500) { break }
	
}

it

if (!is.null(dim(keep.pars))) { 
	output.name = paste('Params_', Sys.Date(), '_of_', it*length(it.pars), '.Rdata', sep = '')
	while(output.name %in% list.files('Model_Results')) {
		output.name = paste('Params_', Sys.Date(), '_of_', it*length(it.pars), sample(letters[c(1:26)], 1), '.Rdata', sep = '')
	}

	save(keep.pars, file = paste('Model_Results/', output.name,  sep = '')) 
}

date()

dim(keep.pars)

table(keep.pars[ , ncol(keep.pars)])


# merge all parameters


files = list.files('Model_Results')
all.pars = NULL
for (i in grep('Params_', files)) {
	load(paste('Model_Results/', files[i], sep = ''))
	all.pars = rbind(all.pars, keep.pars)
	rm(keep.pars)
}

dim(all.pars)
table(all.pars[ , ncol(all.pars)])

save(all.pars, file = 'Model_Results/all_params.Rdata')


cat('\n', date(), '\t-\t', table(all.pars[ , ncol(all.pars)]), '\n\n', file = 'Model_Results/Count_pars.txt', append = TRUE)



