
# Perturbing one parameter at a time, 10% perturbation

date()

rm(list = ls())

library(parallel)


setwd('~/work/Param_Search')
source('Incidence_checks.R')

R.types = c('RR', 'RRRR', 'RRRRRR', 'R16R', 'R48R')
L.types = c('L', '2L', 'Ll')
I.types = c('nI', 'I')



for (rr in R.types) {
	for (ll in L.types) {
		for (ii in I.types) {
			

setwd(paste('~/work/Param_Search/', rr, '_', ll, '_', ii, sep = ''))

if (rr %in% c('RR', 'RRRR', 'RRRRRR')) {
	num.r = nchar(rr)
} else {
	num.r = as.numeric(as.character(substr(rr, 2, 3))) 
}


if (ll == 'L') {
	group_labels = c('E')
	populations = c(1000000)
} else if (ll == '2L') {
	group_labels = c('E', 'e')
	populations = c(500000, 500000)
} else if (ll == 'Ll') {
	group_labels = c('E', 'e')
	populations = c(800000, 200000)
}

starts = c('S' = 0.6, 'I' = 10/1000000, 'II' = 10/1000000, 'R' = round(0.25 / num.r / 3, 3))
starts[1] = 1 - sum(starts[2:3]*3) - (starts[4]*num.r*3)

mod.code = paste(rr, ll, ii, sep = '_')
file.name = paste('Code/ModelCode', mod.code, '_a.txt', sep = '')

source('Code/write_model.R')
source(file.name)

if (rr %in% c('RR', 'RRRR', 'RRRRRR')) {
	load('Model_Results/All_pars.Rdata')
} else {
	load('Model_Results/all_params.Rdata')
}

all.pars = all.pars[which(all.pars[ , ncol(all.pars)] == 7), ]
if (nrow(all.pars) > 3000) { all.pars = all.pars[1:3000, ] }

ch.10 = rep(0, min(30000, nrow(all.pars)))
in.10 = matrix(0, nrow = min(30000, nrow(all.pars)), ncol = ncol(all.pars)-8)

start.time = proc.time() 
for (p in 1:length(ch.10)) {
	
	new.pars = all.pars[p, -c((ncol(all.pars)-7):ncol(all.pars))]
	
	
	for (ip in 3:length(new.pars)) {
	qq = mclapply(c(1:20), function(x) { 
		pert = round(runif(length(new.pars), 0.9, 1.1), 3); pert[1:3] = sample(pert[1:3], 1)
		if (ip == 3) { 
			pert[4:length(pert)] = 1
		} else {
			pert[-ip] = 1
		}
		mod.pars = new.pars * pert
		mod.pars['sigma12'] = min(mod.pars['sigma12'], 1)
		mod.pars['sigma13'] = min(mod.pars['sigma13'], 1)
		mod.pars['sigma23'] = min(mod.pars['sigma23'], 1)
		if (ll != 'L') {
			mod.pars['eta12'] = min(mod.pars['eta12'], 1)
			mod.pars['eta21'] = min(mod.pars['eta21'], 1)
		}
		if (ii == 'I') { mod.pars['impfreq'] = round(mod.pars['impfreq']) }
	

		if ((ll %in% c('2L', 'Ll')) && (ii == 'nI')) {
			out = sirout_groups(beta_H1 = mod.pars['beta_H1'], beta_B = mod.pars['beta_B'], beta_H3 = mod.pars['beta_H3'], 
				h_H1 = 0, h_B = 0, h_H3 = 0, 
				d_H1 = 0, d_B = 0, d_H3 = 0, 
				nu = 1/(5/2), rho = 1/(mod.pars['duration']/num.r), 
				sigma12 = 1 - mod.pars['sigma12'], sigma13 = 1 - mod.pars['sigma13'], sigma23 = 1 - mod.pars['sigma23'],
				eta12 = mod.pars['eta12'], eta21 = mod.pars['eta21'], 
				N0 = 1000000, tf = 365*20)
		} else if ((ll %in% c('2L', 'Ll')) && (ii == 'I')) { 
			out = sirout_groups(beta_H1 = mod.pars['beta_H1'], beta_B = mod.pars['beta_B'], beta_H3 = mod.pars['beta_H3'], 
				h_H1 = 0, h_B = 0, h_H3 = 0, 
				d_H1 = 0, d_B = 0, d_H3 = 0, 
				nu = 1/(5/2), rho = 1/(mod.pars['duration']/num.r), 
				sigma12 = 1 - mod.pars['sigma12'], sigma13 = 1 - mod.pars['sigma13'], sigma23 = 1 - mod.pars['sigma23'],
				eta12 = mod.pars['eta12'], eta21 = mod.pars['eta21'], 
				N0 = 1000000, tf = 365*20, days_import = mod.pars['impfreq'])
		} else if (!(ll %in% c('2L', 'Ll')) && (ii == 'nI')) { 
			out = sirout_groups(beta_H1 = mod.pars['beta_H1'], beta_B = mod.pars['beta_B'], beta_H3 = mod.pars['beta_H3'], 
				h_H1 = 0, h_B = 0, h_H3 = 0, 
				d_H1 = 0, d_B = 0, d_H3 = 0, 
				nu = 1/(5/2), rho = 1/(mod.pars['duration']/num.r), 
				sigma12 = 1 - mod.pars['sigma12'], sigma13 = 1 - mod.pars['sigma13'], sigma23 = 1 - mod.pars['sigma23'],
				N0 = 1000000, tf = 365*20)
		} else {
			out = sirout_groups(beta_H1 = mod.pars['beta_H1'], beta_B = mod.pars['beta_B'], beta_H3 = mod.pars['beta_H3'], 
				h_H1 = 0, h_B = 0, h_H3 = 0, 
				d_H1 = 0, d_B = 0, d_H3 = 0, 
				nu = 1/(5/2), rho = 1/(mod.pars['duration']/num.r), 
				sigma12 = 1 - mod.pars['sigma12'], sigma13 = 1 - mod.pars['sigma13'], sigma23 = 1 - mod.pars['sigma23'],
				N0 = 1000000, tf = 365*20, days_import = mod.pars['impfreq'])
		}
		
		
		names(out) = all.cols

		mod.out = as.data.frame(out)	
		out.inc = get.inc(mod.out, 365*10)
	
		checks = check.par.traj(out.inc)
		return(sum(checks))
		
	}, mc.cores = 10)
	in.10[p, ip] = sum(unlist(qq) >= 7)
	}
	
	
	
	
}
end.time = proc.time(); cat((end.time-start.time)[3], 'seconds \n') #toc()

colnames(in.10) = paste(colnames(all.pars)[1:ncol(in.10)], '10', sep = '_'); colnames(in.10)[3] = 'beta_10'


new.all.pars = cbind(all.pars, in.10[ , -c(1:2)])
write.csv(new.all.pars, 'Model_Results/Allpars_pert_ind10.csv')



} } }
