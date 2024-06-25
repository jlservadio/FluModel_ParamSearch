


check.par.traj = function(inc.dat) {
	
	checks = rep(NA, 7)

	#check 1: sd of daily incidence for each subtype is > 100
	checks[1] = ifelse(sd(inc.dat$inc1) > 100 && sd(inc.dat$inc2) > 100 && sd(inc.dat$inc3) > 100, 1, 0)
	
	#check 2: all pairwise correlation coefficients have absolute value less than 0.7
	cc = cor(inc.dat); diag(cc) = NA
	checks[2] = ifelse(max(abs(cc), na.rm = TRUE) < 0.7, 1, 0)
	
	# check 3: Among the 3 subtype maxima, the lowest max is no less than 80% of the highest max
	mxs = apply(inc.dat, 2, max)
	checks[3] = ifelse(min(mxs) >= 0.8*max(mxs), 1, 0)
	
	# check 4: Maximum values of each subtype in each 5-year half occur at least 100 days apart
	tm.1 = inc.dat[1:round(nrow(inc.dat)/2), ]
	tm.2 = inc.dat[(round(nrow(inc.dat)/2)+1):nrow(inc.dat), ]
	
	mx.t.1 = apply(tm.1, 2, function(x) { which.max(x) } )
	mx.t.2 = apply(tm.2, 2, function(x) { which.max(x) } )
	
	checks[4] = ifelse(min(diff(sort(mx.t.1))) > 100 && min(diff(sort(mx.t.2))) > 100, 1, 0)
	
	# check 5: Resembling data: quantiles relative to the median (weekly data)
	weeks = rep(1, 7); while(length(weeks)<nrow(inc.dat)) { weeks = c(weeks, rep(max(weeks)+1, 7)) }; weeks = weeks[1:nrow(inc.dat)]
	tot.cases = aggregate(apply(inc.dat, 1, sum), by = list('week' = weeks), FUN = sum)$x
	qs = quantile(tot.cases, c(0.1, 0.25, 0.75, 0.9)) / median(tot.cases)
	qs.c = qs / c(0.335, 0.608, 1.601, 2.663)	
	checks[5] = prod(qs.c > 2/3 & qs.c < 4/3)
	if (median(tot.cases) == 0) { checks[5] = 0 }

	
	# check 6: All annual attack rates between 15 and 35%
	yrs = rep(1, 365); while(length(yrs) < nrow(inc.dat)) { yrs = c(yrs, rep(max(yrs)+1, 365)) }; yrs = yrs[1:nrow(inc.dat)]
	ars0 = aggregate(inc.dat, by = list('year' = yrs), FUN = sum)[ , -1] / 1000000
	ars = apply(ars0, 1, sum)
	checks[6] = ifelse(any(ars < 0.15) || any(ars > 0.35), 0, 1)
	
	
	# check 7: Annual attack rates do not monotonically increase/decrease
	checks[7] = ifelse(any(diff(ars) < 0) && any(diff(ars) > 0), 1, 0)
	
	return(checks)	
	
	
}

