



cat('\n\nlibrary(Rcpp)\n\ncppFunction( \" \n', file = file.name)

cat('List sirout_groups(double beta_H1, double beta_B, double beta_H3, 
	double h_H1, double h_B, double h_H3, 
	double d_H1, double d_B, double d_H3, 
	double nu, double rho, 
	double sigma12, double sigma13, double sigma23, 
	double eta12, double eta21, 
	double N0, double tf, 
	int days_import) { 
	
	
	double eta11 = 1.0;
	double eta22 = 1.0;
		
	int grp = 0;
	int imptype = 0; 
	int numsubs = 3; \n\n\t // initial conditions, susceptibles \n\tdouble t = 0;\n', 
	file = file.name, append = TRUE)



###############################
#
# MODEL INPUTS:
#	beta_H1, beta_B, beta_H3: transmission paramters
#	h_H1, h_B, h_H3: hospitalization fractions
#	d_H1, d_B, d_H3: death fractions among hospitalized
#	nu: half of the duration of infection (because there are two I compartments)
#	rho: 1/num.r of the duration of infection (divide duration of immunity by num.r)
#	sigmaxy: Cross immunity between subtypes x and y. 1-z refers to a reduction by z
#	etaab: Mixing parameter between populations a and b. 
#	N0: total model population
#	tf = number of days to run the model for 
#	days_import: Days between sporadic importation of cases
#
###############################




#########################
# Initialize compartments
#########################

	
for (i in 1:length(group_labels)) { 
	cat(paste('\tdouble S_', group_labels[i], ' = ', populations[i], ' * ', starts['S'], '; \n', sep = ''), file = file.name, append = TRUE)
}
for (i in 1:length(group_labels)) { 
	for (j in 1:3) {
		cat(paste('\tdouble I', j, '_', group_labels[i], ' = ', populations[i], ' * ', starts['I'], '; \n', sep = ''), file = file.name, append = TRUE)
	}
}
for (i in 1:length(group_labels)) { 
	for (j in 1:3) {
		cat(paste('\tdouble II', j, '_', group_labels[i], ' = ', populations[i], ' * ', starts['II'], '; \n', sep = ''), file = file.name, append = TRUE)
	}
}
for (i in 1:length(group_labels)) { 
	for (j in 1:3) {
		cat(paste('\tdouble H', j, '_', group_labels[i], ' = ', populations[i], ' * 0.00; \n', sep = ''), file = file.name, append = TRUE)
	}
}

for(r in 1:num.r) {
	for (i in 1:length(group_labels)) {
		for (j in 1:3) {
			cat(paste('\tdouble ', 'R', r, 'R', j, '_', group_labels[i], ' = ', populations[i], ' * ', starts['R'], '; \n', sep = ''), 
				file = file.name, append = TRUE)
		}
	}
}


cat(paste('\n\tdouble D = 0.0; \n\tdouble N = ', sum(populations), '; \n\n', sep = ''), file = file.name, append = TRUE)


for (i in 1:length(group_labels)) {
for (j in 1:3) {
	cat(paste('\tdouble J', j, '_', group_labels[i], ' = 0.00; \n', sep = ''), file = file.name, append = TRUE)
}
}



all.cols = c('t')
for(i in 1:length(group_labels)) { all.cols = c(all.cols, paste('S_', group_labels[i], sep = '')) }
types = c('I', 'J', 'II', 'H')
for (r in 1:num.r) { types = c(types, paste('R', r, 'R', sep = '')) }
for (comp in types) {
for (j in 1:3) {
for (i in 1:length(group_labels)) {
	all.cols = c(all.cols, paste(comp, j, '_', group_labels[i], sep = ''))
}
}
}
all.cols = c(all.cols, 'D', 'N')

cat('\n\n', file = file.name, append = TRUE)
for (i in 1:length(all.cols)) {
	cat(paste('\tstd::vector<double> ', all.cols[i], 'a;\n', sep = ''), file = file.name, append = TRUE)
}



#####################################
# Start loop for moving forward daily
#####################################


cat('\n\n\tdo { \n\n\t\t //int it = int(t); \n\n', file = file.name, append = TRUE)
for (i in 1:length(all.cols)) {
	cat(paste('\t\t', all.cols[i], 'a.push_back(', all.cols[i], ');\n', sep = ''), file = file.name, append = TRUE)
}




########################
#
# Case Importation
#
########################




cat('		
		// ---------------------------
		// Sporadic case introductions
		// ---------------------------
		
		int ttt = int(t); 	
				
		if (ttt%days_import == 0) {
			
			int impgroup = grp %', length(group_labels), '+ 1; 
			grp += 1;
			imptype += 1;
			int stype = imptype%numsubs; \n', file = file.name, append = TRUE)
			

for (i in 1:length(group_labels)) {
	cat(paste('\n\t\t\tif (impgroup == ', i, ') { ', '\n\t\t\t\tS_', group_labels[i], ' -= 1; \n\t\t\t\t', 
		'if (stype == 0) { \n\t\t\t\t\tI1_', group_labels[i], '+= 1; \n\t\t\t\t} ', 
		'else if (stype == 1) { \n\t\t\t\t\tI2_', group_labels[i], '+= 1; \n\t\t\t\t} ', 
		'else if (stype == 2) { \n\t\t\t\t\tI3_', group_labels[i], '+= 1; \n\t\t\t\t} \n\t\t\t}\n', 
		sep = ''), file = file.name, append = TRUE)
}
			
			
			
	
cat('			
		}
	
		
', file = file.name, append = TRUE)





########################
#
# Differential Equations
#
########################




cat('		
		// ----------------------
		// Differential equations
		// ----------------------
		
', file = file.name, append = TRUE)



###
# S
###

cat('\n\t\t //dS \n', file = file.name, append = TRUE)
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dS_', group_labels[i], ' = ', sep = '')
	for (j in 1:length(group_labels)) {
		line = paste(line, '- ', 
			'eta', j, i, '*beta_H1*(I1_', group_labels[j], '+II1_', group_labels[j], ')*S_', group_labels[i], '/N - ', 
			'eta', j, i, '*beta_B*(I2_', group_labels[j], '+II2_', group_labels[j], ')*S_', group_labels[i], '/N - ', 
			'eta', j, i, '*beta_H3*(I3_', group_labels[j], '+II3_', group_labels[j], ')*S_', group_labels[i], '/N ', 
			sep = '')
		if (j == length(group_labels)) { 
			line = paste(line, '\n\t\t\t+ rho*(', 'R', num.r, 'R', '1_', group_labels[i], '+', 'R', num.r, 'R', '2_', group_labels[i], 
				'+', 'R', num.r, 'R', '3_', group_labels[i], '); \n', sep = '') 
		} else { 
			line = paste(line, '\n\t\t\t', sep = '') 
		}	
	}
	cat(line, file = file.name, append = TRUE)
}


###
# I
###

cat('\n\t\t //dI \n', file = file.name, append = TRUE)
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dI1_', group_labels[i], ' = ', sep = '')
	for (j in 1:length(group_labels)) {
		if (j > 1) { line = paste(line, ' + ', sep = '') }
		line = paste(line, 
			'eta', j, i, '*beta_H1*(I1_', group_labels[j], '+II1_', group_labels[j], ')*S_', group_labels[i], '/N', 
			sep = '')
	}
	line = paste(line, '\n\t\t\t- h_H1*nu*I1_', group_labels[i], ' - (1-h_H1)*nu*I1_', 
		group_labels[i], '\n\t\t\t', sep = '')
	for (r in 1:num.r) {
		for (j in 1:length(group_labels)) {
			line = paste(line, '+ ', 'eta', j, i, '*beta_H1*(I1_', group_labels[j], '+II1_', group_labels[j], 
				')*(sigma12*', 'R', r, 'R', '2_', group_labels[i], '+sigma13*', 'R', r, 'R', '3_', group_labels[i], ')/N', sep = '')
			if (j == length(group_labels) && r == num.r) {
				line = paste(line, '; \n', sep = '')
			} else {
				line = paste(line, ' \n\t\t\t', sep = '')
			}
		}
	}
	
	cat(line, file = file.name, append = TRUE)
}
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dI2_', group_labels[i], ' = ', sep = '')
	for (j in 1:length(group_labels)) {
		if (j > 1) { line = paste(line, ' + ', sep = '') }
		line = paste(line, 
			'eta', j, i, '*beta_B*(I2_', group_labels[j], '+II2_', group_labels[j], ')*S_', group_labels[i], '/N', 
				sep = '')
	}
	line = paste(line, '\n\t\t\t- h_B*nu*I2_', group_labels[i], ' - (1-h_B)*nu*I2_', group_labels[i], '\n\t\t\t', sep = '')
	for (r in 1:num.r) {
		for (j in 1:length(group_labels)) {
			line = paste(line, '+ ', 'eta', j, i, '*beta_B*(I2_', group_labels[j], '+II2_', group_labels[j], 
				')*(sigma12*', 'R', r, 'R', '1_', group_labels[i], '+sigma23*', 'R', r, 'R', '3_', group_labels[i], ')/N', sep = '')
			if (j == length(group_labels) && r == num.r) {
				line = paste(line, '; \n', sep = '')
			} else {
				line = paste(line, ' \n\t\t\t', sep = '')
			}
		}
	}
	
	cat(line, file = file.name, append = TRUE)
}
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dI3_', group_labels[i], ' = ', sep = '')
	for (j in 1:length(group_labels)) {
		if (j > 1) { line = paste(line, ' + ', sep = '') }
		line = paste(line, 
			'eta', j, i, '*beta_H3*(I3_', group_labels[j], '+II3_', group_labels[j],')*S_', group_labels[i], '/N', 
			sep = '')
	}
	line = paste(line, '\n\t\t\t- h_H3*nu*I3_', group_labels[i], ' - (1-h_H3)*nu*I3_', group_labels[i], '\n\t\t\t', sep = '')
	for (r in 1:num.r) {
		for (j in 1:length(group_labels)) {
			line = paste(line, '+ ', 'eta', j, i, '*beta_H3*(I3_', group_labels[j], '+II3_', group_labels[j], 
				')*(sigma13*','R', r, 'R', '1_', group_labels[i], '+sigma23*', 'R', r, 'R', '2_', group_labels[i], ')/N', sep = '')
			if (j == length(group_labels) && r == num.r) {
				line = paste(line, '; \n', sep = '')
			} else {
				line = paste(line, ' \n\t\t\t', sep = '')
			}
		}
	}
	
	cat(line, file = file.name, append = TRUE)
}



###
# J
###

cat('\n\t\t //dJ \n', file = file.name, append = TRUE)
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dJ1_', group_labels[i], ' = ', sep = '')
	for (j in 1:length(group_labels)) {
		if (j > 1) { line = paste(line, ' + ', sep = '') }
		line = paste(line, 
			'eta', j, i, '*beta_H1*(I1_', group_labels[j], '+II1_', group_labels[j], ')*S_', group_labels[i], '/N', 
			sep = '')
		for (r in 1:num.r) {
			line = paste(line, ' + ', 'eta', j, i, '*beta_H1*(I1_', group_labels[j], '+II1_', group_labels[j], 
				')*(sigma12*', 'R', r, 'R', '2_', group_labels[i], '+sigma13*', 'R', r, 'R', '3_', group_labels[i], ')/N', sep = '')
			if (j == length(group_labels) && r == num.r) {
				line = paste(line, '; \n', sep = '')
			} else {
				line = paste(line, ' \n\t\t\t', sep = '')
			}
		}
	}
	
	cat(line, file = file.name, append = TRUE)
}
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dJ2_', group_labels[i], ' = ', sep = '')
	for (j in 1:length(group_labels)) {
		if (j > 1) { line = paste(line, ' + ', sep = '') }
		line = paste(line, 
			'eta', j, i, '*beta_B*(I2_', group_labels[j], '+II2_', group_labels[j], ')*S_', group_labels[i], '/N', 
			sep = '')
		for (r in 1:num.r) {
			line = paste(line, '+ ', 'eta', j, i, '*beta_B*(I2_', group_labels[j], '+II2_', group_labels[j], 
				')*(sigma12*', 'R', r, 'R', '1_', group_labels[i], '+sigma23*', 'R', r, 'R', '3_', group_labels[i], ')/N', sep = '')
			if (j == length(group_labels) && r == num.r) {
				line = paste(line, '; \n', sep = '')
			} else {
				line = paste(line, ' \n\t\t\t', sep = '')
			}
		}
	}
	
	cat(line, file = file.name, append = TRUE)
}
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dJ3_', group_labels[i], ' = ', sep = '')
	for (j in 1:length(group_labels)) {
		if (j > 1) { line = paste(line, ' + ', sep = '') }
		line = paste(line, 
			'eta', j, i, '*beta_H3*(I3_', group_labels[j], '+II3_', group_labels[j], ')*S_', group_labels[i], '/N', 
			sep = '')
		for (r in 1:num.r) {
			line = paste(line, '+ ', 'eta', j, i, '*beta_H3*(I3_', group_labels[j], '+II3_', group_labels[j], 
				')*(sigma13*', 'R', r, 'R', '1_', group_labels[i], '+sigma23*', 'R', r, 'R', '2_', group_labels[i], ')/N', sep = '')
			if (j == length(group_labels) && r == num.r) {
				line = paste(line, '; \n', sep = '')
			} else {
				line = paste(line, ' \n\t\t\t', sep = '')
			}
		}
	}
	
	cat(line, file = file.name, append = TRUE)
}







####
# II
####

cat('\n\t\t //dII \n', file = file.name, append = TRUE)
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dII1_', group_labels[i], ' = ', sep = '')
	line = paste(line, '\n\t\t\t(1-h_H1)*nu*I1_', group_labels[i], ' - nu*II1_', group_labels[i], '; \n', sep = '')
	
	cat(line, file = file.name, append = TRUE)
}

for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dII2_', group_labels[i], ' = ', sep = '')
	line = paste(line, '\n\t\t\t(1-h_B)*nu*I2_', group_labels[i], ' - nu*II2_', group_labels[i], '; \n', sep = '')
	
	cat(line, file = file.name, append = TRUE)
}

for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble dII3_', group_labels[i], ' = ', sep = '')
	line = paste(line, '\n\t\t\t(1-h_H3)*nu*I3_', group_labels[i], ' - nu*II3_', group_labels[i], '; \n', sep = '')
	
	cat(line, file = file.name, append = TRUE)
}





###
# H
###

cat('\n\t\t //dH \n', file = file.name, append = TRUE)
subtypes = c('H1', 'B', 'H3'); subtype.labels = c(1:3)
for (s in 1:length(subtypes)) {
	for (i in 1:length(group_labels)) {
		cat(paste('\t\tdouble dH', subtype.labels[s], '_', group_labels[i], ' = ', 
			'h_', subtypes[s], '*nu*I', subtype.labels[s], '_', group_labels[i], ' - ', 
			'd_', subtypes[s], '*nu*H', subtype.labels[s], '_', group_labels[i], ' - ',
			'(1-d_', subtypes[s], ')*nu*H', subtype.labels[s], '_', group_labels[i], '; \n', sep = ''), 
			file = file.name, append = TRUE)
	}
}





###
# R
###

for (r in 1:num.r) {
cat('\n\t\t //d', 'R', r, 'R', ' \n', file = file.name, append = TRUE)
for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble d', 'R', r, 'R', '1_', group_labels[i], ' = ', sep = '')
	if (r == 1) {
		line = paste(line, 'nu*II1_', group_labels[i], ' + (1-d_H1)*nu*H1_', group_labels[i], '\n', sep = '')
	} else {
		line = paste(line, 'rho*', 'R', r-1, 'R', '1_', group_labels[i], '\n', sep = '')
	}
	for (j in 1:length(group_labels)) {
		line = paste(line, '\t\t\t- ', 'eta', j, i, '*sigma12*beta_B*(I2_', group_labels[j], '+II2_', group_labels[j], 
			')*', 'R', r, 'R', '1_', group_labels[i], '/N - ', 
			'eta', j, i, '*sigma13*beta_H3*(I3_', group_labels[j], '+II3_', group_labels[j], 
			')*', 'R', r, 'R', '1_', group_labels[i], '/N', sep = '')
		if (j == length(group_labels)) {
			cat(paste(line, '\n\t\t\t- rho*', 'R', r, 'R', '1_', group_labels[i], '; \n', sep = ''), 
				file = file.name, append = TRUE)
		} else {
			line = paste(line, '\n', sep = '')
		}
	}
}

for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble d', 'R', r, 'R', '2_', group_labels[i], ' = ', sep = '')
	if (r == 1) {
		line = paste(line, 'nu*II2_', group_labels[i], ' + (1-d_B)*nu*H2_', group_labels[i], '\n', sep = '')
	} else {
		line = paste(line, 'rho*', 'R', r-1, 'R', '2_', group_labels[i], '\n', sep = '')
	}
	for (j in 1:length(group_labels)) {
		line = paste(line, '\t\t\t- ', 'eta', j, i, '*sigma12*beta_H1*(I1_', group_labels[j], '+II1_', group_labels[j], 
			')*', 'R', r, 'R', '2_', group_labels[i], '/N - ', 
			'eta', j, i, '*sigma23*beta_H3*(I3_', group_labels[j], '+II3_', group_labels[j], 
			')*', 'R', r, 'R', '2_', group_labels[i], '/N', sep = '')
		if (j == length(group_labels)) {
			cat(paste(line, '\n\t\t\t- rho*', 'R', r, 'R', '2_', group_labels[i], '; \n', sep = ''), 
				file = file.name, append = TRUE)
		} else {
			line = paste(line, '\n', sep = '')
		}
	}
}



for (i in 1:length(group_labels)) {
	line = paste('\t\tdouble d', 'R', r, 'R', '3_', group_labels[i], ' = ', sep = '')
	if (r == 1) {
		line = paste(line, 'nu*II3_', group_labels[i], ' + (1-d_H3)*nu*H3_', group_labels[i], '\n', sep = '')
	} else {
		line = paste(line, 'rho*', 'R', r-1, 'R', '3_', group_labels[i], '\n', sep = '')
	}
	for (j in 1:length(group_labels)) {
		line = paste(line, '\t\t\t- ', 'eta', j, i, '*sigma13*beta_H1*(I1_', group_labels[j], '+II1_', group_labels[j], 
			')*', 'R', r, 'R', '3_', group_labels[i], '/N - ', 
			'eta', j, i, '*sigma23*beta_B*(I2_', group_labels[j], '+II2_', group_labels[j], 
			')*', 'R', r, 'R', '3_', group_labels[i], '/N', sep = '')
		if (j == length(group_labels)) {
			cat(paste(line, '\n\t\t\t- rho*', 'R', r, 'R', '3_', group_labels[i], '; \n', sep = ''), 
				file = file.name, append = TRUE)
		} else {
			line = paste(line, '\n', sep = '')
		}
	}
}

}





###
# D
###

cat('\t\t//dD\n', file = file.name, append = TRUE)
line = '\t\tdouble dD = '
for (i in 1:length(group_labels)) {
	if (i > 1) { line = paste(line, ' + \n\t\t\t', sep = '') }
	for (s in 1:length(subtypes)) {	
		line = paste(line, 'd_', subtypes[s], '*nu*H', subtype.labels[s], '_', group_labels[i], sep = '')
		if (s != length(subtypes)) { line = paste(line, '+ ') }
	}
}
cat(paste(line, '; \n', sep = ''), file = file.name, append = TRUE)





###################
# Aging/birth/death
###################



cat('
		
		// -----------------
		// Birth/death/aging
		// -----------------
	
', file = file.name, append = TRUE)





comp.types = NULL
for (i in 2:length(all.cols)) { comp.types = c(comp.types, strsplit(all.cols[i], '_')[[1]][1]) }
comp.types = unique(comp.types); comp.types = comp.types[which(!(comp.types %in% c('N', 'D')))]
comp.types = comp.types[-grep('J', comp.types)]

death.rate.per.1000 = 14
birth.rate.per.1000 = 14

for (i in 1:length(group_labels)) {
	birth.line = paste('\t\tdS_', group_labels[i], ' += 0.0 ', sep = '')

	for (c in 1:length(comp.types)) {
		
		line = paste('\t\td', comp.types[c], '_', group_labels[i], ' += ', sep = '')

		line = paste(line, ' - ', format(death.rate.per.1000/365, nsmall = 2), '*', comp.types[c], '_', 
			group_labels[i], '/1000', sep = '') # death
		birth.line = paste(birth.line, ' + ', format(death.rate.per.1000/365, nsmall = 2), '*', comp.types[c], '_', 
			group_labels[i], '/1000', sep = '') # birth (should be same as death)
		if (c %in% seq(8, length(comp.types), by = 6)) { birth.line = paste(birth.line, '\n\t\t\t') }
		cat(paste(line, ';\n', sep = ''), file = file.name, append = TRUE)
		
	}


	birth.line = paste('\t\tdS_', group_labels[i], ' += ', populations[i] / sum(populations), ' * ', 
		birth.rate.per.1000, ' * N / 1000 / 365', sep = '')
	cat(paste('\n', birth.line, ';\n\n', sep = ''), file = file.name, append = TRUE)
}






cat('\n\n\t\t//Apply derivatives\n', file = file.name, append = TRUE)
for (i in 2:(length(all.cols)-1)) {
	cat(paste('\t\t', all.cols[i], ' += d', all.cols[i], '; \n', sep = ''), file = file.name, append = TRUE)
}


line = '\n\t\tN = '
for (i in 2:(length(all.cols)-1)) {
	if (substr(all.cols[i], 1, 1) != 'J') { 
		if (substr(all.cols[i], 1, 1) != substr(all.cols[i-1], 1, 1)) { line = paste(line, '\n\t\t\t') }
		line = paste(line, all.cols[i], sep = '')
		if (i == length(all.cols)-1) {
			line = paste(line, '; \n', sep = '')
		} else {
			line = paste(line, '+ ')
		}
	}	
}
line = paste(line, '\n\t\tt++; \n\t} while (t < tf); \n\n\n\tList output(0);\n', sep = '')
cat(line, file = file.name, append = TRUE)

line = '\toutput.push_back(ta); \n'
for (i in 2:length(all.cols)) {
#	line = paste(line, '\toutput.push_back(', all.cols[i], 'a, \"', all.cols[i], '\"); \n', sep = '')
	line = paste(line, '\toutput.push_back(', all.cols[i], 'a); \n', sep = '')
}

cat(paste(line, '\n\treturn(output); \n\n} \n\" \n)\n\n', sep = ''), file = file.name, append = TRUE)


