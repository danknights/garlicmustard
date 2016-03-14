plotGrowth <- function (inputData='../data/0813_Hinf_fc_abx.txt', 
	outputStats='../results/H_influenzae/0813_Hinf_fc_abx_stats.txt', 
	outputPlot = '../results/H_influenzae/0813_Hinf_fc_abx.pdf', 
	bacteria = '/H. influenzae', abx = T, labels=c('Control', 'Ax', 'Am')) {
	raw_data <- read.table(inputData, header = T)
	avg_data <- matrix(data = 0, nrow = 20, ncol = 8)
	stdev <- matrix (data = 0.01, nrow = 20, ncol = 4)
	stats_data <- array(data=0, dim=c(20,4,4))
	for (i in 1:nrow(avg_data)) {
		for (j in 1:8) {
			avg_data[i,j] <- sum(raw_data[i, 1:4+j*8+11])/4
			for (n in 1:length(raw_data[i, 1:4+j*8+11])) {
				if (j%%2==0) {
					stats_data[(i-1) + (n-1)*20 + (j/2-1)*80 + 1] <- stats_data[(i-1) + (n-1)*20 + (j/2-1)*80 + 1] + raw_data[i, 1:4+j*8+11][n]
				} else {
					stats_data[(i-1) + (n-1)*20 + ((j+1)/2-1)*80 + 1] <- stats_data[(i-1) + (n-1)*20 + ((j+1)/2-1)*80 + 1] -raw_data[i, 1:4+j*8+11][n]
				}
			}
		}
		for (k in 1:4) {
			stdev[i,k] <- sd(raw_data[i, 1:4+k*16+3])/sqrt(4) + sd(raw_data[i, 1:4+k*16+11])/sqrt(4)
		}
	}

	# subtracting each negative control to account for differential absorbance of inhibitors

	diff_data <- matrix(data = 0, nrow = 20, ncol = 4)
	for (i in 1:4) {
		diff_data[,i] <- avg_data[,2*i] - avg_data[,2*i-1]
	}

	# setting up areas under graph for statistical analysis

	stats_areas <- matrix (data = 0, nrow = 4, ncol = 4)
	for (i in 1:4) {
		for (j in 1:4) {
			stats_areas[j,i] <- 0.25*(2*Reduce("+",stats_data[0:19 + (j-1)*20 + (i-1)*80 + 1])-stats_data[0 + (j-1)*20 + (i-1)*80 + 1][[1]]-stats_data[19 + (j-1)*20 + (i-1)*80 + 1][[1]])
		}
	}
	stats_results <- file(outputStats)
	p_values <- c()
	area_ratios <- c()
	for (i in 2:ncol(stats_areas)) {
		p_values <- c(p_values,toString(t.test(stats_areas[,1], stats_areas[,i])$p.value))
		area_ratios <- c(area_ratios,1-mean(stats_areas[,i])/mean(stats_areas[,1]))
	}
	writeLines(c(p_values, area_ratios), stats_results)
	close(stats_results)

	# making growth curve graphs

	pdf(outputPlot,width=5,height=5)
	library("RColorBrewer")
	if (abx) {
		colors <- brewer.pal(4, 'Set1')[c(1,3,2,4)]
		lineSet <- c(1, 2, 4)
	} else {
		colors <- brewer.pal(4, 'Set1')
		lineSet <- c(1,2,3,4)
	}
	times <- 0.5*0:19
	plot(times, diff_data[,1], xlab='Time Point (hrs)', ylab = expression('OD'['600']), xlim = c(0,10), ylim=c(0,1.35), type='n', main=bacteria)
	legend ('topleft', labels, fill=colors[lineSet])
	for (i in lineSet) {
		lines(times, diff_data[,i], col = colors[i], type='l')
		segments(times, diff_data[,i]-stdev[,i],times, diff_data[,i]+stdev[,i])
		segments(times-0.05,diff_data[,i]-stdev[,i],times+0.05,diff_data[,i]-stdev[,i])
		segments(times-0.05,diff_data[,i]+stdev[,i],times+0.05,diff_data[,i]+stdev[,i])
	}
	dev.off()
}

print("run plotGrowth(inputFileName, outputStatsFileName, outputPlotFileNmae, bacerialSpecies, evaluating Antibiotics (T or F), Legend labels)")