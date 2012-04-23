# Benchmarking R script
#
#	Stephen Tredger
#

get_times_from_file = function(path) {
	data = read.csv(path, header=TRUE, strip.white=TRUE, stringsAsFactors=FALSE)

	start_time = c(do.call("cbind",data[1]))
	finish_time = c(do.call("cbind",data[2]))
	elapsed_time = (finish_time - start_time)

	times = list(start_time=start_time, finish_time=finish_time, elapsed_time=elapsed_time)
}

path = "/home/stredger/Documents/464/project/CSC-464/phylip-3.69/analysis"
file = "114bp-5sp-nv.txt"

# read in times
times = get_times_from_file(paste(path, file, sep="/"))

# colours for scatterplot points
all_col="red"
#c_col="blue"

# scatterplot
plot(c(1:length(times$elapsed_time)), times$elapsed_time, pch=20, cex=0.5, xlab="run number", 
		ylab="elapsed time", main=paste(file, "scatterplot"), col=all_col)
mean(times$elapsed_time)

# scatterplot with both c and lind... right now looks like crap as all the c values are crunched together...
#plot(lind_times$std_start_time, lind_times$elapsed_time, log="y", pch=20, cex=0.5, xlab="start time", ylab="log elapsed time", main=paste(file_name, "scatterplot"), col=lind_col)
#points(c_times$std_start_time, c_times$elapsed_time, pch=20, cex=0.5, col=c_col)

# 114-5-sc 1.279288e-05
# 114-5-all 1.804522e-02 includes printing of output...
# 114-5-sm 5.03386e-05
# 114-5-ev 1.128742e-05
# 114-5-nv 6.492615e-06
