args = commandArgs(trailingOnly = TRUE)

library(HMMcopy)
library(rtracklayer)

#rfile <- import.wig(args[1])
#gfile <- import.wig(args[2])
#mfile <- import.wig(args[3])
plotSegments2 <- function (correctOutput, segmentOutput,
    ...)
{
    if (is.null(segmentOutput$segs)) {
        warning("Processed segments now found, automatically processing")
        segmentOutput$segs <- processSegments(segments$segs,
            space(correctOutput), start(correctOutput), end(correctOutput),
            correctOutput$copy)
    }
    segs <- segmentOutput$segs
    correctOutput$state <- segmentOutput$state
    cols <- stateCols()
    range <- quantile(correctOutput$copy, na.rm = TRUE, prob = c(0.01,
        0.99))
		
	last_x = 0	
	last_z = 0
	x = vector()
	y = vector()
	color = vector()
	last = vector()
    for(chr in paste0("chr",c(1:22,"X","Y"))){
#	for(chr in paste0("chr",c(1))){
		region_xs = last_x + end(correctOutput[chr])
		x = c(x,region_xs)
		y = c(y,correctOutput[chr]$copy)
		last_x = tail(x,1)
		color = c(color,as.numeric(as.character(correctOutput[chr]$state)))
		last = c(last,last_x)	
		segs[segs$chr == chr, ]$start = segs[segs$chr == chr, ]$start + last_z
		segs[segs$chr == chr, ]$end = segs[segs$chr == chr, ]$end + last_z
		last_z = tail(segs[segs$chr == chr, ]$end,1)
	}
    xtick = vector()
    for(i in c(1:length(last))){
    if(i==1){xtick[i]=last[i]/2}
    else{xtick[i] = (last[i] - last[i-1])/2 + last[i-1]}
    }	
    print(summary(x))
    print(summary(y))
    pdf("copy_number.pdf",width = 16, height = 4)	
    
    par(cex.main = 0.5, cex.lab = 1, cex.axis = 1, mar = c(2.5, 1.5,0.5, 0.5), mgp = c(1, 0.5, 0))
    plot(x, y,col=cols[color],xlim = c(0,tail(last,1)), ylim = c(-3,3),ylab="", xaxt = "n", xaxs = "i", yaxs = "i",...)
    axis(side = 1, at = xtick, labels = paste0("chr",c(1:22,"X","Y")), tick = FALSE,las = 2, hadj = 0.9)
    axis(side = 2, at = seq(-3,3,0.5), labels = FALSE)	
    for (k in 1:nrow(segs)) {
       lines(c(segs$start[k], segs$end[k]), rep(segs$median[k], 2), lwd = 1,
            col = "green")
    }
    abline(v=last,lty = 3)
}
normal_reads <- wigsToRangedData(args[1],args[2],args[3])
normal_copy <- correctReadcount(normal_reads)
default_param <- HMMsegment(normal_copy, getparam = TRUE)

longseg_param <- default_param
longseg_param$e <- 0.999999999999999
longseg_param$strength <- 1e30
longseg_segments <- HMMsegment(normal_copy, longseg_param, verbose = FALSE)
#par(mfrow = c(3, 1))
#par(cex.main = 0.5, cex.lab = 1, cex.axis = 1, mar = c(0.1, 0.1,0.1, 0.1), mgp = c(1, 0.5, 0),yaxs = "i")
plotSegments2(normal_copy,longseg_segments,pch=".")

#for(i in c(1:22))
#{
#	chr = paste0("chr",i)
#	par(mfrow = c(2, 1))
#	par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, mar = c(2, 2, 2, 2), mgp = c(1, 0.5, 0))
#	plotSegments(normal_copy, longseg_segments, chr = chr, pch = ".", ylab = "Copy Number", xlab = "Chromosome Position")
#	cols <- stateCols() # 6 default state colours
#	legend("topleft", c("HOMD", "HETD", "NEUT", "GAIN", "AMPL", "HLAMP"),fill = cols, horiz = TRUE, bty = "n", cex = 0.5)
#	
#}
#dev.off()
if(FALSE){
par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, mar = c(2, 1.5, 0, 0), mgp = c(1, 0.5, 0))
plotSegments(normal_copy, longseg_segments, pch = ".",ylab = "Normal Copy Number", xlab = "Chromosome Position")
for(i in 1:nrow(longseg_segments$mus)) {abline(h = longseg_segments$mus[i ,ncol(longseg_segments$mus)], col = cols[i],lwd = 2, lty = 3)}
abline(v = 7.68e7, lwd = 2, lty = 3)
abline(v = 8.02e7, lwd = 2, lty = 3)
}
if(FALSE){
newmu_param <- longseg_param
newmu_param$mu <- c(-0.5, -0.4, -0.15, 0.1, 0.4, 0.7)
newmu_segments <- HMMsegment(normal_copy, newmu_param, verbose = FALSE)


#####plot4
par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, mar = c(2, 1.5, 0, 0), mgp = c(1, 0.5, 0))
plotSegments(normal_copy, newmu_segments, pch = ".",ylab = "Normal Copy Number", xlab = "Chromosome Position")


par(cex.main = 0.5, cex.lab = 0.5, cex.axis = 0.5, mar = c(2, 1.5, 0, 0), mgp = c(1, 0.5, 0))
newmu_param$m <- newmu_param$mu
realmu_segments <- HMMsegment(normal_copy, newmu_param, verbose = FALSE)


#####plot5
plotSegments(normal_copy, realmu_segments, pch = ".",ylab = "Normal Copy Number", xlab = "Chromosome Position")
}
#}

