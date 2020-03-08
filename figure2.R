#############################################
## Figure 2: Intervention parameterisation ##
#############################################

## Plotting parameters
width  = 33*12/72.27/2
height = width*9/16

d = seq(as.Date("2017-01-01"), as.Date("2017-12-31"), by = "day")
f = numeric(365)
f[format(d,"%b") %in% c("Dec","Jan","Feb","Mar")] = 1
f[format(d,"%b") %in% "Nov"] = seq( 1/31,30/31,+1/31)
f[format(d,"%b") %in% "Apr"] = seq(30/31, 1/31,-1/31)

xat  = which(format(d, "%d") == "01")

postscript(file = "figure2.eps", width = width, height = height,
           horizontal = FALSE, pointsize = 8, paper = "special")

## Plotting parameters
par(ann = FALSE, las = 1, mar = c(2,2.5,0.5,0.5), mgp = c(2.5,1/3,0),
    ps = 8, tcl = -1/3, xaxs = "i")

## Weighting function
plot(f[1:365], type = "n", xaxt = "n", yaxt = "n")

polygon(c( 90,121, 90), c(0,0,1), border = NA, col = "lightgrey")
polygon(c(304,335,335), c(0,0,1), border = NA, col = "lightgrey")
arrows( 95, 0, 116, 0, code = 3, col = "darkgrey", length = 0.05)
arrows(309, 0, 330, 0, code = 3, col = "darkgrey", length = 0.05)
text(105, 0.1, expression(gamma[2]), col = "darkgrey")
text(320, 0.1, expression(gamma[1]), col = "darkgrey")

abline(v = 121, lty = "dashed", col = "darkgrey")
abline(v = 305, lty = "dashed", col = "darkgrey")
text(295, 0.40, expression(alpha), col = "darkgrey")

arrows(310, 0.5, 365, 0.5, code = 1, col = "darkgrey", length = 0.1)
arrows(  1, 0.5, 115, 0.5, code = 2, col = "darkgrey", length = 0.1)
text(31, 0.40, expression(gamma), col = "darkgrey")

lines(f[1:365], col = "black", lwd = 2)

axis(side = 1, at = xat[seq(1,12,2)], labels = month.abb[seq(1,12,2)])
axis(side = 1, at = xat[seq(2,12,2)], labels = month.abb[seq(2,12,2)])
axis(side = 2, at = seq(0,1,0.2), hadj = 1.2)
title(xlab = "Date", line = 1)
title(ylab = expression(paste("Intervention ", lambda[t])), line = 1.5)

dev.off()

