
library(RColorBrewer)

NumberCompartmts <- 10 # as defined by the domain image 

theta.fn <- function(x, VGP){
  # van Genuchten parameters
  theta_s <- VGP[1]
  theta_r <- VGP[2]
  alpha <- VGP[3]
  n <- VGP[4]
  m <- 1-1/n
  # van Genuchten equation
  theta_r + (theta_s - theta_r)*(1 + (alpha*x)^n)^-m
}

K.fn <- function(x, VGP){
  # Mualem-van Genuchten parameters
  theta_s <- VGP[1]
  theta_r <- VGP[2]
  alpha <- VGP[3]
  n <- VGP[4]
  m <- 1-1/n
  K_s <- VGP[5]
  lambda <- VGP[6]
  # effective saturation as function of theta
  S_e <- (theta.fn(x,VGP) - theta_r)/(theta_s - theta_r)
  # Mualem-van Genuchten equation
  K_s*(S_e^lambda*(1-(1-S_e^(1/m))^m)^2)
}

# Measured 20 November 2018 | Rep 2
VGP <- c(0.428523926, 0.009293583, 0.118361388, 1.260107319, 0.238622357, 0.899235338)
 
# give pressure range
h <- c(seq(0.1,1.9,0.1),seq(2,1e5,2))

# generate random conductivities
# set standard deviation

# Coefficient of variation of hydraulic properties
SigmaSq <- 6
Sigma.y <- sqrt(SigmaSq)

# number of conductivities generated
N.Kh <- 10000

# generate log-normal distribution of conductivities
DistRangeKs <- rlnorm(N.Kh, meanlog=log(VGP[5]), sdlog=Sigma.y)
QRangeKs <- quantile(DistRangeKs, seq(0.05,0.95,length.out = NumberCompartmts), type=5)
RangeKs <- sort(QRangeKs,decreasing = T)

DistRangeAlpha <- rlnorm(N.Kh, mean=log(VGP[3]), sd=Sigma.y)
QRangeAlpha <- quantile(DistRangeAlpha, seq(0.05,0.95,length.out = NumberCompartmts), type=5)
RangeAlpha <- sort(QRangeAlpha, decreasing = T)

# geometric mean function
geo.mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# conductivity plot
# setwd("...")
png(paste0('KhSigmaSq',SigmaSq,'_ensemble.png'), height=12, width=12, units='cm', res = 600)
windowsFonts(FontCarlito=windowsFont("Carlito"))
par(mar=c(2.6,3.1,1.4,1),lend=2, cex=1.25,family="FontCarlito")
plot(NULL,log="xy",ylim = c(1e-12,1e0),xlim = c(0.1,1e5),ylab="",xlab="",axes=F,
     main=bquote(sigma['y']^2==.(SigmaSq)), family="FontCarlito")
# axes specifications
lab.at.x <- 10^(seq(-2,6,1))
lab.x <- sapply(seq(-2,6,1),
                function(i) as.expression(bquote(10^.(i))))
axis(side=1,at=lab.at.x, labels=lab.x,padj=-0.4,family="FontCarlito")
at.y <- 10^(seq(-13,2,2))
labels <- sapply(seq(-13,2,2),
                 function(i) as.expression(bquote(10^ .(i))))
axis(2,at=at.y,labels=labels,las=1,hadj=0.75,family="FontCarlito")
mtext("Pressure head [cm]", side = 1, line = 1.6,cex=1.25,family="FontCarlito")
mtext(expression(paste("K(h) [cm s"^"-1","]")), side = 2, line = 2,cex=1.25,family="FontCarlito")
# grid
abline(h=at.y,v=lab.at.x,col="gray80", lty=3)

Kh.set <- NULL
for (i in 1:NumberCompartmts){
  VGP.range <- NULL
  VGP.range[1] <- VGP[1]
  VGP.range[2] <- VGP[2]
  VGP.range[3] <- RangeAlpha[i]
  VGP.range[4] <- VGP[4]
  VGP.range[5] <- RangeKs[i]
  VGP.range[6] <- VGP[6]
  lines(h,K.fn(h,VGP.range),col="gray50",lwd=1.5)

  Kh.x <- K.fn(h,VGP.range)
  Kh.set <- cbind(Kh.set,Kh.x)
}

Kh_geo <- apply(Kh.set, 1, geo.mean)
Kh_arth <- apply(Kh.set, 1, mean)

lines(h,K.fn(h,VGP),col="darkblue", lwd=3, lty=1)
lines(h,Kh_geo,col="red3", lwd=3,lty=6)
lines(h,Kh_arth,col="green3", lwd=3,lty=3)

legend("topright",c(expression('K'['measured']), expression('K'['G']), expression('K'['A'])),
       lwd=3, col = c('darkblue', 'red3','green3'), bty="n",
       lty = c(1,6,3))
dev.off()


#+++++++++++++ To save calculated ensemble of K(h) ++++++++++++++++++++++++++++
RangeK <- data.frame(RangeAlpha,RangeKs)
#setwd("...")
write.table(RangeK, file=paste0("KhEnsembleAlphaKs_SigmaSq",SigmaSq,"_N",NumberCompartmts,".txt"), row.names = F, sep = ";")
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


