# ----------------------------
# --- Weight Optimisation Graphs
# ============================

## Contour Plots ======================================================================
# Grid of data values
varWS <- seq(1500,5500, length.out = 21)
varPW <- seq(0, 30, length.out = 21)
weightoptim <- expand.grid(WS = varWS, PW = varPW)

weightoptim <- weightoptim %>%
  rowwise() %>%
  do(data.frame(
    # Previous WS and PW values
    .,
    # Determine MTOM
    MTOM = ModifiedSecant(
      function(W_dg_SI)
        W_dg_SI - Weight_Estimate(.$WS, .$PW, W_dg_SI, composite = TRUE, iteration = TRUE),
      6000, 0.001,0.01, positive = TRUE
    )
  )) %>%
  do(data.frame(
    ., 
    Weight_Estimate(.$WS, .$PW, .$MTOM, composite = TRUE, iteration = FALSE)[[1]]
  ))
weightoptim <- data.frame(weightoptim)

library(lattice) 
MTOMlattice <- weightoptim %>% 
  select(WS, PW, MTOM) %>% 
  spread(PW, MTOM) 
rownames(MTOMlattice) <- MTOMlattice$WS 
MTOMlattice <- select(MTOMlattice, -WS) 
asdf <- as.matrix(MTOMlattice) 
filled.contour(varWS,varPW,asdf,nlevels=9,col=brewer.pal(9,"YlOrRd")) 


library(MASS)
data.x <- runif(50);
data.v1 <- runif(50);

x <- matrix(rep(0:50/50,51),nrow=51,ncol=51);
y <- t(x);
theta <- 50;
x0 = 0.2; y0 = 0.5;
z <- abs(sqrt(((x-x0)*cos(theta)+(y-y0)*sin(theta))^2 + ((y-y0)*cos(theta)+(x-    x0)*sin(theta))^2)-1);
Palette <- colorRampPalette(c("lightgrey","black"),
                            interpolate="spline" )
Levels <- (0.1*(0:10))
pt.color <- "black";

filled.contour(x=0:50/50,y=0:50/50,z=z,levels=Levels,
               xlim=c(0,1),ylim=c(0,1),
               color.palette=Palette, xlab="X",ylab="V",
               plot.axes={axis(1);axis(2);
                 points(data.x,data.v1,pch=19,col=pt.color, cex=.5);
                 contour(x=0:50/50,y=0:50/50,z=z,levels=Levels,labcex=1.5,
                         col=grey(0.5),lwd=1.5,add=TRUE, labels="  ", method="flattest"
                 );
                 contour(x=0:50/50,y=0:50/50,z=z,levels=Levels,lwd=1.5,labcex=1.5,
                         lty=0,col=grey(0.2),add=TRUE, method="flattest"
                 );
               }
)
