fourway.decomp <- function( dat, sub.index = 1:nrow(acls), trt, mediator, confounders, survdist = "loglogistic") {
  
  m.string <- sprintf("%s ~ %s + %s", mediator, trt, paste(confounders,collapse = " + "))
  model.m <- glm( as.formula( m.string ), family = "binomial", data = dat[sub.index,])
  y.string <- sprintf("Surv(followyear, deceased) ~ %s * %s + %s", trt, mediator, paste(confounders,collapse = " + "))
  model.y <- survreg( as.formula(y.string), data = dat[sub.index,], dist = survdist )
  
  csr <- coef(model.y)
  beta1 <- csr[2]
  beta2 <- csr[3]
  beta3 <- csr[length(csr)]
  beta4 <- csr[-c(1:3,length(csr))]
  
  cmr <- coef(model.m)
  gamma0 <- cmr[1]
  gamma1 <- cmr[2]
  gamma2 <- cmr[3:length(cmr)]
  
  CM <- as.matrix(acls[,confounders]) ## confounder data matrix
  PM0 <- expit( gamma0 + CM %*% gamma2 )
  PM1 <- expit( gamma0 + gamma1 + CM %*% gamma2)
  
  CDE <- -beta1
  INT.ref <- colMeans(-beta3*PM0)
  INT.med <- colMeans(-beta3*( PM1 - PM0 ))
  PIE <- colMeans(-beta2*( PM1 - PM0 ))
  MED <- INT.med + PIE
  INT <- INT.med + INT.ref
  
  TE <- CDE + INT.ref + INT.med + PIE
  
  pCDE <- CDE/TE
  pINT.ref <- INT.ref/TE
  pINT.med <- INT.med/TE
  pPIE <- PIE/TE
  
  pMED <- MED/TE
  pINT <- INT/TE
  
  MCE <- colMeans( log( PM1/PM0 ) ) ## Causal log RR for the mediator
  
  return( c(TE=TE, MCE = MCE, CDE=CDE, INT.ref=INT.ref, INT.med=INT.med, PIE=PIE, MED=MED, INT=INT,
            pCDE=pCDE, pINT.ref = pINT.ref, pINT.med = pINT.med, pPIE = pPIE,
            pMED = pMED, pINT = pINT))
}


plotEff <- function(eff, pval.eff, name.eff, ln, eff.div = eff) {
  
  if(eff != 0) {
    if(eff < 0) { 
      xpoly <- c(0.1,0.1,0.9)
      ypoly <- ln + min(1,abs(eff/eff.div))*poly.height*c(-1,1,0)
    } else { 
      xpoly <- c(0.1, 0.9, 0.9)
      ypoly <- ln + min(1,abs(eff/eff.div))*poly.height*c(0,-1,1)
    }
    
    polycol <- ifelse(pval.eff < 0.05, rgb(1,0,0,0.3), rgb(0.5, 0.5, 0.5, 0.3))
    
    polygon( xpoly, ypoly, border = NA, col = polycol )
    text(1.1, ln, sprintf("%s = %.2f", name.eff, eff), pos = 4)
  }
  
}

plotInt <- function(eff, pval.eff, name.eff, ln, eff.div = eff) {
  
  if(eff != 0) {
    if(eff < 0) { 
      xpoly <- 0.5 + min(1,abs(eff/eff.div))*poly.width*c(-1,1,0)
      ypoly <- ln + c(-2,-2,-0.3)
    } else { 
      xpoly <- 0.5 + min(1,abs(eff/eff.div))*poly.width*c(-1,1,0)
      ypoly <- ln + c(-0.3,-0.3,-2)
    }

    polycol <- ifelse(pval.eff < 0.05, rgb(1,0,0,0.3), rgb(0.5, 0.5, 0.5, 0.3))
    
    polygon( xpoly, ypoly, border = NA, col = polycol )
    text(0.5, ln-2, sprintf("%s = %.3f", name.eff, eff), pos = 1)
  }
  
}

plot.fourway <- function( med.boot, name.mediator) {
  pval <- 1 - pnorm( abs( med.boot$t0/apply(med.boot$t,2,sd) ) )
  M <- data.frame(rbind( med.boot$t0, pval))
  
  colnames(M) <- c("TE", "MCE", "CDE", "INT.REF", "INT.MED", "PIE", "pCDE", "pINT.REF", "pINT.MED", "pPIE", 
                   "pMED", "pINT")
  
  par(mar=c(0,1,3,1))
  
  poly.height <- 0.4
  poly.width <- 0.2
  
  plot(NULL, xlim = c(0,1.5), ylim = c(0.5,9.5), bty="n", xaxt="n", yaxt="n", xlab = "", ylab = "")
  title(sprintf("Mediator = %s", name.mediator))
  
  text(c(0,0,0,0), c(3,7,8,9), "Z", font = 2, cex=1.5)
  text(c(0.5,0.5), c(3,6.2), "M", font = 2, cex = 1.5)
  text(1, 9, "M", font = 2, cex = 1.5)
  text(c(1,1,1), c(3,7,8), "Y", font = 2, cex = 1.5)
  
  abline(h = c(3.5,7.5,8.5), lty = "dashed", col = gray(0.7))
  
  plotEff(M$MCE[1], M$MCE[2], "MCE", 9, eff.div = M$MCE[1])
  plotEff(M$TE[1], M$TE[2], "TE", 8)
  plotEff(M$CDE[1], M$CDE[2], "CDE", 7, eff.div = M$TE[1])
  plotEff(M$PIE[1], M$PIE[2], "PIE", 3, eff.div = M$TE[1])
  
  plotInt(M$INT.REF[1], M$INT.REF[2], "INT.REF", 6.2, eff.div = M$TE[1])
  plotInt(M$INT.MED[1], M$INT.MED[2], "INT.MED", 3, eff.div = M$TE[1])
  
  text( 1, 5, sprintf("%% Interaction = %d%%", round(100*M$pINT[1],0)), cex = 1.2, col = rgb(as.integer(M$pINT[2]<0.05),0,0))
  text( 1, 1.5, sprintf("%% Mediated = %d%%", round(100*M$pMED[1],0)), cex = 1.2, col = rgb(as.integer(M$pMED[2]<0.05),0,0))

}

table.fourway <- function(L, mediators, transf = function(x) x) {
  M <- do.call("rbind", 
               lapply(L, function(B) {
                 allCIs <- t(sapply(1:length(B$t0), 
                                    function(j) boot.ci(B, index=j, type="perc")$perc[4:5] ))
                 
                 CI.orig <- apply( cbind( transf(allCIs[1:8,1]), transf(allCIs[1:8,2]) ), 1, sort )
                 CI.pct <- apply( cbind( 100*allCIs[-c(1:8),1], 100*allCIs[-c(1:8),2] ), 1, sort )
                 
                 return( rbind( c( sprintf("%.2f", transf(B$t0[1:8])),
                            sprintf("%.0f%%", 100*B$t0[-c(1:8)]) ),
                            c( sprintf("(%.2f, %.2f)", CI.orig[1, ], CI.orig[2, ] ),
                            sprintf("(%.0f%%, %.0f%%)", CI.pct[1, ], CI.pct[2, ] ) ) ) ) }) )
  
  
  colnames(M) <- c("TE", "MCE", "CDE", "INT.REF", "INT.MED", "PIE", "MED", "INT", "pCDE", "pINT.REF", "pINT.MED", "pPIE", 
                   "pMED", "pINT")
  rownames(M) <- rep(mediators, each = 2)
  
  return(M)
}

table.med <- function(L, mediators, transf = function(x) x) {
  M <- do.call("rbind", 
               lapply(L, function(B) {
                 allCIs <- t(sapply(1:length(B$t0), 
                                    function(j) boot.ci(B, index=j, type="perc")$perc[4:5] ))
                 
                 return( cbind( B$t0, allCIs )) }))

  effects <- c("TE", "MCE", "CDE", "INT.REF", "INT.MED", "PIE", "MED", "INT", "pCDE", "pINT.REF", "pINT.MED", "pPIE", 
              "pMED", "pINT")
  
  df <- data.frame(estimate = M[,1], ci.lo = M[,2], ci.hi = M[,3])
  df$effect <- rep(effects, 2)
  df$mediator <- rep( mediators, each = length(effects) )
  
  return(df)
}

barchart.fourway <- function(L, mediators, custom.ylim = NULL, leg.show = FALSE) {

  allCIs <- lapply(L, function(B) {
    t(sapply(1:length(B$t0), 
             function(j) boot.ci(B, index=j, type="perc")$perc[4:5] ))[c(1,7,8),] })
  
  
  V <- as.matrix( t( do.call("rbind", lapply(L, function(Li) {
    return( Li$t0[c(1,7,8)] ) } ) )) )
  colnames(V) <- mediators
  rownames(V) <- c("Total Effect", "Mediated Effect", "Interaction Effect")
  
  if(is.null(custom.ylim)) { ylim <- 2*range*(exp(V)) } else { ylim <- exp(custom.ylim) }
  
  pdf(NULL)
  plot(exp(V), type = "n", ylim = ylim)
  axt <- axTicks( 2 )
  dev.off()
  
  plot(V, type = "n", xlim = c(0,(length(mediators)+2)*3), ylim = log(ylim),
       xlab = "",
       ylab = "Relative Rate of Progression to Death",
       xaxt = "n",
       yaxt = "n")
  axis( side = 2, at = log(axt), labels = axt )
  
  abline(h=0)
  abline(v=0.5+c(0:length(mediators))*4, lty = "dotted", col = gray(0.5))
  
  
  barplot( V, beside = TRUE, add = TRUE,
           border = NA,
           col = gray( c(0.5, 0.7, 0.9) ),
           cex.names = 0.9,
           legend.text = leg.show,
           args.legend = list(x="bottomright"),
           yaxt = "n")
  
  bp <- barplot( V, beside = TRUE, add = TRUE,
                   yaxt = "n",
                   plot = FALSE)
  
  sapply( 1:ncol(bp), function(j) {
    segments( bp[,j], allCIs[[j]][,1], bp[,j], allCIs[[j]][,2])
  })
  
}

gg.fourway <- function(res, exposure.name, mediator.names, transf = function(x) { exp(-x) - 1 }, custom.ylim = NULL) {
  
  effs <-unlist(lapply(res, function(B) {
    return( transf( B$t0[c(1,3:6)] ) ) }))
  
  CIs <- do.call("rbind", lapply(res, function(B) {
    transf( t(sapply(1:length(B$t0), 
             function(j) boot.ci(B, index=j, type="perc")$perc[4:5] ))[c(1,3:6),] )}))

  type.names <- c("Total Effect", "Controlled Direct Effect", "Reference Interaction", "Mediated Interaction", "Pure Indirect Effect")

  df <- data.frame( eff = effs, 
                    ci.lo = CIs[,1], 
                    ci.hi = CIs[,2],
                    med = factor( rep(mediator.names, each = length(type.names)), levels = mediator.names ),
                    type = factor( rep(type.names, length(mediator.names)), levels = type.names ),
                    bar.pos = rep( c(1.5,3,4,5,6), length(mediator.names)))
  
  g.out <- list()
  
  j <- 1
  
  for(med.name in mediator.names) {
    sub.df <- subset(df, med == med.name)
    g <- ggplot(data = sub.df, aes(x = bar.pos, y = eff, fill = type))
    g <- g + geom_bar( stat = "identity", position = "dodge", width = 1 ) +
      geom_errorbar( aes(ymin = ci.lo, ymax = ci.hi ), width = 0.3, colour = gray(0.4)) +
      geom_hline( yintercept = 0, size = 1, colour = gray(0.3)) +
      geom_vline( xintercept = 2.25, linetype = "dashed") +
      xlab(NULL) + ylab(NULL) +
      scale_fill_brewer( name = "", palette = "Set2" ) +
      scale_y_continuous( label = function(x) { round(exp(x),1) } ) +
      theme( axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             plot.title = element_text(size = rel(1.2), face="bold"),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             legend.text = element_text(size = rel(1.2)))
    g <- g + ggtitle(sprintf("Exposure: %s\nMediator: %s", 
                    exposure.name,
                    med.name))

    if(!is.null(custom.ylim)) {
      g.out[[j]] <- g + coord_cartesian( ylim = log( custom.ylim ) )
    } else{
      g.out[[j]] <- g
    }
    
    j <- j+1
  }
  
  return(g.out)

}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols),
                     byrow = TRUE)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
