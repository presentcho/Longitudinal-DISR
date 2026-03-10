plot.fun <- function(uu, vv, alpha, beta, hs, ind.inside){
  if (nrow(alpha)==length(uu)){
    alpha2 <- alpha; beta2 <- beta
  } else{
    alpha2 <- matrix(NA, length(uu), ncol(alpha))
    alpha2[ind.inside,] <- as.matrix(alpha)
    beta2 <- matrix(NA, length(uu), ncol(beta))
    beta2[ind.inside,] <- as.matrix(beta)
  }
  df <- as.data.frame(cbind(uu, vv, alpha2, beta2))
  colnames(df) <- c('uu', 'vv', 'alpha0', 'alpha1', 'alpha2',
                    'beta0', 'beta1', 'beta2')
  hs2 <- data.frame(V1 = hs$V1 - 0.085, V2 = hs$V2)
  hs2$V1 <- ifelse(hs2$V1 < 0, hs2$V1 + 0.085, hs2$V1)
  for (i in c(3:ncol(df))){
    p <- ggplot(df, aes(uu, vv)) + geom_raster(aes(fill = df[,i]))+ geom_polygon(hs2, mapping = aes(V1, V2), inherit.aes=F, color = 'black', fill = NA, size = 0.7) +
      scale_fill_gradientn(colours = matlab.like(104), limits = c(-4.5, 5), na.value = 'transparent') +
      theme(panel.background = element_blank(),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_blank(),
            legend.title=element_blank(),
            text=element_text(size=20),
            legend.position = 'none',
            axis.ticks.y = element_blank()) + labs(x = '', y = '') +theme(plot.title = element_text(hjust = 0.5)) #+ ggtitle(colnames(df)[i])
    print(p)
  }
}

