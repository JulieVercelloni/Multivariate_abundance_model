genBaseProbs <- function(n, base, similarity, digits = 8) {
  
  n_levels <- length(base)
  x <- rdirichlet(n, similarity * base) 
  
  #--- ensure that each vector of probabilities sums exactly to 1
  
  x <- round(floor(x*1e8)/1e8, digits)   # round the generated probabilities
  xpart <- x[, 1:(n_levels-1)]           # delete the base prob of the final level
  partsum <- apply(xpart, 1, sum)        # add the values of levels 1 to K-1
  x[, n_levels] <- 1 - partsum           # the base prob of the level K = 1 - sum(1:[K-1])
  
  return(x)
}

make_data <- function(nreef,nyear){
  
  d_study <- genData(nreef, id = "reef")
  d_ind <- genCluster(d_study, cLevelVar = "reef", numIndsVar = 500, level1ID = "id") # 500 observations per reef
  d_ind[, z := 0] # set to zero to be used after 
  
  
  # Get the probability of occurrence for each year 
  basestudy <- genBaseProbs(
    n = nyear,
    base =  c(0.3, 0.05, 0.5, 0.05, 0.1), # prob of occurrence 
    similarity = 50 # variation of prob across years (highest number = more similarities between years)
  )
  
  list_ind_y <- list() 
  abundance_matrix <- list()
  
  for(j in 1:nyear) {
    list_ind_y[[j]] <- lapply(
      X = 1:nreef, 
      function(i) {
        b <- basestudy[j,]
        d_x <- d_ind[reef == i]
        genOrdCat(d_x, adjVar = "z", b, catVar = "ordY")
      })
    
    list_ind_y <- list_ind_y[[j]]
    abundance_matrix[[j]] <- t(sapply(list_ind_y, function(x) x[, prop.table(table(ordY))]))
  }
  
  # transform into proportion matrix 
  abundance_table <<- do.call(rbind,abundance_matrix)%>% data.frame()%>%
    rename("HC" = "X1",
           "SC" = "X2",
           "Algae" = "X3",
           "Sponge" = "X4",
           "Sediment" = "X5") %>%
    dplyr::mutate(year = factor(rep(1:nyear, each=nreef)))%>%
    dplyr::mutate(reef = factor(rep(paste('reef_', 1:`nreef`),2))
    )
  return(abundance_table)
}

ordination_plot <- function(model){

  #alpha is an arguement in lvsplot(); relates to scaling
  alpha <- 0.7
  testcov <- model$lv.median %*% t(model$lv.coefs.median[, 2:3])
  
  #singular value decom
  do.svd <- svd(testcov, model$num.lv, model$num.lv)
  choose.lvs <- scale(do.svd$u * matrix(do.svd$d[1:model$num.lv]^alpha, 
                                        nrow = model$n, ncol = 2, byrow = T), center = T, scale = F)%>%data.frame()
  
  
  colnames(choose.lvs) <- c("LV1","LV2")
  choose.lvs$Zone <- X$Zone
  choose.lvs$Year <- X$year
  choose.lvs$reef <- abundance_table$reef
  
  choose.lv.coefs <- scale(do.svd$v * matrix(do.svd$d[1:model$num.lv]^(1 - alpha), nrow = model$p, ncol = 2, byrow = T), center = T, 
                           scale = F)%>%data.frame()
  colnames(choose.lv.coefs)<-c("LV1","LV2")
  choose.lv.coefs$group <- group
  
  
  plot <- ggplot()+geom_point(data=choose.lvs,aes(x=LV1,y=LV2,col=reef),alpha=0.6,size=2.8)+
    geom_text_repel(data=choose.lv.coefs,aes(x=LV1,y=LV2,label=group),
                    box.padding = unit(0.35, "lines"),
                    point.padding = unit(0.3, "lines"))+facet_wrap(~Year)+
    geom_point(data=choose.lv.coefs,aes(x=LV1,y=LV2),size=2.8,shape=17)+
    geom_vline(xintercept = 0,linetype="dashed")+
    geom_hline(yintercept = 0,linetype="dashed")+theme_bw()+
    theme(axis.text.x = element_text(size=13),legend.position="right",
          legend.title = element_blank(), 
          legend.text = element_text(colour = "black", size = 13), 
          panel.grid.minor = element_blank(),plot.title = element_text(hjust=0,vjust=2,size=18, face="bold"),
          axis.text.y = element_text(size=13),axis.title.y=element_text(size=15),axis.title.x=element_text(size=15),
          strip.text=element_text(size=14),strip.background = element_rect(fill="gray98"))
  return(plot)
}