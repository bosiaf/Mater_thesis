rm (list = ls())

p <- "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170613/"
nr.ep <- 544

##20170607_3
#INF_RATE_CONST <- c(4.0e-9, 7.0e-9, 1.0e-8, 3.0e-8, 6.0e-8, 9.0e-8, 1.5e-7, 4.0e-7, 6.0e-7, 9.0e-7, 1.5e-6, 
#                    4.0e-6, 6.0e-6, 7.0e-6, 9.0e-6, 1.0e-5, 1.2e-5, 1.4e-5)
#FITNESS <- c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55)
#FIT_NON_SNP <- c(-0.1, -0.15, -0.2, -0.3, -0.35, -0.4, -0.45, -0.5, -0.55, -0.6, -0.65)

##20170613
INF_RATE_CONST <- c(6.0e-8, 9.0e-8, 1.5e-7, 4.0e-7, 6.0e-7, 9.0e-7, 1.5e-6, 4.0e-6, 6.0e-6, 7.0e-6, 9.0e-6,
                    1.0e-5, 1.2e-5, 1.4e-5, 1.7e-5, 2.0e-5, 2.3e-5)
FITNESS <- c(0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)
FIT_NON_SNP <- c(-0.35, -0.4, -0.45, -0.5, -0.55, -0.6, -0.65, -0.7, -0.75, -0.8, -0.85, -0.9, -0.95, -1.0, -1.05, -1.1)


plot.Parmap(path = p, just.adimsys = F, inf = INF_RATE_CONST, fit = FITNESS, nr.ep, threshold = 1.5e6)
for (e in 1:nr.ep) Parameters.Test(path = p, epidemics = e, threshold = 1.5e6)
for (e in 1:nr.ep) plot.ViralLoad(path = p, epidemics = e, per_vir = F)
for (e in 1:nr.ep) plot.ViralLoad(path = p, epidemics = e, per_vir = T)
for (e in 1:nr.ep) plot.InfTree(path = p, epidemics = e)



plot.ViralLoad <- function(path = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/"
                           , epidemics = 1, per_vir = T){
  
  if(suppressWarnings(!require(assertthat, quietly = T))) install.packages("assertthat")
  library(assertthat)
  #Do it for every host
  cat(sprintf("EPIDEMICS %d:\n", epidemics))
  nr_hosts <- (length(list.files(paste(path, "Epidemics_", epidemics, "/dyn/", sep = ""))) - 1)/2
  cat(sprintf("Found %g hosts\n", nr_hosts))
  for (host in 0:(nr_hosts - 1) )
  {
    if (.Platform$OS.type == "windows")
    {
      path_to_epi <- paste(path, "Epidemics_", epidemics, "/dyn/", sep = "")
      inputfile <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, "_healthy_cells.dat", sep = "")
      inputfile2 <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, ".dat", sep = "")
      
    }else
    {
      path_to_epi <- paste(path, "Epidemics_", epidemics, "/dyn/", sep = "")
      inputfile <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, "_healthy_cells.dat", sep = "")
      inputfile2 <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, ".dat", sep = "")
      
    }
    if ( !(paste("host_", host, "_healthy_cells.dat", sep = "") %in% list.files(path_to_epi)) ) next
    if( !not_empty(readLines(con = inputfile)) )
    {
      cat(sprintf("Host number %d is empty\n", host))
      next
    }
    fin <- read.table(inputfile, header = T, sep = "\t")
    if (nrow(fin) == 0)
    {
      cat(sprintf("Host number %d has 0 row dimensions\n.", host))
      next
    }
    cRP <- colorRampPalette(c("blue", "green", "yellow", "red"), space = "rgb")
    
    time <- fin[,1]
    viremy <- fin[,4]
    nr.hc <- fin[,2]
    tot.cells <- fin[,3]
    nr.ic <- tot.cells - nr.hc
    #other stuff
    
    normviremy <- viremy/1e7
    normtotc <- tot.cells/1e6
    if (per_vir == T){
      
      if( !not_empty(readLines(con = inputfile2)) )
      {
        bla <- cat(sprintf("Host number %d is empty\n", host))
        next
      }
      fin2 <- read.table(inputfile2, header = T, sep = "\t")
      
      plot.nr <- unique(fin2[which(fin2[,3]>max(fin2[,3])/200), 2])
      my_palette <- cRP(length(plot.nr))
      
      par(mar=c(5, 5, 5, 5))
      
      plot(time, normviremy, axes = F, ylim = c(0, max(normviremy)), 
           xlab = "", ylab = "", type = "n", col = "black",)
      lines(time, normviremy, lty = 1, lwd = 2)
      axis(2, ylim=c(0, max(normviremy)), col = "black", lwd = 2, las = 1)
      mtext(2, text = expression(paste("Viral Load / ", 10^7)), line = 3)
      
      legend("topright", legend = "Total Viral Load", 
             bty = "n", lwd = 2, col = "black")
      axis(1, xlim = c(time[1], time(length(time))))
      mtext(1, text = "Time [day]", line = 3)
      i <- 1
      for (strain in plot.nr)
      {
        t <- fin2[which(fin2[,2] == strain), 1]
        v <- fin2[which(fin2[,2] == strain), 3]/1e7
        lines(t, v, col = my_palette[i])
        i <- i + 1
      }
      dev.copy2pdf(file = paste(path, "Epidemics_", epidemics, "/ViralLoadPlot_host_", host, "_pervir.pdf", sep = ""), 
                   height = 7, width = 10)
      
    }else{
      par(mar=c(5, 5, 5, 5))
      
      plot(time, normviremy, axes = F, ylim = c(0, max(normviremy)), 
           xlab = "", ylab = "", type = "n", col = "black",)
      lines(time, normviremy, lty = 1, lwd = 2)
      axis(2, ylim=c(0, max(normviremy)), col = "black", lwd = 2, las = 1)
      mtext(2, text = expression(paste("Viral Load / ", 10^7)), line = 3)
      
      par(new = T)
      plot(time, normtotc, axes = F, ylim = c(0, max(normtotc)),
           xlab = "", ylab = "", type = "n")
      lines(time, normtotc, lty = 1, lwd = 2, col = "darkgrey")
      axis(4, ylim = c(0, max(normtotc)), col = "black", lwd = 2, las = 1)
      mtext(4, text = expression(paste("Cells count / ", 10^6)), line = 3)
      
      legend("topright", legend = c("Viral Load", "Cell count"), 
             bty = "n", lwd = 2, col = c("black", "darkgrey"))
      axis(1, xlim = c(time[1], time(length(time))))
      mtext(1, text = "Time [day]", line = 3)
      dev.copy2pdf(file = paste(path, "Epidemics_", epidemics, "/ViralLoadPlot_host_", host, ".pdf", sep = ""), 
                   height = 7, width = 10)
    }
    
    
  }
  
}

#plot to show the relative frequencies of every viral strain (y<-#of virions, x<-time)
plot.Density <- function(path = "C:\\Users\\Francesco\\Desktop\\ReconstHIV\\francesco\\simulationtools\\data_out\\", 
                         patient = 1, which = "VA", ...){
  t <- 1
  time <- numeric(0)
  wholefile <- numeric(0)
  if (.Platform$OS.type == "windows"){
    inputfile <- paste(path, "patient", patient, "\\p", patient, "_", which, "_time_", t, ".csv", sep = "")
  }else{
    inputfile <- paste(path, "patient", patient, "/p", patient, "_", which, "_time_", t, ".csv", sep = "")
  }
  while(file.exists(inputfile)){
    fin <- read.csv(inputfile, header = T, sep = " ")
    wholefile <- rbind(wholefile, cbind(fin[,2], t))
    time <- append(time, t)
    t <- t + 1
    if (.Platform$OS.type == "windows"){
      inputfile <- paste(path, "patient", patient, "\\p", patient, "_", which, "_time_", t, ".csv", sep = "")
    }else{
      inputfile <- paste(path, "patient", patient, "/p", patient, "_", which, "_time_", t, ".csv", sep = "")
    }
  }
  fac <- wholefile[,2]
  virions <- wholefile[,1]
  my.palette <- colorRampPalette(c("blue", "orange", "red"), space = "Lab")
  par(bty = "l", las = 1)
  smoothScatter(fac, virions/1000, pch=".", nbin = 256, 
                ylab = expression(paste("Virions / ", 10^3)), xlab = "Time",
                colramp = my.palette, ylim = c(0, quantile(virions, 0.975)/1e3), nrpoints = 0)
}


plot.InfTree <- function(path = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170517_2/", epidemics = 1)
{
  if(suppressWarnings(!require(treescape))) install.packages("treescape")
  if(suppressWarnings(!require(igraph))) install.packages("igraph")
  if(suppressWarnings(!require(ape))) install.packages("ape")
  library("treescape")
  library("igraph")
  library("ape")
  
  if (.Platform$OS.type == "windows")
  {
    path_to_epi <- paste(path, "Epidemics_", epidemics, "/dyn/", sep = "")
    path_out <- paste(path, "Epidemics_", epidemics, "/", sep = "")
    inputfile <- paste(path, "Epidemics_", epidemics, "/dyn/Infection_history.dat", sep = "")
    
  }else
  {
    path_to_epi <- paste(path, "Epidemics_", epidemics, "/dyn/", sep = "")
    path_out <- paste(path, "Epidemics_", epidemics, "/", sep = "")
    inputfile <- paste(path, "Epidemics_", epidemics, "/dyn/Infection_history.dat", sep = "")
  }
  if ( !("Infection_history.dat" %in% list.files(path_to_epi)) ) 
  {
    print ("No file found!")
    return()
  }
  tmax <- as.numeric(readLines(paste(path, "/Epidemics_", epidemics, "/parameters_cluster.dat", sep = ""))[2])
  fin <- read.table(inputfile, header = T, sep = "\t")
  if (nrow(fin) == 0)
  {
    cat(sprintf("Infection history table cannot have 0 row dimension!\n"))
    return()
  }
  tr.mat <- as.matrix(fin[,2:3] + 1)
  fin_corr <- t(t(fin) + c(0,1,1,0))
  fin_corr <- cbind (fin_corr, "Rev.Time" = tmax - fin_corr[,1])
  
  treegraph <- graph_from_edgelist(tr.mat)
  par(oma = c(1,0,0,0))

  plot (treegraph)
  elem_mat <- (dim(fin)[1] + 1)
  chain <- list()
  times <- numeric(elem_mat)
  for (i in 1:elem_mat) chain[[i]] <- as.character(i)
  for (i in 1:elem_mat) times[i] <- 401
  iter <- elem_mat-1
  inf <- fin_corr[,2]
  infd <- fin_corr[,3]
  rev.time <- fin_corr[,5]
  while (iter > 0)
  {
    
    if (times[infd[iter]] != tmax + 1)
    {
      outd <- rev.time[iter] - times[infd[iter]]
      if (times[inf[iter]] != tmax + 1)  
      {
        out <- rev.time[iter] - times[inf[iter]]
        chain[[inf[iter]]] <- paste( "(", chain[[inf[iter]]], ":", out, ",",
                                     chain[[infd[iter]]], ":", outd, ")", sep = "" )
        times[inf[iter]] <- rev.time[iter]
      }else
      {
        times[inf[iter]] <- rev.time[iter]
        chain[[inf[iter]]] <- paste( "(", chain[[inf[iter]]], ":", times[inf[iter]], ",", 
                                     chain[[infd[iter]]], ":", outd, ")", sep = "" )
        
      }
      times[infd[iter]] <- rev.time[iter]
    }else
    {
      times[infd[iter]] <- rev.time[iter]
      if (times[inf[iter]] != tmax + 1)  
      {
        out <- rev.time[iter] - times[inf[iter]]
        chain[[inf[iter]]] <- paste( "(", chain[[inf[iter]]], ":", out, ",",
                                     chain[[infd[iter]]], ":", times[infd[iter]], ")", sep = "" )
        times[inf[iter]] <- rev.time[iter]
      }else
      {
        times[inf[iter]] <- rev.time[iter]
        chain[[inf[iter]]] <- paste( "(", chain[[inf[iter]]], ":", times[inf[iter]], ",", 
                                     chain[[infd[iter]]], ":", times[infd[iter]], ")", sep = "" )
        
      }
    }
    
    iter <- iter - 1  
  }
  
  NewickTree <- paste(chain[[1]], ";", sep = "")
  phy <- read.tree(text = NewickTree)
  phy$root.time <- fin_corr[1]
  phy$tip.label <- as.character(as.numeric(phy$tip.label) - 1)
  plot.phylo(phy)
  axisPhylo(backward = F, root.time = phy$root.time)
  
  dev.copy2pdf(file = paste(path_out, "Infection_History.pdf"))
  
}


Parameters.Test <- function(path = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/"
                            , epidemics = 1, threshold = 1.5e6){
  
  bla <- ""
  if(suppressWarnings(!require(caTools, quietly = T))) install.packages("caTools")
  if(suppressWarnings(!require(assertthat, quietly = T))) install.packages("assertthat")
  
  library(caTools)
  library(assertthat)
  #Do it for every host
  bla <- sprintf(paste(bla, sprintf("EPIDEMICS %d:\n", epidemics), sep = ""))
  
  nr_hosts <- (length(list.files(paste(path, "Epidemics_", epidemics, "/dyn/", sep = ""))) - 1)/2
  if(nr_hosts < 1) 
  {
    S <- 100
    return(100)
  }else{
    bla <- sprintf(paste(bla, "Found %g hosts\n", sep = ""), nr_hosts)
    S <- 0
    same <- T
    temp_hosts <- nr_hosts
    
    #check for abrupt end (non surviving infections -> bad sign!)
    for(host in 0:(nr_hosts -1))
    {
      inputfile <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, "_healthy_cells.dat", sep = "")
      path_to_epi <- paste(path, "Epidemics_", epidemics, "/dyn/", sep = "")
      
      if ( !(paste("host_", host, "_healthy_cells.dat", sep = "") %in% list.files(path_to_epi)) ) next
      if (!not_empty(readLines(con = inputfile)))
      {
        bla <- sprintf(paste(bla, "Host number ", host, " is empty.\n", sep = ""))
        next
      }
      fin <- read.table(inputfile, header = T, sep = "\t")
      if (nrow(fin) == 0)
      {
        bla <- sprintf(paste(bla, "Host number ", host, " has 0 row dimension.\n", sep = ""))
        next
      }
      
      if (nrow(fin) <= 3)
      {
        bla <- sprintf(paste(bla, "Host number ", host, " has less equal than 3 row, skipping it.\n", sep = ""))
        next
      }
      
      if (host == 0) 
      {
        time_end <- fin[dim(fin)[1],1]
      }
      else 
      {
        same <- (fin[dim(fin)[1],1] == time_end) && same
      }
    }
    
    if (!same)
    {
      bla <- sprintf(paste(bla, "Different End of Epidemics for different hosts!\n", sep = ""))
      
      cat(bla, file = paste(path, "/Parameters_Score.dat", sep = ""), append = T)
      return(100)
    }
    
    
    for (host in 0:(nr_hosts - 1) )
    {
      if (.Platform$OS.type == "windows")
      {
        path_to_epi <- paste(path, "Epidemics_", epidemics, "/dyn/", sep = "")
        inputfile <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, "_healthy_cells.dat", sep = "")
        inputfile2 <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, ".dat", sep = "")
        
      }else
      {
        path_to_epi <- paste(path, "Epidemics_", epidemics, "/dyn/", sep = "")
        inputfile <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, "_healthy_cells.dat", sep = "")
        inputfile2 <- paste(path, "Epidemics_", epidemics, "/dyn/host_", host, ".dat", sep = "")
        
      }
      #check existance of file
      if ( !(paste("host_", host, "_healthy_cells.dat", sep = "") %in% list.files(path_to_epi)) ) next
      #check nonemptiness of file
      if (!not_empty(readLines(con = inputfile)))
      {
        bla <- sprintf(paste(bla, "Host number ", host, " is empty.\n", sep = ""))
        next
      }
      #read in file
      fin <- read.table(inputfile, header = T, sep = "\t")

      #check that file has data inside
      if (nrow(fin) == 0)
      {
        bla <- sprintf(paste(bla, "Host number ", host, " has 0 row dimension.\n", sep = ""))
        temp_hosts <- temp_hosts -1
        next
      }
      
      #check to have more than 1 data point to integrate
      if (nrow(fin) == 1)
      {
        bla <- sprintf(paste(bla, "Host number ", host, " has just 1 row dimension!\n", sep = ""))
        temp_hosts <- temp_hosts -1
        next
      }
      
      tmax <- as.numeric(readLines(paste(path, "/Epidemics_", epidemics, "/parameters_cluster.dat", sep = ""))[2])
#      if((fin[dim(fin)[1],1] - fin[1,1]) < 45 && (fin[1,1] < tmax - 45) )
#      {
#        S <- S + 1000
#        next
#      }
      time <- fin[,1] - fin[1,1]
      #check that it is at least 50 days long
      if (length(time) < 50)
      {
        temp_hosts <- temp_hosts -1
        next
      }
      viremy <- fin[,4]
      nr.hc <- fin[,2]
      tot.cells <- fin[,3]
      nr.ic <- tot.cells - nr.hc
      viremy[which (viremy < threshold)] <- threshold
      #just consider the days after the primary infection
      S <- S + log(trapz(time[30:length(time)], viremy[30:length(time)]) - (max(time)-time[30]) * threshold + 1)
      cat("S = ")
      cat(S)
      cat("\thost: ")
      cat(host)
      cat("\n")
    }
    bla <- sprintf ( paste(bla, "Normed integral over threshold is %g \n", sep = ""), S/temp_hosts)
    cat(bla, file = paste(path, "/Parameters_Score.dat", sep = ""), append = T)
    return (S/temp_hosts)
  }
  
}


plot.Parmap <- function(path = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170607_3/"
                         , just.adimsys = F, inf = INF_RATE_CONST, fit = FITNESS, nr.ep, threshold = 1.5e6)
{
  if(suppressWarnings(!require(lattice, quietly = T))) install.packages("lattice")
  if(suppressWarnings(!require(graphics, quietly = T))) install.packages("graphics")
  
  
  if(just.adimsys)
  {
    scores <- numeric(nr.ep)
    for (i in 1:(nr.ep))
    {
      path_to_epi <- paste(path, "Epidemics_", i, sep = "")
      scores[i] <- Parameters.Test(path = path, epidemics = i, threshold = threshold)
      cat(paste("Collecting score from file ", path_to_epi, "\n", sep = ""))
      
    }
    g <- expand.grid(fit,inf)
    g$z <- scores
    scores[which(scores > 30)] <- 30
    
    x.scale <- list(at=fit)
    y.scale <- list(log = 2, at=inf)
    
    levelplot(z~Var1*Var2, g, cuts = 30, scales = list(x = x.scale, y = y.scale), 
              xlab = "Fitness Bonus", ylab = "Infection Rate Constant")
    dev.copy2pdf(file = paste(path, "Parameters_heatmap_just_adimsys.pdf", sep = ""))
  }else{
    scores <- numeric(nr.ep/2)
    scores_adimsys <- numeric(nr.ep/2)
    for (i in 1:(nr.ep/2))
    {
      path_to_epi <- paste(path, "Epidemics_", 2*i - 1, sep = "")
      path_to_epi_adimsys <- paste(path, "Epidemics_", 2*i, sep = "")
      cat(paste("Collecting score from file ", path_to_epi, " and ", path_to_epi_adimsys, "\n", sep = ""))
      scores[i] <- Parameters.Test(path = path, epidemics = 2*i - 1, threshold = threshold)
      scores_adimsys[i] <- Parameters.Test(path = path, epidemics = 2*i, threshold = threshold)
      
    }
    scores[which(scores > 30)] <- 30
    scores_adimsys[which(scores_adimsys > 30)] <- 30
    g <- expand.grid(fit,inf)
    g_adimsys <- expand.grid(fit,inf)
    g$z <- scores
    g_adimsys$z <- scores_adimsys
    
    x.scale <- list(at=fit)
    y.scale <- list(log = 2, at=inf)
    
    pdf(file = paste(path, "/Parameters_heatmap_no_adimsys.pdf", sep = ""), width = 7, height = 10)
    print(levelplot(z~Var1*Var2, g, cuts = 50, scales = list(x = x.scale, y = y.scale), 
              xlab = "Fitness Bonus", ylab = "Infection Rate Constant"))

    dev.off()
    
    pdf(file = paste(path, "/Parameters_heatmap_adimsys.pdf", sep = ""), width = 7, height = 10)
    
    print(levelplot(z~Var1*Var2, g_adimsys, cuts = 50, scales = list(x = x.scale, y = y.scale), 
              xlab = "Fitness Bonus", ylab = "Infection Rate Constant"))

    dev.off()
    
  }
  
  
}



