rm (list = ls())

p <- "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170523_3/"

for (e in 1:8) plot.ViralLoad(path = p, epidemics = e, per_vir = F)
for (e in 1:8) plot.ViralLoad(path = p, epidemics = e, per_vir = T)
for (e in 1:8) plot.InfTree(path = p, epidemics = e)


plot.ViralLoad <- function(path = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/"
                           , epidemics = 1, per_vir = T){
  
  #Do it for every host
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
    fin <- read.table(inputfile, header = T, sep = "\t")
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

