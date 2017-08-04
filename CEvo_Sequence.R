## Analysis of sequences

CreateHistPlot <- function(path = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170630/", 
                           host = 0, epidemics = 1, which = "VA", bin = 2){
  #For each timestep read the according sequence file, plot the file and repeat
  tmax <- as.numeric(readLines(paste(path, "/Epidemics_", epidemics, "/parameters_cluster.dat", sep = ""))[2])
  timesteps <- seq(0, tmax-1, by = bin)
  
  for (t in timesteps)
  {
    inputfile <- paste(path, "/Epidemics_", epidemics,  "/seq/host_", 
                       host, "_seq_time_", t, ".dat", sep = "")
    fin <- read.table(file = inputfile, header = T, sep = "\t")
    image <- array(0, dim = c(200,200,length(timesteps)))
    dataset <- fin[,3]
    lfin <- strsplit(as.character(dataset), "")
    slen <- length(lfin[[1]])
    n <- length(lfin)
    mfin <- matrix(as.numeric(unlist(lfin)), ncol = slen, nrow = n, byrow = T)
    freq <- matrix(numeric(slen*4),nrow = 4, ncol = slen)
    for (i in 1:slen){
      freq[,i] <- table(factor(mfin[,i], levels = c("A","C","G","T")))/n
    }
    #plot
    # nucleotides: red1=A, green2=C, blue3=G, yellow4=T
    filename <- paste(path, "/Epidemics_", epidemics, "/rectplot_", "host_", host, "_time_", t, sep = "", collapse = "")
    filename <- paste(filename, ".png", sep = "", collapse = "")
    png(filename)
    xcoord <- seq(0, slen, length.out= slen + 1) + 0.5
    beginx <- xcoord[-(slen+1)]
    endx <- xcoord[-1]
    plot (xcoord, type = "n", xlab = "Sequence position", ylab = "Frequency", xlim = c(0, slen + 1), ylim = c(0,1))
    mycol <- c(rgb(1,0,0), rgb(0,1,0), rgb(0,0,1), rgb(1,1,0))
    cumfreq <- apply(freq, 2, cumsum)
    for (i in 1:slen){#A = red, C = green, G = blue, T = yellowish
      ybegin<-0
      for (j in 1:4){
        rect(beginx[i], ybegin, endx[i], cumfreq[j,i], col = mycol[j], border = NA)
        ybegin <- cumfreq[j,i]  
      }
      dev.off()
    }
  }
}

SampleAndConvert <- function(path_to_epi = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170703/1000bp_5SNP_06Malus_seq_ImSys005_HighCap/Epidemics_1", 
                             which = c("random","before","after","infection", "transmitted", "tips"), howMuchTime = 30,
                             HowManyStrains = 100)
{
  path_to_seq <- paste(path_to_epi, "/seq/", sep = "")
  path_to_dyn <- paste(path_to_epi, "/dyn/", sep = "")
  
  InfHis <- read.table(file = paste(path_to_dyn, "/Infection_history.dat", sep = ""), header = T)
  n.hosts <- (length(list.files(path_to_dyn))-1)/2
  consequences <- character(n.hosts)
  times <- numeric(n.hosts)
  
  if(which == "random")
  {
    for (h in 0:(n.hosts - 1))
    {
      a <- grep(paste("host_", h , "_", sep = ""), list.files(path_to_seq))
      fin <- read.table(paste(path_to_dyn, "host_", h, ".dat", sep = ""), header = T, sep = "\t")
      tstep <- sample(fin[1,1]:fin[nrow(fin),1], 1)
      p <- fin[which(fin[,1] == tstep), 3]
      seqs <- read.table(paste(path_to_seq, "/host_", h, "_seq_time_", tstep, ".dat", sep = ""), header = T, sep = "\t")
      if (length(which(p == 0)) != 0) 
      {
        seqs <- seqs[-which(p == 0),]
        p <- p[-(which(p == 0))]
      }
      setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
      seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
      seqlist <- sapply(setos, strsplit, split = "")
      for (i in 1:length(seqlist)) seqmat[i,] <- seqlist[[i]]
      s <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
      consequences[h+1] <- paste(s, sep = "", collapse = "")
      times[h+1] <- tstep
    }
  }else if (which == "before")
  {
    to_smp <- rep(T, n.hosts)
    for (i in 1:nrow(InfHis))
    {
      s <- InfHis[i,2]
      if (to_smp[s+1])
      {
        fin <- read.table(paste(path_to_dyn, "host_", s, ".dat", sep = ""), header = T, sep = "\t")
        tmin <- fin[1,1]
        tstep <- max(InfHis[i,1]-howMuchTime, tmin)
        p <- fin[which(fin[,1] == tstep), 3]
        seqs <- read.table(paste(path_to_seq, "/host_", s, "_seq_time_", tstep, ".dat", sep = ""), header = T, sep = "\t")
        if (length(which(p == 0)) != 0) 
        {
          seqs <- seqs[-which(p == 0),]
          p <- p[-(which(p == 0))]
        }
        
        setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
        seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
        seqlist <- sapply(setos, strsplit, split = "")
        for (j in 1:length(seqlist)) seqmat[j,] <- seqlist[[j]]
        se <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
        consequences[s+1] <- paste(se, sep = "", collapse = "")
        to_smp[s+1] <- F
        times[s+1] <- tstep
      }
      
    }
    for (h in (which(to_smp)-1))
    {
      fin <- read.table(paste(path_to_dyn, "host_", h, ".dat", sep = ""), header = T, sep = "\t")
      
      tstep <- fin[nrow(fin),1]
      p <- fin[which(fin[,1] == tstep), 3]
      seqs <- read.table(paste(path_to_seq, "/host_", h, "_seq_time_", tstep, ".dat", sep = ""), header = T, sep = "\t")
      seqs <- seqs[-which(p == 0),]
      p <- p[-which(p == 0)]
      setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
      seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
      seqlist <- sapply(setos, strsplit, split = "")
      for (j in 1:length(seqlist)) seqmat[j,] <- seqlist[[j]]
      se <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
      consequences[h+1] <- paste(se, sep = "", collapse = "")
      times[h+1] <- tstep
    }
  }else if (which == "after")
  {
    to_smp <- rep(T, n.hosts)
    for (i in 1:nrow(InfHis))
    {
      s <- InfHis[i,2]
      if (to_smp[s+1])
      {
        fin <- read.table(paste(path_to_dyn, "host_", s, ".dat", sep = ""), header = T, sep = "\t")
        tmax <- fin[nrow(fin),1]
        tstep <- min(InfHis[i,1]+howMuchTime, tmax)
        p <- fin[which(fin[,1] == tstep), 3]
        seqs <- read.table(paste(path_to_seq, "/host_", s, "_seq_time_", tstep, ".dat", sep = ""), header = T, sep = "\t")
        if (length(which(p == 0)) != 0) 
        {
          seqs <- seqs[-which(p == 0),]
          p <- p[-(which(p == 0))]
        }
        setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
        seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
        seqlist <- sapply(setos, strsplit, split = "")
        for (j in 1:length(seqlist)) seqmat[j,] <- seqlist[[j]]
        se <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
        consequences[s+1] <- paste(se, sep = "", collapse = "")
        to_smp[s+1] <- F
        times[s+1] <- tstep
      }
    }
    for (h in (which(to_smp)-1))
    {
      fin <- read.table(paste(path_to_dyn, "host_", h, ".dat", sep = ""), header = T, sep = "\t")
      
      tstep <- fin[nrow(fin),1]
      p <- fin[which(fin[,1] == tstep), 3]
      seqs <- read.table(paste(path_to_seq, "/host_", h, "_seq_time_", tstep, ".dat", sep = ""), header = T, sep = "\t")
      seqs <- seqs[-which(p == 0),]
      p <- p[-which(p == 0)]
      setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
      seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
      seqlist <- sapply(setos, strsplit, split = "")
      for (j in 1:length(seqlist)) seqmat[j,] <- seqlist[[j]]
      se <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
      consequences[h+1] <- paste(se, sep = "", collapse = "")
      times[h+1] <- tstep
    }
  }else if (which == "infection")
  {
    to_smp <- rep(T, n.hosts)
    for (i in 1:nrow(InfHis))
    {
      s <- InfHis[i,2]
      if (to_smp[s+1])
      {
        fin <- read.table(paste(path_to_dyn, "host_", s, ".dat", sep = ""), header = T, sep = "\t")
        tstep <- InfHis[i,1]
        p <- fin[which(fin[,1] == tstep), 3]
        seqs <- read.table(paste(path_to_seq, "/host_", s, "_seq_time_", tstep, ".dat", sep = ""), header = T, sep = "\t")
        if (length(which(p == 0)) != 0) 
        {
          seqs <- seqs[-which(p == 0),]
          p <- p[-(which(p == 0))]
        }
        setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
        seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
        seqlist <- sapply(setos, strsplit, split = "")
        for (j in 1:length(seqlist)) seqmat[j,] <- seqlist[[j]]
        se <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
        consequences[s+1] <- paste(se, sep = "", collapse = "")
        to_smp[s+1] <- F
        times[s+1] <- tstep
      }
    }
    for (h in (which(to_smp)-1))
    {
      fin <- read.table(paste(path_to_dyn, "host_", h, ".dat", sep = ""), header = T, sep = "\t")
      
      tstep <- fin[nrow(fin),1]
      p <- fin[which(fin[,1] == tstep), 3]
      seqs <- read.table(paste(path_to_seq, "/host_", h, "_seq_time_", tstep, ".dat", sep = ""), header = T, sep = "\t")
      seqs <- seqs[-which(p == 0),]
      p <- p[-which(p == 0)]
      setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
      seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
      seqlist <- sapply(setos, strsplit, split = "")
      for (j in 1:length(seqlist)) seqmat[j,] <- seqlist[[j]]
      se <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
      consequences[h+1] <- paste(se, sep = "", collapse = "")
      times[h+1] <- tstep
    }
  }else if(which == "transmitted")
  {
    for (h in 0:(n.hosts-1))
    {
      tstep <- unlist(strsplit(readLines(paste(path_to_dyn, "host_", h, ".dat", sep = ""))[2], split = "\t"))[1]
      a <- unlist(strsplit(readLines(paste(path_to_seq, "/host_", h, "_seq_time_", tstep, ".dat", sep = ""))[2], split = "\t"))
      consequences[h+1] <- a[3]
      times[h+1] <- tstep
    }
  }else if (which == "tips")
  {
    tmax <- unlist(strsplit(list.files(path_to_seq)[length(list.files(path_to_seq))], split = "[_ \\.]"))
    tmax <- as.numeric(tmax[length(tmax)-1])
    
    for (h in 0:(n.hosts-1))
    {
      fin <- read.table(paste(path_to_dyn, "host_", h, ".dat", sep = ""), header = T, sep = "\t")
      
      p <- fin[which(fin[,1] == tmax), 3]
      seqs <- read.table(paste(path_to_seq, "/host_", h, "_seq_time_", tmax, ".dat", sep = ""), header = T, sep = "\t")
      if (length(which(p == 0)) != 0) 
      {
        seqs <- seqs[-which(p == 0),]
        p <- p[-(which(p == 0))]
      }
      setos <- as.character(sample(seqs[,3], min(nrow(seqs), HowManyStrains), prob = p, replace = F))
      seqmat <- matrix("", ncol = length(strsplit(setos[1], split = "")[[1]]), nrow = length(setos))
      seqlist <- sapply(setos, strsplit, split = "")
      for (j in 1:length(seqlist)) seqmat[j,] <- seqlist[[j]]
      se <- apply(seqmat, 2, function(x) {names(which.max(table(x)))})
      consequences[h+1] <- paste(se, sep = "", collapse = "")
      times[h+1] <- tmax
    }
  }
  if (any(which == c("random", "infection", "tips")))
  {
    target <- paste(path_to_epi, "/Consensus_Sequences_", which, "_", HowManyStrains, "StrainsSampled.fasta", sep = "")
    t <- paste(path_to_epi, "/Consensus_Sequences_sampling_times_", which, ".dat", sep = "")
  }else if (which == "transmitted")
  {
    target <- paste(path_to_epi, "/Consensus_Sequences_", which, "_", "StrainSampled.fasta", sep = "")
    t <- paste(path_to_epi, "/Consensus_Sequences_sampling_times_", which, ".dat", sep = "")
  }else
  {
    target <- paste(path_to_epi, "/Consensus_Sequences_", which, "_", howMuchTime, 
                    "_", HowManyStrains, "StrainsSampled.fasta", sep = "")
    t <- paste(path_to_epi, "/Consensus_Sequences_sampling_times_", which, "_", howMuchTime, ".dat", sep = "")
  }
  for (i in 1:n.hosts)
  {
    cat(file = target, paste(">Host.", i-1, "\n", sep = ""), append = T)
    cat(file = t, paste("Host.", i-1, "\t", times[i], "\n", sep = ""), append = T)
    cat(file = target, consequences[i], append = T)
    cat(file = target, "\n\n", append = T)
  }
}

CompareTrees <- function(path_to_ref = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170711/HighCap_LowEvo/",
                         path_to_samples = "C:/Users/iotn_/Desktop/Toy_Bayesian/20170710_HighCap3/Epi_1/", epidemics = 1, 
                         dist_method = "treeVec")
{
  if(suppressWarnings(!require(treespace))) install.packages("treespace")
  if(suppressWarnings(!require(ape))) install.packages("ape")
  library("ape")
  library("treespace")

  path_to_epi <- paste(path_to_ref, "Epidemics_", epidemics, "/dyn/", sep = "")
  path_out <- paste(path_to_ref, "Epidemics_", epidemics, "/", sep = "")
  inputfile <- paste(path_to_ref, "Epidemics_", epidemics, "/dyn/Infection_history.dat", sep = "")

  if ( !("Infection_history.dat" %in% list.files(path_to_epi)) ) 
  {
    print ("No file found!")
    return()
  }
  tmax <- as.numeric(readLines(paste(path_to_ref, "/Epidemics_", epidemics, "/parameters_cluster.dat", sep = ""))[2])
  fin <- read.table(inputfile, header = T, sep = "\t")
  if (nrow(fin) == 0)
  {
    cat(sprintf("Infection history table cannot have 0 row dimension!\n"))
    return()
  }
  tr.mat <- as.matrix(fin[,2:3] + 1)
  fin_corr <- cbind (fin, "Rev.Time" = tmax - fin[,1])

  elem_mat <- (dim(fin)[1] + 1)
  chain <- list()
  times <- numeric(elem_mat)
  for (i in 1:elem_mat) times[i] <- tmax + 1
  iter <- elem_mat-1
  inf <- fin_corr[,2] + 1
  infd <- fin_corr[,3] + 1
  rev.time <- fin_corr[,5]
  Dict <- character(elem_mat)
  for (i in 0:(elem_mat-1)) Dict[i+1] <- as.character(i)
  Dict <- sort(Dict)
  for (i in 1:elem_mat) chain[[i]] <- match(i-1, Dict)
  #inf <- sapply(inf, match, Dict)
  #infd <- sapply(infd, match, Dict)
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

  NewickTree_ref <- paste(chain[[1]], ";", sep = "")
  NewickTree_ref2 <- paste("tree TREEREF = ", NewickTree_ref, sep = "")
  
  output <- character(0)

  for (d in list.dirs(path_to_samples, recursive = F))
  {
    p <- paste(d, "/", sep = "")
    MCC_Tree <- paste(p, list.files(p)[grep("*MCC.trees", list.files(p))], sep = "")
    Ref_tree <- genSmpldTree(path_to_ref = path_to_ref, path_to_samples = p, epidemics = epidemics)
    if (gsub("//", "/", d) == list.dirs(path_to_samples, recursive = F)[1])
    {
      output <- paste(readLines(MCC_Tree)[1:(length(readLines(MCC_Tree))-2)], collapse = "\n")
    }
    a <- unlist(strsplit(readLines(MCC_Tree)[length(readLines(MCC_Tree))-1], split = " "))
    a[2] <- unlist(strsplit(d, split = "[// /]"))[length(unlist(strsplit(d, split = "[// /]")))]
    TreeForNexus <- paste(a, collapse = " ")
    output <- paste(output, TreeForNexus, Ref_tree, sep = "\n")
  }
  output <- paste(output, NewickTree_ref2, "End;", sep = "\n")
  
  cat(output, file = paste(path_to_samples, "/Multitree_Epi_", epidemics, ".nexus", sep = ""))
  
  par(xpd = FALSE)
  write.nexus(read.nexus(file = paste(path_to_samples, "/Multitree_Epi_", epidemics, ".nexus", sep = "")), 
              file = paste(path_to_samples, "/Multitree_parsed_Epi_", epidemics, ".nexus", sep = ""))
  MultTree <- read.nexus(file = paste(path_to_samples, "/Multitree_parsed_Epi_", epidemics, ".nexus", sep = ""))
  pdf(paste(path_to_samples, "/Multitree.pdf", sep = ""))
  for (i in 1:length(MultTree)) plot(MultTree[i], main = names(MultTree[i]))
  dev.off()
  
  TS <- treespace(MultTree, nf = 2, method = dist_method)
  distMat <- as.matrix(TS$D)
  TS_Table <- as.matrix(TS$pco$li)
  x_lim <- c(min(TS_Table[,1])-5, max(TS_Table[,1])+5)
  y_lim <- c(min(TS_Table[,2])-5, max(TS_Table[,2])+5)
  pdf(paste(path_to_samples, "/MDS.pdf", sep = ""))
  plot(TS_Table, las = 1, xlab = "", ylab = "", xlim = x_lim, ylim = y_lim)
  abline(v = 0, lty = 2, col = "gray")
  abline(h = 0, lty = 2, col = "gray")
  m_pos <- matrix(0, ncol=2, nrow = length(MultTree))
  initials <- character(length(MultTree))
  for (i in 1:length(MultTree)) 
  {
    m_pos[i,] <- apply(TS_Table[-i,], 2, mean)
    m_pos[i,] <- m_pos[i,]/sqrt(m_pos[i,]%*%m_pos[i,])
    if( length(unlist(strsplit(rownames(TS_Table), split = "_")[[i]])) == 2)
    {
      initials[i] <- paste(unlist(strsplit(rownames(TS_Table), split = "")[[i]][1:2]), collapse = "")
    }else if (length(unlist(strsplit(rownames(TS_Table), split = "_")[[i]])) == 3 && unlist(strsplit(rownames(TS_Table), split = "_")[[i]])[3] == "ref")
    {
      initials[i] <- paste(unlist(strsplit(rownames(TS_Table), split = "")[[i]][1:2]), collapse = "")
    }
     else 
    {
      initials[i] <- paste(paste(unlist(strsplit(rownames(TS_Table), split = "")[[i]][1:2]), collapse = ""), 
                           unlist(strsplit(rownames(TS_Table), split = "_")[[i]][2]), sep = "")
    }
    if (i == length(MultTree)) initials[i] <- "RefTree"
    if (unlist(strsplit(rownames(TS_Table), split = "_")[[i]])[length(unlist(strsplit(rownames(TS_Table), split = "_")[[i]]))] == "ref")
    {
      initials[i] <- paste(initials[i], "_Ref", sep = "")
    }
    text(x = m_pos[i,1] * 2 + TS_Table[i,1], y = m_pos[i,2] * 2 + TS_Table[i,2], label = initials[i])
  } 
  dev.off()
  pdf(paste(path_to_samples, "/Dist_to_Ref.pdf", sep = ""))
  dists <- as.matrix(distMat)[length(MultTree),]
  i_ord <- factor(initials[order(dists)], levels = initials[order(dists)], ordered = T)
  plot(i_ord, sort(dists), las = 1, xlab = "Sampling method", ylab = "Kendall Colijn Distance")
  dev.off()
}


genSmpldTree <- function(path_to_ref = "D:/Documents/ETH/Master/4Semester/Master_thesis/New_Proj/Output/Euler/20170711/HighCap_LowEvo/",
                         path_to_samples = "C:/Users/iotn_/Desktop/Toy_Bayesian/20170710_HighCap3/Epi_1/Random_100Strains", epidemics = 1, which = "random")
{
  if(suppressWarnings(!require(treespace))) install.packages("treespace")
  if(suppressWarnings(!require(ape))) install.packages("ape")
  library("ape")
  library("treespace")
  
  path_to_epi <- paste(path_to_ref, "Epidemics_", epidemics, "/dyn/", sep = "")
  path_out <- paste(path_to_ref, "Epidemics_", epidemics, "/", sep = "")
  inputfile <- paste(path_to_ref, "Epidemics_", epidemics, "/dyn/Infection_history.dat", sep = "")
  InfHis <- read.table(inputfile, header = T, sep = "\t")
  #tips = sampling times
  tips <- read.table(paste(path_to_samples, list.files(path_to_samples)[grep("times", list.files(path_to_samples))], sep = ""))
  tips[,1] <- 0:(nrow(tips)-1)
  #sorted to be nexus/BEAST compatible
  tips.sorted <- tips[order(as.character(tips[,1])),]
  #index + 1 where is it on the ordered vector? (maps tips to tips.sorted)
  whereplusone <- sapply(1:nrow(tips), function(x) which(x == order(as.character(tips[,1]))))
  tips.datesort <- tips[order(tips[,2]),]
  #loop over tips in tips.datesort and reorder the whole tree from the previous to the others.
  #the idea is to add a new branch at the sampling time, and eventually reorder the subtree.
  #Subtree is reordered recursively by substituting the next branching with the previous:
  #Example: Moving the sampling of tip 5 from a complex subtree from ending time (699) to time 271
  #                                  271    5    18
  # 373    5    18                   373    18    21
  # 427    5    21                   427    21    23
  # 474    5    23                   474    23    30
  # 555    5    30                   555    30    34
  # 576    30    33        ------>   576    30    36
  # 584    5    34                   584    34    39
  # 604    18    35                  604    18    35
  # 630    33    37                  630    33    37
  # 683    5    39
  NewInfHis <- InfHis[,1:3]
  for (tip in tips.datesort[,1])
  {
    if (tip %in% NewInfHis[,2])
    {
      time <- tips.datesort[which(tips.datesort[,1] == tip),2]
      #if there is at least an infection propagated by the tip
      #elaborate on the tip as in the example
      
      #At first substitute the original infection after sampling time with the time of sampling
      #get sampling time
      to.change <- which(NewInfHis[,2] == tip & NewInfHis[,1] >= time)
      changed <- NewInfHis[to.change,]
      new.inf <- changed[1,]
      new.inf[1] <- time
      for (i in 1:(nrow(changed)-1))
      {
        changed[i,2] <- changed[i,3]
        changed[i,3] <- changed[i+1,3]
      }
      changed <- rbind(new.inf, changed[-nrow(changed),])
      NewInfHis[to.change,] <- changed
      NewInfHis <- NewInfHis[order(NewInfHis[,1]),]
    }
    
  }
  
  #Construct Newick tree
  tmax <- max(tips[,2])
  
  elem_mat <- (dim(NewInfHis)[1] + 1)
  chain <- list()
  times <- tips[,2]
  iter <- elem_mat-1
  inf <- NewInfHis[,2] + 1
  infd <- NewInfHis[,3] + 1
  Dict <- character(elem_mat)
  for (i in 0:(elem_mat-1)) Dict[i+1] <- as.character(i)
  Dict <- sort(Dict)
  for (i in 1:elem_mat) chain[[i]] <- match(i-1, Dict)
  is_leaf <- rep(T,elem_mat)
  #inf <- sapply(inf, match, Dict)
  #infd <- sapply(infd, match, Dict)
  
  for (i in iter:1)
  {
    outd <- times[infd[i]] - NewInfHis[i,1]
    if (outd == 0) outd <- 0.01
    out <- times[inf[i]] - NewInfHis[i,1]
    if (out == 0) out <- 0.01
    times[inf[i]] <- NewInfHis[i,1]
    chain[[inf[i]]] <- paste( "(", chain[[inf[i]]], ":", out, ",",
                                   chain[[infd[i]]], ":", outd, ")", sep = "" )
  }
  
  NewickTree_ref <- paste(chain[[1]], ";", sep = "")
  lab <- unlist(strsplit(path_to_samples, split = "[// /]"))[length(unlist(strsplit(path_to_samples, split = "[// /]")))]
  NewickTree_ref2 <- paste("tree ", lab, "_ref = ", NewickTree_ref, sep = "")
  
  return (NewickTree_ref2)
}
