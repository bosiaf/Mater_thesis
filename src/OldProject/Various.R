a <- file.choose()
b <- read.table(a, header = T, sep = "\t")
head(b)
plot(b$Time, b$Total.virions/10000, type = "l", las = 1, 
     xlab = "Time", ylab = expression(paste("Viral load / mm"^-3, sep = "")))
abline(v = 30, lwd = 2, col = "red")
abline(h = 150, lwd = 2, col = "blue")
dev.copy2pdf(file = paste(pa, "Scorebounds.pdf", sep = ""))

             

#########
#########
#Draw sampling schemes
#########
#########
path_to_samples = "C:/Users/iotn_/Desktop/Toy_Bayesian/20170719/"
InfHis <- file.choose()
InfHis <- "D:\\Documents\\ETH\\Master\\4Semester\\Master_thesis\\New_Proj\\Output\\Euler\\20170719\\20170719_LowAndHigh_2\\Epidemics_1\\dyn\\Infection_history.dat"
Infhis <- read.table(InfHis, header =T, sep = "\t")
textinfhis <- cbind(paste("Host", Infhis[,2], sep = "."), 
                    paste("Host", Infhis[,3], sep = "."))

lstdrs <- list.dirs(path_to_samples, recursive = F)
infhis_unique <- textinfhis[match(unique(textinfhis[,1]), textinfhis[,1]),]
len <- length(readLines(paste(lstdrs[2], "/DistMat.dat", sep = "")))
initials <- as.character(unlist(read.table(file = paste(lstdrs[2], "/Labels.dat", sep = ""))))
fini <- factor(initials[-len], levels = initials[-len])
MultTree <- read.nexus(file = paste(lstdrs[2], list.files(lstdrs[2])[grep("Multitree_parsed_Epi_", list.files(lstdrs[2]))], sep = "/"))
Rand <- list.dirs(lstdrs[2], recursive = F)[6]
smpl_time <- read.table(paste(Rand, list.files(Rand)[7], sep = "/"),header = F, sep = "\t")
smpl_time <- cbind(smpl_time, 800 - smpl_time[,2])
smpl_time[,2] <- smpl_time[,2]-60
tree <- MultTree[[17]]
tot <- cbind(tree$edge, tree$edge.length)
tree$root.edge <- 60
plot(tree, root.edge = T)
axisPhylo(side = 1, backward = F, root.time = 60)

nd_hg <- node.height(tree)
nd_dpt <- node.depth.edgelength(tree)
tot <- cbind(tot, nd_hg[which(nd_dpt != 739)], nd_dpt[which(nd_dpt != 739)])
edges <- c(5, 42, 27, 34, 21, 70, 14, 68, 56, 63, 36, 43, 69, 61, 9, 37, 10, 51, 74,
           17, 24, 65, 66, 60, 49, 50, 30, 73, 32, 31, 16, 12, 7, 47, 59, 19, 53, 23) 

edgelabels(pch = 21, edge = edges, date = smpl_time[,3], cex = 0.7, col = "red", bg = "red")
edgelabels(edge = edges, date = smpl_time[,3], adj = c(1, 1.2),
           text = paste("Host", 0:37, sep = "."), frame = "none",  cex = 0.7)

dev.copy2pdf(file = paste(path_to_samples, "/sampled_tree_random.pdf", sep = ""), height = 7, width = 10)

Af_30 <- list.dirs(lstdrs[2], recursive = F)[1]
smpl_time <- read.table(paste(Af_30, list.files(Af_30)[7], sep = "/"),header = F, sep = "\t")
smpl_time <- cbind(smpl_time, 800 - smpl_time[,2])
smpl_time[,2] <- smpl_time[,2]-60
tree <- MultTree[[17]]
tot <- cbind(tree$edge)
tree$root.edge <- 60
plot(tree, root.edge = T)
axisPhylo(side = 1, backward = F, root.time = 60)

nd_hg <- node.height(tree)
nd_dpt <- node.depth.edgelength(tree)
tot <- cbind(tot, nd_dpt[tot[,1]], nd_dpt[tot[,2]])
edges <- c(3, 40, 26, 34, 21, 71, 14, 68, 55, 63, 36, 44, 69, 61, 9, 37, 11, 52, 
           74, 18, 24, 65, 66, 60, 49, 50, 30, 73, 32, 31, 16, 12, 7, 47, 59, 19, 53, 23) 

edgelabels(pch = 21, edge = edges, date = smpl_time[,3], cex = 0.7, col = "red", bg = "red")
edgelabels(edge = edges, date = smpl_time[,3], adj = c(1, 1.2),
           text = paste("Host", 0:37, sep = "."), frame = "none",  cex = 0.7)

dev.copy2pdf(file = paste(path_to_samples, "/sampled_tree_Af_30.pdf", sep = ""), height = 7, width = 10)



Be_30 <- list.dirs(lstdrs[2], recursive = F)[3]
smpl_time <- read.table(paste(Be_30, list.files(Be_30)[7], sep = "/"),header = F, sep = "\t")
smpl_time <- cbind(smpl_time, 800 - smpl_time[,2])
smpl_time[,2] <- smpl_time[,2]-60
tree <- MultTree[[17]]
tot <- cbind(tree$edge)
tree$root.edge <- 60
plot(tree, root.edge = T)
axisPhylo(side = 1, backward = F, root.time = 60)

nd_hg <- node.height(tree)
nd_dpt <- node.depth.edgelength(tree)
tot <- cbind(tot, nd_dpt[tot[,1]], nd_dpt[tot[,2]])
edges <- c(38, 25, 33, 20, 70, 13, 67, 54, 62, 35, 43, 69, 61, 8, 37, 10, 51, 74,
           17, 24, 64, 66, 60, 48, 50, 28, 73, 32, 31, 16, 12, 7, 47, 59, 19, 53, 23) 

edgelabels(pch = 21, edge = edges, date = smpl_time[-1,3], cex = 0.7, col = "red", bg = "red")
edgelabels(edge = edges, date = smpl_time[-1,3], adj = c(1, 1.2),
           text = paste("Host", 1:37, sep = "."), frame = "none",  cex = 0.7)
points(30, 22.92578, pch = 21, col = "red", bg = "red", cex = 0.7)
text(30, 22.92578, labels = "Host.0", adj = c(1,1.2), cex = 0.7)
dev.copy2pdf(file = paste(path_to_samples, "/sampled_tree_Be_30.pdf", sep = ""), height = 7, width = 10)



In <- list.dirs(lstdrs[2], recursive = F)[5]
smpl_time <- read.table(paste(In, list.files(In)[7], sep = "/"),header = F, sep = "\t")
smpl_time <- cbind(smpl_time, 800 - smpl_time[,2])
smpl_time[,2] <- smpl_time[,2]-60
tree <- MultTree[[17]]
tot <- cbind(tree$edge)
tree$root.edge <- 60
plot(tree, root.edge = T)
axisPhylo(side = 1, backward = F, root.time = 60)

nd_hg <- node.height(tree)
nd_dpt <- node.depth.edgelength(tree)
tot <- cbind(tot, nd_dpt[tot[,1]], nd_dpt[tot[,2]])
nodes <- c(39, 58, 52, 56, 50, 74, 47, 73, 67, 71, 57, 62, 5, 6, 45, 8, 46, 66, 11, 
           49, 14, 72, 16, 17, 65, 19, 54, 21, 22, 23, 25, 26, 27, 28, 29, 30, 31, 32) 

nodelabels(pch = 21, node = nodes, cex = 0.7, col = "red", bg = "red")
nodelabels(node = nodes, adj = c(1, 1.2),
           text = paste("Host", 0:37, sep = "."), frame = "none",  cex = 0.7)
dev.copy2pdf(file = paste(path_to_samples, "/sampled_tree_Infection.pdf", sep = ""), height = 7, width = 10)


tree$root.edge <- 60
plot(tree, root.edge = T)
axisPhylo(side = 1, backward = F, root.time = 60)
tiplabels(pch = 21, col = "red", bg = "red", cex = 0.7)
dev.copy2pdf(file = paste(path_to_samples, "/sampled_tree_Tips.pdf", sep = ""), height = 7, width = 10)
