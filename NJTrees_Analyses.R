######Code for Cardoso et al. (2024)
######Calculating functional diversity metrics using neighbor-joining trees

library(ape)
library(BAT)

######analyses for Fig 1 and demonstration of how the method works

mat = matrix(c(0, 1, 7, 7, 1, 0, 6, 6, 7, 6, 0, 4, 7, 6, 4, 0), nrow = 4)
colnames(mat) = rownames(mat) = c("A", "B", "C", "D")
mat = as.dist(mat)
tree = nj(mat)
plot(tree)
plot(tree, "u")

comm1 = c(1,1,1,0)
comm2 = c(1,1,0,1)
comm = rbind(comm1, comm2)
colnames(comm) = c("A", "B", "C", "D")

alpha(comm, tree)
originality(comm, tree, relative = FALSE)
uniqueness(comm, tree, relative = FALSE)
contribution(comm, tree, F, relative = FALSE)
dispersion(comm, tree, relative = FALSE)
evenness(comm, tree)
beta(comm, tree)

######analyses for the NJ tree in Fig 3
datJR <- read.csv("kiwi.csv")
datJR = standard(log(datJR))
pca = prcomp(datJR, retx = T)
summary(pca)
tree = tree.build(datJR, distance = "euclidean")
plot(tree, type = "fan", show.tip.label = FALSE, edge.width = 4)

comm1 = rep(1, 105)
comm2 = comm1
comm2[7:11] = 0
rich = alpha(comm = rbind(comm1, comm2), tree)
differ = (rich[1]-rich[2])/rich[1]
differ

tree = tree.build(datJR, func = "upgma", distance = "euclidean")
plot(tree, type = "fan", show.tip.label = TRUE, edge.width = 1)
rich = alpha(comm = rbind(comm1, comm2), tree)
differ = (rich[1]-rich[2])/rich[1]
differ
