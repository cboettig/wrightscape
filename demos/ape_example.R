require(ape)

     ### The example in Phylip 3.5c (originally from Lynch 1991)
     cat("((((Homo:0.21,Pongo:0.21):0.28,",
        "Macaca:0.49):0.13,Ateles:0.62):0.38,Galago:1.00);",
        file = "ex.tre", sep = "\n")
     tree.primates <- read.tree("ex.tre")
     X <- c(4.09434, 3.61092, 2.37024, 2.02815, -1.46968)
     Y <- c(4.74493, 3.33220, 3.36730, 2.89037, 2.30259)
     names(X) <- names(Y) <- c("Homo", "Pongo", "Macaca", "Ateles", "Galago")
     pic.X <- pic(X, tree.primates)
     pic.Y <- pic(Y, tree.primates)
     cor.test(pic.X, pic.Y)

require(wrightscape)
cpu=2; nboot=100

primates <- format_data(tree.primates, data.frame(X,Y))
attach(primates)
bm <- brown(data[1],tree)
ou <- hansen(data[1], tree, regime=regimes, 1, 1)
LR1 <- choose_model(list(bm=bm, ou=ou), nboot, cpu)

bm <- brown(data[2],tree)
ou <- hansen(data[2], tree, regime=regimes, 1, 1)
LR2 <- choose_model(list(bm=bm, ou=ou), nboot, cpu)

par(mfrow=c(2,1))
pretty_plot(LR1[[1]])
pretty_plot(LR2[[1]])



