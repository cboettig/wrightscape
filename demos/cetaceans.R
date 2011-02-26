require(wrightscape)
path = "/home/cboettig/Documents/ucdavis/research/phylotrees/data/artiodactyla/"

tree <- read.nexus(paste(path, "artiodactyla_tree.nex", sep=""))
data <-read.csv(paste(path,"artiodactyla_mass.txt", sep=""),sep="\t")
out <- format_data(tree, data, species_names = data[,2], regimes = data[,1])


bm <- brown(log(out$data[3]), out$tree)
ou1 <- hansen(log(out$data[3]), out$tree, regimes=out$noregimes, 1, 1)
ou2 <- hansen(log(out$data[3]), out$tree, regimes=out$regimes, 1, 1)
ws1 <- wrightscape(log(out$data[3]), out$tree, regime=out$noregimes, 1, 1)
ws2 <- wrightscape(log(out$data[3]), out$tree, regime=out$regimes, (ou2@sqrt.alpha)^2, ou2@sigma, theta=ou2@theta[[1]])

labrid_models <- list(bm = bm, ou1=ou1, ou2=ou2)

LR <- choose_model(labrid_models, 100)
png("cetacean.png",width=1500, height=600) 
par(mfrow=c(1,2))
pretty_plot(LR[[1]], main="support for OU over BM")
pretty_plot(LR[[2]], main="support for 2 peaks over 1")
dev.off()


par_boots <- fast_boot(ws2, nboot = 40)
png("cetacean_parboot.png", width=800, height=400)
plot(par_boots)
legend("topright", levels(out$regimes), lty=c(1,2), lwd=2) 
dev.off()


