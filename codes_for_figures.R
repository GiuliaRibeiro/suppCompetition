library(plyr)
library(dplyr)
library(R2jags)
load("arc1.RData")
load("arc2.RData")
load("arc3.RData")
load("pyx1.RData")
load("pyx2.RData")
load("pyx3.RData")
load("compfit1.RData")
load("compfit2.RData")
load("compfit3.RData")
load("compfit4.RData")
load("compfit5.RData")
load("compfit6.RData")
load("compfit7.RData")

## Functions and codes that generated the figures
## in the paper

## Functions to generate the plots shown in the paper
## For logistics
p1 <- function(zoo.obj, col="blue"){
    zoo.obj %>%
        fortify() %>%
        ggplot(aes(Index, obsm)) +
        geom_point(size=2, colour=col) +
        geom_ribbon(aes(ymin=post.low, ymax=post.up), alpha=0.25, fill=col) +
        scale_y_continuous(name="")+
        scale_x_continuous(name="")+
        theme_bw()
    }

## for competition
## Predicted total abundances
p2 <- function(zoo.obj){
    fortify(zoo.obj) %>%
        ggplot(aes(Index, espA)) +
        geom_line(colour="red", lwd=1) +
        geom_ribbon(aes(ymin=espA.low, ymax=espA.up), alpha=0.25, fill="red") +
        geom_line(aes(Index,espP),colour="blue", lwd=1) +
        geom_ribbon(aes(ymin=espP.low, ymax=espP.up), alpha=0.25, fill="blue") +
        scale_y_continuous(name="Expected density of cells (1/cm2)")+
        scale_x_continuous(name="Time (days)")+
        theme_bw()
}

## Predicted counts and credible intervals (including detection)
p3 <- function(zoo.obj){
    fortify(zoo.obj) %>%
        ggplot(aes(Index, obsA)) +
        geom_point(colour="red", size=2) +
        geom_ribbon(aes(ymin=projA.low, ymax=projA.up), alpha=0.25, fill="red") +
        geom_point(aes(Index,obsP),colour="blue", size=2) +
        geom_ribbon(aes(ymin=projP.low, ymax=projP.up), alpha=0.25, fill="blue") +
        scale_y_continuous(name="")+
        scale_x_continuous(name="")+
        theme_bw()
    }

f1 <- function(obj, index){
    x <- as.data.frame(summary(as.mcmc(obj)[,c("rA","rP","kA","kP","aAP", "aPA")], quantile=c(0.025,0.5,0.975))[[2]])
    names(x) <- c("lwr", "median", "upr")
    x$type <- "both"
    x$experiment <- index
    x$variable <- substr(rownames(x),1,1)
    x$sp <- substr(rownames(x),2,2)
    x
}

f2 <- function(obj, index, sp){
    x <- as.data.frame(summary(as.mcmc(obj)[,c("r","k")], quantile=c(0.025,0.5,0.975))[[2]])
    names(x) <- c("lwr", "median", "upr")
    x$type <- "single"
    x$experiment <- index
    x$variable <- rownames(x)
    x$sp <- sp
    x
    }

tmp1 <- f1(fit1, 1)
lista <- list(fit1,fit2, fit3, fit4, fit5, fit6, fit7)
lista2 <- list(fit.A1,fit.A2, fit.A3, fit.P1,fit.P2, fit.P3)
sps <- rep(c("A","P"), each=3)
exps <- rep(1:3, each=2)
for(i in 2:7)
    tmp1 <- rbind(tmp1, f1(lista[[i]], i))
for(i in 1:6)
    tmp1 <- rbind(tmp1, f2(lista2[[i]], exps[i], sps[i]))


################################################################################
## Figures as in the paper
## ----single species obs x pred, fig.height=12----------------------------
pdf("growth.pdf", width=9, height = 9)
pP1 <- p1(summary.df(fit.P1, data=pyx1)) + ggtitle("P. operculata")
pP2 <- p1(summary.df(fit.P2, data=pyx2))
pP3 <- p1(summary.df(fit.P3, data=pyx3))
## Arcella
pA1 <- p1(summary.df(fit.A1, data=arc1), col="red") + ggtitle("A. intermedia")
pA2 <- p1(summary.df(fit.A2, data=arc2), col="red")
pA3 <- p1(summary.df(fit.A3, data=arc3), col="red")
grid.arrange(pP1,  pA1, pP2, pA2, pP3, pA3,
             bottom=textGrob("Time (days)", gp=gpar(fontsize=18)),
             left=textGrob("Number of cells", rot=90, gp=gpar(fontsize=18)))
dev.off()
################################################################################

## ----competition obs x predicted counts, fig.height=13.5-----------------
pdf("competition.pdf", width=9, height = 12)
pe1 <- p3(summary.df2(fit1, experim1)) + ggtitle("Experiment #1") +
    annotate("text", x = c(8,9.5), y = c(70,5),
             label = c("P.operculata","A. intermedia"),
             color=c("blue","red"), size = 4)
pe2 <- p3(summary.df2(fit2, experim2)) + ggtitle("Experiment #2")
pe3 <- p3(summary.df2(fit3, experim3)) + ggtitle("Experiment #3")
pe4 <- p3(summary.df2(fit4, experim4)) + ggtitle("Experiment #4")
pe5 <- p3(summary.df2(fit5, experim5)) + ggtitle("Experiment #5")
pe6 <- p3(summary.df2(fit6, experim6)) + ggtitle("Experiment #6")
pe7 <- p3(summary.df2(fit7, experim7)) + ggtitle("Experiment #7")
grid.arrange(pe1,pe2,pe3,pe4,pe5,pe6,pe7, ncol=2,
             bottom=textGrob("Time (days)", gp=gpar(fontsize=18)),
             left=textGrob("Number of cells", rot=90, gp=gpar(fontsize=18)))
dev.off()


################################################################################
## The figure with parameters estimates and credibility intervals for all experiments
##pdf("fig5_alt.pdf", width = 10.5)
pdf("fig5_alt.pdf")

##par(mar=c(5,6,4,0), mfrow=c(1,3), las=1)
par(fig=c(0.1,0.4,0,1), mar=c(5,1,4,1), new=TRUE, las=1)
## r
espaco <- 0.1
espaco2 <- espaco*3
################################################################################
## r
## Open frame
plot(1:10, 1:10,
     xlim = range(c(tmp1$upr[tmp1$variable=="r"], tmp1$lwr[tmp1$variable=="r"])),
     type="n", axes=FALSE, xlab="Intrinsic growth ratio (r)", ylab="")
## Competition experiments, Arcella
with(subset(tmp1, sp=="A"&variable=="r"&type=="both"),{
    points(median, (1:7)-espaco, col="red", pch=19)
    segments(lwr, (1:7)- espaco, upr, (1:7)-espaco, col="red")
}
)
## Competition experiments, Pyx
with(subset(tmp1, sp=="P"&variable=="r"&type=="both"),{
    points(median, (1:7)+espaco, col="blue", pch=19)
    segments(lwr, (1:7) + espaco, upr, (1:7)+espaco, col="blue")
}
)
## Single experiments, Arc
with(subset(tmp1, sp=="A"&variable=="r"&type=="single"),{
    points(median, 7.6 +(1:3)*espaco2, col="red", pch=19)
    segments(lwr, 7.6 +(1:3)*espaco2, upr, 7.6 +(1:3)*espaco2, col="red")
}
)
## Single experim, Pyx
with(subset(tmp1, sp=="P"&variable=="r"&type=="single"),{
    points(median, 9.1 +(1:3)*espaco2, col="blue", pch=19.1)
    segments(lwr, 9.1 +(1:3)*espaco2, upr, 9.1+(1:3)*espaco2, col="blue")
}
)
## Axes, scales ...
axis(2, at=c(1:7), labels=paste("#",1:7), tick=FALSE)
axis(1)
box()
abline(h=(1:7)+0.5,
       lty=rep(c(2,1),c(6,1)),
       col=rep(c("gray","black"), c(6,1)))
par(las=0)
mtext("Competition", at=4, side=2, line=2.75, cex=1.5)
mtext("Mono-specific", at=9, side=2, line=2.75, cex=1.5)
################################################################################
## K
## Open frame
#par(mar=c(5,2,4,2), las=1)
par(fig=c(0.4,0.7,0,1), new=TRUE, las=1)
plot(1:10, 1:10,
     xlim = range(c(tmp1$upr[tmp1$variable=="k"], tmp1$lwr[tmp1$variable=="k"])),
     type="n", axes=FALSE, xlab="Carrying capacity (K)", ylab="", log="x")
## Competition experiments, Arcella
with(subset(tmp1, sp=="A"&variable=="k"&type=="both"),{
    points(median, (1:7)-espaco, col="red", pch=19)
    segments(lwr, (1:7)- espaco, upr, (1:7)-espaco, col="red")
}
)
## Competition experiments, Pyx
with(subset(tmp1, sp=="P"&variable=="k"&type=="both"),{
    points(median, (1:7)+espaco, col="blue", pch=19)
    segments(lwr, (1:7) + espaco, upr, (1:7)+espaco, col="blue")
}
)
## Single experiments, Arc
with(subset(tmp1, sp=="A"&variable=="k"&type=="single"),{
    points(median, 7.6 +(1:3)*espaco2, col="red", pch=19)
    segments(lwr, 7.6 +(1:3)*espaco2, upr, 7.6 +(1:3)*espaco2, col="red")
}
)
## Single experim, Pyx
with(subset(tmp1, sp=="P"&variable=="k"&type=="single"),{
    points(median, 9.1 +(1:3)*espaco2, col="blue", pch=19.1)
    segments(lwr, 9.1 +(1:3)*espaco2, upr, 9.1+(1:3)*espaco2, col="blue")
}
)
## Axes, scales ...
axis(1)
box()
abline(h=(1:7)+0.5,
       lty=rep(c(2,1),c(6,1)),
       col=rep(c("gray","black"), c(6,1)))
################################################################################
## coefficients
## Open frame
par(fig=c(0.7,1,0,0.78), las=1, new=TRUE)
box()
par(fig=c(0.7,1,0,1), las=1, new=TRUE)
plot(1:10, 1:10,
     xlim = range(c(tmp1$upr[tmp1$variable=="a"], tmp1$lwr[tmp1$variable=="a"])),
     type="n", axes=FALSE, xlab="Competition coefficients", ylab="", log="x")
## Competition experiments, Arcella
with(subset(tmp1, sp=="A"&variable=="a"&type=="both"),{
    points(median, (1:7)-espaco, col="red", pch=19)
    segments(lwr, (1:7)- espaco, upr, (1:7)-espaco, col="red")
}
)
## Competition experiments, Pyx
with(subset(tmp1, sp=="P"&variable=="a"&type=="both"),{
    points(median, (1:7)+espaco, col="blue", pch=19)
    segments(lwr, (1:7) + espaco, upr, (1:7)+espaco, col="blue")
}
)
## Legend
legend(x = 0.05, y = 10, c("A. intermedia", "P. operculata"),
       pch=19, lty=1, col=c("red", "blue"), bty="n")
## Axes, scales ...
axis(1)
abline(h=(1:6)+0.5,
       lty=2,
       col="gray")

dev.off()
################################################################################
