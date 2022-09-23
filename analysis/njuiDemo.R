## This file reproduces the analysis of longitudinal mood compositional data for unemployed workers in New Jersey.
devtools::load_all()

library(fdapace)
library(ggplot2)
library(ggridges)
library(reshape2)
theme_set(theme_bw(base_size=14))
library(plyr)
library(dplyr)

## Tuning parameters
KUse <- 30
bwMu <- 'GCV'
bwCov <- 'GCV'
kern <- 'epan'
npoly <- 1
useDays <- 84
## End Tuning parameters

ToutRange <- c(0, useDays)
pts <- seq(min(ToutRange), useDays, length.out=51)
varnames <- c('bad', 'low', 'mild', 'good')

figpath <- '../figures/NJUI/demo'
dir.create(figpath, showWarnings=FALSE, recursive=TRUE)

## Load longitudinal compositional data (datUse) and baseline income data (income0)
# The columns of datUse are the subject ID (caseid), proportions of time the subject spent in bad, low, mild, and good moods (bad, low, mild, and good), and the day since the start of study when the subject responded to the survey (day)
# The columns of income0 are the subject ID (caseid) and the annual household income before 2008
load('njuiDat.RData')

## Analysis
(n <- length(unique(datUse$caseid)))
ni <- table(datUse$caseid)
table(ni)
allID <- unique(datUse$caseid)

# Pool data and then smooth 
tList <- dlply(datUse, 'caseid', function(a) a$day)
yList <- dlply(datUse, 'caseid', function(a) t(sqrt(as.matrix(a[varnames]))))

# Requires around 30 minutes to finish on a 2018 macbook pro
print(system.time({
  resSp <- RFPCA(yList, tList, list(mfdName='sphere', userBwMu=bwMu, userBwCov=bwCov, nRegGrid=51, ToutRange=ToutRange, maxK=30, npoly=npoly))
}))

flip <- c(mean(resSp$phi[, 3:4, 1]) < 0, # phi_1 corresponds to positive moods
          mean(resSp$phi[, c(1, 4), 2]) > 0, # phi_2 corresponds to moods stability
          mean(resSp$phi[, c(2, 4), 3]) < 0, # phi_3 corresponds to shift to positive moods
          mean(diff(resSp$phi[, 4, 4])) < 0) # phi_4 corresponds to more positive moods over time
flipInd <- which(flip)
if (length(flipInd) > 0) {
  resSp <- flipComponents(resSp, flipInd)
}

xiDat <- data.frame(caseid=as.numeric(as.character(rownames(resSp$xi))), resSp$xi[, 1:4], row.names=NULL) %>%
  inner_join(income0)

plotDat <- xiDat %>% 
  melt(measure.vars=grep('^xi', colnames(xiDat), value=TRUE))

pd2 <- plotDat %>% filter(variable == 'xi2' & !is.na(income)) %>% mutate(income = factor(income, rev(levels(income))))
pXi2 <- ggplot(pd2, aes(x=value, y=income, fill=income)) + geom_point(data=pd2 %>% group_by(income) %>% summarize(value=mean(value)), aes(color=income)) + geom_density_ridges(scale=7, rel_min_height=0.002, alpha=0.15) + xlab(expression(xi[2])) + ylab('Annual Household Income') + theme_ridges(center_axis_labels=TRUE) + theme(legend.position='none')
ggsave(file.path(figpath, 'xi2Income.pdf'), pXi2, width=9, height=5)

ggformat <- list(scale_color_discrete(name='Mood'), 
                 xlab('Days since the start of study'), 
                 ylab('Square-root proportion'))

# Mean
pMu <- ggplot(melt(resSp$muObs, varnames=c('variable', 'day')), aes(x=day, y=value, color=variable)) + geom_line(size=1.5) + ggformat
ggsave(file.path(figpath, sprintf('muSphere.pdf')), pMu, width=6, height=5)

# Eigenfunctions
phiPlotDat <- structure(resSp$phiObs, dimnames=
                          list(day = resSp$obsGridTrunc, 
                               variable = varnames, 
                               K = seq_len(KUse))) %>% 
  melt %>% 
  mutate(Kname = sprintf('phi[%s] ~ ~ %.2g', K, resSp$lam[K] / sum(resSp$lam)))
pPhi <- ggplot(phiPlotDat %>% filter(K <= 4), aes(x=day, y=value, color=variable)) + geom_line(size=1.5) + facet_grid( ~ Kname, scales='free_y', labeller=label_parsed) + ggformat
ggsave(file.path(figpath, sprintf('phi.pdf')), pPhi, width=10, height=5)

# Cumulative FVE
cumsum(resSp$lam) / sum(resSp$lam)

# Individual curves
ids <- c(255836830, 406714417, 548237021, 681292255)
KFit <- 7
mfdSp <- structure(1, class='Sphere')
fitPool <- Fitted(mfdSp, 
                  resSp$muObs[, resSp$obsGrid %in% resSp$obsGridTrunc, drop=FALSE], 
                  resSp$xi, 
                  resSp$phiObs, 
                  KFit)
dimnames(fitPool) <- list(caseid = rownames(resSp$xi), 
                      variable=varnames, 
                      day = resSp$obsGridTrunc)
fitPool <- melt(fitPool, value.name='fit')
plotDat <- right_join(datUse %>% 
                        melt(id.vars=c('caseid', 'day'), value.name='raw') %>% 
                        filter(caseid %in% ids) %>% 
                        mutate(raw = sqrt(raw)), 
                      fitPool %>% 
                        filter(caseid %in% ids))

# Jitter identical original observations for plotting
jitter <- 0.02
jitterIdentical <- function(x) {
  if (all(!duplicated(x))) {
  } else {
    firstDup <- x[which(duplicated(x))[1]]
    ind <- x == firstDup
    nDup <- sum(ind)
    xDup <- x[ind]
    xDup <- xDup + seq(jitter, jitter * nDup, by=jitter) - jitter * (nDup + 1) / 2
    if (min(xDup) < 1e-15) {
      xDup <- xDup - min(xDup)
    } else if (max(xDup) > 1) {
      stop('impossible')
    }
    x[ind] <- xDup
    x <- jitterIdentical(x)
  }
  x
}

plotDat <- ddply(plotDat, c('caseid', 'day'), function(dat) {
  if (!all(is.na(dat$raw))) {
    dat$raw <- jitterIdentical(dat$raw)
  }
  dat
})

pSamp <- ggplot(plotDat, aes(x=day, y=fit, color=variable)) + geom_line(linetype=1) + geom_point(data=plotDat, aes(x=day, y=raw, color=variable)) + geom_line(data=plotDat, aes(x=day, y=raw, color=variable), linetype=3) + facet_wrap(~caseid) + ggformat + theme(strip.background = element_blank(), strip.text.x = element_blank())
ggsave(file.path(figpath, sprintf('samples%d.pdf', length(ids))), pSamp, width=8, height=6)
