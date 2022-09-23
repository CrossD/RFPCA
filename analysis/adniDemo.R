## This file reproduces the ADNI functional connectivity analysis.
devtools::load_all()

library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(fdapace)
theme_set(theme_bw(base_size=14))
theme_update(legend.title=element_blank())
source('frechetTest.R')

figpathBase <- '../figures'
mfd <- structure(1, class=c('SPD', 'AffInv'))
mfdName <- 'affinv'

## Tuning parameters
bwMu <- 0.3
bwCov <- 0.3
kern <- 'gauss'
ToutMax <- 1.1
nRegGrid <- 21
## End Tuning parameters

set.seed(1)
ToutRange <- c(0, ToutMax)

# The brain regions are Buckner's hubs (Buckner 2009)
buck10Hubs <- c(
  "L I/S parietal lobule", 
  "M S frontal", 
  "M prefrontal", 
  "R I/S parietal lobule", 
  "L mid frontal", 
  "Post cingulate/precuneus", 
  "R supramarginal", 
  "L mid temporal", 
  "R mid temporal", 
  "R mid frontal")
buck10Hubs <- factor(buck10Hubs, buck10Hubs)
d <- length(buck10Hubs)

figpath <- file.path(figpathBase, 'ADNI', 'demo')
dir.create(figpath, showWarnings=FALSE, recursive=TRUE)

## Read in data: load ageDat and spdUse
# ageDat contains subject ID (subjectId), Alzheimer's disease status (researchGroup), number of observations per subject (n), age at the first visit (age0), and time since the first visit (newage) which we use as the time-axis.
# spdUse contains the SPD connectivity correlation matrices at each visit, organized in the long format. The columns are subject ID (subjectId), time since the first visit (newage), region 1 (i1), region (2), and the correlation values (value). 
load(file='adniDat.RData')

# Information at the first visit
ageDat0 <- ageDat %>% group_by(subjectId) %>% filter(row_number() == 1) %>% ungroup

## RPACE
fpcaInput <- spdUse %>% 
  reshape2::dcast(subjectId + newage ~ i1 + i2) %>% 
  MakeInput('subjectId', 'newage')
logLy <- fpcaInput$Ly
Lt <- lapply(fpcaInput$Lt, round, digits=2)

# Requires around 4 minutes to finish on a 2018 macbook pro
print(system.time({
  resMani <- RFPCA(logLy, Lt, 
                   list(mfdName=mfdName, 
                        userBwMu=bwMu, 
                        userBwCov=bwCov, 
                        kernel=kern, 
                        ToutRange=ToutRange, 
                        nRegGrid=nRegGrid
                        ))
}))

## Summary statistics: Number of subjects in each group
ageDat0 %>% filter(subjectId %in% names(Lt)) %>% .$researchGroup %>% table

## Lambda and FVE
lam <- resMani$lam
VE <- lam / sum(lam)
FVE <- cumsum(lam) / sum(lam)
cat(sprintf(
'lam: %s
VE: %s
FVE: %s

', 
paste(round(lam[1:20], 3), collapse=', '), 
paste(round(VE[1:20], 3), collapse=', '), 
paste(round(FVE[1:20], 3), collapse=', '))
)

## Visualization of model components
spdTheme <- list(
  geom_raster(),
  scale_fill_gradient2(low=scales::muted('blue'), high=scales::muted('red')), 
  coord_fixed(), 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

a <- resMani$muWork
d <- sqrt(dim(a)[1])
meanDat <- array(a, c(d, d, dim(a)[2]))
dimnames(meanDat) <- list(i1=buck10Hubs, i2=buck10Hubs, t=resMani$workGrid)
meanDat <- reshape2::melt(meanDat) 
allt <- unique(meanDat$t)
pmu <- ggplot(meanDat %>% 
                filter(t %in% allt[seq.int(1, length(allt), length.out=5)]), 
              aes(x=i1, y=i2, fill=value)) + 
  spdTheme + facet_wrap(~t, nrow=1)
ggsave(file.path(figpath, 'mu.pdf'), pmu, height=5, width=10)


# phi
K <- dim(resMani$phi)[3]
phiDat <- array(resMani$phi, c(dim(resMani$phi)[1], d, d, K))
dimnames(phiDat) <- list(t=allt, i1=buck10Hubs, i2=buck10Hubs, K=seq_len(K))
phiDat <- reshape2::melt(phiDat)


# Plot all
pphi3 <- phiDat %>% filter(K <= 3 & t %in% allt[seq.int(1, length(allt), length.out=5)]) %>% ggplot(aes(x=i1, y=i2, fill=value)) + spdTheme + facet_grid(K~t) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab('') + ylab('')
ggsave(file.path(figpath, 'phiAll.pdf'), pphi3, height=10, width=10, limitsize=FALSE)

# Raw and fitted values
fittedObs <- fitted(resMani, grid='obs') %>% 
  reshape2::melt() %>% 
  separate('j', c('i1', 'i2'), sep='_') %>% 
  rename(newage = t, subjectId = i) %>%
  group_by(newage, subjectId) %>% 
  mutate(fitted=cov2cor(matrix(value, d, d))) %>% 
  ungroup %>%
  select(-value)

rawDat <- spdUse %>% 
  inner_join(ageDat) %>% 
  filter(subjectId %in% names(logLy)) %>%
  rename(raw = value)

rawFitMean <- fittedObs %>% 
  inner_join(rawDat %>% 
               select(subjectId, newage, i1, i2, raw) %>%
               mutate(newage = round(newage, 2)), 
             by=c('subjectId', 'i1', 'i2', 'newage'))


## Plot individuals
allSub <- sort(unique(spdUse$subjectId))
plotNames <- allSub[c(10, 55, 16)]
pList <- lapply(plotNames, function(sname) {
  p <- ggplot(rawFitMean %>% filter(subjectId == sname) %>% gather(type, value, raw, fitted) %>% mutate(type = factor(type, c('raw', 'fitted'))), aes(x = i1, y = i2, fill = value)) + spdTheme + facet_grid(type ~ paste0('t = ', newage)) + xlab('') + ylab('')
  p <- p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(file.path(figpath, sprintf('%s.pdf', sname)), p, width=10, height=4)
  invisible(p)
})

# Frechet ANOVA in the bootstrap version.
# Reference: Dubey P, Müller H-G. Fréchet analysis of variance for random objects. Biometrika. 2019 Dec 1;106(4):803–21. 
a <- rawFitMean %>% filter(newage == 0) %>% inner_join(ageDat0 %>% select(subjectId, researchGroup))
fit1 <- a %>% 
  filter(researchGroup == 'AD') %>% 
  select(i1, i2, fitted, subjectId) %>% 
  reshape2::acast(i1 + i2 ~ subjectId, value.var='fitted') %>% 
  plyr::alply(2, matrix, d, d)
fit2 <- a %>% 
  filter(researchGroup == 'CN') %>% 
  select(i1, i2, fitted, subjectId) %>% 
  reshape2::acast(i1 + i2 ~ subjectId, value.var='fitted') %>% 
  plyr::alply(2, matrix, d, d)
raw1 <- a %>% 
  filter(researchGroup == 'AD') %>% 
  select(i1, i2, raw, subjectId) %>% 
  reshape2::acast(i1 + i2 ~ subjectId, value.var='raw') %>% 
  plyr::alply(2, matrix, d, d)
raw2 <- a %>% 
  filter(researchGroup == 'CN') %>% 
  select(i1, i2, raw, subjectId) %>% 
  reshape2::acast(i1 + i2 ~ subjectId, value.var='raw') %>% 
  plyr::alply(2, matrix, d, d)

subjectMean <- plyr::ddply(
  rawFitMean %>% filter(newage <= 1.1), 
  'subjectId', 
  function(dat) {
    X <- reshape2::acast(dat, i1 + i2 ~ newage, value.var='raw')
    mu <- frechetMean(mfd, X)
    mu <- as.numeric(cov2cor(matrix(mu, d, d)))
    data.frame(value=mu, rn = rownames(X))
  }
) %>% 
  inner_join(ageDat0 %>% select(subjectId, researchGroup))

ave1 <- subjectMean %>% 
  filter(researchGroup == 'AD') %>% 
  reshape2::acast(rn ~ subjectId) %>% 
  plyr::alply(2, matrix, d, d)
ave2 <- subjectMean %>% 
  filter(researchGroup == 'CN') %>% 
  reshape2::acast(rn ~ subjectId) %>% 
  plyr::alply(2, matrix, d, d)

R <- 1000 # Number of bootstrap samples
nCores <- 3L

# Requires around 90 minutes on a 2018 macbook pro with 3 cores
print(system.time({
  resRaw <- test(raw1, raw2, R=R, nCores=nCores)
  resAve <- test(ave1, ave2, R=R, nCores=nCores)
  resFit <- test(fit1, fit2, R=R, nCores=nCores)
}))

cat(sprintf('Raw bootstrap p-value: %.3g\n', resRaw$pval.boot))
cat(sprintf('Ave bootstrap p-value: %.3g\n', resAve$pval.boot))
cat(sprintf('Fit bootstrap p-value: %.3g\n', resFit$pval.boot))
