library("flexmix")
data("tribolium", package = "flexmix")
set.seed(1234)
## uses FLXMRglm() from package.
tribMix <- initFlexmix(cbind(Remaining, Total - Remaining) ~ Species, 
                       k = 2, nrep = 5, data = tribolium,
                       model = FLXMRglm(family = "binomial"))
parameters(tribMix); logLik(tribMix)
source("FLXMRglm.R")
set.seed(1234)
## uses FLXMRglm() from source file which allows to specify link.
tribMixNew <- initFlexmix(cbind(Remaining, Total - Remaining) ~ Species, 
                          k = 2, nrep = 5, data = tribolium,
                          model = FLXMRglm(family = "binomial"))
parameters(tribMixNew); logLik(tribMixNew)
set.seed(1234)
tribMixlogit <- initFlexmix(cbind(Remaining, Total - Remaining) ~ Species, 
                            k = 2, nrep = 5, data = tribolium,
                            model = FLXMRglm(family = binomial(link = logit)))
parameters(tribMixlogit); logLik(tribMixlogit)
set.seed(1234)
tribMixprobit <- initFlexmix(cbind(Remaining, Total - Remaining) ~ Species, 
                             k = 2, nrep = 5, data = tribolium,
                             model = FLXMRglm(family = binomial(link = probit)))
parameters(tribMixprobit); logLik(tribMixprobit)
