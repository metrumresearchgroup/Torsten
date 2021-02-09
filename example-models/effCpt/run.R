library("cmdstanr")
library("bayesplot")
library("posterior")
library("ggplot2")

mod <- cmdstan_model("effCpt.stan")
fit <- mod$sample(data="effCpt.data.R", init="effCpt.init.R"), seed=5432, parallel_chains=4)

pars <- c(
"CLHat"   ,  "QHat"     , "V1Hat"   ,  "V2Hat"   ,  "kaHat"   ,  "ke0Hat"   ,
"EC50Hat" ,  "omega[1]" , "omega[2]",  "omega[3]",  "omega[4]",  "omega[5]" ,
"rho[1,2]",  "rho[1,3]" , "rho[1,4]",  "rho[1,5]",  "rho[2,3]",  "rho[2,4]" ,
"rho[2,5]",  "rho[3,4]" , "rho[3,5]",  "rho[4,5]",  "omegaKe0",  "omegaEC50",
"sigma"   ,  "sigmaResp")

subset.pars <- subset_draws(fit$draws(), variable=pars)

## Write summary
write.csv(summarise_draws(subset.pars), file="summary_pars.csv")

## density plot
mcmc_dens_overlay(subset.pars, facet_args=list(ncol=4))+facet_text(size=10)+theme(axis.text=element_text(size=8))

## get data
source("effCpt.data.R")

## 25 subject times 14 observations per dose times 4 doses.
cobs.pred.summary <-
as_draws_df(fit$draws(variables=c("cObsPred"))) %>%
    summarize_all(function(x) {quantile(x,probs=c(0.05,0.5,0.95))}) %>%
    select(starts_with("cObsPred"))

## PK prediction
data <- rbind(cobs.pred.summary,
              unlist(mapply(rep, 1:nSubjects, (end - start + 1))), # subject id
              time) %>% t %>% as_tibble() %>%
        rename(lb=V1, median=V2, ub=V3, subject=V4, time=V5)

## PK obs data
obs.data <- tibble(time=time[iObs], y=cObs,
                   subject=unlist(mapply(rep, 1:nSubjects, (end - start + 1)))[iObs])

ppc.cobs <- function(start.id, end.id) {
    ggplot(subset(data, subject >= start.id & subject <= end.id)) + geom_ribbon(aes(x=time,ymin=lb,ymax=ub),fill="#b3cde0",alpha=0.8) + geom_line(aes(x=time,y=median),color="#005b96") + geom_point(data=subset(obs.data, subject >=start.id & subject <=end.id),aes(x=time,y=y),size=1.0)+
scale_x_continuous(name="time (h)") +
scale_y_continuous(name="plasma drug concentration (ng/mL)") +
    facet_wrap(.~subject)
}

ggsave("ppc_study_1_5mg.pdf", ppc.cobs(1,25))
ggsave("ppc_study_1_10mg.pdf", ppc.cobs(26,50))
ggsave("ppc_study_1_20mg.pdf", ppc.cobs(51,75))
ggsave("ppc_study_1_40mg.pdf", ppc.cobs(76,100))
ggsave("ppc_study_2_20mg.pdf", ppc.cobs(101,140), width=7, height=9)

## PD
resp.pred.summary <-
as_draws_df(fit$draws(variables=c("respObsPred"))) %>%
    summarize_all(function(x) {quantile(x,probs=c(0.05,0.5,0.95))}) %>%
    select(starts_with("respObsPred"))

pd.data <- rbind(resp.pred.summary,
                 unlist(mapply(rep, 1:nSubjects, (end - start + 1))), # subject id
                 time) %>% t %>% as_tibble() %>%
           rename(lb=V1, median=V2, ub=V3, subject=V4, time=V5)

resp.data <- tibble(time=time[iObs], y=respObs,
                    subject=unlist(mapply(rep, 1:nSubjects, (end - start + 1)))[iObs])

ppc.resp <- function(start.id, end.id) {
    ggplot(subset(pd.data, subject >= start.id & subject <= end.id)) + geom_ribbon(aes(x=time,ymin=lb,ymax=ub),fill="#b3cde0",alpha=0.8) + geom_line(aes(x=time,y=median),color="#005b96") + geom_point(data=subset(resp.data, subject >=start.id & subject <=end.id),aes(x=time,y=y),size=1.0)+
scale_x_continuous(name="time (h)") +
scale_y_continuous(name="response") +
    facet_wrap(.~subject)
}

ggsave("ppc_study_1_5mg_resp.pdf", ppc.resp(1,25))
ggsave("ppc_study_1_10mg_resp.pdf", ppc.resp(26,50))
ggsave("ppc_study_1_20mg_resp.pdf", ppc.resp(51,75))
ggsave("ppc_study_1_40mg_resp.pdf", ppc.resp(76,100))
ggsave("ppc_study_2_20mg_resp.pdf", ppc.resp(101,140), width=7,height=9)
