## "summary" comes from combining each benchmark results
## "extended.summary" has mutate results ess/time and ess/(time *
## chain), it can be also read in from "summary.csv"

library("tidyverse")
library("latex2exp")

## effect of # of parallel chains on ESS
ggplot(lp.summary %>% filter(method=="cross-chain"),
        aes(group=as.factor(target_ESS), color=as.factor(target_ESS), x=parallel_chains, y=ess_tail_per_sec_per_chain)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    geom_hline(yintercept=lp.summary[["ess_tail_per_sec_per_chain"]][1],linetype = "dashed") +
    scale_x_continuous(trans='log2') +
    labs(x = TeX("n_{chain}"), y = TeX("ESS_{tail}/(Time_{sampling} x n_{chain})"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_tail_per_time_chain.pdf",width=4, height=3.2)

ggplot(lp.summary %>% filter(method=="cross-chain"),
        aes(group=as.factor(target_ESS), color=as.factor(target_ESS), x=parallel_chains, y=ess_bulk_per_sec_per_chain)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    geom_hline(yintercept=lp.summary[["ess_bulk_per_sec_per_chain"]][1],linetype = "dashed") +
    scale_x_continuous(trans='log2') +
    labs(x = TeX("n_{chain}"), y = TeX("ESS_{bulk}/(Time_{sampling} x n_{chain})"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_bulk_per_time_chain.pdf",width=4, height=3.2)

ggplot(lp.summary %>% filter(method=="cross-chain"),
        aes(group=as.factor(target_ESS), color=as.factor(target_ESS), x=parallel_chains, y=ess_tail_per_sec)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    geom_hline(yintercept=lp.summary[["ess_tail_per_sec"]][1],linetype = "dashed") +
    scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') +
    labs(x = TeX("n_{chain}"), y = TeX("ESS_{tail}/Time_{sampling}"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_tail_per_time.pdf",width=4, height=3.2)

ggplot(lp.summary %>% filter(method=="cross-chain"),
        aes(group=as.factor(target_ESS), color=as.factor(target_ESS), x=parallel_chains, y=ess_bulk_per_sec)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    geom_hline(yintercept=lp.summary[["ess_bulk_per_sec"]][1],linetype = "dashed") +
    scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') +
    labs(x = TeX("n_{chain}"), y = TeX("ESS_{bulk}/Time_{sampling}"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_bulk_per_time.pdf",width=4, height=3.2)

## ESS/total time
ggplot(lp.summary %>% filter(method=="cross-chain"),
        aes(group=as.factor(target_ESS), color=as.factor(target_ESS), x=parallel_chains, y=ess_bulk_per_total_time_per_chain)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    geom_hline(yintercept=lp.summary[["ess_bulk_per_total_time_per_chain"]][1],linetype = "dashed") +
    scale_x_continuous(trans='log2') +
    labs(x = TeX("n_{chain}"), y = TeX("ESS_{bulk}/(Time_{total} x n_{chain})"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_bulk_per_total_time_chain.pdf",width=4, height=3.2)

ggplot(lp.summary %>% filter(method=="cross-chain"),
        aes(group=as.factor(target_ESS), color=as.factor(target_ESS), x=parallel_chains, y=ess_tail_per_total_time_per_chain)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    geom_hline(yintercept=lp.summary[["ess_tail_per_total_time_per_chain"]][1],linetype = "dashed") +
    scale_x_continuous(trans='log2') +
    labs(x = TeX("n_{chain}"), y = TeX("ESS_{tail}/(Time_{total} x n_{chain})"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_tail_per_total_time_chain.pdf",width=4, height=3.2)

## parallel efficiency
## speedup is based on ess/total_time, cross-chain parallel runs / sequential run
lp.speedup.200 <- lp.summary %>% filter(target_ESS==200) %>%
    mutate(ess_bulk_speedup=ess_bulk_per_total_time/(lp.summary$ess_bulk_per_total_time[1]),
           ess_tail_speedup=ess_tail_per_total_time/(lp.summary$ess_tail_per_total_time[1]))
lp.speedup.400 <- lp.summary %>% filter(target_ESS==400) %>%
    mutate(ess_bulk_speedup=ess_bulk_per_total_time/(lp.summary$ess_bulk_per_total_time[1]),
           ess_tail_speedup=ess_tail_per_total_time/(lp.summary$ess_tail_per_total_time[1]))
lp.speedup.600 <- lp.summary %>% filter(target_ESS==600) %>%
    mutate(ess_bulk_speedup=ess_bulk_per_total_time/(lp.summary$ess_bulk_per_total_time[1]),
           ess_tail_speedup=ess_tail_per_total_time/(lp.summary$ess_tail_per_total_time[1]))
lp.speedup <- bind_rows(lp.speedup.200, lp.speedup.400, lp.speedup.600) %>%
    mutate(ess_bulk_efficiency=ess_bulk_speedup/(parallel_chains/4), ess_tail_efficiency=ess_tail_speedup/(parallel_chains/4))

ggplot(lp.speedup, aes(group=as.factor(target_ESS),
                                   color=as.factor(target_ESS), x=parallel_chains,
                                   y=ess_bulk_speedup)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') +
    labs(x = TeX("n_{chain}(=n_{proc})"), y = TeX("Speedup of ESS_{bulk}/(Time_{total})"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_bulk_per_total_time_speedup.pdf",width=4, height=3.2)

ggplot(lp.speedup, aes(group=as.factor(target_ESS),
                                   color=as.factor(target_ESS), x=parallel_chains,
                                   y=ess_tail_speedup)) +
    geom_line(aes(linetype=as.factor(target_ESS))) +
    geom_point(aes(shape=as.factor(target_ESS))) +
    scale_x_continuous(trans='log2') + scale_y_continuous(trans='log2') +
    labs(x = TeX("n_{chain}(=n_{proc})"), y = TeX("Speedup of ESS_{tail}/(Time_{total})"),  color="Target ESS", shape="Target ESS",linetype="Target ESS")
ggsave(filename="ess_tail_per_total_time_speedup.pdf",width=4, height=3.2)
