library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

##Draw Runtime Plots
runtime_k5 <- read.csv("./Results/runtime_k5.txt", header = F, sep = "\t")
names(runtime_k5) <- c("n", "k", "s", "total", "one_round", "opttime", "num_rounds")
head(runtime_k5)

runtime_k5_longer <- runtime_k5[,c("n", "k", "s", "total", "opttime")] %>% pivot_longer(cols = c("total", "opttime"))
head(runtime_k5_longer)


ggplot(data=runtime_k5_longer,
       aes(x=n,y=value/60, color = name, group = name))+
  geom_point(alpha=0.5)+
  geom_smooth(data = runtime_k5_longer[runtime_k5_longer$n > 100,],method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="#d95f02",x=1000,y=20,geom="text",
           label=round(with(runtime_k5_longer[runtime_k5_longer$n > 100 & runtime_k5_longer$name == "total",], lm(log10(value)~log10(n)))[[1]][[2]],1))+
  annotate(color="#1b9e77",x=1000,y=1e-3,geom="text",
            label=round(with(runtime_k5_longer[runtime_k5_longer$n > 100 & runtime_k5_longer$name == "opttime",], lm(log10(value)~log10(n)))[[1]][[2]],1))+
  scale_y_continuous(trans="log10",expression(runtime~(minutes)))+
  scale_color_manual(values =c("#1b9e77", "#d95f02"),name="", labels=c("small problem", "total"))+
  scale_shape(name="")+
  theme(
    legend.position = c(0.3, 0.95),   # place inside: (x,y) in [0,1] relative to plot  # center alignment
    legend.box = "vertical"
  )
# ggsave("./Plots/total_runtime_k5.pdf",width = 3,height =2.5)
ggsave("./Plots/total_runtime_k5_combined.pdf",width = 3,height =3)

ggplot(data=runtime_k5,
       aes(x=n,y=total/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(data = runtime_k5[runtime_k5$n > 100,],method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=1000,y=20,geom="text",
           label=round(with(runtime_k5[runtime_k5$n > 100,], lm(log10(total)~log10(n)))[[1]][[2]],1))+
  scale_y_continuous(trans="log10",expression(total~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
# ggsave("./Plots/total_runtime_k5.pdf",width = 3,height =2.5)
ggsave("./Plots/total_runtime_k5.pdf",width = 2,height =2)

ggplot(data=runtime_k5,
       aes(x=n,y=opttime/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(data =runtime_k5[runtime_k5$n > 100,], method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=1000,y=35e-5,geom="text",
           label=round(with(runtime_k5[runtime_k5$n > 100,], lm(log10(opttime)~log10(n)))[[1]][[2]],1))+
  scale_y_continuous(trans="log10",expression(opt~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/optimization_runtime_k5.pdf",width = 2,height =2)


summary(runtime_k5$num_rounds / runtime_k5$k)
ggplot(data=runtime_k5,
       aes(x=n,y=num_rounds/k, color = "#1b9e77"))+
  # geom_point(alpha=0.5)+
  # stat_summary()+
  geom_boxplot(aes(group = n))+
  geom_smooth(data =runtime_k5[runtime_k5$n > 100,],method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n), breaks = c(50, 100, 200, 500, 1000, 2000))+
  annotate(color="grey20",x=1000,y=3.4,geom="text",
           label=round(with(runtime_k5[runtime_k5$n > 100,], lm(log10(num_rounds/k)~log10(n)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(rounds))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/num_rounds_k5.pdf",width = 3,height =2.5)

##Runtime for k
runtime_n500 <- read.csv("./Results/runtime_n500.txt", header = F, sep = "\t")
names(runtime_n500) <- c("n", "k", "s", "total", "one_round", "opttime", "num_rounds")
head(runtime_n500)

runtime_n500_longer <- runtime_n500[,c("n", "k", "s", "total", "opttime")] %>% pivot_longer(cols = c("total", "opttime"))
head(runtime_n500_longer)


ggplot(data=runtime_n500_longer,
       aes(x=k,y=value/60, color = name, group = name))+
  geom_point(alpha=0.5)+
  geom_smooth(data = runtime_n500_longer[runtime_n500_longer$k >= 20,],method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="#d95f02",x=70,y=350,geom="text",
           label=round(with(runtime_n500_longer[runtime_n500_longer$k >= 20 & runtime_n500_longer$name == "total",], lm(log10(value)~log10(k)))[[1]][[2]],1))+
  annotate(color="#1b9e77",x=70,y=5e-3,geom="text",
           label=round(with(runtime_n500_longer[runtime_n500_longer$k >= 20 & runtime_n500_longer$name == "opttime",], lm(log10(value)~log10(k)))[[1]][[2]],1))+
  scale_y_continuous(trans="log10",expression(runtime~(minutes)))+
  scale_color_manual(values =c("#1b9e77", "#d95f02"),name="", labels=c("small problem", "total"))+
  scale_shape(name="")+
  theme(
    legend.position = "none",   # place inside: (x,y) in [0,1] relative to plot  # center alignment
    legend.box = "vertical"
  )
# ggsave("./Plots/total_runtime_k5.pdf",width = 3,height =2.5)
ggsave("./Plots/total_runtime_n500_combined.pdf",width = 3,height =3)

unique(runtime_n500$k)
ggplot(data=runtime_n500,
       aes(x=k,y=total/60, color = "#d95f02"))+
  geom_point(alpha=0.5)+
  geom_smooth(data = runtime_n500[runtime_n500$k >= 20,], method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=35,y=200,geom="text",
           label=paste0("slope=",round(with(runtime_n500[runtime_n500$k >= 20,], lm(log10(total)~log10(k)))[[1]][[2]],1)))+
  scale_y_continuous(trans="log10",expression(total~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
# ggsave("./Plots/total_runtime_n500.pdf",width = 3,height =2.5)
ggsave("./Plots/total_runtime_n500.pdf",width = 2,height =2)

ggplot(data=runtime_n500,
       aes(x=k,y=opttime/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(data = runtime_n500[runtime_n500$k >= 20,], method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=65,y=2e-3,geom="text",
           label=round(with(runtime_n500[runtime_n500$k >= 20,], lm(log10(opttime)~log10(k)))[[1]][[2]],1))+
  scale_y_continuous(trans="log10",expression(opt~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/optimization_runtime_n500.pdf",width = 2,height =2)

ggplot(data=runtime_n500,
       aes(x=k,y=num_rounds/k, color = "#1b9e77"))+
  # geom_point(alpha=0.5)+
  # stat_summary()+
  geom_boxplot(aes(group = k))+
  geom_smooth(data =runtime_n500,method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=100,y=3,geom="text",
           label=round(with(runtime_n500, lm(log10(num_rounds/k)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(rounds))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/num_rounds_n500.pdf",width = 3,height =2.5)

ggplot(data=runtime_n500,
       aes(x=k,y=num_rounds/(k*(k+1))*2, color = "#1b9e77"))+
  # geom_point(alpha=0.5)+
  geom_boxplot(aes(group = k))+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=11,y=17,geom="text",
           label=round(with(runtime_n100, lm(log10(num_rounds/(k*(k+1))*2)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(rounds))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/num_rounds_n100.pdf",width = 4,height =4)

runtime_n100 <- read.csv("./Results/runtime_n100.txt", header = F, sep = "\t")
names(runtime_n100) <- c("n", "rep", "s", "k", "runtime")
head(runtime_n100)

ggplot(data=runtime_n100,
       aes(x=k,y=runtime/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=10,y=15,geom="text",
           label=round(with(runtime_n100, lm(log10(runtime)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(runtime~(minutes)))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/runtime_n100.pdf",width = 4,height =4)


##real data running time
data <- read.csv("16k_data_placements/decodiphy_runtime.txt", sep = "\t", header=F)
names(data) <- c("sample", "tree", "runtime")
head(data)

data %>% group_by(sample) %>% summarise(n =sum(runtime/3600)) %>% summarise(mean(n))


##draw ROC 
pr_n100 <- read.csv("./Results/pr_n100.txt", header = F, sep = "\t")
names(pr_n100) <- c("n", "rep", "s", "k", "ck", "precision", "recall", "weighted unifrac", "loss", "|p=0|")
head(pr_n100)

min_loss <- pr_n100 %>% group_by(n, rep, s) %>% summarise(min_loss_log = min(log(loss)))

pr_longer <- merge(pr_n100, min_loss) %>% mutate(`1 - wUniFrac` = 1 - `weighted unifrac`, `p=0 ratio` = `|p=0|` / ck, "norm log loss" = log(loss) / min_loss_log) %>% pivot_longer(cols = c("precision", "recall", "1 - wUniFrac", "norm log loss", "p=0 ratio"), names_to = "meassure", values_to = "value")

pr_longer$meassure <- factor(pr_longer$meassure, levels = c("1 - wUniFrac", "precision", "recall", "norm log loss", "p=0 ratio"), labels = c("1 - wUniFrac", "precision", "recall", "norm log loss", "p ≈ 0 (ratio)"))

ggplot(pr_longer, aes(x = as.factor(ck), y = value, group = meassure, color = meassure)) +
  geom_vline(xintercept = 4, color = "grey30", linetype = "dashed")+
  geom_hline(yintercept = 1, color = "grey30", linetype = "dashed")+
  stat_summary()+
  stat_summary(geom = "line")+
  theme_bw() +
  
  labs(
    x = "number of anchors",
    color = "metric"
  ) +
  scale_color_manual(values = c("#1b9e77", "#377eb8", "#96aef4", "#ff7f00", "#d53e4f"))+
  theme(
    legend.position = "right"
  )+
  guides(color = guide_legend(nrow = 5))
 ggsave("./Plots/pr_n100.pdf",width = 5,height =2.5)
 
 
##loss for noise
loss <- read.csv("./Results/loss_n100_k5.txt", sep = "\t", header = F)
loss <- read.csv("./Results/minp_n100_k5.txt", sep = "\t", header = F)
names(loss) <- c("n","k","s","e","round","loss")
head(loss)

ggplot(loss %>% filter(round < 10), aes(x = as.factor(round), y = loss, group = e, color = e)) +
  geom_vline(xintercept = 5, color = "grey30", linetype = "dashed")+
  geom_hline(yintercept = 0.01, color = "grey30", linetype = "dashed")+
  stat_summary()+
  stat_summary(geom = "line")+
  theme_bw()+
  scale_y_log10()+
  labs(
    x = "number of anchors",
    color = "noise (exp scale)"
  ) +
  # scale_color_manual(values = c("#1b9e77", "#377eb8", "#96aef4", "#ff7f00", "#d53e4f"))+
  theme(
    legend.position = "right"
  )+
guides(color = guide_legend(nrow = 5))
ggsave("./Plots/loss_n100_k5.pdf",width = 5,height =2.5)
 
 ##distances
 all_dist <- read.csv("./Demixing_Data/biotrees/all_distances.txt", header = F, sep = "\t")
 names(all_dist) <- c("data", "k", "s", "n", "l", "species", "distance")
head(all_dist) 

true_dist <- all_dist[all_dist$n == 0, c("data", "k", "s", "species", "distance")]
names(true_dist)[5] <- "true_distance"
head(true_dist)

merged_dist <- merge(all_dist[all_dist$n > 0, ], true_dist)
head(merged_dist)

merged_dist$diff <- merged_dist$distance - merged_dist$true_distance

ggplot(data = merged_dist[merged_dist$n == 1,], aes(x=as.factor(l), y = abs(diff), color = as.factor(n), group = as.factor(interaction(l,n))))+
  stat_summary()+
  facet_wrap(~k)

ggplot(data = merged_dist[merged_dist$data == "bees-bosseret" & merged_dist$s == 1 & merged_dist$n == 1,], aes(x=as.factor(l), y = abs(diff), color = as.factor(n), group = as.factor(interaction(l,n))))+
  stat_summary()+
  facet_wrap(~k)

ggplot(merged_dist[merged_dist$n ==2 ,], aes(x = distance, y = true_distance, color = k)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  geom_point()+
  facet_grid(n~l)

all_dist <- read.csv("./Demixing_Data/biotrees/bees-bosseret/all_distances.txt", header = F, sep = "\t")
names(all_dist) <- c("k", "s", "n", "l", "node1", "node2", "diff")
head(all_dist) 

# true_dist <- all_dist[all_dist$n == 0, c("k", "s", "node", "length")]
# names(true_dist)[4] <- "true_length"
# head(true_dist)

# merged_dist <- merge(all_dist[all_dist$n > 0, ], true_dist)
# head(merged_dist)

# merged_dist$diff <- merged_dist$length - merged_dist$true_length

ggplot(data = all_dist, aes(x=as.factor(l), y = abs(diff), color = as.factor(n), group = as.factor(interaction(l,n))))+
  stat_summary()+
  facet_wrap(~k)

ggplot(merged_dist, aes(x = length, y = true_length, color = k)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey30") +
  geom_point()+
  facet_grid(n~l)
