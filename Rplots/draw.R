library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

runtime_bio <- read.csv("./Demixing_Data/biotrees/all_runtime_res.txt", header = F, sep = "\t")
names(runtime_bio) <- c("data", "k", "s", "noise", "seq_length", "n", "est_k", "total", "one_round")
head(runtime_bio)

runtime_bio$total_n <- runtime_bio$k + runtime_bio$n
data_size <- runtime_bio %>% group_by(data) %>% summarise(n = mean(total_n))


ggplot(data=runtime_bio[runtime_bio$noise == 0,],
       aes(x=n,y=total/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=420,y=9,geom="text",
           label=round(with(runtime_bio[runtime_bio$noise == 0,], lm(log10(total)~log10(n)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(total~runtime~(minutes)))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/total_runtime_bio.pdf",width = 3,height =2)

ggplot(data=runtime_bio[runtime_bio$noise == 0,],
       aes(x=n,y=one_round/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=400,y=4,geom="text",
           label=round(with(runtime_bio[runtime_bio$noise == 0,], lm(log10(one_round)~log10(n)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(one~round~runtime~(minutes)))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/one_round_runtime_bio.pdf",width = 4,height =4)

ggplot(data=runtime_bio[runtime_bio$noise == 0,],
       aes(x=k,y=total/60, color = "#d95f02"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=11,y=9,geom="text",
           label=round(with(runtime_bio[runtime_bio$noise == 0,], lm(log10(total)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(total~runtime~(minutes)))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/total_runtime_k_bio.pdf",width = 4,height =4)

ggplot(data=runtime_bio[runtime_bio$noise == 0,],
       aes(x=k,y=one_round/60, color = "#d95f02"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=11,y=0.5,geom="text",
           label=round(with(runtime_bio[runtime_bio$noise == 0,], lm(log10(one_round)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(one~round~runtime~(minutes)))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/one_round_runtime_k_bio.pdf",width = 4,height =4)

##Draw Runtime Plots
runtime_k5 <- read.csv("./Results/runtime_k5.txt", header = F, sep = "\t")
names(runtime_k5) <- c("n", "rep", "s", "k", "n-k", "est_k", "total", "one_round", "opttime", "num_rounds")
head(runtime_k5)

ggplot(data=runtime_k5,
       aes(x=n,y=total/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(data = runtime_k5[runtime_k5$n >= 200,],method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=600,y=30,geom="text",
           label=round(with(runtime_k5[runtime_k5$n >= 200,], lm(log10(total)~log10(n)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(total~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/total_runtime_k5.pdf",width = 3,height =2.5)

ggplot(data=runtime_k5,
       aes(x=n,y=one_round/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=1000,y=8,geom="text",
           label=round(with(runtime_k5, lm(log10(one_round)~log10(n)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(one~round~runtime~(minutes)))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/one_round_runtime_k5.pdf",width = 3,height =2.5)

ggplot(data=runtime_k5,
       aes(x=n,y=opttime/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(data =runtime_k5[runtime_k5$n >= 200,], method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=800,y=8e-4,geom="text",
           label=round(with(runtime_k5[runtime_k5$n >= 200,], lm(log10(opttime)~log10(n)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(optimization~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/optimization_runtime_k5.pdf",width = 3,height =2.5)

ggplot(data=runtime_k5,
       aes(x=n,y=num_rounds, color = "#1b9e77"))+
  # geom_point(alpha=0.5)+
  # stat_summary()+
  geom_boxplot(aes(group = n))+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(n))+
  annotate(color="grey20",x=1200,y=11.3,geom="text",
           label=round(with(runtime_k5, lm(log10(num_rounds)~log10(n)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(rounds))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/num_rounds_k5.pdf",width = 3,height =2.5)

##Runtime for k
runtime_n100 <- read.csv("./Results/runtime_n100.txt", header = F, sep = "\t")
names(runtime_n100) <- c("n", "rep", "s", "k", "n-k", "est_k", "total", "one_round", "opttime", "num_rounds")
head(runtime_n100)

ggplot(data=runtime_n100[runtime_n100$total < 8*60,],
       aes(x=k,y=total/60, color = "#d95f02"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=8,y=10,geom="text",
           label=round(with(runtime_n100[runtime_n100$total < 8*60,], lm(log10(total)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(total~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/total_runtime_n100.pdf",width = 3,height =2.5)

ggplot(data=runtime_n100,
       aes(x=k,y=one_round/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=10,y=0.8,geom="text",
           label=round(with(runtime_n100, lm(log10(one_round)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(one~round~runtime~(minutes)))+
  scale_color_brewer(palette = "Dark2",name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/one_round_runtime_n100.pdf",width = 4,height =4)

ggplot(data=runtime_n100[runtime_n100$opttime < 1e-3 * 60,],
       aes(x=k,y=opttime/60, color = "#1b9e77"))+
  geom_point(alpha=0.5)+
  geom_smooth(method="lm",se=F)+
  theme_classic()+
  scale_x_continuous(trans="log10",name=expression(k))+
  annotate(color="grey20",x=9,y=3.1e-4,geom="text",
           label=round(with(runtime_n100[runtime_n100$opttime < 1e-3 * 60,], lm(log10(opttime)~log10(k)))[[1]][[2]],2))+
  scale_y_continuous(trans="log10",expression(optimization~runtime~(minutes)))+
  scale_color_manual(values =c("#01665e"),name="")+
  scale_shape(name="")+
  theme(legend.position = "none")
ggsave("./Plots/optimization_runtime_n100.pdf",width = 3,height =2.5)

ggplot(data=runtime_n100,
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


##draw ROC 
pr_n100 <- read.csv("./Results/pr_n100.txt", header = F, sep = "\t")
names(pr_n100) <- c("n", "rep", "s", "k", "ck", "precision", "recall", "weighted unifrac", "loss", "|p=0|")
head(pr_n100)

min_loss <- pr_n100 %>% group_by(n, rep, s) %>% summarise(min_loss_log = min(log(loss)))

pr_longer <- merge(pr_n100, min_loss) %>% mutate(`1 - wUniFrac` = 1 - `weighted unifrac`, `p=0 ratio` = `|p=0|` / ck, "norm log loss" = log(loss) / min_loss_log) %>% pivot_longer(cols = c("precision", "recall", "1 - wUniFrac", "norm log loss", "p=0 ratio"), names_to = "meassure", values_to = "value")

pr_longer$meassure <- factor(pr_longer$meassure, levels = c("1 - wUniFrac", "precision", "recall", "norm log loss", "p=0 ratio"))

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

 


