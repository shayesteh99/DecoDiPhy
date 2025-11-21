library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

total_p = read.csv("./16k_sim/all_p.txt", sep="\t", header = F)
names(total_p) <- c("sample", "tree", "p")
head(total_p)

total_count = read.csv("./16k_sim/all_c.txt", sep="\t", header = F)
names(total_count) <- c("sample", "tree", "count")
head(total_count)

all_rounds = read.csv("./16k_sim/all_rounds_data.txt" , sep = "\t", header = F)
names(all_rounds) <- c("sample", "tree", "true_k", "round", "wunifrac", "loss")
head(all_rounds)

all_random_new = read.csv("./16k_sim/all_random_fix_prev.txt" , sep = "\t", header = F)
names(all_random_new) <- c("sample", "tree", "true_k", "round", "not_so_random")
head(all_random_new)

all_random = read.csv("./16k_sim/all_random.txt" , sep = "\t", header = F)
names(all_random) <- c("sample", "tree", "true_k", "round", "random")
head(all_random)

all_lower = read.csv("./16k_sim/all_lower_bound.txt" , sep = "\t", header = F)
names(all_lower) <- c("sample", "tree", "lower")
head(all_lower)

all_n = read.csv("./16k_sim/all_tree_size.txt" , sep = "\t", header = F)
names(all_n) <- c("sample", "tree", "n")
head(all_n)

all_num = read.csv("./16k_sim/all_placement_size.txt" , sep = "\t", header = F)
names(all_num) <- c("sample", "tree", "placements")
head(all_num)

all_rounds <- merge(merge(all_rounds, all_random[,c("sample", "tree", "round", "random")]), total_p)
all_rounds <- merge(all_rounds, all_random_new[,c("sample", "tree", "round", "not_so_random")])
all_rounds <- merge(all_rounds, all_lower)
all_rounds <- merge(all_rounds, all_n)
all_rounds <- merge(all_rounds, total_count)
all_rounds <- merge(all_rounds, all_num)
head(all_rounds)

all_rounds$diff <- all_rounds$round - all_rounds$true_k
all_rounds$var <- all_rounds$lower * all_rounds$count / all_rounds$n
all_rounds$b2 <- all_rounds$loss - (all_rounds$var * (all_rounds$n - (2*all_rounds$round + 1)) / all_rounds$count)
all_rounds$SNR <- all_rounds$count * all_rounds$b2 / (all_rounds$n * all_rounds$var) *all_rounds$p
# all_rounds$MSB <- all_rounds$b2 / all_rounds$n *all_rounds$p

all_rounds$norm_loss <- all_rounds$loss * all_rounds$count / all_rounds$n
all_rounds$thresh <- exp(-2*log(all_rounds$n)/all_rounds$n)
all_rounds$BIC <- (all_rounds$n) * log(all_rounds$loss / (all_rounds$n)) + (3*all_rounds$round + 1) * log(all_rounds$n)
head(all_rounds)
temp <- all_rounds[(all_rounds$loss - all_rounds$lower)* all_rounds$count * all_rounds$round**2 / all_rounds$n < 1,] %>% group_by(sample, tree, true_k) %>% summarise(min_k = min(round))

ggplot(data= temp, aes(x=true_k, y=min_k-true_k))+
  geom_hline(yintercept = 0, color = "grey30", linetype = "dashed")+
  stat_summary()

all_rounds <- all_rounds %>%
  arrange(sample, tree, round) %>%
  group_by(sample, tree) %>%
  mutate(loss_ratio = loss / lag(loss),loss_delta= lag(loss)-loss) %>%
  ungroup()
head(all_rounds)
head(all_rounds[,c(c("sample", "tree", "true_k", "round", "thresh", "loss_ratio"))])
View(all_rounds[,c("sample", "tree", "true_k", "round", "thresh", "loss_ratio")])

ggplot(all_rounds, 
       aes(x = round, y = 1/loss_ratio)) +
  geom_vline(aes(xintercept = 0), color = "grey30", linetype = "dashed")+
  stat_summary()+
  stat_summary(geom = "line")+
  # geom_boxplot(varwidth = FALSE)+
  # geom_boxplot(aes(y=random),color="blue",varwidth = FALSE)+
  # geom_boxplot(aes(y=not_so_random),color="red",varwidth = FALSE)+
  theme_bw()+
  # scale_color_viridis_b(n.breaks=10)+
  # stat_smooth(se=F)+
  # geom_hline(yintercept = 5e-2)+
  facet_wrap(~cut(true_k, c(0,2,3,4,5,10,20,30,40,50)),ncol=5, scales = "free")+
  coord_cartesian(ylim = c(1,5))
  scale_y_log10()
  
  all_rounds %>%
    filter(true_k<15) %>%
    filter(round<2*p*100) %>%
    filter(tree==8) %>%
    filter(sample %in% unique(all_rounds$sample)[20:50])%>%
    arrange(sample,round)%>%
    group_by(sample) %>%
    mutate(wd = loss/((lag(loss,1)+lag(loss,2)+lag(loss,3))/3))%>%
    ungroup()%>%
  ggplot( 
         aes(y = wd, x = round,color=as.factor(true_k),shape=(true_k==round),group=sample)) +
    geom_vline(aes(xintercept = 0), color = "grey30", linetype = "dashed")+
    geom_point()+
    geom_path()+
    #geom_smooth(se=F)+
    #geom_abline(slope = 10^4)+
    #stat_summary(geom = "line")+
    # geom_boxplot(varwidth = FALSE)+
    # geom_boxplot(aes(y=random),color="blue",varwidth = FALSE)+
    # geom_boxplot(aes(y=not_so_random),color="red",varwidth = FALSE)+
    theme_bw()+
    facet_wrap(~sample,scales="free")+
    # scale_color_viridis_b(n.breaks=10)+
    # stat_smooth(se=F)+
    # geom_hline(yintercept = 5e-2)+
    #facet_wrap(~cut(round, c(0,3,6,10,20,30,40,50)),ncol=5, scales = "free")+
    #coord_cartesian(ylim = c(1,5))+
    #scale_y_log10()+
    scale_color_viridis_d()
scale_y_log10()


  all_rounds %>%
    filter(round==4) %>%
    ggplot( 
      aes(x = tree, y = loss/n)) +
    geom_vline(aes(xintercept = 0), color = "grey30", linetype = "dashed")+
    stat_summary()+
    theme_bw()

ggplot(all_rounds, 
       aes(x = cut(diff, c(-50, -20, -10, -5, -1, 1, 5, 10, 20, 50)), y = wunifrac /random , group = cut(diff, c(-50, -20, -10, -5, -1, 1, 5, 10, 20, 50)))) +
  geom_vline(aes(xintercept = which(levels(cut(diff, c(-50, -20, -10, -5, -1, 1, 5, 10, 20, 50), include.lowest = TRUE)) == "(-1,1]")), color = "grey30", linetype = "dashed")+
  # stat_summary()+
  # stat_summary(geom = "line")+
  geom_boxplot(varwidth = FALSE)+
  geom_boxplot(aes(y=random),color="blue",varwidth = FALSE)+
  geom_boxplot(aes(y=not_so_random),color="red",varwidth = FALSE)+
  theme_bw()+
  # scale_color_viridis_b(n.breaks=10)+
  # stat_smooth(se=F)+
  geom_hline(yintercept = 5e-2)+
  facet_wrap(~cut(true_k, 5),ncol=5)+
  scale_y_log10()

ggplot(all_rounds, 
       aes(x = cut(diff, c(-50, -20, -10, -5, -1, 1, 5, 10, 20, 50)), y = SNR , group = cut(diff, c(-50, -20, -10, -5, -1, 1, 5, 10, 20, 50)))) +
  geom_vline(aes(xintercept = which(levels(cut(all_rounds$diff, c(-50, -20, -10, -5, -1, 1, 5, 10, 20, 50), include.lowest = TRUE)) == "(-1,1]")), color = "grey30", linetype = "dashed")+
  # stat_summary()+
  # stat_summary(geom = "line")+
  geom_boxplot(varwidth = FALSE)+
  theme_bw()+
  # scale_color_viridis_b(n.breaks=10)+
  # stat_smooth(se=F)+
  geom_hline(yintercept = 1e+9)+
  facet_wrap(~cut(true_k, 5),ncol=5)+
  scale_y_log10()


info_data <- read.csv("./16k_sim/all_info.txt", sep="\t", header = F)
names(info_data) <- c("sample", "tree", "true_k", "round", "wunifrac", "wunifrac_true", "wunifrac_krepp", "loss", "min_p", "last_p")
head(info_data)

total_p = read.csv("./16k_sim/all_p.txt", sep="\t", header = F)
names(total_p) <- c("sample", "tree", "p")
head(total_p)

total_count = read.csv("./16k_sim/all_c.txt", sep="\t", header = F)
names(total_count) <- c("sample", "tree", "count")
head(total_count)

info_data <- merge(info_data, total_p)
info_data <- merge(info_data, total_count)

krepp_wunifrac <- read.csv("./16k_sim/all_krepp_wunifrac.txt", sep = "\t", header = F)
names(krepp_wunifrac) <- c("sample", "tree", "krepp_wunifrac_to_true")
head(krepp_wunifrac)

mean(krepp_wunifrac$krepp_wunifrac_to_true)

info_data <- merge(info_data, krepp_wunifrac)

info_data %>% group_by(sample, tree) %>% 
  slice_min(order_by = wunifrac_true, n = 1) %>%
  mutate(min_k = round) %>%
  filter(true_k > 1 & round == min_k) %>%
  ungroup() %>% 
  summarise(mean(wunifrac_true))


info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == true_k) %>%
  ungroup() %>% summarise(mean(wunifrac_true))

info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == sample(2:max(round), 1)) %>%
  ungroup() %>% summarise(mean(wunifrac_true))


info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == sample(2:max(round), 1)) %>%
  ggplot(aes(x = wunifrac_true))+
  geom_freqpoly()+
  geom_freqpoly(aes(x = min_error), color = "blue3")+
  geom_vline(aes (xintercept = mean(wunifrac_true)))+
  geom_vline(aes (xintercept = mean(min_error)),color = "blue3")+
  # facet_wrap(~cut(true_k, 10), scales = "free")+
  scale_x_log10()
ggsave("./16k_sim/random.png",width = 4,height = 3)

info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(round == true_k) %>%
  ggplot(aes(x = wunifrac_true))+
  geom_freqpoly()+
  geom_freqpoly(aes(x = min_error), color = "blue4")+
  geom_vline(aes (xintercept = mean(wunifrac_true)))+
  geom_vline(aes (xintercept = mean(min_error)),color = "blue3")+
  # facet_wrap(~cut(true_k, 10), scales = "free")+
  scale_x_log10()
ggsave("./16k_sim/true_k.png",width = 4,height = 3)


info_data <- info_data %>% 
  # filter(loss * p > 1e-5) %>%
  group_by(sample, tree) %>%
  # mutate(knee = max(round))
  mutate(loss_knee = kneedle(x = round, y = log(loss),s=0.5)[1],
         p_knee = kneedle(x = round, y = min_p)[1],
         w_knee = kneedle(x = round, y = wunifrac_krepp)[1])

info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == w_knee) %>%
  ungroup() %>% summarise(mean(wunifrac_true))

info_data_ratio %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(round == knee) %>%
  ggplot(aes(x = wunifrac_true))+
  geom_freqpoly()+
  geom_freqpoly(aes(x = min_error), color = "blue4")+
  geom_vline(aes (xintercept = mean(wunifrac_true)))+
  geom_vline(aes (xintercept = mean(min_error)),color = "blue3")+
  # facet_wrap(~cut(true_k, 10), scales = "free")+
  scale_x_log10()
ggsave("./16k_sim/knee.png",width = 4,height = 3)

info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(round == w_knee) %>%
ggplot(aes(x = true_k, y = w_knee - true_k))+
  geom_hline(yintercept = 0)+
  geom_point(alpha= 0.2)


ggplot(data = info_data[info_data$sample %in% unique(info_data$sample)[15:19] &
                                info_data$p < 10.1& info_data$p > -0.03 
                              # info_data$loss * info_data$p >= 1e-5
                              ,], 
       aes(x = round)) +
  geom_vline(aes(xintercept = true_k), color = "grey30", linetype = "dashed")+
  # geom_vline(aes(xintercept = loss_knee), color = "purple3", linetype = "dashed")+
  # geom_vline(aes(xintercept = p_knee), color = "blue3", linetype = "dashed")+
  # geom_vline(aes(xintercept = w_knee), color = "pink3", linetype = "dashed")+
  geom_hline(aes(yintercept = krepp_wunifrac_to_true))+
  geom_point(aes(y=wunifrac_true), color = "green4")+
  geom_line(aes(y=wunifrac_true), color = "green4")+
  geom_point(aes(y=wunifrac_krepp), color = "red4")+
  geom_line(aes(y=wunifrac_krepp), color = "red4")+
  facet_wrap(~p, scales = "free")+
  # geom_hline(yintercept = 0.05)+
  scale_y_log10()

#read threshold
table(cut(info_data$count, c(0, 500, 1000, 5000, 10000, 50000, 1e5, 1e6, 1e7, 1e10)))

# info_data <- info_data %>% filter(round <= count/1000)
info_data <- merge(info_data, info_data %>% filter(min_p < 1000/count) %>% group_by(sample, tree, true_k) %>% summarise(max_k = min(round)))
# info_data_ratio <- info_data %>% filter(round <= max_k)

#log loss ratio
info_data_ratio <- info_data %>% 
  filter(round <= max_k) %>%
  arrange(sample, tree, round) %>%
  group_by(sample, tree) %>%
  mutate(loss_diff_ratio = (lag(loss) - loss)/lag(loss)) %>%
  ungroup()

info_data_ratio %>% select(sample, tree, true_k, round, wunifrac_true, loss_diff_ratio) %>% 
  filter(loss_diff_ratio < 0.1) %>% group_by(sample, tree) %>% summarise() %>% nrow

thresh <- 0.1
merge(info_data_ratio, info_data_ratio %>% filter(loss_diff_ratio < thresh) %>% group_by(sample, tree) %>% 
  summarise(est_k = min(round))) %>%
  group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(round == est_k) %>%
  ggplot(aes(x = wunifrac_true))+
  geom_freqpoly()+
  geom_freqpoly(aes(x = min_error), color = "blue4")+
  geom_vline(aes (xintercept = mean(wunifrac_true)))+
  geom_vline(aes (xintercept = mean(min_error)),color = "blue3")+
  # facet_wrap(~cut(true_k, 10), scales = "free")+
  scale_x_log10()


merge(info_data_ratio, info_data_ratio %>% filter(loss_diff_ratio < thresh) %>% group_by(sample, tree) %>% 
        summarise(est_k = min(round))) %>%
  group_by(sample, tree) %>%
  # mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == est_k) %>%
  ungroup() %>% summarise(mean(wunifrac_true))

merge(info_data_ratio, info_data_ratio %>% filter(loss_diff_ratio < thresh) %>% group_by(sample, tree) %>% 
        summarise(est_k = min(round))) %>%
  group_by(sample, tree) %>%
  ggplot(aes(x = true_k, y = est_k - true_k))+
  geom_hline(yintercept = 0)+
  geom_point(alpha= 0.2)


merge(info_data_ratio, info_data_ratio %>% filter(loss_diff_ratio < thresh) %>% group_by(sample, tree) %>% 
        summarise(est_k = min(round))) %>%
  group_by(sample, tree) %>%
  # mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == est_k) %>%
  ungroup() %>% summarise(mean(est_k))


merge(info_data_ratio, info_data_ratio %>% filter(loss_diff_ratio < thresh) %>% group_by(sample, tree) %>% 
        summarise(est_k = min(round))) %>%
  group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true))
  
  info_data_ratio %>%
  group_by(sample, tree) %>%
  filter(round == knee) %>%
  ungroup() %>%
  ggplot(aes(x = true_k, y = knee))+
  geom_hline(yintercept = 0)+
  geom_point(alpha=0.3)


merge(info_data_ratio, info_data_ratio %>% filter(loss_diff_ratio < 0.1) %>% group_by(sample, tree) %>% 
  summarise(est_k = min(round))) %>%
  group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(round == est_k) %>%
  ungroup() %>% summarise(mean(wunifrac_true))




ggplot(data = info_data, aes(x = true_k, y =  knee - true_k))+
  geom_hline(yintercept = 0)+
  geom_point()


info_data %>% slice_min(order_by = wunifrac_true, n = 1) %>%
ggplot(aes(x = true_k, y =  round - true_k))+
  geom_hline(yintercept = 0)+
  geom_point()

info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == est_k) %>%
  ungroup() %>% summarise(mean(wunifrac_true))

info_data %>% group_by(sample, tree) %>%
             mutate(min_error = min(wunifrac_true)) %>%
            filter(round == est_k) %>%
ggplot(aes(x = wunifrac_true))+
  geom_freqpoly()+
  geom_freqpoly(aes(x = min_error), color = "blue4")+
  geom_vline(aes (xintercept = mean(wunifrac_true)))+
  geom_vline(aes (xintercept = mean(min_error)),color = "blue3")+
  # facet_wrap(~cut(true_k, 10), scales = "free")+
  scale_x_log10()



info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == true_k) %>%
  ggplot(aes(x = wunifrac_true))+
  geom_freqpoly()+
  geom_freqpoly(aes(x = min_error), color = "blue3")+
  geom_vline(aes (xintercept = mean(wunifrac_true)))+
  # geom_text(aes (label = mean(wunifrac_true), y = 100))+
  geom_vline(aes (xintercept = mean(min_error)),color = "blue3")+
  # facet_wrap(~cut(true_k, 10), scales = "free")+
  scale_x_log10()
ggsave("./16k_sim/true_k.png",width = 4,height = 3)

info_data %>% group_by(sample, tree) %>%
  mutate(min_error = min(wunifrac_true)) %>%
  filter(true_k > 1 & round == true_k) %>%
  ungroup() %>% summarise(mean(wunifrac_true))

# kneedle(y = info_data[info_data$sample ==  unique(info_data$sample)[1] & info_data$tree == 5,]$wunifrac_krepp, 
#         x= info_data[info_data$sample ==  unique(info_data$sample)[1] & info_data$tree == 5,]$round)[1]

ggplot(data = info_data_ratio[info_data_ratio$sample %in% unique(info_data_ratio$sample)[10:14] &
                                info_data_ratio$p < 10.1& info_data_ratio$p > -0.03 
                          # info_data$loss * info_data$p >= 1e-5
                        ,], 
  aes(x = round)) +
  geom_vline(aes(xintercept = true_k), color = "grey30", linetype = "dashed")+
  # geom_vline(aes(xintercept = est_k), color = "purple3", linetype = "dashed")+
  geom_point(aes(y=wunifrac_true), color = "green4")+
  geom_line(aes(y=wunifrac_true), color = "green4")+
  geom_point(aes(y=loss_diff_ratio), color = "red4")+
  geom_line(aes(y=loss_diff_ratio), color = "red4")+
  facet_wrap(~p, scales = "free")+
  # geom_hline(yintercept = 0.05)+
  scale_y_log10()
  
  
ggplot(data = info_data[info_data$sample %in% 
                          unique(info_data$sample)[1:100] & 
                          #info_data$p < 0.1& info_data$p >0.01&
                          info_data$round == 2,]) +
    geom_point(aes(y=true_k,x=wunifrac_krepp,color=p))

##smaller trees
all_rounds_10 = read.csv("./16k_sim/all_rounds_data_10k.txt" , sep = "\t", header = F)
names(all_rounds_10) <- c("tree", "true_k", "s", "dist", "round", "wunifrac", "random", "loss", "min_p", "last_p")
head(all_rounds_10)

all_rounds_10_new = read.csv("./16k_sim/all_rounds_data_10k_new.txt" , sep = "\t", header = F)
names(all_rounds_10_new) <- c("tree", "true_k", "s", "dist", "round", "wunifrac", "not_so_random", "loss")
head(all_rounds_10_new)

all_lower_10 = read.csv("./16k_sim/all_lower_bound_10k.txt" , sep = "\t", header = F)
names(all_lower_10) <- c("tree", "true_k", "s", "dist", "lower")
head(all_lower_10)

all_num_10 = read.csv("./16k_sim/all_placement_size_10k.txt" , sep = "\t", header = F)
names(all_num_10) <- c("tree", "true_k", "s", "dist", "placements")
head(all_num_10)

all_rounds_10 <- merge(all_rounds_10, all_rounds_10_new[,c("tree", "true_k", "s", "dist", "round", "not_so_random")])
all_rounds_10 <- merge(all_rounds_10, all_lower_10)
all_rounds_10 <- merge(all_rounds_10, all_num_10)

all_rounds_10$diff <- all_rounds_10$round - all_rounds_10$true_k
all_rounds_10$p <- 1
all_rounds_10$count <- 100000

all_rounds_10$n = 0
all_rounds_10[all_rounds_10$tree == 15, ]$n = 111
all_rounds_10[all_rounds_10$tree == 17, ]$n = 34
all_rounds_10$n <- all_rounds_10$n - all_rounds_10$true_k

all_rounds_10$var <- all_rounds_10$lower * all_rounds_10$count / all_rounds_10$n
all_rounds_10$b2 <- all_rounds_10$loss - (all_rounds_10$var * (all_rounds_10$n - (2*all_rounds_10$round + 1)) / all_rounds_10$count)
all_rounds_10$SNR <- all_rounds_10$count * all_rounds_10$b2 / (all_rounds_10$n * all_rounds_10$var)

# all_rounds_10$thresh <- exp(-3*log(all_rounds_10$n)/all_rounds_10$n)
head(all_rounds_10)

all_rounds_10 <- all_rounds_10 %>%
  arrange(tree, true_k, s, dist, round) %>%
  group_by(tree, true_k, s, dist) %>%
  mutate(loss_ratio = loss / lag(loss)) %>%
  ungroup()
head(all_rounds_10[,c(c("tree", "true_k", "s", "dist", "round", "thresh", "loss_ratio"))])
# View(all_rounds_10[,c("tree", "true_k", "s", "dist", "round", "thresh", "loss_ratio")])


summary(all_rounds %>% filter(is.na(loss_ratio))%>%select(round))

temp <- all_rounds_10[(all_rounds_10$loss - all_rounds_10$lower)* all_rounds_10$count * all_rounds_10$round**2 / all_rounds_10$n < 1,] %>% group_by(tree, true_k, s, dist) %>% summarise(min_k = min(round))
ggplot(data= temp, aes(x=true_k, y=min_k-true_k))+
  geom_hline(yintercept = 0, color = "grey30", linetype = "dashed")+
  facet_wrap(~tree)+
  stat_summary()

ggplot(all_rounds_10[all_rounds_10$round <= all_rounds_10$true_k + 5,], 
       aes(x = round, y = wunifrac /not_so_random , group = round)) +
  geom_vline(aes(xintercept = true_k), color = "grey30", linetype = "dashed")+
  # stat_summary()+
  # stat_summary(geom = "line")+
  geom_boxplot(varwidth = FALSE)+
  geom_boxplot(aes(y=random),color="blue",varwidth = FALSE)+
  geom_boxplot(aes(y=not_so_random),color="red",varwidth = FALSE)+
  theme_bw()+
  # scale_color_viridis_b(n.breaks=10)+
  # stat_smooth(se=F)+
  geom_hline(yintercept = 1e-3)+
  facet_wrap(~true_k,ncol=5)+
  scale_y_log10()

ggplot(all_rounds_10[all_rounds_10$round <= all_rounds_10$true_k + 5,], 
       aes(x = round, y = SNR, group = round)) +
  geom_vline(aes(xintercept = true_k), color = "grey30", linetype = "dashed")+
  # stat_summary()+
  # stat_summary(geom = "line")+
  geom_boxplot(varwidth = FALSE)+
  # geom_boxplot(aes(y=lower),color="blue",varwidth = FALSE)+
  theme_bw()+
  # scale_color_viridis_b(n.breaks=10)+
  # stat_smooth(se=F)+
  geom_hline(yintercept = 1e+9)+
  facet_wrap(~tree+true_k)+
  scale_y_log10()#+
  coord_cartesian(ylim = c(1e-3, 1e+2))
  
  
# how loss/lower is related to n, k, and p

temp_10 <- all_rounds_10[all_rounds_10$diff == 0,c("n", "true_k", "p", "loss", "lower")]
temp <- all_rounds[all_rounds$diff == 0,c("n", "true_k", "p", "loss", "lower")]

temp <- rbind(temp_10, temp)

head(temp)

ggplot(temp, 
       aes(x = n, y = loss * p/n)) +
  # geom_vline(aes(xintercept = true_k), color = "grey30", linetype = "dashed")+
  # stat_summary()+
  # stat_summary(geom = "line")+
  geom_point(color="grey30", alpha=0.5)+
  # geom_boxplot(aes(y=lower),color="blue",varwidth = FALSE)+
  theme_bw()+
  # scale_color_viridis_b(n.breaks=10)+
  # stat_smooth(se=F)+
  # geom_hline(yintercept = 1e+10)+
  # facet_wrap(~cut(true_k, 5))+
  scale_y_log10()#+
coord_cartesian(ylim = c(1,1e+14))
