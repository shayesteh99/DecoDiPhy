library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)
library(rdirichlet)
library(MCMCpack)

#decodiphy results
unifrac_data <- read.csv("./biotrees/all_unifrac_res.txt", header= F, sep = "\t")
names(unifrac_data) <- c("data", "k", "s", "noise", "scale", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc", "method")
head(unifrac_data)
unique(unifrac_data$method)

unifrac_data %>% filter(method=="f1_exhaustive", scale == 0) %>% summary(jaccard)
unifrac_data %>% filter(scale==0, k>4) %>% group_by(method) %>% summarise(mean(jaccard))

runtime_bio <- read.csv("./biotrees/all_runtime_res.txt", header = F, sep = "\t")
names(runtime_bio) <- c("data", "k", "s", "noise", "seq_length", "n", "est_k", "total", "one_round")
head(runtime_bio)

runtime_bio$total_n <- runtime_bio$k + runtime_bio$n
data_size <- runtime_bio %>% group_by(data) %>% summarise(n = mean(total_n))

# exh_data <- rbind(unifrac_data_q1_r2, unifrac_data_f1_q1_r2, unifrac_data_baseline, unifrac_data_q1_r2_perm_rest, unifrac_data_q1_r2_perm_rest_1e1)
# exh_data$seq_length <- factor(exh_data$seq_length, levels = c("0", "10000", "1000", "100"))
# exh_data$method <- factor(exh_data$method, levels = c("fast_r2_fixed", "fast_r2_perm_rest", "fast_r2_perm_rest_1e-1", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2", "baseline"))
                            # c("f1_exh", "f1", "standard", "fast_r2_fixed", "fast_r2", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest", "baseline"))

combined <- unifrac_data %>% mutate(
  delta_k=(est_k-k),
  "j_comp"=1-jaccard
) %>% pivot_longer(cols=c(j_comp,jaccard, wunifrac, unifrac, bc, delta_k)) 

combined <- merge(combined, data_size)
combined$name <- factor(combined$name, levels = c("j_comp", "bc",  "unifrac", "wunifrac", "delta_k","jaccard"), labels=c("1 - Jaccard", "Bray-Curtis", "UniFrac", "wUniFrac",  "delta k","Jaccard"))
# combined$method <- factor(combined$method, levels = c("f1_exh", "f1", "standard", "f1_reverse_hill", "reverse_hill"))
combined$method <- factor(combined$method, levels = c("f1_exhaustive", "q0_f1", "q0_minp_count", "q1_r2_minp_count", "f1_reverse_hill", "f1_closest_iterative", "f1_closest"), labels = c("exact-exact", "exact-heu", "heu-heu(complete)", "heu-heu(fast-default)", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves"))

unique(combined$data)


ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                       & combined$method %in% c("exact-exact", "exact-heu", "heu-heu(complete)", "heu-heu(fast-default)")
                       # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                       # & combined$scale < 4
                       # & combined$data == "birds-stiller"
                       # & combined$name %in% c("delta_k")
                       ,]
       , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_grid(name~k,scale="free", labeller = labeller(k = function(x) paste0("k = ", x)))+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#b3de69", "#8dd3c7", "#80b1d3", "#fdb462", "#fb8072", "purple", "red4"))+
  scale_x_discrete(labels = c("0" = "0", "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  labs(x = "noise (exp scale)", y = "", fill = "method(k-search)") +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))
ggsave(paste0("./Plots/bio_comp_only_decodiphy.pdf"),width = 10,height = 5.5)

library(cowplot)

a=ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                       & combined$method %in% c("exact-exact", "exact-heu", "heu-heu(complete)", "heu-heu(fast-default)")
                       # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                       # & combined$scale < 4
                       # & combined$data == "birds-stiller"
                        & !combined$name %in% c("delta k","UniFrac","Bray-Curtis","1 - Jaccard","Jaccard")
                       ,]
       , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_wrap(.~ifelse(k<=3,"k<4","k>4"),scale="free", labeller = labeller(k = function(x) paste0("k = ", x)),nrow=1)+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#b3de69", "#8dd3c7", "#80b1d3", "#fdb462", "#fb8072", "purple", "red4"))+
  scale_x_discrete(name = "noise (mean edges)",
                   labels = c("0" = "0", "2" = 1.56, "1" = "0.6", "100" = expression(10^-2))) +
  labs(y = "Weighted UniFrac", fill = "method(k-search)") +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))

b=ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                         & combined$method %in% c("exact-exact", "exact-heu", "heu-heu(complete)", "heu-heu(fast-default)")
                         # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                         # & combined$scale < 4
                         # & combined$data == "birds-stiller"
                         & !combined$name %in% c("delta k","UniFrac","Bray-Curtis","1 - Jaccard","wUniFrac")
                         ,]
         , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_wrap(.~ifelse(k<=3,"k<4","k>4"),scale="free", labeller = labeller(k = function(x) paste0("k = ", x)),nrow=1)+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#b3de69", "#8dd3c7", "#80b1d3", "#fdb462", "#fb8072", "purple", "red4"))+
  scale_x_discrete(name = "noise (mean edges)",
    labels = c("0" = "0", "2" = 1.56, "1" = "0.6", "100" = expression(10^-2))) +
  labs(y = "Jaccard", fill = "method(k-search)") +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))
b


c=ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                         & combined$method %in% c("exact-heu", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves")
                         # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                         # & combined$scale < 4
                         # & combined$data == "birds-stiller"
                         & !combined$name %in% c("delta k","UniFrac","Bray-Curtis","1 - Jaccard","Jaccarda")
                         ,]
         , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_wrap(.~name,scale="free", 
             labeller = labeller(name=c(wUniFrac="Weighted UniFrac",Jaccard="Jaccard")),nrow=1)+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#8dd3c7", "#fdb462", "#fb8072", "#c2a5cf"), labels =c("DecodiPhy", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves"))+
  scale_x_discrete(name = "noise (mean edges)",
                   labels = c("0" = "0", "2" = 1.56, "1" = "0.6", "100" = expression(10^-2))) +
  labs(y = "", fill = "method(k-search)") +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))
c
plot_grid(b,a,c,nrow=1,labels=c("a)","b)","c)"))
ggsave("main-e1.pdf",width=9,height = 2.5)


b=ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                       & combined$method %in% c("exact-heu", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves","exact-exact", "exact-heu", "heu-heu(complete)", "heu-heu(fast-default)")
                       # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                       # & combined$scale < 4
                       # & combined$data == "birds-stiller"
                       & !combined$name %in% c("delta k","UniFrac","Bray-Curtis","1 - Jaccard","Jaccard")
                       ,]
       , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_wrap(.~ifelse(k<=3,"k<4","k>4"),scale="free_x", 
             labeller = labeller(name=c(wUniFrac="Weighted UniFrac",Jaccard="Jaccard")),nrow=1)+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#b3de69", "#8dd3c7", "#80c1f3", "#80b1d3", "#fdb462", "#fb8072", "#c2a5cf"), 
                    #labels =c("DecodiPhy", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves")
                    )+
  scale_x_discrete(name = "noise (mean edges)",
                   labels = c("0" = "0", "2" = 1.56, "1" = "0.6", "100" = expression(10^-2))) +
  labs(y = "Weighted UniFrac", fill = "method(k-search)") +
  theme(legend.position = "bottom",,panel.grid.major.x = element_blank())+
  guides(fill = guide_legend(nrow = 1))

a=ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                         & combined$method %in% c("exact-heu", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves","exact-exact", "exact-heu", "heu-heu(complete)", "heu-heu(fast-default)")
                         # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                         # & combined$scale < 4
                         # & combined$data == "birds-stiller"
                         & !combined$name %in% c("delta k","UniFrac","Bray-Curtis","1 - Jaccard","wUniFrac")
                         ,]
         , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_wrap(.~ifelse(k<=3,"k<4","k>4"),scale="free_x", 
             labeller = labeller(name=c(wUniFrac="Weighted UniFrac",Jaccard="Jaccard")),nrow=1)+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#b3de69", "#8dd3c7", "#80c1f3", "#80b1d3", "#fdb462", "#fb8072", "#c2a5cf"), 
                    #labels =c("DecodiPhy", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves")
  )+
  scale_x_discrete(name = "noise (mean edges)",
                   labels = c("0" = "0", "2" = 1.56, "1" = "0.6", "100" = expression(10^-2))) +
  labs(y = "Jaccard", fill = "method(k-search)") +
  theme(legend.position = "bottom",panel.grid.major.x = element_blank())+
  guides(fill = guide_legend(nrow = 1))
a

plot_grid(a,b,nrow=1,labels=c("a)","b)"))
ggsave("main-e1.pdf",width=8,height = 2.7)


b=ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                         & combined$method %in% c("exact-heu", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves")
                         # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                         # & combined$scale < 4
                         # & combined$data == "birds-stiller"
                         & !combined$name %in% c("delta k","UniFrac","Bray-Curtis","1 - Jaccard","wUniFrac")
                         ,]
         , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_wrap(.~ifelse(k<=3,"k<4","k>4"),scale="free", labeller = labeller(k = function(x) paste0("k = ", x)),nrow=1)+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#8dd3c7", "#fdb462", "#fb8072", "#c2a5cf"), labels =c("DecodiPhy", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves"))+
  scale_x_discrete(name = "noise (mean edges)",
                   labels = c("0" = "0", "2" = 1.56, "1" = "0.6", "100" = expression(10^-2))) +
  labs(y = "Jaccard", fill = "method(k-search)") +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))
b

plot_grid(b,a,nrow=1,labels=c("a)","b)"))
ggsave("main-e1-base.pdf",width=7,height = 2.5)


summary(pmin(as.integer(rexp(n = 10000, rate = 1/2))))

ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k & combined$scale %in% c(0,1,2)
                       & combined$method %in% c("exact-heu", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves")
                       & combined$name != "delta k",
                       # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                       # & combined$scale < 4
                       # & combined$data == "birds-stiller"
                       # & combined$name %in% c("delta_k")
                       ,]
       , aes(x = as.factor(scale), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_grid(name~k,scale="free", labeller = labeller(k = function(x) paste0("k = ", x)))+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#8dd3c7", "#fdb462", "#fb8072", "#c2a5cf"), labels =c("DecodiPhy", "reverse_search", "k-nearest leaves (iterative)", "k-nearest leaves"))+
  scale_x_discrete(labels = c("0" = "0", "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  labs(x = "noise (exp scale)", y = "", fill = "method(correct k)") +
  theme(legend.position = "bottom")+
  guides(fill = guide_legend(nrow = 1))
ggsave(paste0("./Plots/bio_comp_others.pdf"),width = 10,height = 5.5)



ggplot(data = combined[combined$k <= 10 & combined$n /3 > combined$k 
                       & !combined$method %in% c("fast_r2", "baseline", "fast_r2_perm_rest")
                       # & combined$method %in% c("fast_r2_fixed", "fast_r2_perm_last", "fast_r2_perm_min", "fast_r2_perm_rest")
                       # & combined$scale < 4
                       # & combined$data == "birds-stiller"
                       # & combined$name %in% c("delta_k")
                       ,]
       , aes(x = (1/scale), y = value, linetype = as.factor(method), color=as.factor(k), group = interaction(k,method))) +
  stat_summary(geom = "line", position = position_dodge(width = 0.05))+
  stat_summary(position = position_dodge(width = 0.05), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_grid(name~.,scale="free", labeller = labeller(k = function(x) paste0("k = ", x)))+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#b3de69", "#8dd3c7", "#80b1d3", "#fdb462", "#fb8072", "purple"))+
  #scale_x_discrete(labels = c("0" = "0", "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  labs(x = "noise (exp scale)", y = "", fill = "method(k-search)") +
  theme(legend.position = "bottom")

#distances
all_distances <- read.csv("./Demixing_Data/biotrees/all_dist.txt", header = F, sep = "\t")
names(all_distances) <- c("data", "k", "s", "rate", "leaf", "true", "noisy")
head(all_distances)

all_distances <- read.csv("./Demixing_Data/biotrees/all_dist_new.txt", header = F, sep = "\t")
names(all_distances) <- c("data", "k", "s", "noise", "rate", "leaf", "dist")
all_distances <- merge(all_distances, all_distances[all_distances$noise == 0,c("data", "k", "s", "leaf", "dist")], by = c("data", "k", "s", "leaf"))
names(all_distances) <- c("data", "k", "s", "leaf", "noise", "rate", "true", "noisy")
head(all_distances)

ggplot(data = all_distances, aes(x=true, y = noisy))+
  geom_abline(slope = 1, intercept = 0, color = "grey30", linetype = "dashed")+
  geom_point(alpha = 0.5)+
  facet_grid(data+k ~ rate)+
  stat_cor(aes(label = paste(..r.label.., sep = "~`,`~")), method = "pearson")+
  theme_bw()+
  theme(legend.position = "bottom")



