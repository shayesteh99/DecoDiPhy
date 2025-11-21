library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(ggpubr)

our_method = "DecoDiPhy"

est_mash_data <- read.csv("./10k_dataset/all_metrics_est_mash.txt", header = F, sep = "\t")
names(est_mash_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
est_mash_data$method <- "sequence distances"
head(est_mash_data)

# est_mash_f1_data <- read.csv("./10k_dataset/all_metrics_est_mash_f1.txt", header = F, sep = "\t")
# names(est_mash_f1_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# est_mash_f1_data$method <- "est-mash_f1"
# head(est_mash_f1_data)

est_placement_data <- read.csv("./10k_dataset/all_metrics_krepp_placements.txt", header = F, sep = "\t")
names(est_placement_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
est_placement_data$method <- "placements"
head(est_placement_data)

est_multiplacement_data <- read.csv("./10k_dataset/all_metrics_krepp_multiplacements.txt", header = F, sep = "\t")
names(est_multiplacement_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
est_multiplacement_data$method <- "multi-placements"
head(est_multiplacement_data)

# est_wmultiplacement_data <- read.csv("./10k_dataset/all_metrics_krepp_wmultiplacements.txt", header = F, sep = "\t")
# names(est_wmultiplacement_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# est_wmultiplacement_data$method <- "weighted_multiplacement"
# head(est_wmultiplacement_data)

# est_placement_0.1_data <- read.csv("./10k_dataset/all_metrics_krepp_placements_1e-1.txt", header = F, sep = "\t")
# names(est_placement_0.1_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# est_placement_0.1_data$method <- "placement"
# head(est_placement_0.1_data)

# est_placement_fixed_data <- read.csv("./10k_dataset/all_metrics_krepp_placements_f1.txt", header = F, sep = "\t")
# names(est_placement_fixed_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# est_placement_fixed_data$method <- "placement_f1"
# head(est_placement_fixed_data)

# est_multiplacement_fixed_data <- read.csv("./10k_dataset/all_metrics_krepp_multiplacements_f1.txt", header = F, sep = "\t")
# names(est_multiplacement_fixed_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# est_multiplacement_fixed_data$method <- "multiplacement_f1"
# head(est_multiplacement_fixed_data)

# est_woltka_placement_data <- read.csv("./10k_dataset/all_metrics_woltka_placements.txt", header = F, sep = "\t")
# names(est_woltka_placement_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# est_woltka_placement_data$method <- "woltka_placement"
# head(est_woltka_placement_data)

true_mash_data <- read.csv("./10k_dataset/all_metrics_true_mash.txt", header = F, sep = "\t")
names(true_mash_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
true_mash_data$method <- "true distances"
head(true_mash_data)

# krep_multiplacement_data <- read.csv("./10k_dataset/all_krepp_multiplacement_metrics.txt", header = F, sep = "\t")
# names(krep_multiplacement_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# krep_multiplacement_data$method <- "krepp_multi"
# head(krep_multiplacement_data)

# krep_wmultiplacement_data <- read.csv("./10k_dataset/all_krepp_wmultiplacement_metrics.txt", header = F, sep = "\t")
# names(krep_wmultiplacement_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# krep_wmultiplacement_data$method <- "krepp_weighted_multi"
# head(krep_wmultiplacement_data)

# voltka_mash_data <- read.csv("./10k_dataset/all_metrics_voltka_mash.txt", header = F, sep = "\t")
# names(voltka_mash_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# # voltka_data$jaccard <- 0
# # voltka_data <- voltka_data[,c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")]
# voltka_mash_data$method <- "woltka_mash"
# head(voltka_mash_data)

# all_data <- rbind(est_data[, setdiff(names(est_data), c("jaccard"))], true_data[, setdiff(names(true_data), c("jaccard"))], voltka_data)
# all_data_mash <- rbind(est_mash_data, true_mash_data, voltka_mash_data, est_multiplacement_data, est_wmultiplacement_data, est_placement_0.1_data, est_woltka_placement_data, krep_multiplacement_data, krep_wmultiplacement_data)
all_data_mash <- rbind(est_mash_data, true_mash_data, est_multiplacement_data, est_placement_data)

nrow(all_data_mash)
# all_data_mash <- merge(est_mash_data[,c("tree", "k", "s", "dist")], all_data_mash)
# all <- rbind(all_data, all_data_mash)

head(all_data_mash)
# View(all_data_mash[all_data_mash$tree %in% c(15, 17),] %>% group_by(tree, k, s, dist) %>% summarise(count = n()))

# View(all_data_mash[all_data_mash$tree == 15 & all_data_mash$k == 10 & all_data_mash$method %in% c("placement", "est-mash"),] %>% group_by(tree, k, s, dist))
  
all_data_mash$n = 0
all_data_mash[all_data_mash$tree == 7, ]$n = 50
all_data_mash[all_data_mash$tree == 4, ]$n = 89
all_data_mash[all_data_mash$tree == 15, ]$n = 111
all_data_mash[all_data_mash$tree == 17, ]$n = 34
all_data_mash[all_data_mash$tree == 1322, ]$n = 30

all_mash_metrics <- all_data_mash[all_data_mash$tree %in% c(15, 17, 1322) & all_data_mash$n / 3 > all_data_mash$k,] %>%
  #mutate(`delta k`=est_k - k) %>%
  mutate(`delta k`=(est_k - k)/k) %>%
  # mutate(bc=bc/2) %>%
  pivot_longer(cols = c("jaccard", "unifrac", "wunifrac", "delta k", "bc"), names_to = "metric")
all_mash_metrics$tree <- factor(all_mash_metrics$tree, levels = c("15", "17", "1322"), labels = c("T1 (n=111)", "T2 (n=34)", "T3 (n=30)"))
all_mash_metrics$metric <- factor(all_mash_metrics$metric, levels = c("jaccard", "wunifrac","unifrac","delta k", "bc"), labels = c("Jaccard", "wUniFrac","UniFrac", "delta k", "Bray-Curtis"))
all_mash_metrics$method <- factor(all_mash_metrics$method, levels = c("true distances", "multi-placements", "placements","sequence distances"))
                       # , labels = c(paste0("true distances + ", our_method), paste0("krepp + ", our_method), "krepp + Woltka"))
all_mash_metrics$k <- factor(all_mash_metrics$k, levels = c(3,5,10), labels = c("k=3","k=5","k=10"))

ggplot(data = all_mash_metrics[
                                 all_mash_metrics$tree != "subtree 3 (n=30)"
                               # & ! all_mash_metrics$method %in% c("woltka_placement", "placement_f1", "multiplacement_f1", "krepp_weighted_multi", "weighted_multiplacement")
                               # & all_mash_metrics$method %in% c("woltka_mash", "multiplacement", "true_mash", "krepp_multi")
                               ,], aes(x = as.factor(tree), y = (value), 
                                       color = method, fill= method, group = as.factor(method)))+
  stat_summary(geom = "col", width = 0.6, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4, color="black")+
  theme_bw()+
  facet_grid(metric~k, scales = "free")+
  scale_color_manual(values = c("#1b9e77", "#b3de69", "#7570b3", "#d95f02"), guide="none"
                     # labels = c("true distances + DecoDiPhy", "krepp + DecoDiPhy", "krepp + Woltka", "krepp")
  ) +
  scale_fill_manual(values = c("#1b9e77", "#b3de69", "#7570b3", "#d95f02"), name = ""
                     # labels = c("true distances + DecoDiPhy", "krepp + DecoDiPhy", "krepp + Woltka", "krepp")
  ) +
  labs(x = "", y= "")+
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, vjust = 0.6))+
  guides(color = guide_legend(nrow = 1, title = ""))
ggsave(paste0("./Plots/10k_placement_bar_allk.pdf"),width = 7,height = 7)
ggsave(paste0("./Plots/10k_placement_bar.pdf"),width = 12.5/1.2188582011,height = 3)


ggplot(data = all_mash_metrics[all_mash_metrics$metric != "Jaccard" &
                                all_mash_metrics$tree != "subtree 3 (n=30)"
                               # & ! all_mash_metrics$method %in% c("woltka_placement", "placement_f1", "multiplacement_f1", "krepp_weighted_multi", "weighted_multiplacement")
                               # & all_mash_metrics$method %in% c("woltka_mash", "multiplacement", "true_mash", "krepp_multi")
                               ,], aes(x = as.factor(k), y = (value), 
                                    color = method, group = as.factor(method)))+
  # stat_summary(fun.data =median_iqr,position = position_dodge(width=0.4))+
  stat_summary()+
  stat_summary(geom = "line")+
  geom_hline(yintercept = 0, color = "grey50", linetype="dashed")+
  # geom_boxplot(aes(group=interaction(k,method)))+
  theme_bw()+
  facet_grid(metric~tree, scales = "free")+
  scale_color_manual(values = c("#1b9e77", "#b3de69", "#7570b3", "#d95f02"),
                     # labels = c("true distances + DecoDiPhy", "krepp + DecoDiPhy", "krepp + Woltka", "krepp")
                     ) +
  # scale_color_manual(values = c("#1b9e77", "#d95f02", "#80b1d3", "#b3de69", "#7570b3", "red4")) +
  # scale_y_log10()+
  labs(x = "number of query taxa", y= "")+
  theme(legend.position = "bottom")+
  guides(color = guide_legend(nrow = 1, title = ""))
ggsave(paste0("./Plots/10k_res_new_placement.pdf"),width = 6,height = 5)

ggplot(data = all_mash_metrics[
  all_mash_metrics$tree != "subtree 3 (n=30)"
  & ! all_mash_metrics$method %in% c("placement_f1", "multiplacement_f1", "krepp_weighted_multi", "weighted_multiplacement")
  & all_mash_metrics$method %in% c("multiplacement", "true_mash", "krepp_multi", "woltka_placement", "est-mash")
  ,], aes(x = metric, y = (value), 
          fill = method, group = as.factor(method)))+
  # stat_summary(fun.data =median_iqr,position = position_dodge(width=0.4))+
  stat_summary(geom = "errorbar",position = position_dodge(width=0.7),color="black",width=0.5)+
  stat_summary(geom = "bar",position = position_dodge(width=0.7),color="black",width=0.7)+
  #geom_hline(yintercept = 0, color = "grey50", linetype="dashed")+
  # geom_boxplot(aes(group=interaction(k,method)))+
  theme_bw()+
  #facet_grid(metric~tree, scales = "free")+
  scale_fill_manual(values = c("#7570b3","#1b9e77","#b3de69", "#d95f02", "red4"),name="",
                   labels = c(expression(PDD(bold(d))), expression(PDD(hat(bold(d))[1])), 
                              expression(PDD(hat(bold(d))[2])), "Closed Ref.\n(Woltka)", "Read Phylo.\nPlace.(krepp)", "dist")
                   ) +
  # scale_color_manual(values = c("#1b9e77", "#d95f02", "#80b1d3", "#b3de69", "#7570b3", "red4")) +
  # scale_y_log10()+
  labs(x = "Error metric", y= "Error value", color = "method")+
  theme(legend.position = "right")

ggsave("simMetagenomeWide.pdf",width=6,height = 4.25)


all_mash_metrics[
  all_mash_metrics$tree != "subtree 3 (n=30)"
  & ! all_mash_metrics$method %in% c("placement_f1", "multiplacement_f1", "krepp_weighted_multi", "weighted_multiplacement")
  & all_mash_metrics$method %in% c("multiplacement", "true_mash", "krepp_multi", "woltka_placement")
  ,] %>%
  group_by(method,k) %>%
  summarise(e=mean(est_k/k))

all_mash_metrics[all_mash_metrics$tree != "subtree 3 (n=30)",] %>% group_by(tree,method) %>% summarise(n=n())

all_mash_metrics <- all_data_mash[all_data_mash$tree %in% c(7, 4, 15) & all_data_mash$n >= 5*all_data_mash$k,] %>%
  mutate(`delta k`=est_k - k) %>%
  pivot_longer(cols = c("unifrac", "wunifrac", "delta k"), names_to = "metric")
all_mash_metrics$tree <- factor(all_mash_metrics$tree, levels = c("7", "4", "15"), labels = c("small (n=50)", "medium (n=89)", "large (n=111)"))
all_mash_metrics$metric <- factor(all_mash_metrics$metric, levels = c("delta k", "unifrac", "wunifrac"), labels = c("delta k", "UniFrac", "weighted UniFrac"))
all_mash_metrics$method <- factor(all_mash_metrics$method, levels = c("true_mash", "est-mash", "voltka_mash"), labels = c(paste0("true distances + ", our_method), paste0("krepp + ", our_method), "krepp + Woltka"))

ggplot(data = all_mash_metrics, aes(x = as.factor(k), y = (value), 
                                    color = method, group = as.factor(method)))+
  stat_summary()+
  stat_summary(geom = "line")+
  # geom_boxplot(aes(group=interaction(k,method)))+
  theme_bw()+
  facet_grid(metric~tree, scales = "free")+
  scale_color_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
  # scale_y_log10()+
  labs(x = "k", y= "", color = "method")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/10k_res_size.pdf"),width = 6,height = 5)


all_mash_metrics <- all_data_mash[all_data_mash$tree %in% c(15, 17, 1322) & all_data_mash$n / 3 > all_data_mash$k,] %>%
  mutate(jaccard_comp=1-jaccard) %>%
  pivot_longer(cols = c("jaccard_comp", "bc"), names_to = "metric")
all_mash_metrics$tree <- factor(all_mash_metrics$tree, levels = c("15", "17", "1322"), labels = c("subtree 1 (n=111)", "subtree 2 (n=34)", "subtree 3 (n=30)"))
all_mash_metrics$metric <- factor(all_mash_metrics$metric, levels = c("jaccard_comp", "bc"), labels = c("1 - Jaccard", "Bray-Curtis"))
all_mash_metrics$method <- factor(all_mash_metrics$method, levels = c("true_mash", "est-mash", "voltka_mash"), labels = c(paste0("true distances + ", our_method), paste0("krepp + ", our_method), "krepp + Woltka"))

ggplot(data = all_mash_metrics, aes(x = as.factor(k), y = (value), 
                                    color = method, group = as.factor(method)))+
  stat_summary()+
  stat_summary(geom = "line")+
  # geom_boxplot(aes(group=interaction(k,method)))+
  theme_bw()+
  facet_grid(metric~tree, scales = "free")+
  scale_color_manual(values = c("#1b9e77", "#7570b3", "#d95f02")) +
  # scale_y_log10()+
  labs(x = "k", y= "", color = "method")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/10k_supp_bc.pdf"),width = 6,height = 4)





ggplot(data = all_data, aes(x = as.factor(k), y = (est_k - k), color = as.factor(method), group = as.factor(method)))+
  geom_hline(yintercept = 0, color = "gray40", linetype = "dashed")+
  stat_summary()+
  stat_summary(geom = "line")+
  theme_bw()+
  facet_wrap(~tree, nrow = 1)+
  scale_color_brewer(palette = "Dark2") +
  # scale_y_log10()+
  labs(x = "k", y = "estk - true k", color = "method")+
  theme(legend.position = "bottom")
ggsave("./Plots/10k_est_k.pdf",width = 4,height = 3)

ggplot(data = all_data, aes(x = as.factor(k), y = (bc), color = as.factor(method), group = as.factor(method)))+
  stat_summary()+
  stat_summary(geom = "line")+
  theme_bw()+
  facet_wrap(~tree, nrow = 1)+
  scale_color_brewer(palette = "Dark2") +
  # scale_y_log10()+
  labs(x = "k", y = "bray-curtis", color = "method")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/10k_wunifrac.pdf"),width = 4,height = 3)


d <- read.csv("../newdist.txt",)
head(d)

View(d[d$V1 == 1322 & d$V2 == 5 & d$V3 == 1 & d$V4 == 'exp', ])

##draw distances

true_distances <- read.csv("./10k_dataset/true_distances.txt", header = F, sep = "\t")
names(true_distances) <- c("tree", "k", "s", "dist", "leaf", "true_distance")
head(true_distances)
unique(true_distances$tree)

# krepp_placement_distances <- read.csv("./10k_dataset/krepp_placement_distances.txt", header = F, sep = "\t")
# names(krepp_placement_distances) <- c("tree", "k", "s", "dist", "leaf", "distance")
# krepp_placement_distances$method <- "krepp_placement"
# head(krepp_placement_distances)

krepp_multiplacement_distances <- read.csv("./10k_dataset/krepp_multiplacement_distances.txt", header = F, sep = "\t")
names(krepp_multiplacement_distances) <- c("tree", "k", "s", "dist", "leaf", "distance")
krepp_multiplacement_distances$method <- "krepp_multiplacement"
head(krepp_multiplacement_distances)

# woltka_placement_distances <- read.csv("./10k_dataset/woltka_placement_distances.txt", header = F, sep = "\t")
# names(woltka_placement_distances) <- c("tree", "k", "s", "dist", "leaf", "distance")
# woltka_placement_distances$method <- "woltka_placement"
# head(woltka_placement_distances)

krepp_distances <- read.csv("./10k_dataset/krepp_distances.txt", header = F, sep = "\t")
names(krepp_distances) <- c("tree", "k", "s", "dist", "leaf", "distance", "portion", "true_dist")
krepp_distances <- krepp_distances[krepp_distances$true_dist != "TRUE_DISTANCE",]
krepp_distances$method <- "krepp"
# krepp_distances <- krepp_distances[krepp_distances$portion > 0.3,]
head(krepp_distances)
nrow(krepp_distances)

all_krep <- rbind( 
                  krepp_multiplacement_distances,
                  # woltka_placement_distances, 
                  krepp_distances[,c("tree", "k", "s", "dist", "leaf", "distance", "method")])

all_dist <- merge(true_distances, all_krep)
all_dist$distance <- as.numeric(all_dist$distance)
head(all_dist)

nrow(all_dist[all_dist$method == "krepp",])
nrow(all_dist[all_dist$method == "krepp_placement",])
nrow(all_dist[all_dist$method == "krepp_multiplacement",])
nrow(all_dist[all_dist$method == "woltka_placement",])

all_dist %>% group_by(tree,k,s,dist,method) %>%
  summarise(R = cor(true_distance, distance, use="na.or.complete")) %>%
  group_by(method) %>%
  summarise(mean_R = mean(R))

all_dist$tree <- factor(all_dist$tree, levels = c("15", "17"), labels = c("T1", "T2"))
all_dist$method <- factor(all_dist$method, levels = c("krepp_multiplacement", "krepp"), label=c("multi-placements","sequence distances"))
all_dist$k <- factor(all_dist$k, levels = c(3,5,10), labels = c("k=3","k=5","k=10"))
ggplot(data = all_dist, aes(x=true_distance, y = distance, color = method))+
  geom_abline(slope = 1, intercept = 0, color = "grey30", linetype = "dashed")+
  geom_point(alpha = 0.5)+
  facet_grid(s ~ tree+k)+
  stat_cor(aes(label = paste(..r.label.., sep = "~`,`~")), method = "pearson")+
  theme_bw()+
  scale_color_brewer(palette = "Dark2", name="krepp")+
  theme(legend.position = "bottom")+
  labs(y="estimated", x="true")
ggsave("./Plots/corr_distances_10k.pdf",width = 10,height = 8)

ggplot(data = krepp_distances, aes(x=as.numeric(true_dist), y = as.numeric(distance), color = as.numeric(portion)))+
  geom_abline(slope = 1, intercept = 0, color = "grey30", linetype = "dashed")+
  geom_point(alpha = 0.8)+
  facet_grid(s ~ tree+k)+
  stat_cor(aes(label = paste(..r.label.., sep = "~`,`~")), method = "pearson")+
  scale_color_viridis_c()+
  theme_bw()
ggsave("./Plots/krepp_ditances.pdf",width = 10,height = 8)

krepp_distances[as.numeric(krepp_distances$true_dist) - as.numeric(krepp_distances$distance) > 0.1 & krepp_distances$portion > 0.2,]
krepp_distances[is.na(as.numeric(krepp_distances$distance)),]


##results per k
est_placement_per_k <- read.csv("./10k_dataset/all_per_k_metrics_krepp_placements_k30.txt", header = F, sep = "\t")
names(est_placement_per_k) <- c("tree", "k", "s", "dist", "round", "jaccard", "unifrac", "wunifrac")
est_placement_per_k$method <- "placement+decodiphy"
head(est_placement_per_k)

est_woltka_placement_per_k <- read.csv("./10k_dataset/all_per_k_metrics_woltka_placements_k30.txt", header = F, sep = "\t")
names(est_woltka_placement_per_k) <- c("tree", "k", "s", "dist", "round", "jaccard", "unifrac", "wunifrac")
est_woltka_placement_per_k$method <- "woltka placement+decodiphy"
head(est_woltka_placement_per_k)

est_multiplacement_per_k <- read.csv("./10k_dataset/all_per_k_metrics_krepp_multiplacements_k30.txt", header = F, sep = "\t")
names(est_multiplacement_per_k) <- c("tree", "k", "s", "dist", "round", "jaccard", "unifrac", "wunifrac")
est_multiplacement_per_k$method <- "multiplacement+decodiphy"
head(est_multiplacement_per_k)

est_wmultiplacement_per_k <- read.csv("./10k_dataset/all_per_k_metrics_krepp_wmultiplacements_k30.txt", header = F, sep = "\t")
names(est_wmultiplacement_per_k) <- c("tree", "k", "s", "dist", "round", "jaccard", "unifrac", "wunifrac")
est_wmultiplacement_per_k$method <- "weighted_multiplacement+decodiphy"
head(est_wmultiplacement_per_k)

est_mash_per_k <- read.csv("./10k_dataset/all_per_k_metrics_est_mash_k30.txt", header = F, sep = "\t")
names(est_mash_per_k) <- c("tree", "k", "s", "dist", "round", "jaccard", "unifrac", "wunifrac")
est_mash_per_k$method <- "distances+decodiphy"
head(est_mash_per_k)

woltka_per_k <- read.csv("./10k_dataset/all_per_k_voltka_metrics.txt", header = F, sep = "\t")
names(woltka_per_k) <- c("tree", "k", "s", "dist", "round", "jaccard", "unifrac", "wunifrac")
woltka_per_k$method <- "woltka"
head(woltka_per_k)

all_per_k <- rbind(woltka_per_k, est_woltka_placement_per_k, est_placement_per_k, est_multiplacement_per_k, est_wmultiplacement_per_k, est_mash_per_k)

ggplot(data = all_per_k, aes(x = round, y = wunifrac, color = method))+
  geom_vline(aes(xintercept = k), color = "grey30", linetype= "dashed")+
  geom_smooth(aes(group = method))+
  # geom_line(aes(group = interaction(s,dist)))+
  facet_grid(tree~k, scales = "free")

ggplot(data = all_per_k, aes(x = round, y = wunifrac, color = method))+
  geom_vline(aes(xintercept = k), color = "grey30", linetype= "dashed")+
  stat_summary(aes(group = method), geom = "line")+
  # geom_line(aes(group = method))+
  facet_grid(tree~k, scales = "free")+
  coord_cartesian(ylim = c(0, 0.3), xlim = c(0,20))
