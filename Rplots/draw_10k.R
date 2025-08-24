our_method = "DecoDiPhy"

est_data <- read.csv("./10k_dataset/all_metrics_est.txt", header = F, sep = "\t")
names(est_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
est_data$method <- "est"
head(est_data)

est_mash_data <- read.csv("./10k_dataset/all_metrics_est_mash.txt", header = F, sep = "\t")
names(est_mash_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
est_mash_data$method <- "est-mash"
head(est_mash_data)

est_mash_corrected_data <- read.csv("./10k_dataset/all_metrics_est_mash_corrected.txt", header = F, sep = "\t")
names(est_mash_corrected_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
est_mash_corrected_data$method <- "est-mash-corrected"
head(est_mash_corrected_data)

true_mash_data <- read.csv("./10k_dataset/all_metrics_true_mash.txt", header = F, sep = "\t")
names(true_mash_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
true_mash_data$method <- "true_mash"
head(true_mash_data)

true_data <- read.csv("./10k_dataset/all_metrics_true.txt", header = F, sep = "\t")
names(true_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
true_data$method <- "true"
head(true_data)

voltka_data <- read.csv("./10k_dataset/all_metrics_voltka.txt", header = F, sep = "\t")
names(voltka_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# voltka_data$jaccard <- 0
# voltka_data <- voltka_data[,c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")]
voltka_data$method <- "voltka"
head(voltka_data)

voltka_mash_data <- read.csv("./10k_dataset/all_metrics_voltka_mash.txt", header = F, sep = "\t")
names(voltka_mash_data) <- c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# voltka_data$jaccard <- 0
# voltka_data <- voltka_data[,c("tree", "k", "s", "dist", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")]
voltka_mash_data$method <- "voltka_mash"
head(voltka_mash_data)

# all_data <- rbind(est_data[, setdiff(names(est_data), c("jaccard"))], true_data[, setdiff(names(true_data), c("jaccard"))], voltka_data)
all_data <- rbind(est_data, true_data, voltka_data)
all_data_mash <- rbind(est_mash_data, true_mash_data, voltka_mash_data)

all <- rbind(all_data, all_data_mash)

head(all_data_mash)
all_data_mash$n = 0
all_data_mash[all_data_mash$tree == 7, ]$n = 50
all_data_mash[all_data_mash$tree == 4, ]$n = 89
all_data_mash[all_data_mash$tree == 15, ]$n = 111
all_data_mash[all_data_mash$tree == 17, ]$n = 34
all_data_mash[all_data_mash$tree == 1322, ]$n = 30

all_mash_metrics <- all_data_mash[all_data_mash$tree %in% c(15, 17, 1322) & all_data_mash$n / 3 > all_data_mash$k,] %>%
  mutate(`delta k`=est_k - k) %>%
  pivot_longer(cols = c("unifrac", "wunifrac", "delta k"), names_to = "metric")
all_mash_metrics$tree <- factor(all_mash_metrics$tree, levels = c("15", "17", "1322"), labels = c("subtree 1 (n=111)", "subtree 2 (n=34)", "subtree 3 (n=30)"))
all_mash_metrics$metric <- factor(all_mash_metrics$metric, levels = c("delta k", "unifrac", "wunifrac"), labels = c("delta k", "UniFrac", "weighted UniFrac"))
all_mash_metrics$method <- factor(all_mash_metrics$method, levels = c("true_mash", "est-mash", "voltka_mash"), labels = c(paste0("true distances + ", our_method), paste0("krepp + ", our_method), "krepp + Woltka"))

ggplot(data = all_mash_metrics, aes(x = as.factor(k), y = (value), 
                                    color = method, group = as.factor(method)))+
  stat_summary()+
  stat_summary(geom = "line")+
  # geom_boxplot(aes(group=interaction(k,method)))+
  theme_bw()+
  facet_grid(metric~tree, scales = "free")+
  scale_color_brewer(palette = "Dark2") +
  # scale_y_log10()+
  labs(x = "k", y= "", color = "method")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/10k_res.pdf"),width = 6,height = 5)


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
  scale_color_brewer(palette = "Dark2") +
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
  scale_color_brewer(palette = "Dark2") +
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


