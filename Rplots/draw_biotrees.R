library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

all_data <- c("bees-bosseret","mammals-song","birds-jarvis", "tilapia-Ciezarek","1kp", "pancrustacean-Bernot","beetles-Johnson", "fish-Troyer", "hemipteroid-johnson", "mammals-foley")

unifrac_data <- read.csv("./Demixing_Data/biotrees/all_unifrac_res.txt", header= F, sep = "\t")
names(unifrac_data) <- c("data", "k", "s", "noise", "seq_length", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
head(unifrac_data)

temp <- unifrac_data[unifrac_data$noise == 0, ]
temp$noise = 1
unifrac_data <- rbind(unifrac_data, temp)

temp$noise = 2
unifrac_data <- rbind(unifrac_data, temp)

unifrac_data <- unifrac_data[unifrac_data$noise != 0,]

# unifrac_data$error <- interaction(unifrac_data$noise, unifrac_data$seq_length, sep = "_")
# unifrac_data$error <- factor(
#   unifrac_data$error,
#   levels = c("0_0", "1_10000", "2_10000", "1_1000", "2_1000", "1_100", "2_100"),
#   labels = c("Noise=0", "Noise=1(l=10000)", "Noise=2(l=10000)", "Noise=1(l=1000)", "Noise=2(l=1000)", "Noise=1(l=100)", "Noise=2(l=100)")
# )

unifrac_data$seq_length <- factor(unifrac_data$seq_length, levels = c("0", "10000", "1000", "100"))

temp <- unifrac_data[unifrac_data$noise == 1, c("data", "k", "s", "seq_length", "jaccard")] %>% pivot_wider(names_from = seq_length, values_from = jaccard)
View(temp[temp$`0` < temp$`10000`,])

ggplot(data = unifrac_data, aes(x = as.factor(seq_length), y = (jaccard), color = as.factor(noise), group = noise == '1'))+
  stat_summary(fun = mean)+
  stat_summary(geom = "line")+
  theme_bw()+
  facet_wrap(.~k, nrow = 1)+
  scale_color_brewer(palette = "Dark2", labels = c("1" = "correlated", "2" = "independent")) +
  scale_x_discrete(labels = c("0" = 0, "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  # scale_y_log10()+
  labs(x = "noise (1/sequence length)", y = "Jaccard", color = "noise type")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/jaccard_bio.pdf"),width = 8,height = 2.8)

ggplot(data = unifrac_data, aes(x = as.factor(seq_length), y = (unifrac), color = as.factor(noise), group = noise == '1'))+
  stat_summary(fun = mean)+
  stat_summary(geom = "line")+
  theme_bw()+
  facet_wrap(.~k, nrow = 1)+
  scale_color_brewer(palette = "Dark2", labels = c("1" = "correlated", "2" = "independent")) +
  scale_x_discrete(labels = c("0" = 0, "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  # scale_y_log10()+
  labs(x = "noise (1/sequence length)", y = "UniFrac", color = "noise type")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/unifrac_bio.pdf"),width = 8,height = 2.8)

ggplot(data = unifrac_data, aes(x = as.factor(seq_length), y = (wunifrac), color = as.factor(noise), group = noise == '1'))+
  stat_summary(fun = mean)+
  stat_summary(geom = "line")+
  theme_bw()+
  facet_wrap(.~k, nrow = 1)+
  scale_color_brewer(palette = "Dark2", labels = c("1" = "correlated", "2" = "independent")) +
  scale_x_discrete(labels = c("0" = 0, "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  # scale_y_log10()+
  labs(x = "noise (1/sequence length)", y = "weighted UniFrac", color = "noise type")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/wunifrac_bio.pdf"),width = 8,height = 2.8)

ggplot(data = unifrac_data, aes(x = as.factor(seq_length), y = (est_k - k), color = as.factor(noise), group = noise == '1'))+
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed")+
  stat_summary(fun = mean)+
  stat_summary(geom = "line")+
  theme_bw()+
  facet_wrap(.~k, nrow = 1)+
  scale_color_brewer(palette = "Dark2", labels = c("1" = "correlated", "2" = "independent")) +
  scale_x_discrete(labels = c("0" = 0, "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  # scale_y_log10()+
  labs(x = "noise (1/sequence length)", y = "delta k", color = "noise type")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/est_k_bio.pdf"),width = 8,height = 2.8)

ggplot(data = unifrac_data, aes(x = as.factor(seq_length), y = (bc), color = as.factor(noise), group = noise == '1'))+
  stat_summary(fun = mean)+
  stat_summary(geom = "line")+
  theme_bw()+
  facet_wrap(.~k, nrow = 1)+
  scale_color_brewer(palette = "Dark2", labels = c("1" = "correlated", "2" = "independent")) +
  scale_x_discrete(labels = c("0" = 0, "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  # scale_y_log10()+
  labs(x = "noise (1/sequence length)", y = "Bray-Curtis", color = "noise type")+
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/bc_bio.pdf"),width = 8,height = 2.8)

##compare with other baseline methods
unifrac_data <- read.csv("./Demixing_Data/biotrees/all_unifrac_res.txt", header= F, sep = "\t")
names(unifrac_data) <- c("data", "k", "s", "noise", "seq_length", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
unifrac_data$method <- "standard"
head(unifrac_data)

unifrac_data_f1 <- read.csv("./Demixing_Data/biotrees/all_unifrac_res_f1.txt", header= F, sep = "\t")
names(unifrac_data_f1) <- c("data", "k", "s", "noise", "seq_length", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
unifrac_data_f1$method <- "f1"
head(unifrac_data_f1)

unifrac_data_f1_exh <- read.csv("./Demixing_Data/biotrees/all_unifrac_res_f1_exh.txt", header= F, sep = "\t")
names(unifrac_data_f1_exh) <- c("data", "k", "s", "noise", "seq_length", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
unifrac_data_f1_exh$method <- "f1_exh"
head(unifrac_data_f1_exh)

unifrac_data_f1_clos <- read.csv("./Demixing_Data/biotrees/all_unifrac_res_f1_clos.txt", header= F, sep = "\t")
names(unifrac_data_f1_clos) <- c("data", "k", "s", "noise", "seq_length", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
unifrac_data_f1_clos$method <- "f1_clos"
head(unifrac_data_f1_clos)

unifrac_data_f1_clos_iter <- read.csv("./Demixing_Data/biotrees/all_unifrac_res_f1_clos_iter.txt", header= F, sep = "\t")
names(unifrac_data_f1_clos_iter) <- c("data", "k", "s", "noise", "seq_length", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
unifrac_data_f1_clos_iter$method <- "f1_clos_iter"
head(unifrac_data_f1_clos_iter)

runtime_bio <- read.csv("./Demixing_Data/biotrees/all_runtime_res.txt", header = F, sep = "\t")
names(runtime_bio) <- c("data", "k", "s", "noise", "seq_length", "n", "est_k", "total", "one_round")
head(runtime_bio)

runtime_bio$total_n <- runtime_bio$k + runtime_bio$n
data_size <- runtime_bio %>% group_by(data) %>% summarise(n = mean(total_n))

#compare all
exh_data <- rbind(unifrac_data, unifrac_data_f1, unifrac_data_f1_exh, unifrac_data_f1_clos, unifrac_data_f1_clos_iter)
exh_data$seq_length <- factor(exh_data$seq_length, levels = c("0", "10000", "1000", "100"))
# exh_data$method <- factor(exh_data$method, levels = c("f1_exh", "f1", "standard"), labels = c("fixed_k exhaustive", "fixed_k", "standard"))
all_data <- unique(unifrac_data_f1_exh$data)
summary(exh_data$jaccard)

View(unifrac_data_f1_exh[unifrac_data_f1_exh$jaccard < 1 & unifrac_data_f1_exh$noise ==0 ,])
temp <- exh_data[, c("data", "noise", "seq_length", "k", "s", "jaccard", "method")] %>% pivot_wider(names_from = method, values_from = jaccard)
temp[temp$f1 < temp$f1_exh,]
head(temp)

combined <- exh_data %>% mutate(
                            delta_k=(est_k-k),
                            # wunifrac_comp = 1-wunifrac,
                            # unifrac_comp = 1-unifrac,
                            # bc_comp = 1-bc,
                            "j_comp"=1-jaccard
                            ) %>% pivot_longer(cols=c(j_comp, wunifrac, unifrac, bc, delta_k)) 
# %>% filter (!data %in% c("bees-bosseret","mammals-song","birds-jarvis") | k<6)
combined$name <- factor(combined$name, levels = c("j_comp", "bc",  "unifrac", "wunifrac", "delta_k"), labels=c("1 - Jaccard", "Bray-Curtis", "UniFrac", "wUniFrac",  "delta_k"))
combined$method <- factor(combined$method, levels = c("f1_exh", "f1", "standard", "f1_clos_iter", "f1_clos"), labels = c("exact-exact", "exact-heu", "heu-heu(default)", "k-nearest leaves (iterative)", "k-nearest leaves"))
combined <- merge(combined, data_size)

ggplot(data = combined[combined$n /3 > combined$k & combined$name != "delta_k",], aes(x = as.factor(seq_length), y = value, fill = as.factor(method), group = as.factor(method))) +
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6),color="black")+
  stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_grid(name~k,scale="free", labeller = labeller(k = function(x) paste0("k = ", x)))+
  # scale_fill_brewer(palette = "Set2") +
  scale_fill_manual(values = c("#b3de69", "#8dd3c7", "#80b1d3", "#fdb462", "#fb8072"))+
  scale_x_discrete(labels = c("0" = "0", "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  labs(x = "noise (1/sequence length)", y = "", fill = "method(k-search)") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/bio_comp_all.pdf"),width = 10,height = 7)


ggplot(data = combined[combined$n /3 > combined$k & combined$name == "delta_k" & combined$method == "heu-heu(default)",], aes(x = as.factor(seq_length), y = value, color = as.factor(method), group = as.factor(method))) +
  geom_hline(yintercept = 0, color = "grey50", linetype = "dashed")+
  stat_summary()+
  stat_summary(geom = "line")+
  # stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6),color="black")+
  # stat_summary(position = position_dodge(width = 0.6), geom="errorbar", width=0.4)+
  theme_bw() +
  facet_wrap(.~k, nrow = 1, labeller = labeller(k = function(x) paste0("k = ", x)))+
  # scale_fill_brewer(palette = "Set2") +
  scale_color_manual(values = c("#377eb8", "#fdb462", "#fb8072"))+
  scale_x_discrete(labels = c("0" = "0", "10000" = expression(10^-4), "1000" = expression(10^-3), "100" = expression(10^-2))) +
  labs(x = "noise (1/sequence length)", y = "delta k", fill = "method(k-search)") +
  theme(legend.position = "none")
ggsave(paste0("./Plots/bio_comp_delta_k_new.pdf"),width = 6,height = 2)




ggplot(data = exh_data[exh_data$data %in% all_data & exh_data$k %in% c(2,3) ,], aes(x = as.factor(seq_length), y = unifrac, fill = as.factor(method))) +
  # geom_col(position = "dodge") +
  # stat_summary(position = position_dodge(width = 0.2))+
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "unifrac", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/exh_unifrac_bio_comp.pdf"),width = 4.5,height = 3)

ggplot(data = exh_data[exh_data$data %in% all_data & exh_data$k %in% c(2,3) ,], aes(x = as.factor(seq_length), y = wunifrac, fill = as.factor(method))) +
  # geom_col(position = "dodge") +
  # stat_summary(position = position_dodge(width = 0.2))+
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "weighted unifrac", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/exh_wunifrac_bio_comp.pdf"),width = 4.5,height = 3)

ggplot(data = exh_data[exh_data$data %in% all_data & exh_data$k %in% c(2,3) ,], aes(x = as.factor(seq_length), y = wasser, fill = as.factor(method))) +
  # geom_col(position = "dodge") +
  # stat_summary(position = position_dodge(width = 0.2))+
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "wasserstein", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/exh_wasser_bio_comp.pdf"),width = 4.5,height = 3)

ggplot(data = exh_data[exh_data$data %in% all_data & exh_data$k %in% c(2,3) ,], aes(x = as.factor(seq_length), y = bc, fill = as.factor(method))) +
  # geom_col(position = "dodge") +
  # stat_summary(position = position_dodge(width = 0.2))+
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "bray-curtis", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/exh_bc_bio_comp.pdf"),width = 4.5,height = 3)

#compare the baseline methods
all_f1 <- rbind(unifrac_data_f1, unifrac_data_f1_clos, unifrac_data_f1_clos_iter)
all_f1$seq_length <- factor(all_f1$seq_length, levels = c("0", "10000", "1000", "100"))
all_f1$method <- factor(all_f1$method, levels = c("f1", "f1_clos", "f1_clos_iter"), labels = c("fixed_k", "fixed_k closest", "fixed k closest (iterative)"))

ggplot(data = all_f1, aes(x = as.factor(seq_length), y = jaccard, fill = as.factor(method))) +
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k, nrow = 1)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "jaccard", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/all_jaccard_bio_comp.pdf"),width = 7,height = 3)

ggplot(data = all_f1, aes(x = as.factor(seq_length), y = unifrac, fill = as.factor(method))) +
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k, nrow = 1)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "unifrac", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/all_unifrac_bio_comp.pdf"),width = 7,height = 3)

ggplot(data = all_f1, aes(x = as.factor(seq_length), y = wunifrac, fill = as.factor(method))) +
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k, nrow = 1)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "weighted unifrac", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/all_wunifrac_bio_comp.pdf"),width = 7,height = 3)

ggplot(data = all_f1, aes(x = as.factor(seq_length), y = wasser, fill = as.factor(method))) +
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k, nrow = 1)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "wasserstein", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/all_wasser_bio_comp.pdf"),width = 7,height = 3)

ggplot(data = all_f1, aes(x = as.factor(seq_length), y = bc, fill = as.factor(method))) +
  stat_summary(geom = "col", width = 0.5, position = position_dodge(width = 0.6))+
  theme_bw() +
  facet_wrap(~k, nrow = 1)+
  # facet_grid(noise~k, scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(labels = c("0" = expression(Inf), "10000" = expression(10^4), "1000" = expression(10^3), "100" = expression(10^2))) +
  labs(x = "sequence length", y = "bray-curtis", fill = "method") +
  theme(legend.position = "bottom")
ggsave(paste0("./Plots/all_bc_bio_comp.pdf"),width = 7,height = 3)
