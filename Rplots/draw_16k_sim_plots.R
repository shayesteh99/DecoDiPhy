library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

# all_metrics <- read.csv("./16k_sim/n100/all_metrics.txt", sep = "\t", header = F)
# names(all_metrics) <- c("sample", "k", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
# all_metrics$method <- "krepp"
# all_metrics$state <- "after"
# head(all_metrics)

all_metrics <- read.csv("./16k_sim/all_metrics_new.txt", sep = "\t", header = F)
names(all_metrics) <- c("model", "sample", "k", "n", "method", "state", "jaccard", "unifrac", "wunifrac", "est_k", "wasser", "bc")
head(all_metrics)

#all_metrics <- 
all_metrics <- all_metrics %>% 
  mutate(delta_k = (est_k - k)/k
           # log1p(abs(est_k - k))
         , jaccard_inv = 1-jaccard) %>% 
  select(model, sample, k, n, method, state, unifrac, wunifrac, jaccard, bc, delta_k) %>%
  pivot_longer(cols = c("unifrac", "wunifrac", "jaccard", "bc", "delta_k"))
    
all_metrics$full_name <- paste(all_metrics$method, all_metrics$state)
all_metrics$name <- factor(all_metrics$name, levels = c("delta_k", "unifrac", "wunifrac", "jaccard", "bc"), labels = c("delta k", "UniFrac", "weighted UniFrac", "Jaccard", "Bray-Curtis"))
all_metrics$full_name <- factor(all_metrics$full_name, levels = c("krepp_multiplacement after","krepp_multiplacement before","woltka before","woltka after","sylph before"),
                                labels = c("krepp+DecoDiPhy", "krepp", "Woltka", "Woltka+DecoDiPhy","sylph"))
# all_metrics$model <- factor(all_metrics$model, levels = c("simulated", "simulated_half_pruned"), labels = c("all novel", "partly novel"))

all_metrics %>% filter(name=="UniFrac", model=="simulated") %>% group_by(method,state) %>% summarise(mean(value))

f=c();for (m in levels(all_metrics$name)) {
a=ggplot(data =all_metrics[ !all_metrics$model%in%c("simulated_200k_novel","simulated_200k", "simulated_200k_1000r")&
                           all_metrics$state != "after(100)"&
                          ! all_metrics$full_name %in% c("Woltka+DecoDiPhy")
                          &all_metrics$name==m,
                         ],
       aes(x = interaction(n,model), y = value, 
                              color = full_name))+
  geom_hline(yintercept = 0, color = "grey50", linetype="dashed")+
  geom_boxplot(outlier.size = 0.5,outlier.alpha = 0.33)+
  stat_summary(position=position_dodge2(width = 0.75),size=.25)+
  #facet_wrap(~name, scales = "free") +
  scale_x_discrete(labels=c("Low\nHigh","Mid\nHigh","High\nHigh","Low\nLow","Mid\nLow","High\nLow"))+
  # scale_color_manual(values = c("#1f78b4", "#a6cee3", "#ff7f00", "#6a3d9a", "#fdbf6f", "#cab2d6"),name="")+
  scale_color_manual(values = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#fdbf6f", "#cab2d6"),name="")+
  theme_bw()+
  labs(x="", y=m)+
  theme(legend.position = "bottom",panel.grid.major.x   = element_blank(),
        axis.text.x = element_text())
f=c(f,a)
};
f[[1]] = f[[1]]+coord_cartesian(ylim=c(0,50));
f[[3]]=f[[3]]+coord_cartesian(ylim=c(0,1))
# f[[1]]=f[[1]]+labs(y="k relative error")

plot_grid(f[[4]],f[[2]],f[[3]],f[[1]],nrow=1,labels="auto")
plot_grid(f[[2]],f[[1]],nrow=1,labels="auto")
ggsave(paste0("./Plots/16k_simulated_no_bc.pdf"),width = 13.5,height = 4)
ggsave(paste0("./Plots/16k_simulated_unfiltered.pdf"),width = 7,height = 4)


ggplot(data =all_metrics[all_metrics$n == 100,], aes(x = cut(k, c(0, 20, 40, 60, 80)), y = value, 
                               color = full_name, group = full_name))+
  stat_summary()+ 
  stat_summary(geom = "line")+
  geom_hline(yintercept = 0, color = "grey50", linetype="dashed")+
  facet_wrap(~name, scales ="free_y", nrow=1)+
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#fdbf6f", "#ff7f00"), labels = c("krepp+DecoDiPhy", "krepp", "Woltka+DecoDiPhy", "Woltka"), name="")+
  scale_linetype_manual(values = c(1, 5), labels = c("with DecoDiPhy", "alone"), name = "")+
  theme_bw()+
  labs(x="", y="")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0("./Plots/16k_simulated.pdf"),width = 9,height = 2.8)

cut(all_metrics[all_metrics$n == 1000,]$k, 5)
ggplot(data =all_metrics[all_metrics$n == 1000,], aes(x = cut(k, c(0, 100, 150, 200, 250, 300, 400)), y = value, 
                                                     color = full_name, group = full_name))+
  stat_summary()+
  stat_summary(geom = "line")+
  geom_hline(yintercept = 0, color = "grey50", linetype="dashed")+
  facet_wrap(~name, scales ="free_y", nrow=1)+
  scale_color_manual(values = c("#a6cee3", "#1f78b4", "#fdbf6f", "#ff7f00"))+
  scale_linetype_manual(values = c(1, 5), labels = c("with DecoDiPhy", "alone"), name = "")+
  theme_bw()+
  labs(x="number of query taxa", y="")+
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
ggsave(paste0("./Plots/16k_simulated.pdf"),width = 9,height = 2.8)
