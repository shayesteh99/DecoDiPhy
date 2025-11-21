library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(pheatmap)
library(tidyverse)
require(plotly)
# library(webshot2)
library(reticulate)
library(png)
library(grid)
library(ggpattern)


#DA figure
da_data <- read.csv("./16k_data_placements/ancom_Diagnosis_minp_1000_exported/features.txt", header = F)
names(da_data) <- c("feature", "Control", "diff")
da_data$method <- "DecoDiPhy"

method="sourmash"
data <- read.csv(paste0("./16k_data_placements/ancom_Diagnosis_", method, "_exported/features.txt"), header = F)
names(data) <- c("feature", "Control", "diff")
data$method <- method
data <- data %>% filter(feature %in% unique(da_data$feature))
da_data <- rbind(da_data, data)

da_data$IBD <- da_data$Control + da_data$diff
head(da_data)

labels <- read.csv("./16k_data_placements/ancom_Diagnosis_minp_1000_exported/tax_ids.txt", header = T, sep = "\t")
names(labels) <- c("feature", "species", "shortname")
head(labels)

da_data <- merge(da_data, labels[,c("feature", "shortname")])

da_data$method <- factor(da_data$method, levels = c("DecoDiPhy", method))
da_data$shortname <- factor(da_data$shortname, levels = da_data %>%
                                                    filter(method == "DecoDiPhy") %>%
                                                    group_by(shortname) %>%
                                                    arrange(desc(diff)) %>%
                                                    pull(shortname))
da_data <- da_data %>%
  arrange(method == method)

ggplot(data =da_data, aes(x=shortname, y=diff, color = grepl("N",feature), fill = interaction(method,diff>0))) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf, fill = "#e0f3f8", linewidth = 0, alpha=0.1)+
  geom_rect(xmin=-Inf, xmax=Inf, ymin=0, ymax=-Inf, fill = "#FBF0EE", linewidth = 0, alpha=0.05)+
geom_bar(stat="identity", position="identity", width=0.85, linewidth = 0.8,
         alpha=0.6)+
  # geom_label(aes(label=feature))+
  #scale_color_manual(values = c("#b2182b", "#4575b4"), guide="none")+
  scale_color_manual(values = c("black", "black"), guide="none")+
  scale_fill_manual(values = c("#fb9a99","#b2182b", "#74add1", "#4575b4"), labels = c("", "", "DecoDiPhy", method), name = "")+
  labs(x="features", y="lfc diff", color ="group")+
  theme_bw()+
theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1,size=12))+
  annotate("text", x = 45, y = 3, label = "IBD", size = 6, fontface = "bold", color = "#4575b4") +
  annotate("text", x = 15, y = -3, label = "Control", size = 6, fontface = "bold", color = "#b2182b")+
  guides(fill = guide_legend(nrow = 2, override.aes = list(
    colour = c("black", "black", "black", "black"),
    linewidth = 0.8)))+
  theme(legend.position = c(0.057, 0.19))+
  labs(x="")
ggsave("./Plots/da_figure.pdf",width=14,height = 4.5)
ggsave("./Plots/da_figure.png",width=14,height = 4.5, dpi = 400)

ggsave(paste0("./Plots/da_figure_", method, ".pdf"),width=14,height = 4.5)



#PCOA Plot
decodiphy_pcoa <- read.csv("16k_data_placements/weighted_unifrac_krepp_pcoa_export/ordination_processed.txt", row.names = 1, header=F, sep="\t")
decodiphy_pcoa <- read.csv("16k_data_placements/unweighted_unifrac_minp_1000_pcoa_export/ordination_processed.txt", row.names = 1, header=F, sep="\t")
decodiphy_pcoa <- read.csv("16k_data_placements/weighted_unifrac_krepp_pcoa_top_export/ordination_processed.txt", row.names = 1, header=F, sep="\t")
decodiphy_pcoa <- read.csv("16k_data_placements/unweighted_unifrac_minp_1000_pcoa_top_export/ordination_processed.txt", row.names = 1, header=F, sep="\t")
decodiphy_pcoa <- read.csv("16k_data_placements/weighted_unifrac_minp_1000_top_pcoa_export/ordination_processed.txt", row.names = 1, header=F, sep="\t")
# decodiphy_pcoa <- read.csv("EMP/weighted_unifrac_minp_1000_top_pcoa_export/ordination_processed.txt", row.names = 1, header=F, sep="\t")
head(decodiphy_pcoa)

metadata <- read.table("16k_data_placements/metadata_diagnosis.tsv", header = TRUE, row.names = 1, sep="\t", check.names = FALSE)
# metadata <- read.table("EMP/metadata.tsv", header = TRUE, row.names = 1, sep="\t", check.names = FALSE)
head(metadata)

read.csv("16k_data_placements/test_samples.txt") -> test

nrow(decodiphy_pcoa[rownames(metadata), ])

merged_df <- merge(
  metadata %>% select(Diagnosis_combined, Diagnosis, Cohort),
  # metadata %>% select(EMPO2),
  decodiphy_pcoa %>% select(V2, V3, V4),
  by = "row.names"
)

merged_df$Diagnosis <- factor(merged_df$Diagnosis, levels = c("Control", "UC", "CD"))
merged_df <- merged_df %>%
  mutate(symbol_group = interaction(Diagnosis, Cohort, sep = "_"))

merged_df$Diagnosis <- factor(merged_df$Diagnosis, levels = c("Control", "UC", "CD"))
merged_df$data <- "IBD"

ggplot(merged_df[,], 
       aes(x=V2, y=V3,  color=Diagnosis,shape=Row.names %in% test$SampleID))+
  geom_point(size=1.5,alpha=0.8)+
  scale_color_manual(values = c("#b2182b", "#4575b4","#74add1"), 
                     # labels = c("Non Saline", "Saline"), 
                     name = "")+
  scale_shape_manual(values = c(1,16,18,20), name="",labels=c("train","test"))+
  theme_bw()+
  facet_wrap(~data)+
  labs(y="PC2", x="PC1")
ggsave("./Plots/IBD_pcoa.pdf",width=3.9,height = 2.7)



p <- plot_ly(
  data = merged_df,
  x = ~V3,
  y = ~V2,
  z = ~V4,
  color = ~Diagnosis,
  symbol = ~symbol_group,
  symbols = c("circle", "circle", "circle", "diamond", "diamond", "diamond"),
  colors = c("#b2182b", "#74add1", "#4575b4"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 8, opacity = 0.8)
) %>%
  plotly::layout(
    scene = list(camera = list(eye = list(x = 1.8, y = 0.6, z = 0.4)),
     xaxis = list(title = list(text = "PC1", font = list(size = 20)),
                  tickfont = list(size = 16)),
     yaxis = list(title = list(text = "PC2", font = list(size = 20)),
                  tickfont = list(size = 16)),
     zaxis = list(title = list(text = "PC3", font = list(size = 20)),
                  tickfont = list(size = 16))
      ),
    legend = list(font = list(size = 16)))
p


decodiphy_pcoa <- read.csv("EMP/weighted_unifrac_minp_1000_top_pcoa_export/ordination_processed.txt", row.names = 1, header=F, sep="\t")
head(decodiphy_pcoa)
metadata <- read.table("EMP/metadata.tsv", header = TRUE, row.names = 1, sep="\t", check.names = FALSE)
head(metadata)

read.csv("./EMP/test_samples_EMPO2.txt") -> test_EMP

merged_df <- merge(
  # metadata %>% select(Diagnosis_combined, Diagnosis, Cohort),
  metadata %>% select(EMPO2, EMPO3),
  decodiphy_pcoa %>% select(V2, V3, V4),
  by = "row.names"
)

unique(merged_df$EMPO3)

merged_df <- merged_df %>%
  mutate(EMPO3 = recode(EMPO3,
            "Solid (n-s)" = "Solid",
            "Solid (s)" = "Solid",
            "Aqueous (n-s)" = "Aqueous",
            "Aqueous (s)" = "Aqueous"
))
merged_df$data <- "EMP"

ggplot(merged_df, aes(x=V2, y=V3, color=EMPO2, fill=interaction(EMPO2, EMPO3,Row.names %in% test_EMP$SampleID), shape=interaction(EMPO3,Row.names %in% test_EMP$SampleID)))+
  geom_point(alpha=0.8, size=1)+
  scale_color_manual(values = c("#b35806", "#542788"), labels = c("Non Saline", "Saline"), name = "")+
  scale_shape_manual(values = c(2,6,24,25), name="", labels=c("Aqueous(train)", "Solid(train)", "Aqueous(test)", "Solid(test)"))+
  # scale_shape_manual(values = c(23,21,18,16), name="", labels=c("Aqueous(train)", "Solid(train)", "Aqueous(test)", "Solid(test)"))+
  # scale_size_manual(values = c(2,2,2,2), name="", labels=c("Aqueous(train)", "Solid(train)", "Aqueous(test)", "Solid(test)"))+
  scale_fill_manual(values = c("white","white","white","white","#b35806","#542788","#b35806", "#542788"), name="", guide="none",
                    # labels=c("Aqueous(train)", "Solid(train)", "Aqueous(test)", "Solid(test)")
                    )+
  theme_bw()+
  facet_wrap(~data)+
  labs(y="PC2", x="PC1")+
  guides(
    shape = guide_legend(
      override.aes = list(
        fill = c("white", "white", "grey10", "grey10")
      )
    ))
ggsave("./Plots/EMP_pcoa.pdf",width=4.3,height = 2.7)


p <- plot_ly(
  data = merged_df,
  x = ~V3,
  y = ~V2,
  z = ~V4,
  color = ~EMPO2,
  symbol = ~EMPO3,
  symbols = c("diamond", "diamond", "circle", "circle", "diamond"),
  colors = c("#b2182b", "#74add1", "#4575b4"),
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 8, opacity = 0.8)
) %>%
  plotly::layout(
    scene = list(camera = list(eye = list(x = 1.8, y = 0.6, z = 0.4)),
                 xaxis = list(title = list(text = "PC1", font = list(size = 20)),
                              tickfont = list(size = 16)),
                 yaxis = list(title = list(text = "PC2", font = list(size = 20)),
                              tickfont = list(size = 16)),
                 zaxis = list(title = list(text = "PC3", font = list(size = 20)),
                              tickfont = list(size = 16))
    ),
    legend = list(font = list(size = 16)))
p


plotly::export(p, file = "./Plots/decodiphy_IBD_pcoa.pdf")
plotly::save_image(p, "./Plots/decodiphy_IBD_pcoa.png", width = 2000, height = 2000)

html_file <- "./Plots/decodiphy_IBD_pcoa.html"
htmlwidgets::saveWidget(p, file = html_file)

webshot2::webshot(
  url = html_file,
  file = "./Plots/decodiphy_IBD_pcoa.png",
  vwidth = 1600,
  vheight = 1600,
  cliprect = c(-300,-100,1500,1500)
)

#get percentage
feat_table <- read.csv("./16k_data_placements/top_feature-table_minp_1000/feature-table.txt", sep = "\t")
head(feat_table)

feat_table <- feat_table %>%
  pivot_longer(
    cols = -X.OTU.ID,
    names_to = "Sample ID",
    values_to = "value"
  )

metadata <- read.table("16k_data_placements/metadata_diagnosis.tsv", header = TRUE, sep="\t", check.names = FALSE)

plt_data <- merge(feat_table, metadata %>% select("Sample ID", "Diagnosis_combined")) %>%
  group_by(X.OTU.ID, Diagnosis_combined) %>%
  summarise(n = mean(value > 0) * 100) %>%
  mutate(fill_color = case_when(
    Diagnosis_combined == "IBD" ~ scales::col_numeric(
      palette = c("white", "#2166ac"), domain = c(0, 100)
    )(n),
    Diagnosis_combined == "Control" ~ scales::col_numeric(
      palette = c("white", "#b2182b"), domain = c(0, 100)
    )(n)
  ))

plt_data <- merge(plt_data, labels[,c("feature", "shortname")], by.x ="X.OTU.ID", by.y = "feature")

# da_data$method <- factor(da_data$method, levels = c("DecoDiPhy", "krepp"))
plt_data$shortname <- factor(plt_data$shortname, levels = da_data %>%
                              filter(method == "DecoDiPhy") %>%
                              group_by(shortname) %>%
                              arrange(desc(diff)) %>%
                              pull(shortname))

ggplot(plt_data, aes(x = shortname, y = Diagnosis_combined, fill = fill_color)) +
  geom_tile(color = "black", size=0.3) +
  geom_text(aes(label = sprintf("%.0f%%", n)), size = 3) +
  scale_fill_identity() +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(fill = "% nonzero")
ggsave("./Plots/IBD_percentage.pdf",width=17,height = 2.7)


#height distributions

IBD_h <- read.csv("16k_data_placements/IBD_top_features_height.txt", header=F, sep ="\t")
names(IBD_h) <- c("feature", "height")
IBD_h$data <- "IBD(n=57)"
head(IBD_h)

EMP_h <- read.csv("EMP/EMP_top_features_height.txt", header=F, sep ="\t")
names(EMP_h) <- c("feature", "height")
EMP_h$data <- "EMP(n=157)"
head(EMP_h)

h_data <- rbind(IBD_h, EMP_h)
h_data$data <- factor(h_data$data, levels = c("IBD(n=57)","EMP(n=157)"))
ggplot(h_data, aes(x = data, y=height, color=data, fill=data, group = data))+
  # geom_violin(draw_quantiles = c(0.25,0.5,0.75), linewidth = 0.5)+
  # geom_point(aes(shape = data),
  #            position = position_jitter(width = 0),
  #            size = 2, alpha = 0.8) +
  geom_boxplot(alpha=0.7)+
  scale_color_manual(values =c("#8c510a","#01665e"), guide="none")+
  scale_fill_manual(values =c("#f6e8c3","#c7eae5"), guide="none")+
  theme_bw()+
  labs(x="",y="depth")
ggsave("./Plots/hight_dist.pdf",width=2.3,height = 3)

ggplot(h_data, aes(x=cut(height,c(-1,0,0.05,0.1,0.5,0.8)), color=data, fill=data, group = data))+
  # geom_violin(draw_quantiles = c(0.25,0.5,0.75), linewidth = 0.5)+
  # geom_point(aes(shape = data),
  #            position = position_jitter(width = 0),
  #            size = 2, alpha = 0.8) +
  # stat_ecdf()+
  stat_count(position = position_dodge())+
  scale_color_manual(values =c("#8c510a","#01665e"))+
  scale_fill_manual(values =c("#f6e8c3","#c7eae5"), guide="none")+
  theme_bw()+#scale_y_continuous(labels = percent)+
  theme(axis.title.y = element_blank(),legend.position = c(0.15,0.75),
        legend.title = element_blank())+
  labs(y="",x="depth")
ggsave("./Plots/hight_dist_hist.pdf",width=2.3,height = 3)

  

# bar plots for separation

#IBD
IBD_data <- read.csv("./16k_data_placements/f_results.txt", header = F, sep = "\t")
names(IBD_data) <- c("method", "data", "metric", "features", "f", "p_val")
head(IBD_data)
unique(IBD_data$method)


rf_data <- data.frame(method = c("minp_1000", "krepp", "woltka", "sylph", "sourmash", "kraken", "braken", "kraken_g", "braken_g"),
                   features=c(0.8,0.8,0.83,0.85,0.78,0.8,0.8,0.82,0.85),
                   features_top=c(0.83,0.86,0.71,0.89,0.75,0.82,0.8,0.82,0.8))
rf_data <- rf_data%>% pivot_longer(cols = c(features,features_top), names_to = "features", values_to = "rf")
IBD_data <- merge(IBD_data, rf_data)

num_data <- data.frame(method = c("minp_1000", "krepp", "woltka", "sylph", "sourmash", "kraken", "braken", "kraken_g", "braken_g"),
                       features=c(4715,30940,1533,973,1059,9393,12293,2948,2948),
                       features_top=c(48,1070,22,76,47,3711,6213,2892,2890))

num_data <- num_data%>% pivot_longer(cols = c(features,features_top), names_to = "features", values_to = "num")
IBD_data <- merge(IBD_data, num_data)
IBD_data$dataset <- "IBD"
head(IBD_data)

#EMP
EMP_data <- read.csv("./EMP/f_results_EMPO2.txt", header = F, sep = "\t")
names(EMP_data) <- c("method", "data", "metric", "features", "f", "p_val")
head(EMP_data)
unique(EMP_data$method)

rf_data <- data.frame(method = c("minp_1000", "krepp", "woltka", "sylph", "sourmash", "kraken", "braken"),
                      features=c(0.92,0.97,0.98,0.84,0.98,0.95,0.94),
                      features_top=c(0.93,0.95,0.97,0,1,0.95,0.95))
rf_data <- rf_data%>% pivot_longer(cols = c(features,features_top), names_to = "features", values_to = "rf")
EMP_data <- merge(EMP_data, rf_data)
EMP_data$dataset <- "EMP"
head(EMP_data)

num_data <- data.frame(method = c("minp_1000", "krepp", "woltka", "sylph", "sourmash", "kraken", "braken", "kraken_g", "braken_g", "minp_1000_woltka"),
                       features=c(6269,31120,14286,565,1855,12467,12467,2948,2948,2304),
                       features_top=c(67,12124,3151,0,9,11383,11789,2892,2890,0))

num_data <- num_data%>% pivot_longer(cols = c(features,features_top), names_to = "features", values_to = "num")
EMP_data <- merge(EMP_data, num_data)

data <- rbind(IBD_data, EMP_data)

data$method <- factor(data$method, levels = c("minp_1000", "krepp", "sourmash", "sylph", "woltka", "kraken", "braken"),
                      labels = c("krepp+DecoDiPhy", "krepp", "sourmash", "sylph", "Woltka", "Kraken", "Bracken"))

unique(data$metric)
data$metric <- factor(data$metric, levels=c("unweighted_unifrac", "weighted_unifrac", "BC"), 
                      labels = c("UniFrac", "wUniFrac", "Bray-Curtis"))

data%>%
  #filter(data == "data_test",features=="features_top",method!="sourmash") %>% 
  filter(method %in% c("sourmash", "sylph", "krepp+DecoDiPhy", "krepp", "Woltka", "Kraken", "Bracken")) %>%
  filter(!grepl("_g",method))%>%
  filter(data == "data_test",features=="features_top",method!="sylph" | dataset=="IBD") %>% 
  #select(method, features,metric,p_val,dataset) %>% 
  #mutate(p_val=log10(p_val)) %>%
  #pivot_wider(names_from = dataset, values_from = (p_val)) %>%
  ggplot(aes(x=num, y=p_val, color = method,shape=dataset))+
  geom_rect(xmin=-Inf, xmax=log10(100), ymin=-Inf, ymax=log10(0.01), fill = "#E6E6CE",
            linewidth = .2, alpha=0.05,color="grey40",linetype=2)+
  # scale_color_brewer(palette = "Dark2", name = "")+
  scale_color_manual(values=c("#1b9e77", "#d95f02", "#66a61e", "#e7298a", "#7570b3", "#e6ab02", "#a6761d"),name="")+
  facet_grid(.~metric)+
  scale_x_log10()+
  geom_point(size=3,alpha=0.99)+
  geom_line(aes(group=method),linetype=1, linewidth = 0.8,alpha=0.75)+
  scale_y_log10()+
  theme_bw()+
  scale_shape_manual(values = c(17,16), name="")+
  labs(y="p-value", x= "number of features")+
  theme(
    legend.spacing.y = unit(-10, "pt"),
    legend.spacing.x = unit(2, "pt")
  )
ggsave("./Plots/pval_scatter.pdf",width=9,height = 2.7)


data$method <- factor(data$method, levels =  c("Bracken","Kraken","Woltka","sylph","sourmash","krepp","krepp+DecoDiPhy"))
 data%>%
  filter(data == "data_test",features=="features_top",method!="sylph" | dataset=="IBD") %>% 
  filter(!grepl("_g",method))%>%
  filter(method %in% c("sourmash", "sylph", "krepp+DecoDiPhy", "krepp", "Woltka", "Kraken", "Bracken")) %>%
  # mutate(method = fct_relevel(method,"sourmash", "sylph","krepp+DecoDiPhy","krepp","Woltka","Kraken","Bracken")) %>%
  #select(method, features,metric,p_val,dataset) %>% 
  mutate(p_val=-log10(p_val)) %>%
  #pivot_wider(names_from = dataset, values_from = (p_val)) %>%
  ggplot(aes(y=method, x=dataset, label=num,fill = p_val,color=p_val< -log10(0.05)))+
  facet_wrap(.~metric,nrow=1,scales="free_x")+
  geom_tile(linewidth = 1)+
  theme_classic()+
  ylab("")+
  #scale_fill_viridis_b( breaks=c(0,1,2,2.5,3))+
  scale_fill_viridis_b( breaks=c(0,1,-log10(0.05),2,-log10(0.005),3),labels=function(x) 10^(-x),
                        name="p-val")+
  geom_text()+
  scale_color_manual(name="p-val",values=c("#1b7837","#f46d43"),
                     labels=c("<0.05",">0.05")
                     )+
  labs(x="")
ggsave("./Plots/real_tile_all.pdf",width=6,height = 4.5)
ggsave("./Plots/real_tile_top.pdf",width=6,height = 4.5)
 
data %>% filter(data=="data", features=="features_top", method=="krepp+DecoDiPhy")


ggplot(data %>% filter(data=="data_test"), aes(x = num, y = f, color = method)) +
  geom_point(size = 2, show.legend = T )+
  geom_line(aes(group = interaction(method,metric)))+
  # scale_color_manual(values = c("#238b45", "#cc4c02", "#6a51a3"), name = "") +
  facet_wrap(~metric)+
  scale_x_log10()
ggsave("./Plots/IBD_fplot_arrows.pdf",width=8.5,height = 2.5)


#pie chart
df <- data.frame(
  category = c("Known species association","Known genus association", "Gut (no prior association)", "Internal", "No evidence in human gut"),
  value = c(12,11,9,13,12)
)

df$category <- factor(df$category, levels = c("Known species association","Known genus association", "Gut (no prior association)", "Internal", "No evidence in human gut"))


ggplot(df, aes(x = "", y = value, fill = category, pattern = category)) +
  geom_col_pattern(
    width = 1,
    color = "black",
    pattern_fill = "black",
    pattern_colour = "black",
    pattern_angle = 45,
    pattern_density = 0.1,
    pattern_key_scale_factor = 0.5,
    pattern_spacing = 0.02,
    pattern_size = 0.1,
    pattern_linetype = 1,
    pattern_alpha = 0.5
  ) +
  # geom_text(
  #   aes(y = label_pos,  label = pct_label),
  #   color = "black",
  #   size = 5
  # ) +
  coord_polar(theta = "y", start = 1.57) +
  scale_fill_manual(values = c("#080808",  "#565656", "#BBBBBB", "white", "#FFFFFF"), name = "") +
  scale_pattern_manual(values = c("none", "none", "none", "stripe", "none"), name = "") +
  theme_void() +
  theme(legend.position = "right")
ggsave("./Plots/IBD_pie.pdf",width=5,height = 3)

#average y
data <- read.csv("./16k_data_placements/all_y.txt", sep="\t", header =F)
names(data) <- c("sample", "tree", "y")
data <- data %>% mutate(Sample.ID = paste0("SRR",sample))
head(data)

meta <- read.table("16k_data_placements/metadata_diagnosis.tsv", header=TRUE, sep="\t", row.names=NULL)
head(meta)

data <- merge(data, meta%>%select(Sample.ID,Diagnosis_combined))

ggplot(data, aes(x=as.factor(tree), y=y, color= Diagnosis_combined, group=interaction(tree,Diagnosis_combined)))+
  geom_boxplot()
ggsave("./Plots/y_plot_ibd.pdf", width = 5, height = 3)


#rf
data <- data.frame(method = c("krepp+DecoDiPhy", "krepp", "Woltka", "sylph", "sourmash"), 
                      num=c(4715,30940,1533,973,1059),
                      rf=c(0.8,0.8,0.83,0.85,0.78))
data$features <- "all features"

d <- data.frame(method = c("krepp+DecoDiPhy", "krepp", "Woltka", "sylph", "sourmash"), 
                   num=c(48,1070,22,76,47),
                    rf=c(0.83,0.86,0.71,0.89,0.75))
d$features <- "DA features"

data <- rbind(data, d)
# data <- data %>% pivot_longer(cols = c("wunifrac", "unifrac", "bc"), names_to = "metric", values_to = "value")
# head(data)

data <- data %>% mutate(full_method = paste0(method, "(n=", num, ")"))

# data$method <- factor(data$method, levels = c("krepp+DecoDiPhy", "krepp", "Woltka", "sylph", "sourmash"))
# data$metric <- factor(data$metric, levels = c("unifrac", "wunifrac", "bc"), labels=c("UniFrac", "weighted UniFrac", "Bray Curtis"))

data$full_method <- factor(data$full_method, levels = data %>% filter(metric=="UniFrac")
                                                          %>%  arrange(desc(full_method))
                                                          %>% pull(full_method))

arrows_df <- data %>%
  pivot_wider(
    id_cols = c(method, metric),
    names_from = features,
    values_from = c(num, value),
    names_sep = "."
  )

ggplot(data %>% filter(metric!="Bray Curtis"), aes(x = num, y = rf, color = method)) +
  geom_point(data = data %>% filter(features=="all features",metric!="Bray Curtis") ,size = 2, show.legend = FALSE ) +
  geom_segment(
    data = arrows_df%>% filter(metric!="Bray Curtis"),
    aes(
      x = `num.all features`, y = `value.all features`,
      xend = `num.DA features`, yend = `value.DA features`,
      color = method,
      linetype = metric
    ),
    linewidth = 1,
    show.legend = TRUE
  ) +
  geom_curve(data = arrows_df%>% filter(metric!="Bray Curtis"),
    aes(
      x = `num.all features`, y = `value.all features`,
      xend = `num.DA features` *0.98, yend = `value.DA features`*1.002,
      color = method
    ),
    curvature = 0,
    arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
    linewidth = 0,
    linetype = "solid"
  )+
  # geom_segment(
  #   data = arrows_df,
  #   aes(
  #     x = `num.all features`, y = `value.all features`,
  #     xend = `num.DA features` *0.98, yend = `value.DA features`*1.002,
  #     color = method
  #   ),
  #   arrow = arrow(length = unit(0.3, "cm"), type = "closed", ends = "last"),
  #   lineend = "round",
  #   linewidth = 0.0
  # ) +
  scale_color_manual(values = c("#238b45", "#cc4c02", "#6a51a3"), name = "") +
  scale_linetype_manual(values = c(1,2,3), name = "") +
  scale_x_log10() +
  theme_minimal() +
  theme(legend.position = "right")+
  labs(x = "number of features", y = "pseudo-F")
  
ggsave("./Plots/IBD_fplot_arrows.pdf",width=6.5,height = 3.5)


ggplot(data %>% filter(features=="all features"), aes(x = as.factor(metric), y=value, fill=full_method))+
  stat_summary(geom = "col",
               width = 0.5,
               position = position_dodge(width = 0.6),
               color = "black",
               linewidth = 0.6)+
  scale_fill_manual(values = c("#238b45", "#cc4c02", "#6a51a3"), name="")+
  facet_wrap(~features, ncol=2)+
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing.y = unit(0.8, "lines"),
    axis.text.x = element_text(angle=15, hjust=0.5, vjust=1,size = 12, color = "grey20"),
    axis.text.y = element_text(size = 11, color = "grey30"),
    strip.text = element_text(size = 13, color = "grey20"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    axis.title.y = element_text(size = 13, margin = margin(r = 10))
  )+
  labs(x="", y="pseudo-F")
ggsave("./Plots/IBD_fplot.pdf",width=6.5,height = 3.5)

#computing correlation between samples

method="minp_1000"
method="krepp"
method="sylph"
method="woltka"
metric="weighted_unifrac"
all_feat <- read.csv(paste0("./16k_data_placements/", metric,"_", method,"/distance-matrix.tsv"), row.names=NULL, sep="\t")
head(all_feat)

all_feat <- all_feat %>% pivot_longer(cols = -X, names_to = "id", values_to = "dist") %>% mutate(pair=paste(X,id)) %>% select(pair, dist)

da_feat <- read.csv(paste0("./16k_data_placements/", metric,"_", method,"_top/distance-matrix.tsv"), row.names=NULL, sep="\t")
head(da_feat)
da_feat <- da_feat %>% pivot_longer(cols = -X, names_to = "id", values_to = "dist") %>% mutate(pair=paste(X,id)) %>% select(pair, dist)
merged <- merge(all_feat, da_feat, by="pair")

cor(merged$dist.x, merged$dist.y, method = "pearson")

#k distributions
p_info <- read.csv("./16k_data_placements/info/p_info.txt", sep = "\t", header = F)
names(p_info) <- c("sample", "tree", "count", "p")
head(p_info)

k_results <- read.csv("./16k_data_placements/info/k_results_diff_0.1+minp_5000.txt", sep = "\t", header = F)
names(k_results) <- c("sample", "tree", "k")
k_results$method <- "diff_0.1+minp_5000"
head(k_results)

for (m in c("minp_1000", "diff_0.4", "diff_0.1+minp", "diff_0.05+minp_1000", "diff_0.01+minp_1000")) {
  k_data <- read.csv(paste0("./16k_data_placements/info/k_results_", m, ".txt"), sep = "\t", header = F)
  names(k_data) <- c("sample", "tree", "k")
  k_data$method <- m
  head(k_data)
  k_results <- rbind(k_results, k_data)
}

meta <- read.table("16k_data_placements/metadata_diagnosis.tsv", header=TRUE, sep="\t", row.names=NULL)
names(meta)[1] <- "sample"
meta$sample <- sub("SRR", "", meta$sample)
head(meta)

merge(meta, k_results, by="sample") %>% grou

merge(meta, k_results, by="sample") %>% group_by(Diagnosis_combined, method) %>% summarise(n())

p_info %>% mutate(est_k = as.integer((count*p)/1000))

merge(merge(meta, k_results, by="sample"), p_info)%>% group_by(method, Diagnosis_combined) %>% summarise(n=sum(k)/n())
merge(k_results, p_info)%>% group_by(method) %>% summarise(n=sum(k))

ggplot(data = merge(k_results, p_info), aes(x = cut(p, 10), y = k, group = cut(p, 10)))+
  geom_boxplot()+
  facet_wrap(~method)

#number of reads distribution
read_counts <- read.csv("./16k_data_placements/read_count_data.txt", sep = "\t", header = F)
names(read_counts) <- c("sample", "tree", "count")
head(read_counts)

table(cut(read_counts$count, c(0, 500, 1000, 5000, 10000, 50000, 1e5, 1e6, 1e7, 1e10)))

krepp_ancom <- read.csv("./16k_data_placements/ancom_Diagnosis_krepp_exported/features.txt", header = F)
names(krepp_ancom) <- c("feature", "int", "diff")
krepp_ancom$control <- exp(krepp_ancom$int)
krepp_ancom$IBD <- exp(krepp_ancom$int + krepp_ancom$diff)
head(krepp_ancom)

loss_ancom <- read.csv("./16k_data_placements/ancom_Diagnosis_loss_exported/features.txt", header = F)
names(loss_ancom) <- c("feature", "int", "diff")
loss_ancom$control <- exp(loss_ancom$int)
loss_ancom$IBD <- exp(loss_ancom$int + loss_ancom$diff)
head(loss_ancom)

values <- c(krepp_ancom$control, krepp_ancom$IBD)
group <- rep(c("Control", "IBD"), each = nrow(krepp_ancom))

#R2
group_means <- tapply(values, group, mean)
grand_mean <- mean(values)

SS_between <- sum(table(group) * (group_means - grand_mean)^2)
SS_total   <- sum((values - grand_mean)^2)

R2 <- SS_between / SS_total
R2

#F_statistics
n_groups <- length(unique(group))
n_total <- length(values)
MS_between <- SS_between / (n_groups - 1)
MS_within  <- (SS_total - SS_between) / (n_total - n_groups)
F_statistics <- MS_between / MS_within
F_statistics


#heat map
dist_mat <- read.table(
  "16k_data_placements/weighted_unifrac_diff_0.35/distance-matrix.tsv",
  # "16k_data_placements/weighted_unifrac_krepp/distance-matrix.tsv",
  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
dist_mat <- as.matrix(dist_mat)

meta <- read.table("16k_data/metadata_diagnosis.tsv", header=TRUE, sep="\t", row.names=1)
head(meta)

ordered_samples <- rownames(meta)[order(meta$Diagnosis_combined)]
dist_mat <- dist_mat[ordered_samples, ordered_samples]
meta <- meta[ordered_samples, , drop=FALSE]

pheatmap(
  dist_mat,
  cluster_rows = T,
  cluster_cols = T,
  annotation_col = meta["Diagnosis_combined"],
  annotation_row = meta["Diagnosis_combined"],
  show_rownames = FALSE,
  show_colnames = FALSE
)
pheatmap(
  dist_mat,
  clustering_method = "none",
  cluster_rows ="euclidean",
  cluster_cols = "euclidean",
  annotation_col = meta,
  annotation_row = meta,
  show_rownames = FALSE,
  show_colnames = FALSE
)

