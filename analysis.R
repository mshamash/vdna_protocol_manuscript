library(tidyverse)
library(phyloseq)
library(vegan)
library(flextable)
library(cowplot)


### DATA CLEANING
phage_cov <- read_tsv("phage-coverage.tsv")
phage_cov_filt <- phage_cov %>% filter(meandepth >= 1)
phage_cov_filt.breadth <- phage_cov_filt %>% mutate(breadth = covbases/endpos)
phage_cov_filt.breadth <- phage_cov_filt.breadth %>% filter(breadth >= 0.75)

phage_cov.summary <- phage_cov_filt.breadth %>%
  group_by(sample, contig) %>%
  summarize(phage_reads = numreads, phage_cov = meandepth, phage_len = endpos)

phage_cov.totals <- phage_cov.summary %>% 
  group_by(sample) %>% 
  summarize(sample_total_reads_mapped = sum(phage_reads), sample_total_cov = sum(phage_cov))

phage_cov.joined <- inner_join(phage_cov.summary, phage_cov.totals, by = c("sample" = "sample"))

phage_cov.joined$relabund_cov <- (phage_cov.joined$phage_cov / phage_cov.joined$sample_total_cov)

phage_cov.joined.filt <- phage_cov.joined %>%
  ungroup()

phage_metadata <- read_csv("vlp-metadata.csv")
phage_cov.joined.filt <- left_join(phage_cov.joined.filt, phage_metadata, by = c("sample" = "SampleName"))


### PHYLOSEQ OBJECT
votu.table <- as.data.frame(pivot_wider(phage_cov.joined.filt[, c("sample", "contig", "phage_cov")], id_cols = "contig", names_from = "sample", values_from = "phage_cov"))
rownames(votu.table) <- votu.table$contig
votu.table <- votu.table[, -1]
votu.table <- floor(votu.table)
votu.table[is.na(votu.table)] <- 0

OTU <- otu_table(votu.table, taxa_are_rows = T, )

metadata <- read.csv("vlp-metadata.csv")
rownames(metadata) <- metadata$SampleName
metadata$Environment <- factor(metadata$Environment)
metadata$ExtractionMethod <- factor(metadata$ExtractionMethod)
metadata$SampleID <- factor(metadata$SampleID)
METADATA <- sample_data(metadata)

ps <- phyloseq(OTU, METADATA)


### BETA-DIVERSITY ANALYSIS
ord.nmds.bray <- ordinate(ps, method="NMDS", distance="bray")
nmds.bray.plot <- plot_ordination(ps, ord.nmds.bray, color="Environment", shape = "ExtractionMethod") +
  geom_point(size = 5) +
  scale_color_manual(values = c("#CC79A7", "#F0E442", "#0072B2")) +
  theme_linedraw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

nmds.bray.plot
ord.nmds.bray$stress

sampledf <- data.frame(sample_data(ps))
ps_BC <- vegdist(t(otu_table(ps)), method = "bray")

adonis2(ps_BC ~ ExtractionMethod + Environment + SampleID, sampledf)


### TABLE WITH PERMANOVA RESULTS
permanova.table <- data.frame(
  `Effect` = c("ExtractionMethod", "Environment", "Sample"),
  `R2` = c("0.00185", "0.36245", "0.61344"),
  `P.value` = c("0.426", "0.001", "0.001"))

permanova.flextable <- flextable(permanova.table) %>% 
  add_header_row(colwidths = c(3),
                 values = c("Bray-Curtis Distance ~ ExtractionMethod + Environment + Sample")) %>% 
  theme_vanilla() %>% 
  add_footer_lines("PERMANOVA (adonis2) with 999 permutations") %>% 
  autofit(add_w = 0.8) %>% 
  bold(i = 3, j = 3) %>% 
  bold(i = 2, j = 3) %>% 
  flextable::compose(i = 2, j = 2, part = "header", value = as_paragraph("R", as_sup("2"))) %>% 
  flextable::compose(i = 2, j = 3, part = "header", value = as_paragraph("P-value"))

permanova.flextable

permanova.flextable.grob <- gen_grob(permanova.flextable, fit = "fixed", just = "center")


### ALPHA-DIVERSITY ANALYSIS
ps_richness <- estimate_richness(ps, measures=c("Observed", "Shannon"))
rownames(ps_richness) <- rownames(ps_richness) %>% str_replace_all("[.]", "-")
ps_richness$SampleName <- rownames(ps_richness)
ps_richness <- left_join(ps_richness, metadata)

richness.plot <- ggplot(ps_richness, 
       aes(x = Environment, y = Observed, fill = ExtractionMethod, group = interaction(Environment, ExtractionMethod), shape = ExtractionMethod)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 2) +
  scale_fill_brewer(palette = "Paired", name = "") +
  scale_shape(name = "") +
  ylab("Observed richness") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

richness.plot

shannon.plot <- ggplot(ps_richness, 
       aes(x = Environment, y = Shannon, fill = ExtractionMethod, group = interaction(Environment, ExtractionMethod), shape = ExtractionMethod)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitterdodge(), size = 2) +
  scale_fill_brewer(palette = "Paired", name = "") +
  scale_shape(name = "") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

shannon.plot

### ALPHA-DIVERSITY ANALYSIS (STATS)
# SHANNON
wilcox.test((ps_richness %>% subset(Environment == "HumanFecal" & ExtractionMethod == "PC"))$Shannon, 
       (ps_richness %>% subset(Environment == "HumanFecal" & ExtractionMethod == "KIT"))$Shannon,
       paired = TRUE)$p.value

wilcox.test((ps_richness %>% subset(Environment == "MouseFecal" & ExtractionMethod == "PC"))$Shannon, 
       (ps_richness %>% subset(Environment == "MouseFecal" & ExtractionMethod == "KIT"))$Shannon,
       paired = TRUE)$p.value

wilcox.test((ps_richness %>% subset(Environment == "Soil" & ExtractionMethod == "PC"))$Shannon, 
       (ps_richness %>% subset(Environment == "Soil" & ExtractionMethod == "KIT"))$Shannon,
       paired = TRUE)$p.value

# OBSERVED RICHNESS
wilcox.test((ps_richness %>% subset(Environment == "HumanFecal" & ExtractionMethod == "PC"))$Observed, 
            (ps_richness %>% subset(Environment == "HumanFecal" & ExtractionMethod == "KIT"))$Observed,
            paired = TRUE)$p.value
wilcox.test((ps_richness %>% subset(Environment == "MouseFecal" & ExtractionMethod == "PC"))$Observed, 
            (ps_richness %>% subset(Environment == "MouseFecal" & ExtractionMethod == "KIT"))$Observed,
            paired = TRUE)$p.value
wilcox.test((ps_richness %>% subset(Environment == "Soil" & ExtractionMethod == "PC"))$Observed, 
            (ps_richness %>% subset(Environment == "Soil" & ExtractionMethod == "KIT"))$Observed,
            paired = TRUE)$p.value

### CONTIG OVERLAP ANALYSIS
ps.melted <- ps %>% psmelt() %>% subset(Abundance > 0)

soil.samples <- ps.melted %>% subset(Environment == "Soil")
soil.samples.contigs.PC <- (soil.samples %>% subset(ExtractionMethod == "PC"))$OTU
soil.samples.contigs.KIT <- (soil.samples %>% subset(ExtractionMethod == "KIT"))$OTU

mouse.samples <- ps.melted %>% subset(Environment == "MouseFecal") %>% group_by(ExtractionMethod)
mouse.samples.contigs.PC <- (mouse.samples %>% subset(ExtractionMethod == "PC"))$OTU
mouse.samples.contigs.KIT <- (mouse.samples %>% subset(ExtractionMethod == "KIT"))$OTU

human.samples <- ps.melted %>% subset(Environment == "HumanFecal") %>% group_by(ExtractionMethod)
human.samples.contigs.PC <- (human.samples %>% subset(ExtractionMethod == "PC"))$OTU
human.samples.contigs.KIT <- (human.samples %>% subset(ExtractionMethod == "KIT"))$OTU


contig.summary <- data.frame(
  sample = c("Soil", "HumanFecal", "MouseFecal"),
  shared = c(length(intersect(soil.samples.contigs.PC, soil.samples.contigs.KIT)),
             length(intersect(human.samples.contigs.PC, human.samples.contigs.KIT)),
             length(intersect(mouse.samples.contigs.PC, mouse.samples.contigs.KIT))),
  unique.PC = c(length(setdiff(soil.samples.contigs.PC, soil.samples.contigs.KIT)),
                length(setdiff(human.samples.contigs.PC, human.samples.contigs.KIT)),
                length(setdiff(mouse.samples.contigs.PC, mouse.samples.contigs.KIT))),
  unique.KIT = c(length(setdiff(soil.samples.contigs.KIT, soil.samples.contigs.PC)),
                 length(setdiff(human.samples.contigs.KIT, human.samples.contigs.PC)),
                 length(setdiff(mouse.samples.contigs.KIT, mouse.samples.contigs.PC)))) %>% 
  pivot_longer(cols = c(shared, unique.PC, unique.KIT))

contig.summary$name <- factor(contig.summary$name, levels = c("unique.PC", "unique.KIT", "shared"))

bar.plot <- contig.summary %>%
  ggplot(aes(fill = name, y = value, x = sample)) +
  geom_bar(position="dodge", stat="identity") + # can change dodge to stack
  geom_text(aes(label=value), position=position_dodge(width=0.9), vjust=-0.4) + 
  xlab("Environment") +
  ylab("Number of vOTUs") +
  ylim(0,925) +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"), name = "",  labels = c("Unique PC", "Unique KIT", "Shared")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

bar.plot

### GENERATE FINAL FIGURES
legend.B <- get_plot_component(
  richness.plot + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom"),
  'guide-box-bottom', return_all = TRUE
)

panel.B <- plot_grid(richness.plot, shannon.plot)
panel.B <- plot_grid(panel.B, legend.B, rel_heights = c(1, 0.14), ncol = 1)

fig.2 <- plot_grid(bar.plot, panel.B, nmds.bray.plot, permanova.flextable.grob, 
          labels = c("A", "B", "C"), 
          label_size = 15,
          rel_heights = c(1, 0.9))

fig.2

ggsave2("Figure2.pdf", fig.2, width = 12, height = 8.4)


### ANALYSIS OF UNIQUE VOTUS
phage_cov.joined.filt %>%
  filter(contig %in% c(setdiff(soil.samples.contigs.PC, soil.samples.contigs.KIT), 
                       setdiff(human.samples.contigs.PC, human.samples.contigs.KIT),
                       setdiff(mouse.samples.contigs.PC, mouse.samples.contigs.KIT),
                       setdiff(soil.samples.contigs.KIT, soil.samples.contigs.PC), 
                       setdiff(human.samples.contigs.KIT, human.samples.contigs.PC),
                       setdiff(mouse.samples.contigs.KIT, mouse.samples.contigs.PC))) %>% 
  summarize(mean_abund = mean(relabund_cov), sd_abund = sd(relabund_cov))

phage_cov.joined.filt %>% 
  filter(Environment == "HumanFecal") %>% 
  filter(contig %in% c(setdiff(human.samples.contigs.PC, human.samples.contigs.KIT), setdiff(human.samples.contigs.KIT, human.samples.contigs.PC))) %>% 
  group_by(sample) %>% 
  summarize(sum_relabund = sum(relabund_cov)) %>% 
  summarize(mean = mean(sum_relabund), sd = sd(sum_relabund))

phage_cov.joined.filt %>% 
  filter(Environment == "MouseFecal") %>% 
  filter(contig %in% c(setdiff(mouse.samples.contigs.PC, mouse.samples.contigs.KIT), setdiff(mouse.samples.contigs.KIT, mouse.samples.contigs.PC))) %>% 
  group_by(sample) %>% 
  summarize(sum_relabund = sum(relabund_cov)) %>% 
  summarize(mean = mean(sum_relabund), sd = sd(sum_relabund))

phage_cov.joined.filt %>% 
  filter(Environment == "Soil") %>% 
  filter(contig %in% c(setdiff(soil.samples.contigs.PC, soil.samples.contigs.KIT), setdiff(soil.samples.contigs.KIT, soil.samples.contigs.PC))) %>% 
  group_by(sample) %>% 
  summarize(sum_relabund = sum(relabund_cov)) %>% 
  summarize(mean = mean(sum_relabund), sd = sd(sum_relabund))


### IPHOP HOST PREDICTION ANALYSIS
iphop_out <- read_csv("Host_prediction_to_genus_m90.csv")

length(unique(iphop_out$Virus)) 
length(unique(phage_cov$contig)) # 623 unique vOTUs with host assignment, out of 2,141 total vOTUs

iphop_out.dedup <- iphop_out %>% group_by(Virus) %>% arrange(desc(`Confidence score`)) %>% slice(1)

iphop_out.dedup <- iphop_out.dedup[, c("Virus", "Host genus")]

iphop_out.split <- iphop_out.dedup %>% 
  separate(col = `Host genus`, into = c("domain", "phylum", "class", "order", "family", "genus"), sep = ";")

iphop_out.split$domain <- gsub("d__", "", iphop_out.split$domain)
iphop_out.split$phylum <- gsub("p__", "", iphop_out.split$phylum)
iphop_out.split$class <- gsub("c__", "", iphop_out.split$class)
iphop_out.split$order <- gsub("o__", "", iphop_out.split$order)
iphop_out.split$family <- gsub("f__", "", iphop_out.split$family)
iphop_out.split$genus <- gsub("g__", "", iphop_out.split$genus)

iphop_out.split$phylum <- sapply(strsplit(iphop_out.split$phylum, "_"), `[`, 1)
iphop_out.split$genus <- sapply(strsplit(iphop_out.split$genus, "_"), `[`, 1)

phage_cov.joined.host <- left_join(phage_cov.joined.filt, iphop_out.split, by = c("contig" = "Virus")) %>% 
  filter(!is.na(family))

phage_cov.joined.host %>% 
  group_by(sample, order, phylum, ExtractionMethod, Environment) %>% 
  summarize (total_order_abund = sum(relabund_cov)) %>% view()
