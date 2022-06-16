#############################################################
##          SELECT RUNS USING A THRESHOLD ON ESS
##
## Creation Date: 18-06-2020
## Last Update: -06-2020
## Maylis Layan
#############################################################

rm(list = ls())
library(ape)
library(treeio)
library(ggtree)
library(viridis)
library(ggimage)
library(tidyverse)
library(lubridate)
library(ggpubr)

setwd("X:/Real_data/rabies/")
source("X:/R_Functions/plot_results.R")

isle_colors = c("cornflowerblue", "gold", "orange", "red", "forestgreen", "purple")
names(isle_colors) = c("Catanduanes", "Cebu", "Luzon", "Mindanao", "NegrosOriental", "OrientalMindoro")

mrst = max(read.table("data/tohma_dates.txt", header = T)$Date)

model_col = c("#00979c", "mediumorchid4", "mediumorchid2", "#883000")
names(model_col) = c("DTA", "BASTA", "BASTA DEMES", "MASCOT")

#############################################################
# Compare root height
#############################################################
ageroots = data.frame()
for (f in list.files("figures/", "ageroot_", full.names = T)) {
  ageroots = bind_rows(ageroots, read.table(f, header = T, sep = "\t"))
}

ageroots = ageroots %>%
  mutate(model = factor(model, c("dta", "basta", "mascot"), c("DTA", "BASTA", "MASCOT")))

ageroots_ss = ageroots %>%
  group_by(model) %>%
  summarise(
    median = median(ageRoot),
    hpd_97_5 = hpd(ageRoot, conf = 0.95, b = "upper"),
    hpd_2_5 = hpd(ageRoot, conf = 0.95, b = "lower")
    )

ggplot(ageroots, aes(x = model, y = ageRoot, group = model, fill = model), col = "grey20") +
  geom_violin(alpha = 0.5) +
  geom_pointrange(data=ageroots_ss, aes(x = model, y = median, ymin = hpd_97_5, ymax = hpd_2_5), col = "grey20") +
  scale_fill_manual(values=model_col) +
  theme_light() +
  theme(legend.position = "none") +
  labs(x = "", y = "Root age", fill = "")
ggsave("figures/ageroot.pdf", width = 3, height = 3)


#############################################################
# Compare genetic parameters
#############################################################
geneticParams = data.frame()
for (f in list.files("figures/", "geneticparams_", full.names = T)) {
  geneticParams = bind_rows(geneticParams, read.table(f, header = T, sep = "\t"))
}

geneticParams = geneticParams %>%
  mutate(model = factor(model, c("dta", "basta", "mascot"), c("DTA", "BASTA", "MASCOT")))

geneticParams_ss = geneticParams %>%
  group_by(model, parameter) %>%
  summarise(
    median = median(value),
    hpd_97_5 = hpd(value, conf = 0.95, b = "upper"),
    hpd_2_5 = hpd(value, conf = 0.95, b = "lower")
  ) %>%
  filter(parameter %in% c("gammaShape", "kappa", "ucldMean", "ucldStdev"))

geneticParams %>%
  filter(parameter %in% c("gammaShape", "kappa", "ucldMean", "ucldStdev")) %>%
  ggplot(., aes(x = model, y = value, group = model, fill = model), col = "grey20") +
  facet_wrap(.~parameter, ncol = 2, scales = "free_y") +
  geom_violin(alpha = 0.5, scale = "width") +
  geom_pointrange(data=geneticParams_ss, aes(x = model, y = median, ymin = hpd_97_5, ymax = hpd_2_5), col = "grey20") +
  scale_fill_manual(values=model_col) +
  theme_light() +
  theme(legend.position = "none") +
  labs(x = "", y = "Median (95% HPD)", fill = "")
ggsave("figures/geneticparams.pdf", width = 8, height = 8)

#############################################################
# Compare lineages identified 
#############################################################
getLineageSupport = function(tree, name, t1, t2) {
  require(tidytree)
  require(ape)
  t = MRCA(tree, t1, t2)
  p_support = as.numeric(tree@data$posterior[tree@data$node == t])
  t = tree@phylo$tip.label[offspring(tree, t)]
  t = sort(t[!is.na(t)])
  t_tree = drop.tip(tree, tree@phylo$tip.label[!tree@phylo$tip.label %in% t])
  write.beast(t_tree, paste0(name, ".nex"))
  list(p = p_support, taxa = t)
}

lineages = data.frame(
  name = c("NegrosOriental", "Cebu", "Mindanao", "Luzon"),
  t1 = c("AB563950|Negros_Oriental|2009.17534246575", "AB563914|Cebu|2008.85792349727", 
         "AB563807|Mindanao|2005.26301369863", "AB683587|Luzon|2004.61202185792"),
  t2 = c("AB563949|Negros_Oriental|2009.11506849315", "AB563893|Cebu|2008.6912568306",
         "AB563921|Mindanao|2008.94262295082", "AB733977|Luzon|2009.8")
  
)

mrsd = strftime(lubridate::date_decimal(mrst),  format="%Y-%m-%d")
lab = c("A", "B", "C")
names(lab) = c("dta", "basta", 'mascot')

taxa = list()

for (m in c("basta", "dta", "mascot")) {
  tree = read.beast(paste0(m, "/tohma_", m, ".mcc.tree"))
  colnames(tree@data) = gsub("typeTrait|traits","location", colnames(tree@data))
  colnames(tree@data) = gsub("\\.states","", colnames(tree@data))
  
  tree@data$location = gsub("_| ", "", tree@data$location)
  tree@data$location.set = sapply(tree@data$location.set, function(x) gsub("_| ", "", x))
  tree@data$location.prob = as.numeric(tree@data$location.prob) 
  
  # Piechart of root location
  ancestorlocs = mapply(
    function(x,y) {
      nI = length(x)
      nT = length(isle_colors)
      x = as.numeric(x)
      names(x) = y
      if (nI < nT) {
        x = c(x, rep(0, nT-nI))
        names(x)[(nI+1):nT] = names(isle_colors)[!names(isle_colors) %in% y]
      }
      
      as.data.frame(t(x))
    },
    x=tree@data$location.set.prob,
    y=tree@data$location.set,
    SIMPLIFY=FALSE
  )
  ancestorlocs = do.call("rbind", ancestorlocs) 
  ancestorlocs$node = tree@data$node
  ancestorlocs[ancestorlocs$node != tree@phylo$Nnode+2, names(isle_colors)] = NA 
  ancestorlocs = ggtree::nodepie(ancestorlocs, cols = 1:6, color = isle_colors)
  #lapply(tree@data$continent.states.set, function(x) ifelse(length(x) < 6, c(x, rep(0, 6-length(x))), x))
  
  # Verify the existence of lineages
  p_support = apply(lineages, 1, function(x) getLineageSupport(tree, paste0("figures/", x[["name"]], "_", m), x[["t1"]], x[["t2"]]))
  taxa[[m]] = sapply(p_support, function(x) x$taxa)
  lineage_nodes = apply(lineages[, c(2,3)], 1, function(x) MRCA(tree, x[1], x[2]))
  
  cat(paste(m, ":", paste(sapply(p_support, function(x) x$p), collapse = " "), "\n"))
  
  # tree
  rootcol = isle_colors[tree@data$location[tree@data$node == tree@phylo$Nnode+2]]
  tree@data = tree@data[order(as.numeric(tree@data$node)), ]
  
  treeplot = ggtree(tree, aes(col = location, size = location.prob), mrsd = mrsd, ladderize = T, right=T) +
    theme_tree2() +
    # geom_tiplab(size=2, align=TRUE, linesize=.5, aes(col = "black")) +
    geom_tippoint(size = 2) +
    geom_nodepoint(size = 2) +
    geom_nodelab(aes(label=round(as.numeric(posterior), 2), subset = (node %in% lineage_nodes), col = "black"), 
                 nudge_x = 0.8, size = 3) + # nudge_y = 1.5, size=3) + #vjust=-.5, ) +
    geom_rootedge(rootedge = 1, col = rootcol, size = 0.8) +
    geom_inset(ancestorlocs, width=.1, height=.1, vjust=-13) +
    scale_color_manual(values = isle_colors) +
    scale_size_continuous(range = c(0.05, 0.7)) +
    ggplot2::labs(col = "") +
    guides(size = "none")
  treeplot = ggarrange(treeplot, labels = lab[m])

  ggsave(paste0("figures/mcc_", m, ".pdf"), treeplot, height = 7, width = 7)
  ggsave(paste0("figures/mcc_", m, ".png"), treeplot, height = 7, width = 7)
   
  # Node posterior support
  node_posterior = tree@data %>%
    mutate(node = as.numeric(node), posterior = as.numeric(posterior)) %>%
    filter(node > tree@phylo$Nnode+1) %>%
    ggplot(., aes(x = factor(trunc(posterior*10)/10))) +
      geom_bar() +
      theme_light() +
    labs(x = "Posterior", y = "Count")
  ggsave(paste0("figures/posterior_", m, ".pdf"), node_posterior, height = 5, width = 5)
  
  # Node location posterior support
  node_location_post = tree@data %>%
    mutate(node = as.numeric(node)) %>%
    filter(node > tree@phylo$Nnode+1) %>%
    ggplot(., aes(x = location.prob)) +
    geom_histogram(binwidth = 0.1, colour = "grey50", fill = "white") +
    theme_light() +
    scale_x_continuous(breaks = seq(0,1,0.1), limits = c(-0.05,1.05)) + #, limits = c(-0.1,1.1)) +
    labs(x = "Location posterior", y = "Count")
  ggsave(paste0("figures/posterior_location_", m, ".pdf"), node_posterior, height = 5, width = 5)
  
}

#############################################################
# Verify that the lineages are complete 
#############################################################
for (i in 1:4) print(sum(!Reduce(`&`,Map(`==`, list(taxa$basta[[i]], taxa$mascot[[i]], taxa$dta[[i]]), list(taxa$basta[[i]])))))

#############################################################
# Plot markov jumps
#############################################################
# Load and extract markov jump history travel history
allPossibleMigs = expand.grid(names(isle_colors), names(isle_colors))
allPossibleMigs = apply(allPossibleMigs, 1, paste, collapse = "_")

outputs = data.frame()
for (f in list.files("figures/", "output_", full.names = T)) {
  outputs = bind_rows(outputs, read.table(f, header = T, sep = "\t"))
}


outputs %>%
  filter(grepl("bssvs_nMigration", parameter)) %>%
  mutate(
    source = gsub("bssvs_nMigration_|_[a-zA-Z]+$", "", parameter), 
    destination = gsub("^.*_", "", parameter),
    X50 = ifelse(BF < 3, NA, X50),
    X2_5_hpd = ifelse(BF < 3, NA, X2_5_hpd),
    X97_5_hpd = ifelse(BF < 3, NA, X97_5_hpd),
    supp = ifelse(BF < 3, "NS", "")
  ) %>%
  mutate(model = factor(model, c("dta", "basta", "basta_demes", "mascot"), names(model_col))) %>%
  ggplot(., aes(x = model, y = X50, ymax = X97_5_hpd, ymin = X2_5_hpd, col = model, label = supp)) +
  geom_pointrange(size = 1, fatten = 1.2) +
  geom_text(aes(y = 0), size = 3) +
  facet_grid(source~destination, switch = "y", scales = "free_y") +
  scale_color_manual(values = model_col) +
  theme_light() +
  theme(axis.text.x = element_blank()) +
  labs(x = "", y = "Total migration counts (median, 95% HPD)", col = "")

ggsave("figures/total_migration_counts.pdf", height = 7, width = 9)


M = outputs %>%
  filter(grepl("bssvs_nMigration", parameter))
M = max(M$X50, na.rm = T)

for (m in c("dta", "basta", "basta_demes", "mascot")) {
  outputs %>%
    filter(grepl("bssvs_nMigration", parameter), model == m) %>%
    mutate(
      source = gsub("bssvs_nMigration_|_[a-zA-Z]+$", "", parameter), 
      destination = gsub("^.*_", "", parameter),
      median = ifelse(BF < 3, NA, X50)
    ) %>%
    ggplot(., aes(x = source, y = destination, fill = median)) +
    # facet_wrap(. ~model, ncol = 2, nrow = 2) +
    geom_tile(height = 0.95, width = 0.95) +
    scale_y_discrete(position = "right", limits=rev) +
    scale_fill_gradient(low = "#56B1F7", high = "#132B43", na.value = "grey90", limits = c(0,M)) +
    theme_light() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.grid = element_blank()) +
    labs(x = "", y = "", col = "")
  ggsave(paste0("figures/migration_heatmap_", m, ".pdf"), height = 3, width = 4)
}




