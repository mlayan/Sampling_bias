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
library(grid)
library(gtable)
library(ggpubr)

Sys.setlocale("LC_TIME", "English")

setwd("sars_cov_2/")
source("../R/plot_results.R")
source("../R/smoothing_functions.R")

colModels <- c("grey80", "#00979c", "mediumorchid4", "#883000", "#fa7850", "gold")
names(colModels) <- c("CTMC-TRAVEL", "CTMC", "BASTA", "MASCOT", "MASCOT-WID", "MASCOT-WHO")

continent_colors <- c("gold", "cornflowerblue", "orange", "red", "forestgreen", "purple",
                     "grey90", "grey70", "grey50", "grey30")
names(continent_colors) <- c("Africa", "Americas", "Asia", "China", "Europe", "Oceania", 
                            "Asia+Africa", "Asia+China", "China+Africa", "Asia+China+Africa")

Oceania <- c("Australia", "NewZealand")      

Europe <- c("Belgium","Finland", "France", "Germany", "Italy", "Luxembourg", "Netherlands",
           "Portugal", "Spain", "Sweden", "Switzerland", "UK")

Americas <- c("Brazil","Canada", "Mexico", "USA")

Asia <- c("Cambodia", "India", "Iran", "Japan", "Nepal", "Singapore", "SouthKorea",
         "Taiwan", "Thailand", "Vietnam")

China <- c("ChinaAnhui", "ChinaBeijing", "ChinaChongqing", "ChinaFujian", "ChinaGuangdong",
          "ChinaHenan", "ChinaHunan", "ChinaHongKong", "ChinaHubei", "ChinaJiangsu","ChinaJiangxi", 
          "ChinaShandong", "ChinaSichuan", "ChinaYunnan", "ChinaZhejiang")

Africa <- "Nigeria"                     


#############################################################
# Convert travel history mcc tree 
#############################################################
traveltree = read.beast("travel/SC2_282_TravelHistory.MCC.tre")
travel_lineages = c()

# Recode sampled locations
for (k in names(continent_colors)[1:6]) {
  traveltree@data$sampleLoc[traveltree@data$sampleLoc %in% get(k)] = k
}

# Recode sampled locations set
traveltree@data$sampleLoc.set = sapply(traveltree@data$sampleLoc.set, function(x) {
  
  for (k in names(continent_colors)[1:6]) {
    x[x %in% get(k)] = k
  }
  x
})

# Lineage B1
b1 = MRCA(traveltree, 
          "Germany/BavPat1/2020|EPI_ISL_406862|Germany|Bavaria|Munich|NA|2020-01-28", 
          "Netherlands/Helmond_1363548/2020|EPI_ISL_413574|Netherlands|Helmond|NA|NA|2020-02-29")
travel_lineages = c(travel_lineages, "B1" = b1)
b1 = traveltree@phylo$tip.label[offspring(traveltree, b1)]
b1 = b1[!is.na(b1)]

# Lineage B4
b4 = MRCA(traveltree, 
          "Shandong/IVDC-SD-001/2020|EPI_ISL_408482|China|Shandong|Qingdao|NA|2020-01-19", 
          "Australia/NSW07/2020|EPI_ISL_413214|Australia|New_South_Wales|Sydney|NA|2020-02-29")
travel_lineages = c(travel_lineages, "B4" = b4)
b4 = traveltree@phylo$tip.label[offspring(traveltree, b4)]
b4 = b4[!is.na(b4)]

# Lineage A
a = MRCA(traveltree, 
     "Sichuan/IVDC-SC-001/2020|EPI_ISL_408484|China|Sichuan|Chengdu|NA|2020-01-15", 
     "China/WF0028/2020|EPI_ISL_413791|China|Shandong|NA|NA|2020-02")
travel_lineages = c(travel_lineages, "A" = a)
a = traveltree@phylo$tip.label[offspring(traveltree, a)]
a = a[!is.na(a)]

# Lineage A1
a1 = MRCA(traveltree, 
     "USA/WA12-UW8/2020|EPI_ISL_413563|USA|Washington|NA|NA|2020-03-03", 
     "USA/WA9-UW6/2020|EPI_ISL_413487|USA|Washington|NA|NA|2020-03-01")
travel_lineages = c(travel_lineages, "A1" = a1)
a1 = traveltree@phylo$tip.label[offspring(traveltree, a1)]
a1 = a1[!is.na(a1)]

# Posterior probabilities travel mcc tree
posterior_prob_travel = traveltree@data %>%
  select(node, sampleLoc, sampleLoc.prob, posterior) %>%
  filter(node %in% travel_lineages) %>%
  mutate(lineage = recode(node, !!!setNames(names(travel_lineages), travel_lineages)),
         sampleLoc.prob = as.numeric(sampleLoc.prob), model = "travel") %>%
  dplyr::rename("continent" = "sampleLoc", "continent.prob" = "sampleLoc.prob")


#############################################################
# Plot MCC Trees
#############################################################
# Most recent sampled tip
mrsd = read.table("data/seq_date.txt", header = T)$date
mrsd = max(mrsd)
mrsd = strftime(lubridate::date_decimal(mrsd),  format="%Y-%m-%d")
lab = c("A", "B", "C", "D", "E", "F")
names(lab) = c("dta", "basta_majmode", "basta_minmode", "mascot", "glm_wid", "glm_who")

posterior_prob = data.frame()

for (m in c( "dta", "basta_majmode", "basta_minmode", "mascot", "glm_wid", "glm_who")) {
  
  tree = read.beast(paste0(gsub("_m.*$", "", m), "/sars_cov_2_", m, ".mcc.tree"))  
  colnames(tree@data) = gsub("location","continent", colnames(tree@data))
  colnames(tree@data) = gsub("typeTrait","continent", colnames(tree@data))
  colnames(tree@data) = gsub("\\.states","", colnames(tree@data))
  
  tree@data$continent.prob = as.numeric(tree@data$continent.prob) 
  
  # Piechart of root location
  ancestorlocs = mapply(
    function(x,y) {
      nI = length(x)
      nT = length(continent_colors)
      x = as.numeric(x)
      names(x) = y
      if (nI < nT) {
        x = c(x, rep(0, nT-nI))
        names(x)[(nI+1):nT] = names(continent_colors)[!names(continent_colors) %in% y]
      }
      
      as.data.frame(t(x))
    },
    x=tree@data$continent.set.prob,
    y=tree@data$continent.set,
    SIMPLIFY=FALSE
  )
  ancestorlocs = do.call("rbind", ancestorlocs) 
  ancestorlocs$node = tree@data$node
  ancestorlocs[ancestorlocs$node != tree@phylo$Nnode+2, names(continent_colors)] = NA 
  ancestorlocs = ggtree::nodepie(ancestorlocs, cols = 1:6, color = continent_colors)

  
  # Verify the existence of lineages
  # Lineage B1
  p_support =c()
  lineage_nodes = c()
  t_b1 = MRCA(tree, b1)
  lineage_nodes = c(lineage_nodes, t_b1)
  p_support = as.numeric(tree@data$posterior[tree@data$node == t_b1])
  t_b1 = tree@phylo$tip.label[offspring(tree, t_b1)]
  t_b1 = sort(t_b1[!is.na(t_b1)])
  
  # Lineage B4
  t_b4 = MRCA(tree, b4) 
  lineage_nodes = c(lineage_nodes, t_b4)
  p_support = c(p_support, as.numeric(tree@data$posterior[tree@data$node == t_b4]))
  t_b4 = tree@phylo$tip.label[offspring(tree, t_b4)]
  t_b4 = sort(t_b4[!is.na(t_b4)])
  
  # Lineage A
  t_a = MRCA(tree, a)
  lineage_nodes = c(lineage_nodes, t_a)
  p_support = c(p_support, as.numeric(tree@data$posterior[tree@data$node == t_a]))
  t_a = tree@phylo$tip.label[offspring(tree, t_a)]
  t_a = sort(t_a[!is.na(t_a)])

  # Lineage A1
  t_a1 = MRCA(tree, a1)
  lineage_nodes = c(lineage_nodes, t_a1)
  p_support = c(p_support, as.numeric(tree@data$posterior[tree@data$node == t_a1]))
  t_a1 = tree@phylo$tip.label[offspring(tree, t_a1)]
  t_a1 = sort(t_a1[!is.na(t_a1)])
  
  names(p_support) = c("B1", "B4", "A", "A1")
  out = paste0(m, ": ")
  if (!identical(t_b1, sort(b1))) out = paste0(out, " lineage B1")
  if (!identical(t_b4, sort(b4))) out = paste0(out, " lineage B4")
  if (!identical(t_a, sort(a))) out = paste0(out, " lineage A")
  if (!identical(t_a1, sort(a1))) out = paste0(out, " lineage A1")
  cat(paste0(out, "\n"))
  print(p_support)
  if (sum(!t_b1 %in% b1)) {cat("B1 "); lineage_nodes = lineage_nodes[names(p_support) != "B1"]; p_support = p_support[names(p_support) != "B1"]}
  if (sum(!t_b4 %in% b4)) {cat("B4 "); lineage_nodes = lineage_nodes[names(p_support) != "B4"]; p_support = p_support[names(p_support) != "B4"]}
  if (sum(!t_a %in% a)) {cat("A "); lineage_nodes = lineage_nodes[names(p_support) != "A"]; p_support = p_support[names(p_support) != "A"]}
  if (sum(!t_a1 %in% a1)) {cat("A1 "); lineage_nodes = lineage_nodes[names(p_support) != "A1"]; p_support = p_support[names(p_support) != "A1"]}
  cat("\n")
  
  # Save posterior probabilities
  tmp_posterior_prob = tree@data[tree@data$node %in% lineage_nodes, c("node", "continent.prob", "continent", "posterior")]
  tmp_posterior_prob$model = m
  tmp_posterior_prob$lineage = sapply(tmp_posterior_prob$node, function(x) names(p_support)[lineage_nodes == x])
  posterior_prob = bind_rows(posterior_prob, tmp_posterior_prob)
  
  # Plot tree with piechart displaying root location probabilities
  rootcol = continent_colors[tree@data$continent[tree@data$node == tree@phylo$Nnode+2]]
  tree@data = tree@data[order(as.numeric(tree@data$node)), ]
  tree@data$lineages = apply(tree@data, 1, function(x) {
    if (x[["node"]] %in% lineage_nodes) {
      paste0(names(p_support)[lineage_nodes == x[["node"]]], ": ", round(as.numeric(x[["posterior"]]), 1))
    } else {
      NA
    }
  })
  treeplot = ggtree(tree, aes(col = continent, size = continent.prob), mrsd = mrsd, right=T, ladderize = T) +
    theme_tree2() +
    geom_tippoint(size = 2) +
    geom_nodepoint(size = 2) +
    geom_nodelab(aes(label=lineages, subset = (node %in% lineage_nodes), col = "black"),
                 nudge_x = 0.005, size = 3) +
    geom_rootedge(rootedge = 0.01, col = rootcol, size = 0.8) +
    geom_inset(ancestorlocs, width=.1, height=.1, vjust=-13, hjust = 0.01) +
    scale_color_manual(values = continent_colors) +
    scale_size_continuous(range = c(0.05, 0.8)) +
    ggplot2::labs(col = "") +
    guides(size = "none")
  treeplot = ggarrange(treeplot, labels = lab[m])
  ggsave(paste0("figures/mcc_", m, ".pdf"), treeplot, height = 8, width = 10)
}

#############################################################
# Plot markov jumps
#############################################################
# Get markov jumps from CTMC-TRAVEL
possibleTransitions = expand.grid(c("Africa", "Americas", "Asia", "China", "Europe", "Oceania"),
                                  c("Africa", "Americas", "Asia", "China", "Europe", "Oceania"))
possibleTransitions = apply(possibleTransitions[possibleTransitions[,1] != possibleTransitions[,2],], 1, paste, collapse="_")
possibleTransitions = as.character(possibleTransitions)

dist_travel = read.table("travel/covid10Mar_282.Ed.jumpHistory.log", header = T, sep = "\t")
dist_travel = dist_travel[round(0.1*nrow(dist_travel)):nrow(dist_travel), ]
dist_travel$completeHistory_1 = gsub(" [0-9]+", "", dist_travel$completeHistory_1)
dist_travel = lapply(dist_travel$completeHistory_1, function(x) {
  allTransitions = strsplit(x, "},{", fixed = T)
  allTransitions = sapply(allTransitions, function(x) gsub('.*[0-9],|\\}}', "", x))
  sources = sapply(allTransitions, function(x) strsplit(x, ",", fixed = T)[[1]][1])
  destinations = sapply(allTransitions, function(x) strsplit(x, ",", fixed = T)[[1]][2])
  transitions = data.frame(source = sources, destination = destinations)
  for (loc in c("Africa", "Americas", "Asia", "China", "Europe", "Oceania")) {
    transitions$source[transitions$source %in% get(loc)] = loc
    transitions$destination[transitions$destination %in% get(loc)] = loc
  }
  
  out = transitions %>%
    group_by(source, destination) %>%
    summarise(value = n(), .groups = "keep") %>%
    mutate(parameter = paste0(source, "_", destination)) %>%
    ungroup() %>%
    select(value, parameter)
  
  if (any(!possibleTransitions %in% out$parameter)) {
    out_temp = data.frame(
      parameter = possibleTransitions[!possibleTransitions %in% out$parameter],
      value = 0
      )
    out = bind_rows(out, out_temp)
  }
  
  return(out)
})
dist_travel = do.call("rbind", dist_travel) %>%
  mutate(model = "travel")


sumstats_travel = dist_travel %>%
  group_by(parameter) %>%
  summarise(
    X50	= median(value), 
    X97_5_hpd = hpd(value, 0.95, "upper"),
    X2_5_hpd = hpd(value, 0.95, "lower")
  ) %>%
  mutate(model = "travel")

# Load summary statistics
sumstats = data.frame()
for (f in list.files("figures/", pattern = "output_", full.names = T)) {
  sumstats_temp = read.table(f, header=T, sep="\t") %>%
    filter( (grepl("glm", model) & grepl("nMigration", parameter)) | (!grepl("glm", model) & grepl("bssvs_nMigration", parameter))) 
  if (!grepl("glm", f)) sumstats_temp = mutate(sumstats_temp, supp = ifelse(BF < 3, "NS", ""))
  sumstats = bind_rows(sumstats, sumstats_temp)
}
sumstats = bind_rows(sumstats, sumstats_travel) %>%
  mutate(parameter = gsub(".*nMigration_", "", parameter))

# Load distributions data 
dist = data.frame()
for (f in list.files("figures/", pattern = "distribution_", full.names = T)) {
  md = gsub(".*distribution_|.txt$", "", f)
  dist_to_remove = sumstats$parameter[sumstats$model == md & sumstats$BF < 3]
  dist_to_remove = gsub("bssvs_forwards_", "", dist_to_remove)
  
  dist_temp = read.table(f, header = T, sep = "\t") %>%
    mutate(parameter = gsub("nMigration_", "", parameter)) %>%
    filter(!parameter %in% dist_to_remove)
  dist = bind_rows(dist, dist_temp)
} 

dist = bind_rows(dist, dist_travel) %>%
  group_by(model, parameter, value) %>%
  summarise(n = n())

for (mod in unique(dist$model)) {
  for (p in unique(dist$parameter[dist$model == mod])) {
    dist$n[dist$model == mod & dist$parameter == p] = dist$n[dist$model == mod & dist$parameter == p] / sum(dist$n[dist$model == mod & dist$parameter == p])
  }
}
dist$side = "right"

# Create left side of the distributions
dist_left = dist
dist_left$n = -dist$n
dist_left$side = "left"
dist = bind_rows(dist, dist_left)

# Change model names
dist = dist %>%
  mutate(model = recode(model, "basta" = "BASTA", "dta" ="CTMC", "glm_who" = "MASCOT-WHO",
                        "glm_wid" = "MASCOT-WID", "mascot" = "MASCOT", "travel" = "CTMC-TRAVEL")) %>%
  mutate(model = factor(model, names(colModels)))

sumstats = sumstats %>%
  mutate(model = recode(model, "basta" = "BASTA", "dta" ="CTMC", "glm_who" = "MASCOT-WHO",
                        "glm_wid" = "MASCOT-WID", "mascot" = "MASCOT", "travel" = "CTMC-TRAVEL")) %>%
  mutate(model = factor(model, names(colModels)))
  
# Create matrix panels
allPlots = list()
i = 1
for (s in c("Africa", "Americas", "Asia", "China", "Europe", "Oceania")) {
  for (d in c("Africa", "Americas", "Asia", "China", "Europe", "Oceania")) {
    
    if (s == d) {
      allPlots[[i]] = NULL
    } else {
      
      sumstats_sub = sumstats %>%
        filter(parameter == paste0(s, "_", d)) %>%
        mutate(X50 = ifelse(supp %in% "NS", NA, X50), 
               X2_5_hpd = ifelse(supp %in% "NS", NA, X2_5_hpd),
               X97_5_hpd = ifelse(supp %in% "NS", NA, X97_5_hpd))
      
      notSupported = as.character(sumstats_sub$model[sumstats_sub$supp %in% "NS"])
      
      dist_sub = dist %>%
        filter(parameter == paste0(s, "_", d)) %>%
        mutate(n = ifelse(model %in% notSupported, NA, n))
      
      allPlots[[i]] = ggplot(dist_sub, aes(x = value, y = n, fill = model)) +
        coord_flip() +
        geom_bar(stat = "identity", alpha = 0.8) +
        geom_pointrange(data = sumstats_sub, aes(x = X50, xmin = X2_5_hpd, xmax = X97_5_hpd, y = 0), 
                        col = "grey30", alpha = 0.9, fatten = 1.5) +
        geom_text(data = sumstats_sub, aes(x = 0, y = 0, label = supp)) +
        facet_grid(cols = vars(model)) +
        scale_fill_manual(values = colModels) +
        theme_light() +
        scale_y_continuous(sec.axis = dup_axis()) +
        scale_x_continuous(sec.axis = dup_axis()) +
        theme(axis.title = element_blank(), 
              axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(),
              axis.ticks.y.left = element_line(color = "black"),
              axis.ticks.y.right = element_blank(),
              axis.text.y.right = element_blank(),
              panel.grid = element_blank(),
              strip.text = element_blank(),
              panel.spacing.x=unit(0, "lines"),
              panel.border = element_blank(),
              axis.line = element_line(color = "black"),
              axis.line.x.bottom = element_line(color = "black"),
              axis.line.x.top = element_line(color = "black"),
              axis.line.y.left = element_line(color = "black"),
              axis.line.y.right = element_line(color = "black"),
              plot.margin=unit(c(0,0,0,0),"mm"),
              legend.position = "none"
              ) + 
        labs(fill = "") +
        rremove("ylab") +
        rremove("xlab") 
      
    }
    
    i = i+1
  }
}

# Create total matrix
arranged = ggarrange(plotlist = allPlots, ncol = 6, nrow = 6, align = "hv",
                     common.legend = T, legend = "bottom")
ggsave("figures/total_migration_counts.pdf", arranged, width = 12, height = 12)

#############################################################
# Plot mutation rate for all runs 
#############################################################
mutationRate = data.frame()
for (f in list.files("figures/", pattern = "posterior_.*.txt", full.names = T)) {
  mutationRate_temp = read.table(f, header = T, sep = "\t")
  model = gsub(".*posterior_|.txt", "", f)
  mutationRate_temp$model = model
  mutationRate_temp = mutationRate_temp[, c("model", 'kappa','gammaShape','pInv','clockRate')]
  mutationRate = bind_rows(mutationRate, mutationRate_temp)
}

g1 = mutationRate %>%
  mutate(model = factor(model, c("dta", "basta", "mascot", "glm_wid", "glm_who"), names(colModels)[-1])) %>%
  ggplot(., aes(x = clockRate, fill = model)) +
  geom_density(alpha = 0.4) +
  theme_light() +
  scale_fill_manual(values = colModels[-1]) +
  labs(x = "",  y ="", fill = "") +
  rremove("xlab")


g2 = mutationRate %>%
  filter(model == "basta") %>%
  mutate(model = "BASTA") %>%
  ggplot(., aes(x = clockRate, fill = model)) +
  geom_density(alpha = 0.5) +
  theme_light() +
  geom_vline(xintercept = 11e-4) +
  geom_vline(xintercept = 10e-4) +
  scale_fill_manual(values = colModels["BASTA"]) +
  xlim(c(min(mutationRate$clockRate), max(mutationRate$clockRate))) +
  labs(x = "Evolutionary rate (substitution/site/year)",  y ="", fill = "")
  
arranged = ggarrange(g1, g2, nrow = 2, labels = c("A", "B"), align = "hv", common.legend = T, legend = "bottom")
ggsave("figures/posterior_mutationRate_arranged.png", height = 5, width = 5)
