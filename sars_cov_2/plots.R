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

setwd("X:/Real_data/sars_cov_2/")
source("X:/R_Functions/plot_results.R")
source("X:/R_Functions/smoothing_functions.R")

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
# WHO case count data - Asian countries
#############################################################
# Asian countries with negative incidence over the study period
problem <- data.frame()
for (f in list.files("data/who_asia/", pattern = "plot-", full.names = T)) {
  temp_incidence <- read.csv(f, header = T) %>%
    filter(source == "WHO", date < as.Date("2020-03-08")) %>%
    mutate(country = gsub('.*plot-|.csv', "", f), date = as.Date(date))
  if (any(diff(temp_incidence$cases) < 0)) {
    print(gsub('.*plot-|.csv', "", f))
    temp_incidence$incidence <- c(diff(temp_incidence$cases), NA)
    problem <- bind_rows(problem, temp_incidence)
  }
}

ggplot(problem, aes(x = date, y = cases, col = country)) +
  geom_line() +
  theme_light() +
  labs(x = "", y = "Cumulative cases")
ggsave("figures/cumulative_asia.png", width = 5, height = 4)

ggplot(problem, aes(x = date, y = incidence, col = country)) +
  geom_line() +
  theme_light()

#############################################################
# Smooth incidence data and write it to parameterize GLM 
#############################################################
# Most recent sampled tip 
mrst <- max(read.table("data/seq_date.txt", header = T)$date)
mrst <- as.Date(lubridate::date_decimal(mrst))

# Load WHO data----------------------------------------------
incidence <- list()
for (f in list.files("data/", pattern = "plot-", full.names = T)) {
  # Get WHO data
  loc <- gsub('.*plot-|.csv', "", f)
  temp_data <- read.csv(f, header = T) %>%
    filter(source == "WHO") %>%
    mutate(country = loc, date = as.Date(date))
  
  # Compute incidence from the cumulative number of cases
  temp_data$incidence <- c(diff(temp_data$cases), diff(temp_data$cases)[nrow(temp_data)-1])
  incidence[[loc]] <- temp_data
}

# Specific cases of Asia and America
incidence$asia$incidence <- incidence$asia$incidence - incidence$china$incidence
incidence$north_america$incidence <- incidence$north_america$incidence + incidence$south_america$incidence
incidence <- incidence[names(incidence) != "south_america"]

# Smooth incidence
incidence <- lapply(incidence, function(x) {
  x$smooth <- ma(x$incidence)
  x[x$date<=mrst, ]
})

# Plot
g1 = do.call("rbind", incidence) %>%
  mutate(
    country = recode(country, !!!c("europe" = "Europe", "asia" = "Asia",
                                   "africa" = "Africa", "china" = "China", 
                                   "north_america" = "Americas", "oceania" = "Oceania"))
    ) %>%
  ggplot(., aes(x = date, y = smooth, col = country)) +
  geom_line() +
  scale_color_manual(values=continent_colors[1:6]) +
  theme_light() +
  labs(x = "", y = "New cases", col ="")
  # scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #               labels = scales::trans_format("log10", scales::math_format(10^.x)))
#ggsave("figures/moving_average_who.png", height = 3, width = 4)

# Save data for Mascot 
do.call("rbind", incidence) %>%
  select(date, country, smooth) %>%
  filter(!is.na(smooth)) %>%
  mutate(
    date = as.numeric(as.Date(mrst) - date),
    smooth = ifelse(smooth == 0, 0.1, smooth),
    country = recode(country, !!!c("europe" = "Europe", "asia" = "Asia",
                                   "africa" = "Africa", "china" = "China", 
                                   "north_america" = "Americas", "oceania" = "Oceania"))
    ) %>%
  arrange(date, country) %>%
  pivot_wider(names_from = date, values_from = smooth) %>%
  rename(continent = country) %>%
  write.table("data/incidence_who.txt", col.names = T, row.names = F, sep = ",", quote = F)

# Load WID data--------------------------------------------
incidence <- read.csv("data/owid-covid-data.csv", header = T) %>%
  select(continent, location, date, new_cases) %>%
  filter(
    location %in% c("China", "Africa", "Asia", "Europe", "North America", "South America", "Oceania"),
    date < as.Date("2020-03-31")
    ) 

incidenceL <- lapply(unique(incidence$location), function(x) incidence[incidence$location == x, ])
names(incidenceL) <- unique(incidence$location)

# Modify asia and america
incidenceL$Asia$new_cases = incidenceL$Asia$new_cases - incidenceL$China$new_cases

toComplete = data.frame(
  continent = "", 
  location = "South America", 
  new_cases = 0, 
  date = seq(as.Date("2020-01-22"), as.Date("2020-02-21"), by = "day")) 
toComplete$date <- as.character(toComplete$date)
incidenceL$`South America` = rbind(toComplete, incidenceL$`South America`)
incidenceL$`North America`$new_cases = incidenceL$`North America`$new_cases + incidenceL$`South America`$new_cases
incidenceL = incidenceL[names(incidenceL) !="South America"]

# Smooth incidence
incidenceL <- lapply(incidenceL, function(x) {
  x$smooth <- ma(x$new_cases)
  x[x$date<=mrst, ]
})

# Plot
g2 = do.call("rbind", incidenceL) %>%
  ggplot(., aes(x = as.Date(date), y = smooth, col = location)) +
  geom_line() +
  scale_color_manual(values = continent_colors[1:6]) +
  xlim(c(as.Date("2020-01-04"), mrst)) + 
  theme_light() +
  labs(x = "", y = "New cases", col = "")
#ggsave("figures/moving_average_wid.png", height = 4, width = 7)


# Save data for mascot
do.call("rbind", incidenceL) %>%
  select(date, location, smooth) %>%
  filter(!is.na(smooth)) %>%
  mutate(
    date = as.numeric(as.Date(mrst) - as.Date(date)),
    smooth = ifelse(smooth == 0, 0.1, smooth),
    location = recode(location, "North America" = "Americas")
  ) %>%
  arrange(date, location) %>%
  pivot_wider(names_from = date, values_from = smooth) %>%
  rename(continent = location) %>%
  replace(is.na(.), 0.1) %>%
  write.table("data/incidence_wid.txt", col.names = T, row.names = F, sep = ",", quote = F)

# Supplementary figure with moving averaged incidence data
G = ggarrange(g1, g2, ncol = 1, nrow = 2, labels = c("A", "B"), 
              align = "hv", common.legend = T, legend = "bottom")
ggsave("figures/supplementary_figure_moving_average.png", G, height = 6, width = 5)

#############################################################
# Plot smoothed incidence data
#############################################################
incidence_data <- data.frame()

for (f in list.files("data", "incidence_day*", full.names = T)) {
  temp_data <- read.table(f, header = T, sep = ",") %>%
    pivot_longer(!continent, values_to = "count", names_to = "day") %>%
    mutate(day = as.Date(lubridate::date_decimal(2020.1748633879781)) - as.integer(gsub("X", "", day)),
           data = ifelse(grepl("bis", f), "WID", "WHO")) 
  incidence_data = bind_rows(incidence_data, temp_data)
}

incidence_data %>%
  mutate(count = ifelse(count == 0.1, 0, count)) %>%
  ggplot(., aes(x = day, y = count, col = continent)) +
  geom_line() +
  facet_grid(data~., scales = "free_y") +
  theme_light() +
  scale_y_continuous(trans=scales::pseudo_log_trans(base = 10),
                     breaks = c(0, 10^seq(1,5)),
                     labels = c(0, format(10^seq(1,5), scientific = T))) +
  scale_color_manual(values = continent_colors[1:6]) +
  labs(x = "", y = "Daily case counts", col = "")
ggsave("figures/incidence_data.png", height = 5, width = 6)

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

# Write mcc tree 
# write.beast(traveltree, "travel/sars_cov_2_travel.mcc.tree")

# Lineage B1
b1 = MRCA(traveltree, 
          "Germany/BavPat1/2020|EPI_ISL_406862|Germany|Bavaria|Munich|NA|2020-01-28", 
          "Netherlands/Helmond_1363548/2020|EPI_ISL_413574|Netherlands|Helmond|NA|NA|2020-02-29")
travel_lineages = c(travel_lineages, "B1" = b1)
b1 = traveltree@phylo$tip.label[offspring(traveltree, b1)]
b1 = b1[!is.na(b1)]

b1_tree = drop.tip(traveltree, traveltree@phylo$tip.label[!traveltree@phylo$tip.label %in% b1])
write.beast(b1_tree, paste0("figures/b1_travel.nex"))

# Lineage B4
b4 = MRCA(traveltree, 
          "Shandong/IVDC-SD-001/2020|EPI_ISL_408482|China|Shandong|Qingdao|NA|2020-01-19", 
          "Australia/NSW07/2020|EPI_ISL_413214|Australia|New_South_Wales|Sydney|NA|2020-02-29")
travel_lineages = c(travel_lineages, "B4" = b4)
b4 = traveltree@phylo$tip.label[offspring(traveltree, b4)]
b4 = b4[!is.na(b4)]

b4_tree = drop.tip(traveltree, traveltree@phylo$tip.label[!traveltree@phylo$tip.label %in% b4])
write.beast(b4_tree, paste0("figures/b4_travel.nex"))

# Lineage A
a = MRCA(traveltree, 
     "Sichuan/IVDC-SC-001/2020|EPI_ISL_408484|China|Sichuan|Chengdu|NA|2020-01-15", 
     "China/WF0028/2020|EPI_ISL_413791|China|Shandong|NA|NA|2020-02")
travel_lineages = c(travel_lineages, "A" = a)
a = traveltree@phylo$tip.label[offspring(traveltree, a)]
a = a[!is.na(a)]
a_tree = drop.tip(traveltree, traveltree@phylo$tip.label[!traveltree@phylo$tip.label %in% a])
write.beast(a_tree, paste0("figures/a_travel.nex"))

# Lineage A1
a1 = MRCA(traveltree, 
     "USA/WA12-UW8/2020|EPI_ISL_413563|USA|Washington|NA|NA|2020-03-03", 
     "USA/WA9-UW6/2020|EPI_ISL_413487|USA|Washington|NA|NA|2020-03-01")
travel_lineages = c(travel_lineages, "A1" = a1)
a1 = traveltree@phylo$tip.label[offspring(traveltree, a1)]
a1 = a1[!is.na(a1)]
a1_tree = drop.tip(traveltree, traveltree@phylo$tip.label[!traveltree@phylo$tip.label %in% a1])
write.beast(a1_tree, paste0("figures/a1_travel.nex"))

# Posterior probabilities travel mcc tree
posterior_prob_travel = traveltree@data %>%
  select(node, sampleLoc, sampleLoc.prob, posterior) %>%
  filter(node %in% travel_lineages) %>%
  mutate(lineage = recode(node, !!!setNames(names(travel_lineages), travel_lineages)),
         sampleLoc.prob = as.numeric(sampleLoc.prob), model = "travel") %>%
  dplyr::rename("continent" = "sampleLoc", "continent.prob" = "sampleLoc.prob")

#############################################################
# Compare root height
#############################################################
ageroots = data.frame()
for (f in list.files("figures/", "ageroot_", full.names = T)) {
  ageroots = bind_rows(ageroots, read.table(f, header = T, sep = "\t"))
}

ageroots = ageroots %>%
  mutate(model = factor(model, 
                        c("dta", "basta", "mascot", "glm_wid", "glm_who"), 
                        c("CTMC", "BASTA", "MASCOT", "MASCOT-WID", "MASCOT-WHO")))

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
  scale_fill_manual(values=colModels) +
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
  mutate(model = factor(model, c("dta", "basta", "mascot", "glm_wid", "glm_who"), 
                        c("CTMC", "BASTA", "MASCOT", "MASCOT-WID", "MASCOT-WHO")))

geneticParams_ss = geneticParams %>%
  group_by(model, parameter) %>%
  summarise(
    median = median(value),
    hpd_97_5 = hpd(value, conf = 0.95, b = "upper"),
    hpd_2_5 = hpd(value, conf = 0.95, b = "lower")
  )

ggplot(geneticParams, aes(x = model, y = value, group = model, fill = model), col = "grey20") +
  facet_wrap(.~parameter, ncol = 3, scales = "free_y") +
  geom_violin(alpha = 0.5) +
  geom_pointrange(data=geneticParams_ss, aes(x = model, y = median, ymin = hpd_97_5, ymax = hpd_2_5), col = "grey20") +
  scale_fill_manual(values=colModels) +
  theme_light() +
  theme(legend.position = "none") +
  labs(x = "", y = "Median (95% HPD)", fill = "")
ggsave("figures/geneticparams.pdf", width = 8, height = 8)

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
  t_b1_tree = drop.tip(tree, tree@phylo$tip.label[!tree@phylo$tip.label %in% t_b1])
  write.beast(t_b1_tree, paste0("figures/b1_", m, ".nex"))
  
  # Lineage B4
  t_b4 = MRCA(tree, b4) 
  lineage_nodes = c(lineage_nodes, t_b4)
  p_support = c(p_support, as.numeric(tree@data$posterior[tree@data$node == t_b4]))
  t_b4 = tree@phylo$tip.label[offspring(tree, t_b4)]
  t_b4 = sort(t_b4[!is.na(t_b4)])
  if (m == "dta") {
    t_b4_tree = drop.tip(tree, tree@phylo$tip.label[!tree@phylo$tip.label %in% c(t_b4, 
                                                                                 "Wuhan/HBCDC-HB-05/2020|EPI_ISL_412981|China|Hubei|Wuhan|NA|2020-01-18")])
  } else {
    t_b4_tree = drop.tip(tree, tree@phylo$tip.label[!tree@phylo$tip.label %in% t_b4])
  }
  write.beast(t_b4_tree, paste0("figures/b4_", m, ".nex"))
  
  # Lineage A
  t_a = MRCA(tree, a)
  lineage_nodes = c(lineage_nodes, t_a)
  p_support = c(p_support, as.numeric(tree@data$posterior[tree@data$node == t_a]))
  t_a = tree@phylo$tip.label[offspring(tree, t_a)]
  t_a = sort(t_a[!is.na(t_a)])
  t_a_tree = drop.tip(tree, tree@phylo$tip.label[!tree@phylo$tip.label %in% t_a])
  write.beast(t_a_tree, paste0("figures/a_", m, ".nex"))

  # Lineage A1
  t_a1 = MRCA(tree, a1)
  lineage_nodes = c(lineage_nodes, t_a1)
  p_support = c(p_support, as.numeric(tree@data$posterior[tree@data$node == t_a1]))
  t_a1 = tree@phylo$tip.label[offspring(tree, t_a1)]
  t_a1 = sort(t_a1[!is.na(t_a1)])
  t_a1_tree = drop.tip(tree, tree@phylo$tip.label[!tree@phylo$tip.label %in% t_a1])
  write.beast(t_a1_tree, paste0("figures/a1_", m, ".nex"))
  
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
  ggsave(paste0("figures/mcc_", m, ".png"), treeplot, height = 8, width = 10)

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


# Lineage posterior support
posterior_prob %>%
  select(model, posterior, lineage) %>%
  mutate(posterior = round(as.numeric(posterior), 3),
         model = recode(model, "dta" = "CTMC", "basta_majmode" = "BASTA - 1st mode", "basta_minmode" = "BASTA - 2nd mode", 
                        "mascot" = "MASCOT", "glm_wid" = "MASCOT-WID", "glm_who" = "MASCOT-WHO")) %>%
  pivot_wider(names_from = lineage, values_from = posterior) %>%
  write.table(., "tables/posterior_support.csv", sep = ",", quote = F, col.names = T, row.names = F)

# Lineage location posterior support
posterior_prob %>%
  bind_rows(., posterior_prob_travel) %>%
  mutate(continent.prob = paste0(continent, " (", round(as.numeric(continent.prob), 3), ")"),
         model = recode(model, "travel" = "CTMC-TRAVEL", "dta" = "CTMC", "basta_majmode" = "BASTA - 1st mode", 
                        "basta_minmode" = "BASTA - 2nd mode", "mascot" = "MASCOT", "glm_wid" = "MASCOT-WID", 
                        "glm_who" = "MASCOT-WHO")) %>%
  select(model, continent.prob, lineage) %>%
  pivot_wider(names_from = lineage, values_from = continent.prob) %>%
  write.table(., "tables/posterior_location_support.csv", sep = ",", quote = F, col.names = T, row.names = F)

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
ggsave("figures/total_migration_counts.png", arranged, width = 10, height = 10)


# Look at the bimodality (no bimodality for CTMC and TRAVEL-CTMC)
sapply(which(!sumstats$supp %in% "NS" & grepl("MASCOT-WID", sumstats$model)), function(x) {
  p = sumstats$parameter[x]
  m = as.character(sumstats$model[x])
  
  out = dist %>%
    filter(parameter == p, model == m, side == "right") %>%
    ggplot(., aes(x = value, y = n)) +
    geom_bar(stat = "identity") +
    theme_light() +
    labs(x="", y="", title=paste(m, p))
  
  print(out)
})


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
ggsave("figures/posterior_mutationRate_all.png", g1, height = 3, width = 6)


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
ggsave("figures/posterior_mutationRate_basta.png", g2, height = 3, width = 6)
  
arranged = ggarrange(g1, g2, nrow = 2, labels = c("A", "B"), align = "hv", common.legend = T, legend = "bottom")
ggsave("figures/posterior_mutationRate_arranged.png", height = 5, width = 5)


# pdf("figures/migfromChina.pdf", height = 3, width = 6)
# invisible(sapply(16:20, function(x) print(all_plots[[x]])))
# dev.off()
# 
# pdf("figures/migtoChina.pdf", height = 3, width = 6)
# invisible(sapply(c(3,8,13,24,29), function(x) print(all_plots[[x]])))
# dev.off()
# 
# # Plots from China
# fromChina_s = sumstats$source == "China" & sumstats$destination %in% c("Africa", "Asia", "Oceania")
# fromChina_d = distributions$source == "China" & distributions$destination %in% c("Africa", "Asia", "Oceania")
# 
# ggplot() +
#   geom_bar(data = droplevels(distributions[fromChina_d, ]), aes(x = value, fill = model), position = position_dodge()) +
#   coord_flip() +
#   facet_wrap(~ model, nrow = 1)
# 
# 
#   geom_pointrange(data = droplevels(sumstats[fromChina_s, ]),
#                   aes(x = model, y = X50, ymax = X97_5_hpd, ymin = X2_5_hpd, col = model), size = 1, fatten = 1.2) +
#   geom_text(data = droplevels(sumstats[fromChina_s, ]), aes(x = model, label = supp,  y = height), size = 3) +
#   facet_grid(destination~., scales = "free_y") +
#   scale_color_manual(values = colModels) +
#   # scale_y_continuous(position = "right") +
#   # scale_y_continuous(trans = scales::pseudo_log_trans()) +
#   theme_light() +
#   theme(legend.position = "none") +
#   labs(x = "", y = "Total migration counts (median, 95% HPD)", col = "")
# 
# 
# 
# # Plots to China
# toChina_s = sumstats$source %in% c("Africa", "Asia", "Oceania") & sumstats$destination == "China"
# toChina_d = distributions$source %in% c("Africa", "Asia", "Oceania") & distributions$destination == "China"
# 
# ggplot() +
#   geom_violin(data = droplevels(distributions[toChina_d, ]), aes(x = model, y = value)) +
#   geom_pointrange(data = droplevels(sumstats[toChina_s, ]),
#                   aes(x = model, y = X50, ymax = X97_5_hpd, ymin = X2_5_hpd, col = model), size = 1, fatten = 1.2) +
#   geom_text(data = droplevels(sumstats[toChina_s, ]), aes(x = model, label = supp,  y = height), size = 3) +
#   facet_grid(source~., scales = "free_y") +
#   scale_color_manual(values = colModels) +
#   # scale_y_continuous(position = "right") +
#   # scale_y_continuous(trans = scales::pseudo_log_trans()) +
#   theme_light() +
#   theme(legend.position = "none") +
#   labs(x = "", y = "Total migration counts (median, 95% HPD)", col = "")
# 
# 
# 
# 
# 
# gg = ggplot() +
#   geom_violin(data = distributions, aes(x = model, y = value)) +
#   geom_pointrange(data = sumstats,
#                   aes(x = model, y = X50, ymax = X97_5_hpd, ymin = X2_5_hpd, col = model), size = 1, fatten = 1.2) +
#   geom_text(data = sumstats, aes(x = model, label = supp,  y = height), size = 3) +
#   facet_grid(source~destination, scales = "free_y") +
#   scale_color_manual(values = colModels) +
#   # scale_y_continuous(position = "right") +
#   # scale_y_continuous(trans = scales::pseudo_log_trans()) +
#   theme_light() +
#   theme(axis.text.x = element_blank()) +
#   labs(x = "", y = "Total migration counts (median, 95% HPD)", col = "")
# 
# grob <- ggplotGrob(gg)
# idx <- which(grob$layout$name %in% paste0("panel-", 1:6, "-", 1:6))
# for (i in idx) grob$grobs[[i]] <- nullGrob()
# 
# idx <- which(grob$layout$name %in% "axis-b-6")
# grob$layout[idx, c("t", "b")] <- grob$layout[idx, c("t", "b")] - 2
# idx <- which(grob$layout$name %in% "axis-l-1")
# grob$layout[idx, c("l", "r")] <- grob$layout[idx, c("l", "r")] + 2
# grid.draw(grob)

ggsave("figures/total_migration_counts.png", height = 6, width = 10)


#### Plot for Guy
log = read.table("X:/Real_data/sars_cov_2/basta/sars_cov_2_basta3.log.txt", sep = "\t", header = T)
log = log[log$Sample > 5000000, ]

lognames = c("Africa", "Americas", "Asia", "China", "Europe", "Oceania")
names(lognames) = 0:5
mignames = colnames(log)[grepl("Flag", colnames(log))]

qk = 5 / (6*5)
prior_or = qk / (1-qk)

par(mfrow = c(6,5), mar = c(2,2,2,2))
for (m in mignames) {
  N = gsub("[a-zA-Z\\.]+_", "", m)
  for (l in names(lognames)) {
    N = gsub(l, lognames[l], N)
  }
  y = log[[gsub("Flag", "", m)]][log[[m]] == 1]
  pk = mean(log[[m]] == 1) 
  if (pk == 1) pk = 1 - 1/nrow(log)
  post_or = pk / (1 - pk)
  bf = round(post_or / prior_or, 1)
  
  plot(density(y), main = paste(N, bf), xlab = "")
}

par(mfrow = c(6,5), mar = c(2,2,2,2))
for (m in mignames) {
  N = gsub("[a-zA-Z\\.]+_", "", m)
  for (l in names(lognames)) {
    N = gsub(l, lognames[l], N)
  }
  cn = gsub(".*Flag_", "treePrior.Count", m)
  cn = gsub("_", "to", cn)
  y = log[[cn]][log[[m]] == 1]
  pk = mean(log[[m]] == 1) 
  if (pk == 1) pk = 1 - 1/nrow(log)
  post_or = pk / (1 - pk)
  bf = round(post_or / prior_or, 1)
  
  hist(y, main = paste(N, bf), xlab = "")
}


