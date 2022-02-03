library(RevGadgets)
setwd("~/Documents/bom_phy/analyses/dec")

## plot params

#analysis <- "single"
analysis <- "multitrees"

# file <- paste0("output/dec_", analysis, ".model.log")
# params <- readTrace(file)
# plotTrace(params, vars = c("extirpation_rate", "rate_bg"))

## plot tree 

file <- paste0("output/dec_", analysis, ".ase.tre")

source("scripts/makeStateLabels.R")
# labels that correspond to each region/ possible combination of regions
labs <- makeStateLabels("data/bomarea_ranges.nex", 
                        paste0("output/dec_", analysis, "_state_labels.txt"))
dec_example  <- processAncStates(file , state_labels = labs)

ncol <- length(dec_example@state_labels)


# A, C, E, N, S
colors <- c("#94C08B", "#1E88E5", "#004D40", "#D81B60", "#FFC107",
            "#fa78fa", "#8214a0", "#fae6be", "darkblue", "#a0fa82", "red4","darkorange", "turquoise3",
            "slategray1", "thistle1")
#scales::show_col(colors)

state_labels <- dec_example@state_labels
state_labels <- state_labels[c(1,2,7,11,14,
                               3,5,6,8,10,12,13,15,
                               4,9)]
names(colors) <- state_labels

# standard DEC plot
pdf(paste0("figures/plot_MAP_", analysis, ".pdf"), height = 20, width = 16)
plotAncStatesMAP(t = dec_example, node_color = colors, 
                 tip_labels_states_offset = 1,
                 node_size = c(1,4),
                 tip_states_size = 2,
                 tip_labels_states_size = 3.5,
                 tip_labels_size = 4, tip_labels_states = T,
                 cladogenetic = TRUE, tip_labels_offset = 4, 
                 timeline = T) +
  ggplot2::theme(legend.position = c(.18,.75),
                 legend.background = 
                   ggplot2::element_rect(fill="transparent"),
                 legend.text = element_text(size = 15),
                 #legend.title = element_blank(),
                 text = element_text(size = 20))
dev.off()

pdf(paste0("figures/plot_MAP_andes_", analysis, ".pdf"), height = 20, width = 16)
plotAncStatesMAP(t = dec_example, node_color = colors, 
                 tip_labels_states_offset = 1,
                 node_size = c(1,4),
                 tip_states_size = 2,
                 tip_labels_states_size = 3.5,
                 tip_labels_size = 4, tip_labels_states = T,
                 cladogenetic = TRUE, tip_labels_offset = 4, 
                 timeline = T) +
  ggplot2::theme(legend.position = c(.15,.77)) + 
  ggplot2::geom_vline(xintercept = -3.5, linetype = "longdash") +
  ggplot2::annotate(geom = "text", 
                    x = -3.5, y = 93, 
                    label = "Isthmus Rises",
                    hjust = 1,
                    size = 8) +
  
  # andean uplift peaks 
  ggplot2::geom_vline(xintercept = -23, linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -12, linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -4.5, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -23, y = 97, 
                    label = "Andean Pulses",
                    hjust = 1, 
                    size = 8) +
  
  # theme adjustments
  ggplot2::theme(legend.position = c(.18,.75),
                 legend.background = 
                   ggplot2::element_rect(fill="transparent"),
                 legend.text = element_text(size = 15),
                 #legend.title = element_blank(),
                 text = element_text(size = 20))
  
dev.off()

pdf(paste0("figures/plot_pies_", analysis, ".pdf"),height = 20, width = 20)
plotAncStatesPie(t = dec_example, pie_colors = colors, 
                 tip_pie_nudge_x = 0.5, tip_labels_states_offset = 0.25,
                 tip_labels_size = 3, node_pie_size = 0.5,
                 tip_pie_size = 0.4, tip_labels_states = T,
                 cladogenetic = TRUE, tip_labels_offset = 0.75, 
                 timeline = T) +

  # isthmus 
  ggplot2::geom_vline(xintercept = -3.5, linetype = "longdash") +
  ggplot2::annotate(geom = "text", 
                    x = -3.5, y = 97, 
                    label = "Isthmus Rises",
                    hjust = 1,
                    size = 10) +
  
  # andean uplift peaks 
  ggplot2::geom_vline(xintercept = -23, linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -12, linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -4.5, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -23, y = 97, 
                    label = "Andean Pulses",
                    hjust = 1, 
                    size = 10) +
  
  # theme adjustments
  ggplot2::theme(legend.position = c(0.15, 0.75),
                 legend.background = 
                   ggplot2::element_rect(fill="transparent"),
                 legend.text = element_text(size = 20),
                 legend.title = element_blank(),
                 legend.key = element_blank(),
                 text = element_text(size = 30))
dev.off()

## zoom in 
tips <- dec_example@phylo$tip.label
non_bom <- c(tips[grep("Alstroemeria",tips)],
             tips[grep("Luzuriaga",tips)],
             tips[grep("Drymophila",tips)])
dec_example_bom <- treeio::drop.tip(dec_example, tip = non_bom)

pdf(paste0("figures/plot_pies_zoomed_", analysis, ".pdf"),height = 20, width = 20)
plotAncStatesPie(t = dec_example_bom, 
                 pie_colors = colors, 
                 tip_labels_size = 3,
                 cladogenetic = TRUE, 
                 tip_labels_offset = 0.25, 
                 timeline = T) +
  # isthmus 
  ggplot2::geom_vline(xintercept = -3.5, linetype = "longdash") +
  ggplot2::annotate(geom = "text", 
                    x = -3.5, y = 80, 
                    label = "Isthmus Rises",
                    hjust = 1,
                    size = 10) +
  
  # andean uplift peaks 
  ggplot2::geom_vline(xintercept = -23, linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -12, linetype = "dotted") +
  ggplot2::geom_vline(xintercept = -4.5, linetype = "dotted") +
  ggplot2::annotate(geom = "text", 
                    x = -4.5, y = 83, 
                    label = "Andean Pulses",
                    hjust = 1, 
                    size = 10) +
  ggplot2::theme(legend.position = c(0.1, 0.75),
                 legend.background = 
                   ggplot2::element_rect(fill="transparent"),
                 legend.key=element_blank(),
                 text = element_text(size = 30))
dev.off()


#pdf("figures/plot_maps_restricted_02_bom.pdf",height = 20, width = 10)
#plotAncStatesMAP(t = dec_example_bom, node_color = colors, tip_labels_size = 3,
#                 node_labels_as = "state", node_labels_centered = T,
#                 node_labels_offset = 0,
#                 cladogenetic = TRUE, tip_labels_offset = 0, timeline = T) +
#  ggplot2::theme(legend.position = c(0.1, 0.75))
#dev.off()


#pdf("plot_tips.pdf",height = 20, width = 10)
#plotAncStatesMAP(t = dec_example, node_color = colors, 
#                 tip_labels_states = T, tip_labels_states_offset = 0,
#                 node_size = c(0.005, 0.005), tip_labels_states_size = 3,
#                 tip_states_size = 5,
#                 cladogenetic = TRUE, tip_labels_offset = 0.01, timeline = F) +
#  ggplot2::theme(legend.position = c(0.1, 0.75))
#dev.off()
#