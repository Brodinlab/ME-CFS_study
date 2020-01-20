require(visNetwork)
library(igraph)
devtools::install_github("datastorm-open/visNetwork")
shiny::runApp(system.file("shiny", package = "visNetwork"))

# Load edge file (to/from)
v_edges <- read.csv('visnetwork_edges_filter.csv', sep = ';')
# Load node file (all nodes plus relative info.)
v_nodes <- read.csv('all_AT_nodes.csv', sep = ';') 

# Custom build color gradient 
rbPal <- colorRampPalette(c('orange', '#BCB4C1', 'purple'))
myBreaks <- c(seq(min(v_nodes$log2.8AT.), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(v_nodes$log2.8AT.)/paletteLength, max(v_nodes$log2.16AT.), length.out=floor(paletteLength/2)))

v_nodes$Colors <- rbPal(50)[as.numeric(cut(v_nodes$log2.16AT., breaks = myBreaks))]


edges <- as.data.frame(v_edges)
nodes <- as.data.frame(v_nodes, color = list(background = v_nodes$Colors))

# Write and load new dataframe with colors based on log2 expression
write.csv(as.data.frame(nodes), file='AT0_nodes_fgrad.csv')
v_nodes0AT <- read.csv('AT0_nodes_fgrad.csv', sep = ';')

# Build network one at a time for each timepoint
visNetwork(v_nodes0AT, edges) %>% 
  visGroups(groupname = "Transcriptional regulators", shape = "diamond") %>% 
  visGroups(groupname = "Regulatory targets", shape = "dot") %>% 
  visLegend() %>% 
  visIgraphLayout(layout = "layout_nicely", smooth = TRUE) 

