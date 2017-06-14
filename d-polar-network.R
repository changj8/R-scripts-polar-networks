library(ggraph)
library(igraph)
library(scales)

nodes <- read.csv("d-nodes.csv")
links <- read.csv("d-edges.csv")

net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)

polar.layout <- function(radii, angles) {
  cbind(radii*cos(angles), -radii*sin(angles))        
}

mRNA <- which(V(net)$type == "mRNA")
miRNA <- which(V(net)$type == "miRNA")
radii <- ifelse(V(net)$type == "miRNA", 1, 2)
angles <- rep.int(0, vcount(net))
angles[mRNA] <- (1:length(mRNA)-1) * 2 * pi / length(mRNA)
angles[miRNA] <- (1:length(miRNA)-1) * 2 * pi / length(miRNA)
layout <- polar.layout(radii, angles)

# Modify label conditions
V(net)$label <- ifelse(V(net)$type=='miRNA' | V(net)$name=='MTOR' |
                         V(net)$name=='RHOA' | V(net)$name=='CDC42'|
                         V(net)$name=='RAC1' | V(net)$name=='RICTOR'|
                         V(net)$name=='RPTR', V(net)$name, NA)

# Vertex size based on attributes
Vsize <- ifelse(V(net)$type=='miRNA', 10, 8)

# color edges of graph based on their source node color
edge.start <- ends(net, es=E(net), names=F)[,1]
edge.col <- V(net)$color[edge.start]

# Get the labels aligned consistently around the edge of the circle for any n of nodes.
# start = offset from 12 o'clock in radians
# direction = 1 for clockwise; -1 for anti-clockwise.

radian.rescale <- function(x, start=0, direction=-1, net) {
  if(V(net)$type == 'mRNA'){
  c.rotate <- function(x) (x + start) %% (2 * pi) * direction
  c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
}

Vdist <- ifelse(V(net)$type == 'mRNA', 0, ifelse(V(net)$name == 'RICTOR', 0, 0))
V.label.size <- ifelse(V(net)$type == 'mRNA', 0.8, 1.2)

lab.locs <- radian.rescale(x=1:91, direction=-1, start=0, net)

# Edge width
ew <- E(net)$weight

#tiff(filename='GFPdown-miRNA-150dpi.tiff', 14*150, 14*150, res=150)

plot(net, layout=layout, edge.arrow.size=0, edge.curved=FALSE, vertex.size=Vsize, vertex.color=V(net)$color,
     edge.color=edge.col, vertex.frame.color=NA, vertex.label.family='Arial', vertex.label.dist=Vdist,
     vertex.label.degree=lab.locs, edge.width=ew, vertex.label.font=2, vertex.label.cex=V.label.size
     )

#dev.off()