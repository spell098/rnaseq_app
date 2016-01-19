updateCytoscapeNetwork = function(network,nodes,edges){
  network$nodeData = merge(network$nodeData,nodes,by = "symbol",all.x=TRUE)
  network$edgeData = merge(network$edgeData,edges,by.y = c("src","dest"),by.x = c("fromNode","toNode"),all.x=TRUE)
}