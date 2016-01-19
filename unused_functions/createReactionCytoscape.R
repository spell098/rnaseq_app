createReactionCytoscape = function(network,nodes,edges,topTable3){
  network$nodeData = merge(topTable3,nodes,by = "symbol",all.y=TRUE)
  network$edgeData = merge(network$edgeData,edges,by.y = c("src","dest"),by.x = c("fromNode","toNode"),all.y=TRUE)
  return(network)
}