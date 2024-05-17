#NB: didn't put much error checking in yet, so need to put in sensible filtering values where applicable

load.breaks = function(prefix, chromosome=NULL) {
	file.name = paste0(prefix, chromosome, ".breaks")
	if (!file.exists(file.name)) {
		if (!is.null(chromosome) && chromosome == 23) return(load.breaks(prefix, "X"))
		else return(NULL)
	} else return(read.table(file.name, header=T))
}

#bounds are inclusive, ie. should include SNPs with BP equal to the boundary points
#size currently reflects just the number of SNPs per block after the MAF filtering applied
make.blocks = function(breaks) {
  N = dim(breaks)[1]
  blocks = data.frame(
    start=breaks$POSITION[-N],
    stop=breaks$POSITION[-1]-1,
    size.filt=breaks$INDEX_FILT[-1] - breaks$INDEX_FILT[-N],
    size.all=breaks$INDEX_ALL[-1] - breaks$INDEX_ALL[-N],
    stringsAsFactors=F
  )
  return(blocks)
}


#only meant to be applied to a complete, loaded breaks data.frame once
#min.size is for filtering on the number of SNPs per block after the MAF filtering, min.size.all is on SNPs per block total
filter.breaks = function(breaks, max.blocks=NULL, max.metric=1, min.size=0, min.size.all=0) {
  tree = build.tree(breaks)
  discard = breaks$RANK[breaks$METRIC_MIN > max.metric]

  checked = rep(F, length(tree))  
  for (n in 1:length(checked)) {
    if (!checked[n]) {
      curr = tree[[n]]
      if (!is.null(curr$break.rank)) {
        children = tree[curr$children]
        if (!check.nodesize(children[[1]], min.size.all, min.size) || !check.nodesize(children[[2]], min.size.all, min.size)) discard = c(discard, curr$break.rank)
        if (curr$break.rank %in% discard) {
          cascade = curr$children
          while (length(cascade) > 0) {
            checked[cascade[1]] = T
            curr = tree[[cascade[1]]]; cascade = cascade[-1]
            if (!is.null(curr$break.rank)) {
              discard = c(discard, curr$break.rank)
              cascade = c(cascade, curr$children)
            }
          }
        }
      }
    }
  }
  
  breaks = breaks[!(breaks$RANK %in% discard),]
  if (!is.null(max.blocks) && dim(breaks)[1] - 1 > max.blocks) {
    ranks = sort(breaks$RANK[breaks$RANK>0])
    breaks = breaks[breaks$RANK < ranks[max.blocks],]    
  }

  return(breaks)
}






####################
# helper functions #
####################

#create binary tree of blocks successively split into smaller blocks
#storing BP range as inclusive for end point, not inclusive for index and filtered index
build.tree = function(breaks) {
  bounds = breaks[breaks$RANK==0,] #should be two
  breaks = breaks[breaks$RANK > 0,]
  breaks = breaks[order(breaks$RANK),]

  tree = list(make.node(1, bounds$POSITION, bounds$INDEX_ALL, bounds$INDEX_FILT, TRUE)); 
  for (i in 1:dim(breaks)[1]) {
    curr.pos = breaks$POSITION[i] 
    split.id = find.node(tree, 1, curr.pos)
    if (is.null(split.id)) stop("NULL id")
    parent = tree[[split.id]]
    
    child.ids = length(tree) + 1:2
    curr.index = breaks$INDEX_ALL[i]; curr.filt = breaks$INDEX_FILT[i]
    tree[[child.ids[1]]] = make.node(child.ids[1], c(parent$position[1], curr.pos), c(parent$index[1], curr.index), c(parent$index.filt[1], curr.filt), TRUE)
    tree[[child.ids[2]]] = make.node(child.ids[1], c(curr.pos, parent$position[2]), c(curr.index, parent$index[2]), c(curr.filt, parent$index.filt[2]), TRUE)
    
    tree[[split.id]]$children = child.ids    
    tree[[split.id]]$break.rank = breaks$RANK[i]
  }
  return(tree)
}

make.node = function(id, position, index, index.filt, shift.end) {
  if (shift.end) {position[2] = position[2] - 1}
  list(id=id, position=position, index=index, index.filt=index.filt, break.rank=NULL, children=NULL)
}
check.nodesize = function(node, index.size, filt.size) {return( (node$index[2] - node$index[1]) >= index.size && (node$index.filt[2] - node$index.filt[1]) >= filt.size )}
in.node = function(node, position) {return(position >= node$position[1] && position <= node$position[2])}
find.node = function(tree, node.id, position) {
  children = tree[[node.id]]$children
  if (!is.null(children)) {
    for (i in 1:length(children)) {if (in.node(tree[[children[i]]], position)) return(find.node(tree, children[i], position))}  
    return(NULL)
  } else {
    return(node.id)
  }
}





############################

# I used this for making the plot I sent before, you could use this to make similar plots if you want to have a look at the splits for a particular chromosome

#orig = read.table("~/tmp/partition/eur_chr22.bed", header=T)
#N.orig = dim(orig)[1]
#breaks.orig = c(orig$start[1], round((orig$stop[-N.orig] + orig$start[-1])/2), orig$stop[N.orig])
#
#
#ld = load.breaks("chr22_max1000")
#ld.sub47 = filter.breaks(ld, total=47)
#ld.filt0 = load.breaks("chr22_max1000_filter100")
#
#x.lim = range(c(breaks.orig, ld$POSITION))
#
#png("vsOrig.png", 1500, 500, type="cairo")
#plot(breaks.orig, rep(1,length(breaks.orig)), pch=19, type="b", lwd=2, col="red", xlim=x.lim, ylim=c(1,4))
#points(ld.sub47$POSITION, rep(2,dim(ld.sub47)[1]), pch=19, type="b", lwd=2, col="blue")
#points(ld$POSITION, rep(3,dim(ld)[1]), pch=19, type="b", lwd=2)
#points(ld.filt0$POSITION, rep(4,dim(ld.filt0)[1]), pch=19, type="b", lwd=2, col="green")
#
#abline(v=breaks.orig, col=2, lty=3)
#add = ld.filt0$POSITION[!ld.filt0$POSITION %in% ld$POSITION][c(2,3,5)]
#points(add, rep(4, length(add)), pch=19, cex=1.5)
#
#dev.off()
#
#

