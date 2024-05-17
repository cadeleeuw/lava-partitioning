## script expects a folder containing all .breaks files, by chromosome; must contain file for chr 1-22, chr X is optional (can code as X or 23 in file name)
## files should have file format [DIR]/[PREFIX][CHR].breaks


flags = commandArgs(T)

## input arguments from command line, should be provided in specified order (only input.dir argument is mandatory)
input.dir = flags[1]  			                                             ## name of folder containing the .breaks files
file.prefix = ifelse(length(flags) > 1, flags[2], "chr")                 ## prefix of .breaks files, expects "chr[CHR].breaks" by default
min.size = ifelse(length(flags) > 2, as.numeric(flags[3]), 0)            
max.metric = ifelse(length(flags) > 3, as.numeric(flags[4]), 100)
max.blocks = ifelse(length(flags) > 4, as.numeric(flags[5]), 1e7)

out.name = input.dir
if (min.size > 0) out.name = paste0(out.name, "_size", min.size)
if (max.metric < 100) out.name = paste0(out.name, "_max", max.metric)
if (max.blocks < 1e7) out.name = paste0(out.name, "_cap", max.blocks)

source("ldblock.r")

out = NULL
breaks = list(); counts = c(); 
for (chr in 1:23) {
	curr = load.breaks(paste0(input.dir, "/", file.prefix), chromosome=chr)
	if (!is.null(curr)) {
		breaks[[chr]] = curr
		counts[chr] = nrow(filter.breaks(breaks[[chr]], min.size=min.size, max.metric=max.metric/100)) - 1
	} else if (chr != 23) stop(paste0("missing input file for chromosome ", chr))
}
chr.range = seq_along(breaks)

if (sum(counts) > max.blocks) {
	frac = counts / sum(counts) * max.blocks
	chr.max = floor(frac)
	residual = frac - chr.max
	if (sum(chr.max) < max.blocks) {
		add = residual >= sort(residual, decreasing=T)[max.blocks - sum(chr.max)]
		chr.max[add] = chr.max[add] + 1
	}
} else {
	chr.max = rep(max.blocks, length(chr.range))
}

for (chr in chr.range) {
  curr.breaks = filter.breaks(breaks[[chr]], max.blocks=chr.max[[chr]], min.size=min.size, max.metric=max.metric/100)
  curr = cbind(chr=chr, make.blocks(curr.breaks), stringsAsFactors=F)
  if (is.null(out)) {
    out = curr
  } else {
    out = rbind(out, curr)
  }
}

write.table(out, file=paste0(out.name, ".blocks"), row.names=F, sep="\t", quote=F)

dim(out)
summary(out$size.filt)
summary(out$size.all)
table(out$chr)
