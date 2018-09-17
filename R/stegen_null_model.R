## bbs.sp.site is a matrix with species as columns and local samples within a region as rows. The matrix contains species abundances.

## James adapting null model for general use with a community matrix

null.stegen<-function(mat){

null.alpha.comp = numeric()

for (randomize in 1:1000) {  
	null.dist = mat;
	for (species in 1:ncol(null.dist)) {
		tot.abund = sum(null.dist[,species]);
		null.dist[,species] = 0;
		for (individual in 1:tot.abund) {
			sampled.site = sample(c(1:nrow(mat)),1);
			null.dist[sampled.site,species] = null.dist[sampled.site,species] + 1;
		};
 
	};
	null.dist = ceiling(null.dist/max(null.dist));
	null.alpha =  sum(null.dist)/nrow(mat); #null.alpha;
	null.alpha.comp = c(null.alpha.comp,null.alpha);
} ## end randomize loop

return(list(null.dist, null.alpha.comp))

}

bbs.sp.site = ceiling(bbs.sp.site/max(bbs.sp.site)); 
mean.alpha = sum(bbs.sp.site)/nrow(bbs.sp.site); #mean.alpha;
gamma = ncol(bbs.sp.site); #gamma;
z.score = ((1-mean.alpha/gamma) - mean(1 - null.alpha.comp/gamma)) / sd(1 - null.alpha.comp/gamma); z.score; ## null departure turnover
obs.turnover = (1-mean.alpha/gamma); obs.turnover; ## raw turnover