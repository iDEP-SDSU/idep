	  if( sum(!is.na( x2$genomeSpan) ) >= minGenes && length(unique(x2$genomeSpan) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(log(genomeSpan)~Set, data=x2[which(!is.na(x2$genomeSpan)),] )$p.value
	   sig = paste("Genome span (coding)  T-test P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ***" )
		p3 <- ggplot(x2, aes(genomeSpan, fill= Set, colour = Set) )+
			geom_density(alpha = 0.1) + 
			scale_x_log10() +
			ggtitle(sig) +
			theme(plot.title = element_text(hjust = 0.5)) +
			labs(x = "Genome span (bp)")		   
	   }		  


	  if( sum(!is.na( x2$FiveUTR) ) >= minGenes && length(unique(x2$FiveUTR) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(log(FiveUTR)~Set, data=x2[which(!is.na(x2$FiveUTR)),] )$p.value
	   sig = paste("5' UTR length (coding genes only)\n T-test P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ***" )
		p4 <- ggplot(x2, aes(FiveUTR, fill= Set, colour = Set) )+
			geom_density(alpha = 0.1) + 
			scale_x_log10() +
			ggtitle(sig) +
			theme(plot.title = element_text(hjust = 0.5)) +
			labs(x = "5' UTR length(bp)")
	}	   
	
	if( sum(!is.na( x2$ThreeUTR) ) >= minGenes && length(unique(x2$ThreeUTR) ) > 2 && length(which(x2$Set == "List") ) > minGenes ) {
       Pval = t.test(log(ThreeUTR)~Set, data=x2[which(!is.na(x2$ThreeUTR)),] )$p.value
	   sig = paste("3' UTR length (coding)  T-test P=",formatC(Pval, digits=2, format="G"),sep="")
	   if( Pval <PvalGeneInfo)  sig = paste(sig," ****************" )
		p5 <- ggplot(x2, aes(ThreeUTR, fill= Set, colour = Set) )+
			geom_density(alpha = 0.1) + 
			scale_x_log10() +
			ggtitle(sig) +
			theme(plot.title = element_text(hjust = 0.5)) +
			labs(x = "3' UTR length(bp)")
	  }	
		  
	  if( sum(!is.na( x2$percentage_gc_content) ) >= minGenes && 
		  length(unique(x2$percentage_gc_content) ) > 2 && 
		  length(which(x2$Set == "List") ) > minGenes ) {
		   Pval = t.test(percentage_gc_content~Set, data=x2[which(!is.na(x2$percentage_gc_content)),] )$p.value
		   sig = paste("%GC content (coding genes only)\n T-test P=",formatC(Pval, digits=2, format="G"),sep="")
		   if( Pval <PvalGeneInfo)  sig = paste(sig," ***" )
			p6 <- ggplot(x2, aes(percentage_gc_content, fill= Set, colour = Set) )+
				geom_density(alpha = 0.1) + 
				scale_x_log10() +
				ggtitle(sig) +
				theme(plot.title = element_text(hjust = 0.5)) +
				labs(x = "GC content")		   
	   }		  
		  