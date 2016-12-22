#match timescale names to stages
stage_match <- read.csv("https://raw.githubusercontent.com/mclapham/shark_extinction/master/stage_match.csv")
time_int <- read.csv("https://paleobiodb.org/data1.2/intervals/list.txt?scale=1&scale_level=5")

#read occurrences
taxa <- c("Brachiopoda,Bivalvia,Gastropoda,Porifera,Cnidaria,Echinodermata,Bryozoa,Trilobita,Malacostraca,Ostracoda")

chond <- read.csv("https://paleobiodb.org/data1.2/occs/list.txt?base_name=Chondrichthyes^Xenacanthimorpha&envtype=marine&interval=Asselian-Bathonian&taxon_reso=genus&pres=regular&show=class")
actino <- read.csv("https://paleobiodb.org/data1.2/occs/list.txt?base_name=Actinopterygii,Actinistia&envtype=marine&interval=Asselian-Bathonian&taxon_reso=genus&pres=regular&show=class")
occ_url <- RCurl::getURL(paste("https://paleobiodb.org/data1.2/occs/list.txt?base_name=",taxa,"&envtype=marine&interval=Asselian-Bathonian&taxon_reso=genus&pres=regular&show=class",sep=""), ssl.verifypeer=F, timeout=300)
invert <- read.csv(text=occ_url)

#function to select occurrences resolved to a single stage
time.cleaning <- function(var_name) {
  #selects time intervals where early interval is in list of stages
  var_sub <- subset(var_name, var_name$early_interval %in% stage_match$int_name)
  
  #replaces early and late interval names with international stage
  var_sub$early_interval <- stage_match$stage_name[match(var_sub$early_interval, stage_match$int_name)]
  var_sub$late_interval <- stage_match$stage_name[match(var_sub$late_interval, stage_match$int_name)]
  
  #finds occurrences in a single interval or in two intervals
  var_single <- subset(var_sub, is.na(var_sub$late_interval) == T)
  var_mult <- subset(var_sub, is.na(var_sub$late_interval) == F)
  
  #finds occurrences where early interval international stage is same as late interval international stage
  var_mult_ident <- subset(var_mult, as.character(var_mult$early_interval) == as.character(var_mult$late_interval))
  
  #combines good data
  var_final <- rbind(var_single, var_mult_ident)
  
  droplevels(var_final)
  
}

#applies time interval cleaning to occurrences
chond_final <- time.cleaning(chond)
actino_final <- time.cleaning(actino)
invert_final <- time.cleaning(invert)

#selects stages that are represented in the occurrence data
time_int_final <- subset(time_int, time_int$interval_name %in% levels(chond_final$early_interval))

#converts occurrence list into taxon/time occurrence matrix
taxon.matrix <- function(occs) {
  #counts number of genus occurrences per time interval
  genus_matrix_raw <- sapply(split(occs, occs$early_interval), function(x) table(x$genus))
  
  #orders genus matrix so that time intervals are in chronologic order
  genus_matrix <- genus_matrix_raw[,match(time_int_final$interval_name, colnames(genus_matrix_raw))]
  
  #reverses order so it's from oldest to youngest
  genus_matrix[,rev(seq(ncol(genus_matrix)))]
}

chond_matrix <- taxon.matrix(chond_final)
chond_matrix <- chond_matrix[which(nchar(rownames(chond_matrix)) > 0),]

actino_matrix <- taxon.matrix(actino_final)
actino_matrix <- actino_matrix[which(nchar(rownames(actino_matrix)) > 0),]

invert_matrix <- taxon.matrix(invert_final)
invert_matrix <- invert_matrix[which(nchar(rownames(invert_matrix)) > 0),]

#Small sample-corrected Akaike Information Criterion calculation
aicc<-function(ml,K,n)
{	
  -2*ml+2*K+(2*K*(K+1))/(n-K-1)

}

#Akaike weights calculation
aicw<-function(x)
{
  exp(-0.5*x)/sum(exp(-0.5*x))

}

#Log likelihoods with boundary-crosser measures
#boundary-crosser extinction q=-log(Nbt/Nb)
aic.comp <- function(var1, var2) {
  #all (sharks and inverts together)
  ll_all <- ((var1$Nbt + var2$Nbt) * log((var1$Nbt + var2$Nbt) / (var1$Nb + var2$Nb))) + ((var1$NbL + var2$NbL) * log((var1$NbL + var2$NbL) / (var1$Nb + var2$Nb)))
  #shark likelihood calculation
  ll_var1 <- (var1$Nbt * log(var1$Nbt / var1$Nb)) + (var1$NbL * log(var1$NbL / var1$Nb))
  #invert likelihood calculation
  ll_var2 <- (var2$Nbt * log(var2$Nbt / var2$Nb)) + (var2$NbL * log(var2$NbL / var2$Nb))
  
  #AIC analysis (boundary-crosser)
  BC_onerate_aic <- aicc(ll_all, 1, var1$Nb + var2$Nb)
  BC_tworate_aic <- aicc(ll_var1 + ll_var2, 2, var1$Nb + var2$Nb)
  
  #Akaike weights
  apply(data.frame(BC_onerate_aic, BC_tworate_aic), 1, aicw)
  
}


#Subsampled boundary-crosser extinction
subsamp.bc <- function(taxon_table, size) {
  taxon_subsamp <- by(taxon_table, taxon_table$early_interval, function(x) x[sample(nrow(x)),][1:size,])
  taxon_matrix_subsamp <- sapply(taxon_subsamp, function(x) table(x$genus))
  
  taxon_matrix_subsamp <- taxon_matrix_subsamp[,match(time_int_final$interval_name, colnames(taxon_matrix_subsamp))]
  
  taxon_matrix_subsamp <- taxon_matrix_subsamp[,rev(seq(ncol(taxon_matrix_subsamp)))]
  
  bc.ext(taxon_matrix_subsamp)
  
}

subsamp.summ <- function(taxon_table, n) {
  stest <- replicate(100, subsamp.bc(taxon_table, n))
  
  Nb <- apply(matrix(unlist(stest[1,]), ncol=100), 1, mean)
  Nbt <- apply(matrix(unlist(stest[2,]), ncol=100), 1, mean)
  NbL <- apply(matrix(unlist(stest[3,]), ncol=100), 1, mean)
  
  raw_ext <- data.frame(Nb, Nbt, NbL)
  
  occ_ct <- table(taxon_table$early_interval)
  small_bin <- which(rev(occ_ct[match(time_int_final$interval_name, names(occ_ct))]) < n)
  small_bin <- subset(small_bin, small_bin < 22)
  
  raw_ext[small_bin,] <- NA
  
  raw_ext
}

chond_subsamp_ext <- subsamp.summ(chond_final, 20)
actino_subsamp_ext <- subsamp.summ(actino_final, 20)
invert_subsamp_ext <- subsamp.summ(invert_final, 1000)

#likelihood comparison
chond_invert_sub <- aic.comp(chond_subsamp_ext, invert_subsamp_ext)
actino_invert_sub <- aic.comp(actino_subsamp_ext, invert_subsamp_ext)

ages <- time_int_final$min_ma[rev(seq(nrow(time_int_final)))]

#BC extinction rates for each group
ext_sum_sub <- data.frame(age=ages[seq(length(ages)-2)], chond_ext= -log(chond_subsamp_ext$Nbt/chond_subsamp_ext$Nb), actino_ext= -log(actino_subsamp_ext$Nbt/actino_subsamp_ext$Nb), invert_ext= -log(invert_subsamp_ext$Nbt/invert_subsamp_ext$Nb))


pdf("extinction_comp_sub.pdf", width=8)
#Plot preparation
colorpal <- RColorBrewer::brewer.pal(9, "PRGn")
bin_midpt <- colMeans(rbind(time_int_final$max_ma, time_int_final$min_ma))
extval <- c(ext_sum_sub$chond_ext, ext_sum_sub$actino_ext, ext_sum_sub$invert_ext)
maxval <- max(extval[is.finite(extval)])
par(mgp=c(2,0.75,0))
par(mfrow=c(2,1))
par(mar=c(3, 3, 2, 1))

#Chondrichthyes comparison plot
chond_invert_col <- colorRampPalette(colorpal)(100)[cut(chond_invert_sub[1,], seq(0,1,by=0.01))]

plot(ext_sum_sub$age, ext_sum_sub$invert_ext, type="n", xlab="Age (Ma)", ylab="Extinction rate", ylim=c(-0.07, maxval), xlim=rev(range(ext_sum_sub$age)))

abline(v=ext_sum_sub$age[c(9, 16)], lty=2)

points(ext_sum_sub$age, ext_sum_sub$invert_ext, type="o", pch=16, cex=0.75, col="gray")

lines(ext_sum_sub$age, ext_sum_sub$chond_ext, lwd=2)
points(ext_sum_sub$age, ext_sum_sub$chond_ext, pch=21, cex=1.5, bg=chond_invert_col)

gradient_image <- as.raster(matrix(colorRampPalette(colorpal)(100), ncol=1))
rasterImage(gradient_image, 290, maxval-0.35, 285, maxval-0.15)
text(285, c(maxval-0.35, maxval-0.25, maxval-0.15), labels=c(0, 0.5, 1), pos=4, cex=0.75)
text(287.5, maxval-0.05, "Two-rate model support", cex=0.75)

rect(time_int_final$max_ma, -0.09, time_int_final$min_ma, -0.03, col=paste(time_int_final$color))
text(bin_midpt, -0.06, strtrim(time_int_final$interval_name,2),cex=0.6)

legend("topright", c("Invertebrates", "Chondrichthyes"), col=c("gray", "black"), lwd=c(1,2), bty="n", cex=0.75)

mtext("A", adj=0)

#Actinopterygian comparison plot
actino_invert_col <- colorRampPalette(colorpal)(100)[cut(actino_invert_sub[1,], seq(0,1,by=0.01))]

plot(ext_sum_sub$age, ext_sum_sub$invert_ext, type="n", xlab="Age (Ma)", ylab="Extinction rate", ylim=c(-0.07, maxval), xlim=rev(range(ext_sum_sub$age)))

abline(v=ext_sum_sub$age[c(9, 16)], lty=2)

points(ext_sum_sub$age, ext_sum_sub$invert_ext, type="o", pch=16, cex=0.75, col="gray")

lines(ext_sum_sub$age, ext_sum_sub$actino_ext, lwd=2)
points(ext_sum_sub$age, ext_sum_sub$actino_ext, pch=21, cex=1.5, bg=actino_invert_col)

gradient_image <- as.raster(matrix(colorRampPalette(colorpal)(100), ncol=1))
rasterImage(gradient_image, 290, maxval-0.35, 285, maxval-0.15)
text(285, c(maxval-0.35, maxval-0.25, maxval-0.15), labels=c(0, 0.5, 1), pos=4, cex=0.75)
text(287.5, maxval-0.05, "Two-rate model support", cex=0.75)

rect(time_int_final$max_ma, -0.09, time_int_final$min_ma, -0.03, col=paste(time_int_final$color))
text(bin_midpt, -0.06, strtrim(time_int_final$interval_name,2),cex=0.6)


legend("topright", c("Invertebrates", "Bony fishes"), col=c("gray", "black"), lwd=c(1,2), bty="n", cex=0.75)

mtext("B", adj=0)
dev.off()
