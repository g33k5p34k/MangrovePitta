#load packages into environment
library(adegenet) # for dapc
library(vcfR) # for reading in genetic data
library(tidyverse) # for manipulating and plotting data
library(LEA) # For sNMF
library(factoextra)
library(rnaturalearth) #for mapping
library(SNPfiltR)
library(cowplot)

#SNPfiltR analysis for phylogenetic dataset

setwd("E:/Pittas/ddRAD/4_Stacksout/HP_BWP_MP_R0_noMAF/")

vcfR <- read.vcfR("populations.snps.vcf")
popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 3,11))

vcfR

#visualize distributions
hard_filter(vcfR=vcfR)
#> no depth cutoff provided, exploratory visualization will be generated.
#>
#hard filter to minimum depth of 5, and minimum genotype quality of 30
vcfR<-hard_filter(vcfR=vcfR, depth = 5, gq = 30)

vcfR<-filter_allele_balance(vcfR)

max_depth(vcfR)

vcfR<-max_depth(vcfR, maxdepth = 150)

missing_by_sample(vcfR=vcfR, popmap = popmap)

#verify that missing data is not driving clustering patterns among the retained samples
miss<-assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds = .8, clustering = FALSE)
#> cutoff is specified, filtered vcfR object will be returned

#if there are still problematic samples, drop them using the following syntax
vcfR_noadmix <- vcfR[,colnames(vcfR@gt) != "MP-MSB-50624-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-MSB-50633-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-MSB-50629-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-N1436-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-N1439-SG"]

#if there are still problematic samples, drop them using the following syntax
#vcfR <- vcfR[,colnames(vcfR@gt) != "MP-MSB-50624-SG"]
#popmap<-popmap[popmap$id %in% colnames(vcfR@gt),]

miss<-assess_missing_data_tsne(vcfR, popmap, clustering = FALSE)

vcfR_noadmix<-missing_by_snp(vcfR_noadmix, cutoff = .8)

#test effect of minor allele count filtering
min_mac(vcfR_noadmix, min.mac = 3)

vcfR.mac<-min_mac(vcfR_noadmix, min.mac = 3)

vcfR::write.vcf(vcfR.mac, file = "E:/Pittas/ddRAD/4_Stacksout/HP_BWP_MP_R0_noMAF/HP_BWP_MP_R80_MACfiltered_noadmix.vcf.gz")

vcfR::write.vcf(vcfR, file = "E:/Pittas/ddRAD/4_Stacksout/HP_BWP_MP_R0_noMAF/BWP_MP_R80.vcf.gz")

vcfR.thin<-distance_thin(vcfR.mac, min.distance = 500)

vcfR::write.vcf(vcfR.thin, file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R80_LDfilt.vcf.gz")

######################################################################

#SNPfiltR analysis for fastsimcoal2

setwd("E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R0_noMAF/")

vcfR <- read.vcfR("populations.snps.vcf")
popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 3,11))

vcfR

#visualize distributions
hard_filter(vcfR=vcfR)
#> no depth cutoff provided, exploratory visualization will be generated.
#>
#hard filter to minimum depth of 5, and minimum genotype quality of 30
vcfR<-hard_filter(vcfR=vcfR, depth = 5, gq = 30)

#vcfR<-filter_allele_balance(vcfR)

max_depth(vcfR)

vcfR<-max_depth(vcfR, maxdepth = 150)

missing_by_sample(vcfR=vcfR, popmap = popmap)

#verify that missing data is not driving clustering patterns among the retained samples
miss<-assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds = .8, clustering = FALSE)
#> cutoff is specified, filtered vcfR object will be returned

#if there are still problematic samples, drop them using the following syntax
vcfR_noadmix <- vcfR[,colnames(vcfR@gt) != "MP-MSB-50624-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-MSB-50633-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-MSB-50629-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-N1436-SG"]
vcfR_noadmix <- vcfR_noadmix[,colnames(vcfR_noadmix@gt) != "BWP-N1439-SG"]

#if there are still problematic samples, drop them using the following syntax
#vcfR <- vcfR[,colnames(vcfR@gt) != "MP-MSB-50624-SG"]
#popmap<-popmap[popmap$id %in% colnames(vcfR@gt),]

miss<-assess_missing_data_tsne(vcfR, popmap, clustering = FALSE)

vcfR_noadmix<-missing_by_snp(vcfR_noadmix, cutoff = .8)

#test effect of minor allele count filtering
#min_mac(vcfR_noadmix, min.mac = 3)

#vcfR.mac<-min_mac(vcfR_noadmix, min.mac = 3)

#vcfR::write.vcf(vcfR.mac, file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R0_noMAF/HP_BWP_MP_R80_MACfiltered_noadmix.vcf.gz")

vcfR::write.vcf(vcfR_noadmix, file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R0_noMAF/BWP_MP_R80_fastsimcoal.vcf.gz")

vcfR.thin<-distance_thin(vcfR_noadmix, min.distance = 500)

vcfR::write.vcf(vcfR.thin, file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R0_noMAF/BWP_MP_R80_LDfilt_fastsimcoal.vcf.gz")

######################################################################

#SNPfiltR analysis for popgen dataset

setwd("E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R0_noMAF/")

vcfR <- read.vcfR("populations.snps.vcf")
popmap<-data.frame(id=colnames(vcfR@gt)[2:length(colnames(vcfR@gt))],pop=substr(colnames(vcfR@gt)[2:length(colnames(vcfR@gt))], 3,11))

vcfR

#visualize distributions
hard_filter(vcfR=vcfR)
#> no depth cutoff provided, exploratory visualization will be generated.
#>
#hard filter to minimum depth of 5, and minimum genotype quality of 30
vcfR<-hard_filter(vcfR=vcfR, depth = 5, gq = 30)

vcfR<-filter_allele_balance(vcfR)

max_depth(vcfR)

vcfR<-max_depth(vcfR, maxdepth = 150)

missing_by_sample(vcfR=vcfR, popmap = popmap)

#verify that missing data is not driving clustering patterns among the retained samples
miss<-assess_missing_data_pca(vcfR=vcfR, popmap = popmap, thresholds = .8, clustering = FALSE)
#> cutoff is specified, filtered vcfR object will be returned


#popmap<-popmap[popmap$id %in% colnames(vcfR@gt),]

miss<-assess_missing_data_tsne(vcfR, popmap, clustering = FALSE)

vcfR<-missing_by_snp(vcfR, cutoff = .8)

#test effect of minor allele count filtering
min_mac(vcfR, min.mac = 3)

vcfR.mac<-min_mac(vcfR, min.mac = 3)

missing_by_sample(vcfR.mac)

vcfR::write.vcf(vcfR.mac, file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R80_MACfiltered_MACfiltered_popgen.vcf.gz")

vcfR::write.vcf(vcfR, file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R80.vcf.gz")

vcfR.thin<-distance_thin(vcfR.mac, min.distance = 500)

vcfR::write.vcf(vcfR.thin, file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R80_LDfilt.vcf.gz")

setwd("E:/Pittas/ddRAD/4_Stacksout/HP_BWP_MP_R0_noMAF")

#read VCF file
pitta_R80_MACfiltered_vcf <- read.vcfR("E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R80_MACfiltered_popgen.vcf")

#convert VCF file to genlight
pitta_R80_MACfiltered_genlight <- vcfR2genlight(pitta_R80_MACfiltered_vcf)

#fill in missing data
pitta_R80_MACfiltered_nomissing_genlight <- tab(pitta_R80_MACfiltered_genlight,NA.method="mean")

#read popmap file
popmap <- read.csv("E:/Pittas/ddRAD/4_Stacksout/bwp_mp_popmap.txt",sep="\t",header=FALSE)

#run PCA
pitta_R80_MACfiltered_pca <- dudi.pca(pitta_R80_MACfiltered_nomissing_genlight)

#plot PCA

fviz_pca_ind(pitta_R80_MACfiltered_pca,axes=c(1,2))
pitta_pca <- fviz_pca_ind(pitta_R80_MACfiltered_pca,geom="point", col.ind = popmap[,2], pointsize = 5, mean.point=F,alpha=0.7)+
  scale_shape_manual(values=c(19,19,19,19))+
  scale_color_manual(values=c("#57d9e9","#2c7bb6","#fdae61","#d7191c"))+
  theme(legend.text = element_text(size=10),
        legend.position=c(0.65,0.15),
        legend.title=element_blank(),
        plot.title = element_text(face="bold"),
        panel.grid = element_blank()
        )+
  xlab("PC1 (26.6%)")+
  ylab("PC2 (4.3%)")+
  ggtitle(NULL)

pitta_pca_final <- ggdraw()+
  draw_plot(pitta_pca,0,0,1,1)+
  draw_plot(fviz_eig(pitta_R80_MACfiltered_pca)+ylab("%Variance")+xlab("Eigenvector")+ggtitle(NULL)+theme(panel.grid=element_bla),0.5,0.35,0.3,0.35)

#convert vcf to geno
vcf2geno("E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R80_MACfiltered_popgen.vcf")

#run SNMF
pitta_R80_MACfiltered_snmf <- snmf(input.file = "E:/Pittas/ddRAD/4_Stacksout/BWP_MP_R80_MACfiltered_popgen.geno",
                              K = 1:10,
                              entropy = TRUE,
                              repetitions = 10,
                              project = "new")
#evaluate cross entropy
best <- which.min(cross.entropy(pitta_R80_MACfiltered_snmf,K=2))
names(popmap) <- c("ID","pop")
popmap$plot_order <- c(1:47)


#plot SNMF when K=2
q_mat <- LEA::Q(pitta_R80_MACfiltered_snmf, K = 2, run = best) 
colnames(q_mat) <- paste0("P", 1:2)
q_df <- q_mat %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = popmap$ID,
         region = popmap$pop,
         order = popmap$plot_order)

q_df_long <- q_df %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q")

q_df_pittas <- q_df_long %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_palette <- c("#2c7bb6","#d7191c")

pitta_snmf_k2 <- q_df_pittas %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette,guide="none") +
  #scale_fill_viridis_d() +
  labs(fill = "Region") +
  theme_minimal() +
  ggtitle("K=2")+
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face="bold",hjust=0.5)
  )

#plot SNMF when K=3
best3 <- which.min(cross.entropy(pitta_R80_MACfiltered_snmf,K=3))
q_mat_k3 <- LEA::Q(pitta_R80_MACfiltered_snmf, K = 3, run = best3) 
colnames(q_mat_k3) <- paste0("P", 1:3)
q_df_k3 <- q_mat_k3 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = popmap$ID,
         region = popmap$pop,
         order = popmap$plot_order)

q_df_long_k3 <- q_df_k3 %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q")

q_df_pittas_k3 <- q_df_long_k3 %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_palette <- c("#fdae61","#d7191c","#2c7bb6")

pitta_snmf_k3 <- q_df_pittas_k3 %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette,guide="none") +
  #scale_fill_viridis_d() +
  labs(fill = "Region") +
  theme_minimal() +
  ggtitle("K=3")+
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face="bold",hjust=0.5)
  )

#plot SNMF when K=4
best4 <- which.min(cross.entropy(pitta_R80_MACfiltered_snmf,K=4))
q_mat_k4 <- LEA::Q(pitta_R80_MACfiltered_snmf, K = 4, run = best4) 
colnames(q_mat_k4) <- paste0("P", 1:4)
q_df_k4 <- q_mat_k4 %>% 
  as_tibble() %>% 
  # add the pops data for plotting
  mutate(individual = popmap$ID,
         region = popmap$pop,
         order = popmap$plot_order)

q_df_long_k4 <- q_df_k4 %>% 
  # transform the data to a "long" format so proportions can be plotted
  pivot_longer(cols = starts_with("P"), names_to = "pop", values_to = "q")

q_df_pittas_k4 <- q_df_long_k4 %>% 
  # arrange the data set by the plot order indicated in Prates et al.
  arrange(order) %>% 
  # this ensures that the factor levels for the individuals follow the ordering we just did. This is necessary for plotting
  mutate(individual = forcats::fct_inorder(factor(individual)))

q_palette <- c("#2c7bb6","#d7191c","gold","#fdae61")

pitta_snmf_k4 <- q_df_pittas_k4 %>% 
  ggplot() +
  geom_col(aes(x = individual, y = q, fill = pop)) +
  scale_fill_manual(values = q_palette,guide="none") +
  #scale_fill_viridis_d() +
  labs(fill = "Region") +
  theme_minimal() +
  ggtitle("K=4")+
  # some formatting details to make it pretty
  theme(panel.spacing.x = unit(0, "lines"),
        axis.line = element_blank(),
        axis.text = element_blank(),
        strip.background = element_rect(fill = "transparent", color = "black"),
        panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.title = element_text(face="bold",hjust=0.5)
  )

# #old plotting script for SNMF
# my.colours <- c("#2c7bb6","#d7191c","#fdae61","gold")
# bp_k2 <- barchart(pitta_R80_MACfiltered_snmf,
#                K = 2,
#                run=best,
#                xlab="Individuals",
#                ylab="Ancestry Proportions",
#                main="SNMF Ancestry Matrix, K=2",
#                border=NA,
#                col = my.colours,
#                space=0,
#                sort.by.Q = FALSE)
# 
# axis(1, at = 1:length(bp$order),
#      labels = bp$order, las = 3,
#      cex.axis = 0.6)
# 
# #evaluate cross entropy
# best <- which.min(cross.entropy(pitta_R80_MACfiltered_snmf,K=4))
# 
# #plot SNMF
# my.colours <- c("lightblue","mediumorchid1","olivedrab","gold")
# bp <- barchart(pitta_R80_MACfiltered_snmf,
#                K = 4,
#                run=best,
#                xlab="Individuals",
#                ylab="Ancestry Proportions",
#                main="SNMF Ancestry Matrix, K=4",
#                border=NA,
#                col = my.colours,
#                space=0,
#                sort.by.Q = FALSE)
# 
# axis(1, at = 1:length(bp$order),
#      labels = bp$order, las = 3,
#      cex.axis = 0.6)

#gwas analysis

#10kb sliding window analysis in vcftools
#vcftools --vcf populations.snps.vcf --weir-fst-pop bwp_mig_population.txt --weir-fst-pop bwp_res_population.txt --fst-window-size 10000 --out BWP_res_vs_mig_10kb

# setwd("E:/Pittas/ddRAD/7_GWAS")
# 
# bwp_mp_R80_MAC3_fst <- read_tsv("E:/Pittas/ddRAD/7_GWAS/BWPvMP/bwp_vs_mp_100kb.windowed.weir.fst")
# 
# #remove negative values
# bwp_mp_R80_MAC3_fst$WEIGHTED_FST[which(bwp_mp_R80_MAC3_fst$WEIGHTED_FST<0)]=0
# bwp_mp_R80_MAC3_fst$MEAN_FST[which(bwp_mp_R80_MAC3_fst$MEAN_FST<0)]=0
# 
# ggplot(bwp_mp_R80_MAC3_fst, aes(CHROM,WEIGHTED_FST))+geom_point()
# 
# #check FST quantiles
# quantile(bwp_mp_R80_MAC3_fst$WEIGHTED_FST, c(0.975,0.995), na.rm=T)
# 
# threshold <- quantile(bwp_mp_R80_MAC3_fst$WEIGHTED_FST, 0.995, na.rm=T)
# 
# bwp_mp_R80_MAC3_fst <- bwp_mp_R80_MAC3_fst %>% mutate(outlier = ifelse(WEIGHTED_FST > threshold, "outlier", "background"))
# 
# bwp_mp_R80_MAC3_fst %>% group_by(outlier) %>% tally()
# 
# #replot with coloured outlier loci
# ggplot(bwp_mp_R80_MAC3_fst, aes(CHROM,log(WEIGHTED_FST), colour = outlier)) + geom_point()

# #100kb sliding window
# 
# bwp_R80_MAC3_100kb_fst <- read_tsv("./BWPvMP/bwp_vs_mp_100kb.txt.windowed.weir.fst")
# 
# ggplot(bwp_R80_MAC3_100kb_fst, aes(CHROM,MEAN_FST))+geom_point()

##plotting genome-wide pixy FST results

setwd("E:/Pittas/ddRAD/7_GWAS/")

bwp_res_mig_pixy_10kb <- read.csv("ResvMig_pixy_10kb.csv", sep =",", header=TRUE)

manhattan(bwp_res_mig_pixy_10kb,chr="CHR",bp="window_pos_1",p="avg_wc_fst",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("blue4","orange3"),snp="CHR",main="Blue-winged Pitta migrants vs residents (10kb windows)",ylim =c(0,1))

bwp_res_mig_pixy_100kb <- read.csv("ResvMig_pixy_100kb.csv", sep =",", header=TRUE)

manhattan(bwp_res_mig_pixy_100kb,chr="CHR",bp="window_pos_1",p="avg_wc_fst",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("blue4","orange3"),snp="CHR",ylim=c(0,1))

bwp_mp_pixy_10kb <- read.csv("bwp_mp_pixy_10kb.csv", sep =",", header=TRUE)

manhattan(bwp_mp_pixy_10kb,chr="CHR",bp="window_pos_1",p="avg_wc_fst",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("blue4","orange3"),snp="CHR",main="Blue-winged Pitta vs Mangrove Pitta (10kb windows)",ylim =c(0,1))

bwp_mp_pixy_100kb <- read.csv("bwp_mp_pixy_100kb.csv", sep =",", header=TRUE)

manhattan(bwp_mp_pixy_100kb,chr="CHR",bp="window_pos_1",p="avg_wc_fst",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("blue4","orange3"),snp="CHR",ylim =c(0,1))

bwp_mp_pixy_500kb <- read.csv("bwp_mp_pixy_500kb.csv", sep =",", header=TRUE)

manhattan(bwp_mp_pixy_500kb,chr="CHR",bp="window_pos_1",p="avg_wc_fst",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("blue4","orange3"),snp="CHR",main="Blue-winged Pitta vs Mangrove Pitta (500kb windows)",ylim =c(0,1))

## ok time to make the damned figure

snmf_plots <- ggarrange(pitta_snmf_k2,pitta_snmf_k3,pitta_snmf_k4,ncol=1)
pca_snmf <- ggarrange(pitta_pca_final,snmf_plots,nrow=1)
manhattan_plots <- ggarrange(manhattan(bwp_res_mig_pixy_100kb,chr="CHR",bp="window_pos_1",p="avg_wc_fst",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("blue4","orange3"),snp="CHR",main="Blue-winged Pitta migrants vs residents (100kb windows)",ylim=c(0,1)),
                             manhattan(bwp_mp_pixy_100kb,chr="CHR",bp="window_pos_1",p="avg_wc_fst",logp=FALSE,ylab="Weir and Cockerham Fst",col=c("blue4","orange3"),snp="CHR",main="Blue-winged Pitta vs Mangrove Pitta (100kb windows)",ylim =c(0,1)),nrow=2)
