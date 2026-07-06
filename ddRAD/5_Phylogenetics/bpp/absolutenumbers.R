#convert stacks phylip to a complete concatenated file

library(phylotools)

setwd("E:/Pittas/ddRAD/5_Phylogenetics/bpp/HP_BWP_MP_R80_MAC3_noadmix_clean")

#before the next step, remember to change the sample names (Stacks automatically truncates to 10 characters) and replace tabs with space

bpp_R80 <- read.phylip("populations.all.phylip")

dat2phylip(bpp_R80,"HP_BWP_MP_R80_MAC3_noadmix_clean_bpp.phy")

#now use Jonathan Chang's python script to split concatenated phylip file into separate partitions.

#mkdir phylip

#cd ./phylip

#python ../unconcatenate_phylip.py ../HP_BWP_MP_R80_MAC3_noadmix_clean_bpp.phy ../populations.all.partitions.phylip

#Sample random phylip loci

#there are 5084 loci in the final shortlist

#ls | shuf -n 1000 | cat > ../HP_BWP_MP_R80_MAC3_bpp_subset12000_4.txt

#{ xargs cat < ../HP_BWP_MP_R80_MAC3_bpp_subset1000_2.txt  ; } > ../HP_BWP_MP_R80_MAC3_bpp_subset1000_2.phy

setwd("E:/Pittas/ddRAD/4_Stacksout/HP_BWP_MP_R0_noMAF/")

bpp_HP_BWP_MP_R80 <- read.phylip("populations.all.phylip")

dat2phylip(bpp_HP_BWP_MP_R80,"HP_BWP_MP_R80_filtered_bpp.phy")

test <- read.phylip("populations.all.phylip")

dat2phylip(test,outfile="test.phy")

mcmc <- read.csv("E:/Pittas/ddRAD/5_Phylogenetics/bpp/HP_BWP_MP_R80_MAC3_noadmix_clean_subset1000_1/HP_BWP_MP_R80_filtered_subset1000_1_mcmc.txt",sep="\t")

pitta.time <- bppr::msc2time.r(mcmc, u.mean = 4.6e-9, u.sd = 1.3e-9, g.mean = 3.142589033, g.sd = 0.47946489)

pitta.time

apply(pitta.time,2,mean)

coda::HPDinterval(coda::as.mcmc(pitta.time))

mcmc_2 <- read.csv("E:/Pittas/ddRAD/5_Phylogenetics/bpp/HP_BWP_MP_R80_MAC3_noadmix_clean_subset1000_2/HP_BWP_MP_R80_filtered_subset1000_2_mcmc.txt",sep="\t")

pitta.time_2 <- bppr::msc2time.r(mcmc_2, u.mean = 4.6e-9, u.sd = 1.3e-9, g.mean = 3.142589033, g.sd = 0.47946489)

apply(pitta.time_2,2,mean)

coda::HPDinterval(coda::as.mcmc(pitta.time_2))

mcmc_test <- read.csv("mcmc.txt", sep = "\t")

pitta.test.time <- bppr::msc2time.r(mcmc_test,u.mean = 4.6e-9, u.sd = 1.3e-9, g.mean = 3.142589033, g.sd = 0.47946489)

apply(pitta.test.time,2,mean)

coda::HPDinterval((coda::as.mcmc(pitta.test.time)))

mcmc_4 <- read.csv("./bpp/subset_500_4_MSCM/mcmc.txt",sep="\t")

pitta.time_4 <- bppr::msc2time.r(mcmc_4, u.mean = 4.6e-9, u.sd = 1.3e-9, g.mean = 3.142589033, g.sd = 0.47946489)

apply(pitta.time_4,2,mean)

coda::HPDinterval(coda::as.mcmc(pitta.time_4))

mcmc_5 <- read.csv("./bpp/subset_500_5_MSCM/mcmc.txt",sep="\t")

pitta.time_5 <- bppr::msc2time.r(mcmc_5, u.mean = 4.6e-9, u.sd = 1.3e-9, g.mean = 3.142589033, g.sd = 0.47946489)

apply(pitta.time_5,2,mean)

coda::HPDinterval(coda::as.mcmc(pitta.time_5))

# Ne_1BR       Ne_2BM       Ne_3MP       Ne_4HP Ne_5BRBMMPHP   Ne_6BRBMMP     Ne_7BRBM  t_5BRBMMPHP    t_6BRBMMP      t_7BRBM     M_BR..BM 
# 6.943544e+03 1.425323e+05 2.786900e+04 1.244822e+05 4.542221e+05 2.771442e+05 1.195096e+04 3.308486e+06 5.135590e+05 4.953130e+05 1.687854e-02 
# M_BM..BR            u            g         rate 
# 1.831001e-01 4.600082e-09 3.141878e+00 1.499047e-09 
# > coda::HPDinterval(coda::as.mcmc(pitta.time))
# lower        upper
# Ne_1BR       3.243154e+03 1.159814e+04
# Ne_2BM       7.069844e+04 2.305597e+05
# Ne_3MP       1.379253e+04 4.530134e+04
# Ne_4HP       6.154156e+04 2.012851e+05
# Ne_5BRBMMPHP 2.254373e+05 7.364180e+05
# Ne_6BRBMMP   1.375186e+05 4.491467e+05
# Ne_7BRBM     1.458594e+03 3.167982e+04
# t_5BRBMMPHP  1.451866e+06 5.603415e+06
# t_6BRBMMP    2.217055e+05 8.738937e+05
# t_7BRBM      2.113837e+05 8.440441e+05
# M_BR..BM     5.130000e-04 3.872400e-02
# M_BM..BR     1.536870e-01 2.136130e-01
# u            2.201018e-09 7.173327e-09
# g            2.238378e+00 4.102366e+00
# rate         6.408766e-10 2.475030e-09
# attr(,"Probability")
# [1] 0.9499997
# > mcmc_2 <- read.csv("E:/Pittas/ddRAD/5_Phylogenetics/bpp/HP_BWP_MP_R80_MAC3_noadmix_clean_subset1000_2/HP_BWP_MP_R80_filtered_subset1000_2_mcmc.txt",sep="\t")
# > pitta.time_2 <- bppr::msc2time.r(mcmc_2, u.mean = 4.6e-9, u.sd = 1.3e-9, g.mean = 3.142589033, g.sd = 0.47946489)
# > apply(pitta.time_2,2,mean)
# Ne_1BR       Ne_2BM       Ne_3MP       Ne_4HP Ne_5BRBMMPHP   Ne_6BRBMMP     Ne_7BRBM  t_5BRBMMPHP    t_6BRBMMP      t_7BRBM     M_BR..BM 
# 4.936997e+03 1.514767e+05 3.062922e+04 1.179791e+05 4.536613e+05 3.068738e+05 1.059510e+04 3.244111e+06 5.324902e+05 5.113155e+05 1.182861e-02 
# M_BM..BR            u            g         rate 
# 1.804488e-01 4.601978e-09 3.143854e+00 1.498790e-09 
# > coda::HPDinterval(coda::as.mcmc(pitta.time_2))
# lower        upper
# Ne_1BR       2.238927e+03 8.277735e+03
# Ne_2BM       7.495064e+04 2.442902e+05
# Ne_3MP       1.516690e+04 4.970507e+04
# Ne_4HP       5.846911e+04 1.906042e+05
# Ne_5BRBMMPHP 2.240875e+05 7.333424e+05
# Ne_6BRBMMP   1.527258e+05 4.977874e+05
# Ne_7BRBM     1.714838e+03 2.484148e+04
# t_5BRBMMPHP  1.420571e+06 5.483606e+06
# t_6BRBMMP    2.326484e+05 9.078100e+05
# t_7BRBM      2.217293e+05 8.723042e+05
# M_BR..BM     2.080000e-04 2.815400e-02
# M_BM..BR     1.520260e-01 2.116110e-01
# u            2.227351e-09 7.198472e-09
# g            2.219433e+00 4.085050e+00
# rate         6.314233e-10 2.467525e-09
# attr(,"Probability")
# [1] 0.9499997