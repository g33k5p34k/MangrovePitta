    seed = -1 
    seqfile = HP_BWP_MP_R80_MAC3__noadmix_clean_bpp_subset1000_2.phy 
    Imapfile = HP_MP_BWP_imap.txt 
    outfile = HP_BWP_MP_R80_filtered_subset1000_2_out.txt 
    mcmcfile = HP_BWP_MP_R80_filtered_subset1000_2_mcmc.txt 

    speciesdelimitation = 0
    speciestree = 0 
    species&tree = 4  BR  BM  MP  HP  
                      4  27  9  3  
                    (((BR, BM), MP), HP); 
    phase =  0 0 0 0
    cleandata = 0
	model = GTR
    usedata = 1
    nloci = 1000  
    thetaprior = 3  0.00033  e
    tauprior = invgamma 3  0.015 
    finetune = 1: 0.02 0.02 0.02 0.02 0.02 0.02 0.02 
    print = 1 0 0 0 
    burnin = 20000 
    sampfreq = 2 
    nsample = 1000000 
    checkpoint = 8000 2000
    migprior = 2 200 
    migration = 2
                BR BM
                BM BR
