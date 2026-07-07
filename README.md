# Mangrove and Blue-winged Pitta Phylogeography
Code repository for Mangrove and Blue-winged Pitta phylogeography analyses

The top-level directory contains three folders:
1. HP_cucullata_Genome: This folder contains scripts used to rescaffold the Pitta sordida cucullata genome (Ericson et al. 2019) based on the Chiroxiphia lanceolata genome, using the RagTag v2.1.0 software (Alonge et al, 2022)
2. ddRAD: This folder contains the scripts used to analyse the ddRAD dataset, including scripts used to generate palaeogeographical models of Southeast Asia. This folder contains the following subfolders:
   - 0_Process_Radtags: Contains script for demultiplexing raw sequence data
   - 2_Reference_Mapping: Scripts for aligning raw reads to the rescaffolded Pitta sordida cucullata genome
   - 3_Ref_Map: Scripts for running the STACKS gstacks pipeline for SNP calling
   - 4_Stacksout: Scripts for filtering and generating SNP matrices. The main script in this folder is SNPfiltr_analysis.R
   - 5_Phylogenetics: Scripts for running IQ-TREE and BPP to reconstruct the phlylogenetic history of the Mangrove and Blue-winged Pitta complex
   - 8_Biogeography: Scripts for reconstructing the ancestral range of the Mangrove Pitta
   - 9_fastsimcoal: Scripts for running fastsimcoal2
3. mtDNA: This folder contains scripts used to generate the time-calibrated BEAST tree based on 13 mitochondrial genes.
4. nichemodelling: Contains the core scripts for the palaeoecological analyses reconstructing the historical ranges of Mangrove and Blue-winged Pittas, including SLiM scripts
