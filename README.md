# Functional-analysis-of-a-gene-locus
These data and code are related to a Developmental Cell 2023 manuscript: "Functional analysis of a gene locus in response to non-canonical combinations of transcription factors". The data include spatiotemporal expression levels of gap proteins in gap mutants and the code includes analyses reported in the manuscript and generates its figures.

# Gap Data:
Files in this folder include protein fluorescent intensities for the 4 main gap gene protein products  along the anterior-posterior (AP) axis of individual drosophila embryos  (see the manuscript text for details). There are 4 data files, each related to one gap mutant, as indicated by the genotype nomenclature in the file name. 

# Eve Data:
Files in this folder include similar AP profiles of protein fluorescent intensities for the Eve protein and a single gap gene protein, for the manipulated gap gene in the given experiment. Here too there are 4 data files, each related to one gap mutant, as indicated by the genotype nomenclature in the file name (see the manuscript text for details). 

# Data files contain the following variables:
“Age” this is the age of an embryo in minutes counted from the estimated onset of nuclear cycle 14 of the drosophila embryo.
“ix” this is the serial index of the embryo in the imaging process.
“Hb” this is the raw fluorescent levels for the antibody labeling the protein Hunchback
“Kr” this is the raw fluorescent levels for the antibody labeling the protein Krupple
“Gt” this is the raw fluorescent levels for the antibody labeling the protein Giant
“Kni” this is the raw fluorescent levels for the antibody labeling the protein Knirps
“Eve” this is the raw fluorescent levels for the antibody labeling the protein Eve


# Code:
This Matlab code is organized according to the paper's figures. The main files are those named "FigXnSx" (for figures 1-5 and their accompanying supplementary figures). The only exception is the chi-square analyses (from figures 5 and S5),that is given in the 6th file and includes an extra panel comparable to this analysis in Petkova et al. 2019.
To run these files, please make sure the code and data folders are added to your Matlab path. Files must run according to figure number order (from 1-6), since the first files are generating processed data that is saved and loaded for generating the following analyzes and/or figures. Other files are being called for, as mentioned in the headers of the relevant compartments.
Before  running the first file, it is necesary to merge the mat files "gt+gt-line4gapPART1" and "gt+gt-line4gapPART2", and save it as "gt+gt-line4gap" under the GapData folder (the original file was to large to be included here).


For further information or assistance, please email me at netta.haroush@gmail.com

