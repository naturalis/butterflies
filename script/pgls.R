#!/usr/bin/env Rscript
library(caper)
library(optparse)

# This script performs a phylogenetically-corrected generalized least squares
# regression. The inputs are a phylogenetic tree in NEXUS format, and a spreadsheet
# with data in tab-separated format. The data and tree are reconciled with one
# another by string matching the tip labels in the tree with the values in the
# column "species" in the table. The analysis attempts to predict the following
# vulnerability indicators:
# red_list   = Red List states, recoded such that NT and upward = 1
# endemicity = binary trait, European endemics = 1
# range_size = after Van Swaay et al., 2010
# ssi        = habitat specificity (SSI, Julliard et al., 2006)
# is_natural = affinity for natural habitats, after Van Swaay et al. (2006)
# By constructing models from pc_b1, pc_b2, pc_b3, pc_c1 and pc_c2

# process command line arguments
option_list = list(
  
  # input tree
  make_option(
    c("-t", "--tree"), 
    type="character", 
    default=NULL, # probably grafted.pruned.nex
    help="input tree file in NEXUS format",
    metavar="character"
  ),
  
  # input trait file
  make_option(
    c("-d", "--data"),
    type="character",
    default=NULL, # probably Tables_EU_traits.tsv
    help="input data file in TSV format",
    metavar="character"
  )
);
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# read the files
data <- read.table(opt$data, sep="\t", header=TRUE);
tree <- read.nexus(opt$tree);

# merge into a comparative data object
cdat <- comparative.data(
  phy=tree,
  data=data,
  names.col="species", # column to join on
  vcv=TRUE, # create variance/covariance matrix, needed for pgls
  na.omit=FALSE, # do NOT omit species with missing data
  warn.dropped=TRUE # warn about any species not shared by tree & data
);

################################################################################
# attempt to predict the Red List
model.pgls.red_list<-pgls(
  red_list~pc_b1+pc_b2+pc_b3+pc_c1+pc_c2, # model
  data=cdat,
  lambda="ML"  
);

# report results
summary(model.pgls.red_list);

################################################################################
# attempt to predict endemicity
model.pgls.endemicity<-pgls(
  endemicity~pc_b1+pc_b2+pc_b3+pc_c1+pc_c2, # model
  data=cdat,
  lambda="ML"  
);

# report results
summary(model.pgls.endemicity);

################################################################################
# attempt to predict range size
model.pgls.range_size<-pgls(
  range_size~pc_b1+pc_b2+pc_b3+pc_c1+pc_c2, # model
  data=cdat,
  lambda="ML"  
);

# report results
summary(model.pgls.range_size);

################################################################################
# attempt to predict range size
model.pgls.ssi<-pgls(
  ssi~pc_b1+pc_b2+pc_b3+pc_c1+pc_c2, # model
  data=cdat,
  lambda="ML"  
);

# report results
summary(model.pgls.ssi);

################################################################################
# attempt to predict affinity for natural habitats
model.pgls.is_natural<-pgls(
  is_natural~pc_b1+pc_b2+pc_b3+pc_c1+pc_c2, # model
  data=cdat,
  lambda="ML"  
);

# report results
summary(model.pgls.is_natural);


