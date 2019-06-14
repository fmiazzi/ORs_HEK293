## Settins & packages ----
# >sessionInfo()
#R version 3.5.2 (2018-12-20)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: macOS Mojave 10.14.3

# Matrix products: default
# BLAS: /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
# 
# locale:
#   [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
#
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
#
# loaded via a namespace (and not attached):
#   [1] compiler_3.5.2 tools_3.5.2    yaml_2.2.0    

# Set working directory
my_wd <- "/Users/basilio8/Documents/Postdoc/InsectORexpression/HighThroughput/Data";
setwd(my_wd);

# Loaded packages
library(scales);
library(graphics)
library(drc); # for curve fitting
library(ggplot2); # for graphs
library(scales); # for custom graph scales
library(FSA); # for Dunnet's post-hoc test
library(DescTools); # Another Dunnet's post-hoc test
library(RVAideMemoire); # For median tests
library(gridExtra); # For custom grids
library(ggthemes); # Extra graph themes for ggplot
library(grid);
library(extrafont);
library(lmtest);
library(stringr);
library(data.table);
library(rcompanion);
library(multcompView);

### -------------------------
### --- General Functions ---
### -------------------------
# Loads as an element of a list each .csv file within the folder 'file.path' in the working directory
read_high.through <- function(file.path) {
  table1 <- read.table(file.path, header = FALSE, skip = 3);
  table2 <- table1[,-(1:4)]; # removes the times and the AVERAGE ROI
  table2 <- table2[-1,c(rep(FALSE,2),TRUE)]; # Keeps only the Ratio values
  colnames(table2) <- paste0(rep("ROI", ncol(table2)), seq(1,ncol(table2))); # Set ROI names as columns names
  rownames(table2) <- seq(1:nrow(table2)); # Set frame number (not frame time) as row names
  output <- apply(table2, 2, function(x) as.numeric(x)); # Making sure that all values are correctly coded as numeric
}

# Loads data
load_data <- function(file.path) {
  files <- dir(file.path, pattern = '*.txt', full.names = TRUE); # Gets the file paths
  names <- dir(file.path, pattern = '*.txt', full.names = FALSE); # Gets the file names
  names.noext <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(names)); # Removes the extension from the file names
  output <- lapply(files, read_high.through); # Imports all the .csv files as elements of a list
  output <- lapply(output, as.data.frame); # Transforms each csv file in a data.frame
  names(output) <- names.noext; # Names each element using the file name w/out extension
  output; # Gives the renamed list as output
}

# High level function for data import
import_data <- function(path.to.data) {
  folders <- list.dirs(path = path.to.data, full.names = TRUE, recursive = FALSE); # Gets the list of all folders
  folder.names <- list.dirs(path = path.to.data, full.names = FALSE, recursive = FALSE); # Gets the names of all folders
  output <- lapply(folders, load_data); # Loads all data in folders
  names(output) <- folder.names; # Names the list objects according to the folder names (unique experiment identifier)
  output;
}

# Metadata import
import_metadata <- function(path.to.metadata) {
  metadata.names <- list.files(path = path.to.metadata, pattern = '*.csv', full.names = TRUE, recursive = FALSE); # Get list of metadata file names
  output <- lapply(metadata.names, read.csv); # Imports all the .csv files as elements of a list
  names(output) <- gsub(".csv","",
                        list.files(path.to.metadata,full.names = FALSE),
                        fixed = TRUE); # Renames list entries with the unique experiment identifier
  output;
}

# Calculates the local maximum value for VUAA1
loc_maxVUAA1 <- function(x){
  output <- x[min(which(diff(sign(diff(x[30:40])))==-2)+1) + 29];
}

# Filtering of input data; removes:
# - ROIs with instable baselevel (standard deviation > 0.05)
# - ROIs with mean beaslevel > 2
# - ROIs with monotonic VUAA1 response (no local maximum -> impossible to calculate comparison to odor response)
# - ROIs non responding to VUAA1 (after VUAA1 > 1.5*baseline)
data_filtering <- function(DataFrame) {
  # Filtering ROIs with sd > 0.05
  data.sd <- apply(as.data.frame(DataFrame[1:10,]), 2, function(x) sd(x, na.rm = TRUE)); #baselevel st.dev calculation
  sd.todelete <- data.sd > 0.05; # st.dev threshold
  data.filtered <- DataFrame;
  data.filtered[,which(sd.todelete == TRUE)] <-  NA; # delete ROis with sd > threshold
  # Filtering ROIs with baselevel mean higher than 2
  data.mean <- apply(DataFrame[1:10,], 2, function(x) mean(x, na.rm = TRUE)); # baselevel mean calculation
  mean.todelete <- data.mean > 2;
  data.filtered[,which(mean.todelete == TRUE)] <-  NA;
  # Filtering ROIs with monotonic VUAA1 increase
  locMax.todelete <- apply(DataFrame, 2, loc_maxVUAA1); # Finds the value of first VUAA1 local maximum
  data.filtered[, which(is.na(locMax.todelete) == TRUE)] <- NA;
  # Filtering ROIs with VUAA1.response < 1.5* baselevels
  VUAA1loc.max <- apply(DataFrame, 2, loc_maxVUAA1);
  VUAA1loc.max[VUAA1loc.max[]==Inf] <- NA;
  VUAA1.int <- VUAA1loc.max - data.mean;
  VUAA1.tokeep <- as.logical(VUAA1loc.max > data.mean*1.5);
  data.filtered[,which(VUAA1.tokeep != TRUE)] <-  NA;
  # Deleting filtered ROIs
  output <- data.filtered[,-(which(colSums(is.na(data.filtered))>0))];
}

ts_resp.calc <- function(ListData, ListMetadata, VariableFactor, FactorLevels) {
  factor.tested <- unique(as.vector(ListMetadata[[VariableFactor]])); # All concentrations tested (for dose response curve)
  output <- vector(mode = "list", length = length(factor.tested)); # Predefines output list
  names(output) <- c(factor.tested);
  time.length <- nrow(ListData[1][[1]]);
  for (i in factor.tested) {
    factor.wells <- ListMetadata$well.name[which(ListMetadata[[VariableFactor]] == i)]; # Gets all wells with the same concentration tested
    data.factor <- list(ListData[c(as.character(factor.wells))]); # Extracts all data from selected wells
    data.factor <- lapply(data.factor, cbind.data.frame); # Makes a data frame with the seleceted data
    data.factor <- data.factor[sapply(data.factor, function(x) length(x)) > 0];
    f0 <- lapply(data.factor, function(x) as.data.frame(apply(x[1:10,], 2, function(x) rep(mean(x), time.length))));
    dF.F0 <- mapply(function(x,y) (x-y)*(100/y), data.factor, f0, SIMPLIFY = FALSE);
    if (length(dF.F0) > 0) {
      output[which(factor.tested == i)] <- lapply(dF.F0, function(x) apply(x, 1, median)); # Calculates the median for each row
    } else {
      output[which(factor.tested == i)] <- as.vector(rep(NA, times = time.length), mode = "double");
      }
    }
  output <- list(output); # Important! Otherwise the names do not match! 
}

odor.df_mkr <- function(DataList, VariableFactor) {
  odor.df <- cbind.data.frame(DataList);
  names.odor <- names(DataList);
  length.odor <- lapply(DataList, length);
  names.elements <- lapply(DataList, names);
  col.names <- paste(rep(names.odor, unlist(length.odor, use.names = FALSE)), unlist(names.elements, use.names = FALSE), sep = ".");
  names(odor.df) <- col.names;
  conc.tested <- unique(unlist(names.elements, use.names = FALSE));
  odor.conc.list <- vector("list", length(conc.tested));
  names(odor.conc.list) <- conc.tested;
  for (i in conc.tested) {
    odor.conc.list.i <- odor.df[,which(substrRight(col.names, 5) == i)];
    names(odor.conc.list.i) <- substr(names(odor.conc.list.i), 1, 7);
    odor.conc.list[[i]] <- odor.conc.list.i;
  }
  odor.mean.list <- lapply(odor.conc.list, function(x) apply(x, 1, mean));
  odor.sd.list <- lapply(odor.conc.list, function(x) apply(x, 1, sd));
  time.list <- lapply(odor.conc.list, function(x) seq(0, (nrow(x)-1))*5);
  odor.conc.list <- rep(conc.tested, each = length(odor.mean.list[[1]]));
  output <- cbind.data.frame(odor.mean = unlist(odor.mean.list, recursive = TRUE, use.names = FALSE),
                             odor.sd = unlist(odor.sd.list, recursive = TRUE, use.names = FALSE),
                             time.sec = unlist(time.list, recursive = TRUE, use.names = FALSE),
                             odor.conc = odor.conc.list);
  colnames(output) <- c("odor.mean", "odor.sd", "time.sec", "odor.conc");
  output$odor.conc <- factor(output$odor.conc, levels = VariableFactor)
  output;
}

# Calculates the intensity of VUAA1 response (internal positive control)
VUAA1_resp.calc <- function(ListData, ListMetadata, VariableFactor, FactorLevels, facets.pent.ac.labels) {
  factor.tested <- unique(as.vector(ListMetadata[[VariableFactor]])); # Gets all unique tested concentrations
  output.list <- vector(mode = "list", length = length(factor.tested)); # Preallocates the output
  names(output.list) <- c(factor.tested);
  for (i in factor.tested) {
    factor.wells <- ListMetadata$well.name[which(ListMetadata[[VariableFactor]] == i)]; # Gets all wells tested with the same odor concentration
    factor.data <- list(ListData[c(as.character(factor.wells))]); # Gets all ROIs for that well
    well.mean <- lapply(factor.data, function(x) lapply(x, function(x) apply(as.data.frame(x), 1, mean))); # Calulates the mean for each row per well
    # Calculates the difference between the peak VUAA1 response and the baselevel before VUAA1 application
    well.VUAA1.resp <- lapply(well.mean, function(x) lapply(x, function(x) (mean(x[32:33]) - mean(x[25:30]))));
    output.list[which(factor.tested == i)] <- well.VUAA1.resp; # Assigns the result to the correspondent output.list entry
  }
  VUAA1.resp <- unlist(output.list, use.names = FALSE); # Unlisted VUAA1 response intensities
  well.name <- unlist(lapply(output.list, function(x) (names(x))), use.names = FALSE);
  well.letters <- as.character(substr(well.name, 1, 1));
  plate.half <- str_replace_all(well.letters, "B|C|D", "1");
  plate.half <- as.factor(str_replace_all(plate.half, "E|F|G", "2"));
  factor.str <- as.character(rep(factor.tested, unlist(lapply(output.list, function(x) length(x)), use.names = FALSE)));
  factor.plate <- paste(factor.str, plate.half, sep = ".");
  file.name <- rep(unique(ListMetadata$file.name), length(VUAA1.resp));
  output <- cbind.data.frame(file.name, VUAA1.resp, well.name, plate.half, factor.str, factor.plate);
  names(output) <- c("file.name", "VUAA1.resp", "well.name", "plate.half", VariableFactor, "factor.plate");
  #  output[[factor.str]] <- factor(output[[factor.str]], levels = FactorLevels);
  output <- list(output);
}

odor.fits_calc <- function(DataFrame) {
  out.l3 <- drm(data = DataFrame, mean.odor.norm ~ odor.mol10, fct = L.3(), na.action = na.omit); # The fit of the actual data with the drm package
  # original one ->  output <- expand.grid(conc=exp(seq(log10(1.00e-7), log10(1.00e-03), length=100))); # Set up the interval of simulated data (x-axis)
  output <- expand.grid(conc=seq(-8, -3, length=100)); # Set up the interval of simulated data (x-axis)
  out.predict <- predict(out.l3, newdata=output, interval="confidence"); # Get the prediction values based on the fitted model
  output$p <- out.predict[,1]; # Appends the mean predicted values
  output$pmin <- out.predict[,2]; # Appends the minimum prediced values
  output$pmax <- out.predict[,3]; # Appends the maximum predicted values
  output$codes <- rep(unique(DataFrame$codes), times = 100);
  output$xvals <- exp(seq(log(1.00e-8), log(1.00e-03), length=100));
  output;
}

# Publication-ready theme for ggplot (Modified from Koundinya Desiraju)
# https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=7, base_family="ArialMT") {
  #theme_Publication <- function(base_size=14) {
  #  (theme_foundation(base_size=base_size)
  (theme_foundation(base_size=base_size, base_family=base_family)
   + theme(plot.title = element_text(face = "bold",
                                     size = rel(1.2), hjust = 0.5),
           text = element_text(),
           panel.background = element_rect(colour = NA),
           plot.background = element_rect(colour = NA),
           panel.border = element_rect(colour = NA),
           axis.title = element_text(size = rel(1)),
           axis.title.y = element_text(angle=90,vjust =0),
           axis.title.x = element_text(vjust = -0.2),
           axis.text.x = element_text(size = 7),
           axis.text.y = element_text(size = 7),
           axis.text = element_text(), 
           axis.line = element_line(colour="black"),
           axis.ticks = element_line(),
           panel.grid.minor = element_blank(),
           panel.grid.major = element_blank(),
           #            panel.grid.major = element_line(colour="#f0f0f0"),
           legend.key = element_rect(colour = NA),
           legend.position = "bottom",
           legend.direction = "horizontal",
           legend.key.size= unit(0.2, "cm"),
           legend.margin = margin(t = 0, unit='cm'),
           legend.title = element_text(face="italic"),
           plot.margin=unit(c(1,1,1,1),"mm"),
           strip.background=element_rect(colour="#00FFFFFF",fill="#00FFFFFF"),
           strip.text = element_text(size = 7, face="bold")
   ))
}

### ------------------------------
### --- Dose-Response Function ---
### ------------------------------

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

# Calculates the treatment grand mean for each independent experiment (unique experiment identifier)
means.calc_per.treat <- function(ListData, ListMetadata) {
  conc.tested <- unique(as.vector(ListMetadata$odor.concentration)); # All concentrations tested (for dose response curve)
  output <- vector(mode = "list", length = length(conc.tested)); # Predefines output list
  names(output) <- c(conc.tested);
  for (i in conc.tested) {
    conc.wells <- ListMetadata$well.name[which(ListMetadata$odor.concentration == i)]; # Gets all wells with the same concentration tested
    data.conc <- list(ListData[c(as.character(conc.wells))]); # Extracts all data from selected wells
    data.conc2 <- lapply(data.conc, cbind.data.frame); # Makes a data frame with the seleceted data
    output[which(conc.tested == i)] <- lapply(data.conc2, function(x) apply(x, 1, mean)); # Calculates the mean for each row
  }
  output <- list(output); # Important! Otherwise the names do not match! 
}

# 
norm_odor.resp.hOr47a <- function(ListData, ListMetadata) {
  file.name.unique <- unique(ListMetadata$file.name);
  conc.tested <- unique(as.vector(ListMetadata$odor.concentration));
  output <- data.frame(file.name = rep(file.name.unique, length(conc.tested)),
                       mean.odor.norm = numeric(length(conc.tested)),
                       odor.conc.num = conc.tested,
                       odor.mol10 = log10(conc.tested),
                       odor.conc = factor(conc.tested, levels = c("1e-08", "1e-07", "3e-07", "1e-06", "3e-06",
                                                                  "1e-05", "3e-05", "1e-04", "3e-04", "0.001")));
  for (i in conc.tested) {
    well.subset <- as.vector(ListMetadata[which(ListMetadata$odor.concentration == i), "well.name"]);
    resp.subset <- cbind.data.frame(ListData[well.subset]);
    VUAA1.resp <- apply(resp.subset, 2, function(x) x[min(which(diff(sign(diff(x[30:40])))==-2)+1) + 29] - mean(x[25:30]));
    odor.resp <- apply(resp.subset, 2, function(x) x[12] - mean(x[1:10]));
    mean.odor.norm <- mean((odor.resp/VUAA1.resp)*100);
    output[which(output$odor.conc.num == i),"mean.odor.norm"] <- mean.odor.norm;
  }
  output <- list(output);
}

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

norm_odor.resp.hOr56a <- function(ListData, ListMetadata) {
  file.name.unique <- unique(ListMetadata$file.name);
  conc.tested <- unique(as.vector(ListMetadata$odor.concentration));
  output <- data.frame(file.name = rep(file.name.unique, length(conc.tested)),
                       mean.odor.norm = numeric(length(conc.tested)),
                       odor.conc.num = conc.tested,
                       odor.mol10 = log10(conc.tested),
                       odor.conc = factor(conc.tested, levels = c("1e-08", "3e-08", "1e-07", "3e-07", "1e-06",
                                                                  "3e-06", "1e-05", "3e-05")));
  for (i in conc.tested) {
    well.subset <- as.vector(ListMetadata[which(ListMetadata$odor.concentration == i), "well.name"]);
    resp.subset <- cbind.data.frame(ListData[well.subset]);
    VUAA1.resp <- apply(resp.subset, 2, function(x) x[min(which(diff(sign(diff(x[30:40])))==-2)+1) + 29] - mean(x[25:30]));
    odor.resp <- apply(resp.subset, 2, function(x) x[12] - mean(x[1:10]));
    mean.odor.norm <- mean((odor.resp/VUAA1.resp)*100);
    output[which(output$odor.conc.num == i),"mean.odor.norm"] <- mean.odor.norm;
  }
  output <- list(output);
}

odor.fits_calc.hOr56a <- function(DataFrame) {
  out.l3 <- drm(data = DataFrame, mean.odor.norm ~ odor.mol10, fct = L.3(), na.action = na.omit); # The fit of the actual data with the drm package
  # original one ->  output <- expand.grid(conc=exp(seq(log10(1.00e-7), log10(1.00e-03), length=100))); # Set up the interval of simulated data (x-axis)
  output <- expand.grid(conc=seq(-8, -5, length=100)); # Set up the interval of simulated data (x-axis)
  out.predict <- predict(out.l3, newdata=output, interval="confidence"); # Get the prediction values based on the fitted model
  output$p <- out.predict[,1]; # Appends the mean predicted values
  output$pmin <- out.predict[,2]; # Appends the minimum prediced values
  output$pmax <- out.predict[,3]; # Appends the maximum predicted values
  output$codes <- rep(unique(DataFrame$codes), times = 100);
  output$xvals <- exp(seq(log(1.00e-8), log(3.00e-05), length=100));
  output;
}

means.calc_per.treat <- function(ListData, ListMetadata) {
  conc.tested <- unique(as.vector(ListMetadata$odor.concentration)); # All concentrations tested (for dose response curve)
  output <- vector(mode = "list", length = length(conc.tested)); # Predefines output list
  names(output) <- c(conc.tested);
  for (i in conc.tested) {
    conc.wells <- ListMetadata$well.name[which(ListMetadata$odor.concentration == i)]; # Gets all wells with the same concentration tested
    data.conc <- list(ListData[c(as.character(conc.wells))]); # Extracts all data from selected wells
    data.conc2 <- lapply(data.conc, cbind.data.frame); # Makes a data frame with the seleceted data
    output[which(conc.tested == i)] <- lapply(data.conc2, function(x) apply(x, 1, mean)); # Calculates the mean for each row
  }
  output <- list(output); # Important! Otherwise the names do not match! 
}

### ----------------------------------------------
### --- Odor Specificity Functions & Constants ---
### ----------------------------------------------

odor_labeller <- function(variable,value){
  return(odor.labels[value])
}

# This cannot be used for concentration data evaluation!
norm_resp <- function(ListData, ListMetadata, VariableFactor, FactorLevels) {
  file.name.unique <- unique(ListMetadata$file.name);
  factor.tested <- unique(as.vector(ListMetadata[[VariableFactor]]));
  output <- data.frame(file.name = rep(file.name.unique, length(factor.tested)),
                       mean.norm = numeric(length(factor.tested)),
                       factor.name = factor(factor.tested, levels = FactorLevels));
  for (i in factor.tested) {
    well.subset <- as.vector(ListMetadata[which(ListMetadata[[VariableFactor]] == i), "well.name"]);
    resp.subset <- cbind.data.frame(ListData[well.subset]);
    VUAA1.resp <- apply(resp.subset, 2, function(x) x[min(which(diff(sign(diff(x[30:40])))==-2)+1) + 29] - mean(x[25:30]));
    odor.resp <- apply(resp.subset, 2, function(x) x[12] - mean(x[1:10]));
    mean.norm <- mean((odor.resp/VUAA1.resp)*100);
    output[which(output$factor.name == i),"mean.norm"] <- mean.norm;
  }
  output <- list(output);
}

specificity.df_mkr <- function(DataList, VariableLevels) {
  odor.tested.list <- vector("list", length(VariableLevels));
  names(odor.tested.list) <- VariableLevels;
  odor.df <- cbind.data.frame(unlist(DataList, use.names = TRUE, recursive = FALSE));
  col.names <- names(odor.df);
  for (i in VariableLevels) {
    odor.tested.list.i <- odor.df[,which(substr(col.names, 9, 100) == i)];
    names(odor.tested.list.i) <- substr(names(odor.tested.list.i), 1, 7);
    odor.tested.list[[i]] <- odor.tested.list.i;
  }
  odor.tested.list;
  odor.mean.list <- lapply(odor.tested.list, function(x) apply(x, 1, function(x) mean(x, na.rm = TRUE)));
  odor.sd.list <- lapply(odor.tested.list, function(x) apply(x, 1, function(x) sd(x, na.rm = TRUE)));
  time.list <- lapply(odor.tested.list, function(x) seq(0, (nrow(x)-1))*5);
  odor.names.list <- rep(VariableLevels, each = length(odor.mean.list[[1]]));
  output <- cbind.data.frame(odor.mean = unlist(odor.mean.list, recursive = TRUE, use.names = FALSE),
                             odor.sd = unlist(odor.sd.list, recursive = TRUE, use.names = FALSE),
                             time.sec = unlist(time.list, recursive = TRUE, use.names = FALSE),
                             odor.name = odor.names.list);
  colnames(output) <- c("odor.mean", "odor.sd", "time.sec", "odor.name");
  output$odor.name <- factor(output$odor.name, levels = VariableLevels)
  output;
}

# Filtering of input data; removes:
# - ROIs with instable baselevel (standard deviation > 0.05)
# - ROIs with mean beaslevel > 2
# - ROIs with monotonic VUAA1 response (no local maximum -> impossible to calculate comparison to odor response)
# - ROIs non responding to VUAA1 (after VUAA1 > 1.5*baseline)
data_filtering.untr <- function(DataFrame) {
  # Filtering ROIs with sd > 0.05
  data.sd <- apply(as.data.frame(DataFrame[1:10,]), 2, function(x) sd(x, na.rm = TRUE)); #baselevel st.dev calculation
  sd.todelete <- data.sd > 0.05; # st.dev threshold
  data.filtered <- DataFrame;
  data.filtered[,which(sd.todelete == TRUE)] <-  NA; # delete ROis with sd > threshold
  # Filtering ROIs with baselevel mean higher than 2
  data.mean <- apply(DataFrame[1:10,], 2, function(x) mean(x, na.rm = TRUE)); # baselevel mean calculation
  mean.todelete <- data.mean > 2;
  data.filtered[,which(mean.todelete == TRUE)] <-  NA;
  # Deleting filtered ROIs
  output <- data.filtered[,-(which(colSums(is.na(data.filtered))>0))];
}

norm_resp.untr <- function(ListData, ListMetadata, VariableFactor, FactorLevels) {
  file.name.unique <- unique(ListMetadata$file.name);
  factor.tested <- unique(as.vector(ListMetadata[[VariableFactor]]));
  output <- data.frame(file.name = rep(file.name.unique, length(factor.tested)),
                       mean.norm = numeric(length(factor.tested)),
                       factor.name = factor(factor.tested, levels = FactorLevels));
  for (i in factor.tested) {
    well.subset <- as.vector(ListMetadata[which(ListMetadata[[VariableFactor]] == i), "well.name"]);
    resp.subset <- cbind.data.frame(ListData[well.subset]);
#    VUAA1.resp <- apply(resp.subset, 2, function(x) x[min(which(diff(sign(diff(x[30:40])))==-2)+1) + 29] - mean(x[25:30]));
    odor.resp <- apply(resp.subset, 2, function(x) x[12] - mean(x[1:10]));
    mean.norm <- mean(odor.resp);
    output[which(output$factor.name == i),"mean.norm"] <- mean.norm;
  }
  output <- list(output);
}

## Data Analysis ----

# -----------------------------------------
# --- DoseResponse data analysis hOr47a ---
# -----------------------------------------

# Loading data
list.of.data.hOr47a <- import_data("Input/DoseResponse/hOr47a"); # Load data files as single list
list.of.metadata.hOr47a <- import_metadata("Metadata/DoseResponse/hOr47a"); # Load metadata as single list

# Filter out non responding and spurious ROIs
rois.filtered.hOr47a <- lapply(list.of.data.hOr47a, function(x) lapply(x, data_filtering));

# -----------------------------
# --- Normalization Example ---
# -----------------------------

example.norm <- norm_resp.example(example.df)

norm_resp.example <- function(DataFrame) {
  roi.names <- names(DataFrame);
  VUAA1.resp <- apply(DataFrame, 2, function(x) x[min(which(diff(sign(diff(x[30:40])))==-2)+1) + 29] - mean(x[25:30]));
  VUAA1.df <- data.frame(matrix(rep(VUAA1.resp, each = 25), nrow = 25));
  odor.resp <- apply(DataFrame, 2, function(x) x[1:25] - mean(x[1:10]));
  norm.resp <- (odor.resp/VUAA1.df)*100;
  output <- data.frame(norm.resp);
  names(output) <- roi.names;
  output;
}

example.palette <- c(
  "C10.ROI1" = "#56B4E9", 
  "10.ROI21" = "#CC79A7");

example.df <- data.frame("C10.ROI1" = rois.filtered.hOr47a[["181113b"]][["C10"]][["ROI1"]],
                         "C10.ROI21" = rois.filtered.hOr47a[["181113b"]][["C10"]][["ROI5"]]);

example.plot.nonorm <- ggplot(data = example.df) + # Feeds the data and the variables
  geom_line(aes(x = seq(0,195,5), y = example.df$C10.ROI1, colour = "#56B4E9")) + # The line displays the mean values aes(colour= example.palette)
  geom_line(aes(x = seq(0,195,5), y = example.df$C10.ROI21, colour = "#CC79A7")) +
#  scale_colour_manual(values = example.palette) +
  # The ribbon displays the mean ± st.dev values
#  facet_wrap( ~ trace.type) + #, labeller = labeller(odor.conc = facets.pent.ac.labels)) + # Grouping according to the stimulus concentrations (groups)
  labs(x = "Time (s)", y = "R (340/380 nm)") +
  #  scale_fill_gradient(name = "Concentration", trans = "log") + # Sets a log scale for the coloring of the different treatments
  scale_x_continuous(breaks=seq(0,200,50)) +
  scale_y_continuous(breaks=seq(0,5,0.5)) +
  coord_cartesian(x = c(0, 200), y = c(0.8,2)) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
example.plot.nonorm;
ggsave("Output/example.pdf", plot = example.plot.nonorm, device = "pdf",  width = 35, height = 35, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

example.plot.norm <- ggplot(data = example.norm) + # Feeds the data and the variables
  geom_line(aes(x = seq(0,120,5), y = example.norm$C10.ROI1, colour = "#56B4E9")) + # The line displays the mean values aes(colour= example.palette)
  geom_line(aes(x = seq(0,120,5), y = example.norm$C10.ROI21, colour = "#CC79A7")) +
  #  scale_colour_manual(values = example.palette) +
  # The ribbon displays the mean ± st.dev values
  #  facet_wrap( ~ trace.type) + #, labeller = labeller(odor.conc = facets.pent.ac.labels)) + # Grouping according to the stimulus concentrations (groups)
  labs(x = "Time (s)", y = expression(paste(rho," (%)"))) +
  #  scale_fill_gradient(name = "Concentration", trans = "log") + # Sets a log scale for the coloring of the different treatments
  scale_x_continuous(breaks=seq(0,200,40)) +
  scale_y_continuous(breaks=seq(0,100,20)) +
  coord_cartesian( y = c(0,100), x = c(0, 120)) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
example.plot.norm;
ggsave("Output/example.norm.pdf", plot = example.plot.norm, device = "pdf",  width = 35, height = 35, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

# -----------------------------

# Time-series of pentyl.acetate

conc.levels.hOr47a <- c(1e-08, 3e-07, 1e-06, 3e-06, 1e-05, 3e-05, 1e-04, 3e-04, 1e-03);

pent.ac.ts <- mapply(function(x,y) ts_resp.calc(x, y, "odor.concentration", conc.levels.hOr47a), rois.filtered.hOr47a, list.of.metadata.hOr47a);

pent.ac.df <- odor.df_mkr(pent.ac.ts, conc.levels.hOr47a);

facets.pent.ac.labels <- c("1e-08" = "DMSO",
                        "3e-07" = "300 nM",
                        "1e-06" = "1 µM",
                        "3e-06" = "3 µM",
                        "1e-05" = "10 µM",
                        "3e-05" = "30 µM",
                        "1e-04" = "100 µM",
                        "3e-04" = "300 µM",
                        "0.001" = "1 mM");

conc.palette.pent.ac <- c(
  "1e-08" = "#000000", 
  "3e-07" = "#e6cc4a", 
  "1e-06" = "#ffc04c",
  "3e-06" = "#e9b061", 
  "1e-05" = "#ffa500", 
  "3e-05" = "#D55E00", 
  "1e-04" = "#e96c61", 
  "3e-04" = "#df2e1f", 
  "0.001" = "#CC79A7");

pent.ac.ts.plot <- ggplot(data = pent.ac.df, aes(x = time.sec, y = odor.mean)) + # Feeds the data and the variables
  geom_line(aes(colour= odor.conc)) + # The line displays the mean values
  scale_colour_manual(values = conc.palette.pent.ac) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax=(odor.mean + odor.sd), ymin=(odor.mean - odor.sd), fill= odor.conc), alpha = 0.4) +
  scale_fill_manual(values = conc.palette.pent.ac) +
  facet_wrap( ~ odor.conc, ncol=3, labeller = labeller(odor.conc = facets.pent.ac.labels)) + # Grouping according to the stimulus concentrations (groups)
  labs(x = "Time (s)", y = "R (340/380 nm)") +
  #  scale_fill_gradient(name = "Concentration", trans = "log") + # Sets a log scale for the coloring of the different treatments
  scale_x_continuous(breaks=seq(0,200,50)) +
#  scale_y_continuous(breaks=seq(0,1.75,0.25)) +
  coord_cartesian(x = c(0, 200)) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
pent.ac.ts.plot;
ggsave("Output/pent.ac.DosResp.pdf", plot = pent.ac.ts.plot, device = "pdf",  width = 90, height = 90, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

# Rescaling of responses
odor.norm.hOr47a <- mapply(norm_odor.resp.hOr47a, rois.filtered.hOr47a, list.of.metadata.hOr47a);
odor.norm.df.hOr47a <- rbindlist(odor.norm.hOr47a, fill=TRUE);

odor.fits.hOr47a <- odor.fits_calc(odor.norm.df.hOr47a);

odor.norm.df.plot.hOr47a <- ggplot(data = odor.norm.df.hOr47a, aes(x = odor.conc.num, y = mean.odor.norm)) + # Feeds the data and the variables
  geom_point(color = "grey", size = 0.5) +
  geom_ribbon(data = odor.fits.hOr47a, aes(x = xvals, y = p, ymin = pmin, ymax = pmax, fill = "#ffa500"), alpha = 0.4) + # Displays the confidence interval of the fit
  geom_line(data = odor.fits.hOr47a, aes(x = xvals, y = p),colour = "#ffa500", size = 1) + # Displays the fit of the data
  annotation_logticks(base = 10, sides = "b", scaled = TRUE, # Sets the custom ticks for the log10 scale of the x-axis
                      short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"),
                      colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL) +
  scale_y_continuous(breaks=seq(-10,100,10)) +
  coord_cartesian(x = c(1e-07, 3e-03), y = c(-5,40)) +
  #  scale_y_continuous(breaks=seq(0,500,100), limits = c(-80, 600)) + # Sets the y-axis
  labs(x = "[Pent. Ac.](M)", y = expression(paste(rho," (%)"))) + # Sets the axis labels and title of the graph
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), # Sets the numbering format of the x-axis
                labels = trans_format("log10", math_format(10^.x))) +
  theme_Publication() +
  theme(legend.position="none");
odor.norm.df.plot.hOr47a;
ggsave("Output/odor.norm.df.plot.hOr47a.pdf", plot = odor.norm.df.plot.hOr47a, device = "pdf",  width =35 , height = 35, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

# EC50 calculation

odor.logfit.hOr47a <- drm(mean.odor.norm ~ odor.mol10, data = odor.norm.df.hOr47a, fct = L.3(names=c("Slope", "Upper Limit","EC50")), na.action = na.omit);
coeftest(odor.logfit.hOr47a);

#  t test of coefficients:
#                           Estimate Std. Error  t value  Pr(>|t|)    
#  Slope:(Intercept)       -2.892865   0.422412  -6.8484 1.036e-08 ***
#  Upper Limit:(Intercept) 30.803021   1.283097  24.0068 < 2.2e-16 ***
#  EC50:(Intercept)        -4.610570   0.063377 -72.7488 < 2.2e-16 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# EC50 value in Molar -> 10^(-4.610570) 95CI (-4.737866, -4.483275) = 2.451489e-05 (1.828664e-05, 3.286435e-05) -> 24.5 (18.3, 32.9) µM

coefci(odor.logfit.hOr47a);
#                             2.5 %    97.5 %
# Slope:(Intercept)       -3.741304 -2.044425
# Upper Limit:(Intercept) 28.225846 33.380197
# EC50:(Intercept)        -4.737866 -4.483275

plot(fitted(odor.logfit.hOr47a), residuals(odor.logfit.hOr47a));
plot(fitted(odor.logfit.hOr47a), residuals(odor.logfit.hOr47a, typeRes = "standard"));
mselect(odor.logfit.hOr47a, list(L.4(), W1.3(), W1.4(), W2.4(), baro5()), linreg = TRUE);

#          logLik       IC Lack of fit  Res var
# L.3   -142.8622 293.7245  0.02316008 13.61838
# L.4   -142.8305 295.6609  0.01302495 13.87964
# Cubic -152.7362 315.4724          NA 20.17058
# Quad  -161.7853 331.5706          NA 27.81284
# Lin   -174.4583 354.9166          NA 43.98815
# W1.3         NA       NA          NA       NA
# W1.4         NA       NA          NA       NA
# W2.4         NA       NA          NA       NA
# baro5        NA       NA          NA       NA

# -----------------------------------------
# --- DoseResponse data analysis hOr56a ---
# -----------------------------------------

list.of.data.hOr56a <- import_data("Input/DoseResponse/hOr56a"); # Load data files as single list
list.of.metadata.hOr56a <- import_metadata("Metadata/DoseResponse/hOr56a"); # Load metadata as single list

# Filter out non responding and spurious ROIs
rois.filtered.hOr56a <- lapply(list.of.data.hOr56a, function(x) lapply(x, data_filtering));

# Time-series of pentyl.acetate

conc.levels.hOr56a <- c(1e-8, 3e-08, 1e-07, 3e-07, 1e-06, 3e-06, 1e-05, 3e-05);

geosmin.ts <- mapply(function(x,y) ts_resp.calc(x, y, "odor.concentration", conc.levels.hOr56a), rois.filtered.hOr56a, list.of.metadata.hOr56a);

geosmin.df <- odor.df_mkr(geosmin.ts, conc.levels.hOr56a);

facets.geosmin.labels <- c("1e-08" = "DMSO",
                           "3e-08" = "30 nM",
                           "1e-07" = "100 nM",
                           "3e-07" = "300 nM",
                           "1e-06" = "1 µM",
                           "3e-06" = "3 µM",
                           "1e-05" = "10 µM",
                           "3e-05" = "30 µM");

conc.palette.geosmin <- c(
  "1e-08" = "#000000", 
  "3e-08" = "#e6cc4a", 
  "1e-07" = "#ffc04c",
  "3e-07" = "#e9b061", 
  "1e-06" = "#ffa500", 
  "3e-06" = "#D55E00", 
  "1e-05" = "#e96c61", 
  "3e-05" = "#df2e1f");

geosmin.ts.plot <- ggplot(data = geosmin.df, aes(x = time.sec, y = odor.mean)) + # Feeds the data and the variables
  geom_line(aes(colour= odor.conc)) + # The line displays the mean values
  scale_colour_manual(values = conc.palette.geosmin) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax=(odor.mean + odor.sd), ymin=(odor.mean - odor.sd), fill= odor.conc), alpha = 0.4) +
  scale_fill_manual(values = conc.palette.geosmin) +
  facet_wrap( ~ odor.conc, ncol=4, labeller = labeller(odor.conc = facets.geosmin.labels)) + # Grouping according to the stimulus concentrations (groups)
  labs(x = "Time (s)", y = "R (340/380 nm)") +
  #  scale_fill_gradient(name = "Concentration", trans = "log") + # Sets a log scale for the coloring of the different treatments
  scale_x_continuous(breaks=seq(0,200,50)) +
  #  scale_y_continuous(breaks=seq(0,1.75,0.25)) +
  coord_cartesian(x = c(0, 200)) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
geosmin.ts.plot;
ggsave("Output/geosmin.DosResp.pdf", plot = geosmin.ts.plot, device = "pdf",  width = 100, height = 100, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

# Rescaling of responses
odor.norm.hOr56a <- mapply(norm_odor.resp.hOr56a, rois.filtered.hOr56a, list.of.metadata.hOr56a);
odor.norm.df.hOr56a <- rbindlist(odor.norm.hOr56a, fill=TRUE);

odor.fits.hOr56a <- odor.fits_calc.hOr56a(odor.norm.df.hOr56a);

odor.norm.df.plot.hOr56a <- ggplot(data = odor.norm.df.hOr56a, aes(x = odor.conc.num, y = mean.odor.norm)) + # Feeds the data and the variables
  geom_point(color = "grey", size = 0.5) +
  geom_ribbon(data = odor.fits.hOr56a, aes(x = xvals, y = p, ymin = pmin, ymax = pmax), fill = "#FF00FF", alpha = 0.4) + # Displays the confidence interval of the fit
  geom_line(data = odor.fits.hOr56a, aes(x = xvals, y = p),colour = "#FF00FF", size = 1) + # Displays the fit of the data
  annotation_logticks(base = 10, sides = "b", scaled = TRUE, # Sets the custom ticks for the log10 scale of the x-axis
                       short = unit(0.05, "cm"), mid = unit(0.1, "cm"), long = unit(0.15, "cm"),
                       colour = "black", size = 0.5, linetype = 1, alpha = 1, color = NULL) +
  scale_y_continuous(breaks=seq(-10,100,10)) +
  coord_cartesian(x = c(1e-8, 1e-04), y = c(-5,25)) +
  #  scale_y_continuous(breaks=seq(0,500,100), limits = c(-80, 600)) + # Sets the y-axis
  labs(x = "[Geosmin](M)", y = expression(paste(rho," (%)"))) + # Sets the axis labels and title of the graph
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x), # Sets the numbering format of the x-axis
                labels = trans_format("log10", math_format(10^.x))) +
  theme_Publication() +
  theme(legend.position="none");
  odor.norm.df.plot.hOr56a;
ggsave("Output/odor.norm.df.plot.hOr56a.pdf", plot = odor.norm.df.plot.hOr56a, device = "pdf",  width =35 , height = 35, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

# EC50 calculation

odor.logfit.hOr56a <- drm(mean.odor.norm ~ odor.mol10, data = odor.norm.df.hOr56a, fct = L.3(names=c("Slope", "Upper Limit","EC50")), na.action = na.omit);
coeftest(odor.logfit.hOr56a);

# t test of coefficients:
#  
#    Estimate Std. Error  t value  Pr(>|t|)    
#    Slope:(Intercept)       -2.46059    0.58201  -4.2277 9.363e-05 ***
#    Upper Limit:(Intercept) 11.93271    0.91978  12.9734 < 2.2e-16 ***
#    EC50:(Intercept)        -6.14331    0.12844 -47.8284 < 2.2e-16 ***
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# EC50 value in Molar -> 10^(-6.14331) 95CI (-6.400936, -5.885680) = 7.189356e-07 (3.972501e-07, 1.301128e-06) -> 0.72 (0.40, 1.30) µM

coefci(odor.logfit.hOr56a);
#                              2.5 %    97.5 %
#  Slope:(Intercept)       -3.627966 -1.293218
#  Upper Limit:(Intercept) 10.087855 13.777561
#  EC50:(Intercept)        -6.400936 -5.885680

plot(fitted(odor.logfit.hOr56a), residuals(odor.logfit.hOr56a));
plot(fitted(odor.logfit.hOr56a), residuals(odor.logfit.hOr56a, typeRes = "standard"));
mselect(odor.logfit.hOr56a, list(L.4(), W1.3(), W1.4(), W2.4(), baro5()), linreg = TRUE);

#           logLik       IC Lack of fit  Res var
#  L.3   -129.6747 267.3495   0.9858509 6.349802
#  L.4   -129.4363 268.8725   0.9943092 6.417027
#  Cubic -129.6755 269.3510          NA 6.472092
#  Lin   -132.7301 271.4602          NA 6.950763
#  Quad  -132.6803 273.3606          NA 7.069326
#  W1.3         NA       NA          NA       NA
#  W1.4         NA       NA          NA       NA
#  W2.4         NA       NA          NA       NA
#  baro5        NA       NA          NA       NA

# 

# --------------------------------------
# --- Odor Specificity Data Analysis ---
# --------------------------------------

specificity.levels <- c("dmso", "3-octanol", "2-heptanone", "isobutyl.acetate", "hexyl.acetate", "3-methylthio-1-propanol",
                        "3-methylthio-1-propanol.B", "3-methylthio-1-propanol.A", "propyl.acetate", "propyl.acetate.A",
                        "propyl.acetate.B", "methyl.hexanoate", "butyl.acetate", "pentyl.acetate");

# .A and .B denote odor samples from different suppliers and with different characteristics:
# 3-methylthio-1-propanol and 3-methylthio-1-propanol.B -> Sigma-Aldrich, Cat. Nr. 318396, 98%
# 3-methylthio-1-propanol.A ->Sigma-Aldrich, Cat. Nr. W341509, ≥98% synthetic, FG grade.
# propyl.acetate" and propyl.acetate.B -> Sigma-Aldrich, Cat. Nr. 133108, 99%
# propyl.acetate.A -> Fluka, Cat. Nr. 40858, analytical standard grade

odor.labels <- c(
  "dmso" = "DMSO",
  "3-octanol" = "3-Octanol",
  "2-heptanone" = "2-Heptanone",
  "isobutyl.acetate" =  "Isobutyl Acetate",
  "hexyl.acetate" = "Hexyl Acetate",
  "3-methylthio-1-propanol" = "3-Methylthio-1-propanol",
  "3-methylthio-1-propanol.B" = "3-Methylthio-1-propanol",
  "3-methylthio-1-propanol.A" = "3-Methylthio-1-propanol",
  "propyl.acetate" = "Propyl Acetate",
  "propyl.acetate.B" = "Propyl Acetate",
  "propyl.acetate.A" = "Propyl Acetate",
  "methyl.hexanoate" = "Methyl Hexanoate",
  "butyl.acetate" = "Butyl Acetate",
  "pentyl.acetate" = "Pentyl Acetate"
);

# Color blind palette
odor.palette <- c(
  "dmso" = "#000000", # black
  "3-octanol" = "#994F00", # brown
  "2-heptanone" = "#E1BE6A",# light brown
  "isobutyl.acetate" =  "#56B4E9", # blue
  "hexyl.acetate" = "#56B4E9", # cyan
  "3-methylthio-1-propanol" = "#009E73", # green
  "3-methylthio-1-propanol.B" = "#009E73", # green
  "3-methylthio-1-propanol.A" = "#009E73", # green
  "propyl.acetate" = "#FFD700", # ocra
  "propyl.acetate.B" = "#FFD700", # ocra
  "propyl.acetate.A" = "#FFD700", # ocra
  "methyl.hexanoate" = "#ffa500", # orange
  "butyl.acetate" = "#D55E00", # vermillion
  "pentyl.acetate" = "#CC79A7" # reddish purple
);

specificity.data.hOr47a <- import_data("Input/ResponseProfile/hOr47a"); # Load data files as single list
specificity.metadata.hOr47a <- import_metadata("Metadata/ResponseProfile/hOr47a"); # Load metadata as single list

# Filter out non responding and spurious ROIs
specificity.filtered.hOr47a <- lapply(specificity.data.hOr47a, function(x) lapply(x, data_filtering));

specificity.ts.hOr47a <- mapply(function(x,y) ts_resp.calc(x, y, "odor", specificity.levels), specificity.filtered.hOr47a, specificity.metadata.hOr47a);
specificity.df <- specificity.df_mkr(specificity.ts.hOr47a, specificity.levels);

specificity.ts.plot <- ggplot(data = specificity.df, aes(x = time.sec, y = odor.mean)) + # Feeds the data and the variables
  geom_line(aes(colour= odor.name)) + # The line displays the mean values
  scale_colour_manual(values = odor.palette) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax=(odor.mean + odor.sd), ymin=(odor.mean - odor.sd), fill = odor.name), alpha = 0.4) +
  scale_fill_manual(values = odor.palette) +
  facet_wrap( ~ odor.name, ncol=2, labeller = labeller(odor.name = odor.labels)) + # Grouping according to the stimulus concentrations (groups)
  labs(x = "Time (s)", y ="R (340/380 nm)") +
  #  scale_fill_gradient(name = "Concentration", trans = "log") + # Sets a log scale for the coloring of the different treatments
  scale_x_continuous(breaks=seq(0,200,50)) +
  #  scale_y_continuous(breaks=seq(0,1.75,0.25)) +
  coord_cartesian(x = c(0, 200)) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
specificity.ts.plot;
ggsave("Output/specificity.hOr47a.pdf", plot = specificity.ts.plot, device = "pdf",  width = 80, height = 100, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

# specificity.shapiro <- by(specificity.df$pr.vals, VUAA1.pr.stat$pr.treat, function(x) shapiro.test(x));

rows.to.delete.specificity.A <- which(specificity.df[,"odor.name"] == "3-methylthio-1-propanol" |
                                      specificity.df[,"odor.name"] == "3-methylthio-1-propanol.A" |
                                      specificity.df[,"odor.name"] == "propyl.acetate" |
                                      specificity.df[,"odor.name"] == "propyl.acetate.A");
specificity.ts.plot.b <-specificity.df[-rows.to.delete.specificity.A,];

specificity.ts.plot.profile <- ggplot(data = specificity.ts.plot.b, aes(x = time.sec, y = odor.mean)) + # Feeds the data and the variables
  geom_line(aes(colour= odor.name)) + # The line displays the mean values
  scale_colour_manual(values = odor.palette) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax=(odor.mean + odor.sd), ymin=(odor.mean - odor.sd), fill = odor.name), alpha = 0.4) +
  scale_fill_manual(values = odor.palette) +
  facet_wrap( ~ odor.name, nrow = 2, labeller = labeller(odor.name = odor.labels)) + # Grouping according to the stimulus concentrations (groups)
  labs(x = "Time (s)", y = "R (340/380 nm)") +
  #  scale_fill_gradient(name = "Concentration", trans = "log") + # Sets a log scale for the coloring of the different treatments
  scale_x_continuous(breaks=seq(0,200,50)) +
  #  scale_y_continuous(breaks=seq(0,1.75,0.25)) +
  coord_cartesian(x = c(0, 200)) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
specificity.ts.plot.profile;
ggsave("Output/specificity.hOr47a.pdf", plot = specificity.ts.plot.profile, device = "pdf",  width = 180, height = 80, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

# Rescaling of responses
specificity.norm <- mapply(function(x, y) norm_resp(x, y, "odor", specificity.levels), specificity.filtered.hOr47a, specificity.metadata.hOr47a);
specificity.norm.df <- rbindlist(specificity.norm, fill=TRUE);

rows.to.delete.ts.plot <- which(specificity.norm.df[,"factor.name"] == "3-methylthio-1-propanol" |
                                  specificity.norm.df[,"factor.name"] == "3-methylthio-1-propanol.A" |
                                  specificity.norm.df[,"factor.name"] == "propyl.acetate" |
                                  specificity.norm.df[,"factor.name"] == "propyl.acetate.A");
specificity.norm.df.b <- specificity.norm.df[-rows.to.delete.ts.plot,];

# Barlett's test for homogeneity of variances
bartlett.test(mean.norm ~ factor.name, specificity.norm.df.b);
# data:  mean.norm by factor.name
# Bartlett's K-squared = 39.792, df = 9, p-value = 8.288e-06

specificity.shapiro <- by(specificity.norm.df.b$mean.norm, specificity.norm.df.b$factor.name, function(x) shapiro.test(x));
specificity.shapiro[["dmso"]]; # W = 0.76707, p-value = 0.008545
# Non parametric tests are needed as the dmso ctrl is not normally distributed!

kruskal.test(mean.norm ~ factor.name, specificity.norm.df.b);
# Kruskal-Wallis rank sum test
# data:  mean.norm by factor.name
# Kruskal-Wallis chi-squared = 50.148, df = 9, p-value = 1.01e-07

DunnettTest(specificity.norm.df.b$mean.norm, specificity.norm.df.b$factor.name, control = "dmso", conf.level = 0.95);
# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# $dmso
#                                     diff     lwr.ci   upr.ci    pval    
# 3-octanol-dmso                  8.212857  0.3622065 16.063507  0.0358 *  
# 2-heptanone-dmso               12.713744  4.8630934 20.564394  0.0003 ***
# isobutyl.acetate-dmso           8.596880  1.3099773 15.883782  0.0129 *  
# hexyl.acetate-dmso             10.260871  2.9739684 17.547774  0.0019 ** 
# 3-methylthio-1-propanol.B-dmso  1.014555 -6.8360955  8.865205  1.0000    
# propyl.acetate.B-dmso           8.207194  0.3565436 16.057844  0.0360 *  
# methyl.hexanoate-dmso          20.120323 12.8334200 27.407225 2.7e-09 ***
# butyl.acetate-dmso             18.952557 11.6656548 26.239460 5.6e-08 ***
# pentyl.acetate-dmso            34.658441 28.6558223 40.661059 < 2e-16 ***
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

specificity.odor.plot <- ggplot(specificity.norm.df.b, aes(x = factor.name, y = mean.norm)) + # Feeds the data and the variables
  geom_jitter(data = specificity.norm.df.b, mapping = aes(x = factor.name, y = mean.norm), color = "808080", size = 1, inherit.aes = FALSE) +
  geom_boxplot(aes(fill = factor(factor.name)), outlier.shape = NA, alpha = 0.4) +
  scale_fill_manual(values = odor.palette) +
  scale_x_discrete(labels=odor.labels) +
  theme(legend.title = element_blank(), legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
  labs(y = expression(paste(rho," (%)"))) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1));
specificity.odor.plot;
ggsave("Output/specificity.odor.pdf", plot = specificity.odor.plot, device = "pdf",  width = 50, height = 50, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

odor.labels.MP <- c( # 3-methylthio-1-propanol
  "dmso" = "DMSO",
  "3-methylthio-1-propanol.A" = "98%",
  "3-methylthio-1-propanol.B" = "\u2265 98% (FG)"
);

odor.labels.PA <- c(
  "dmso" = "DMSO",
  "propyl.acetate.A" = "99%",
  "propyl.acetate.B" = "\u2265 99% (AG)"
);

rows.to.keep.purity.MP <- which(specificity.norm.df[,factor.name] == "dmso" |
                          specificity.norm.df[,factor.name] == "3-methylthio-1-propanol.A" |
                          specificity.norm.df[,factor.name] == "3-methylthio-1-propanol.B");
specificity.norm.df.purity.MP <- specificity.norm.df[rows.to.keep.purity.MP,];
specificity.norm.df.purity.MP$factor.name <- factor(specificity.norm.df.purity.MP$factor.name, levels = c("dmso",
                                                                                                          "3-methylthio-1-propanol.A",
                                                                                                          "3-methylthio-1-propanol.B"));
kruskal.test(mean.norm ~ factor.name, specificity.norm.df.purity.MP);
# Kruskal-Wallis rank sum test
# data:  mean.norm by factor.name
# Kruskal-Wallis chi-squared = 8.6808, df = 2, p-value = 0.01303

DunnettTest(specificity.norm.df.purity.MP$mean.norm, specificity.norm.df.purity.MP$factor.name, control = "dmso", conf.level = 0.95);
# Dunnett's test for comparing several treatments with a control :  
# 95% family-wise confidence level
# $dmso
#                                     diff    lwr.ci    upr.ci    pval    
# 3-methylthio-1-propanol.A-dmso 12.343902  5.891584 18.796220 0.00061 ***
# 3-methylthio-1-propanol.B-dmso  1.014555 -5.437764  7.466873 0.90431   
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

rows.to.keep.purity.PA <- which(specificity.norm.df[,factor.name] == "dmso"|
                                  specificity.norm.df[,factor.name] == "propyl.acetate.A" |
                                  specificity.norm.df[,factor.name] == "propyl.acetate.B");
specificity.norm.df.purity.PA <- specificity.norm.df[rows.to.keep.purity.PA,];
specificity.norm.df.purity.PA$factor.name <- factor(specificity.norm.df.purity.PA$factor.name, levels = c("dmso",
                                                                                                          "propyl.acetate.A",
                                                                                                          "propyl.acetate.B"));
kruskal.test(mean.norm ~ factor.name, specificity.norm.df.purity.PA);
# Kruskal-Wallis rank sum test
# data:  mean.norm by factor.name
# Kruskal-Wallis chi-squared = 10.088, df = 2, p-value = 0.006447

DunnettTest(specificity.norm.df.purity.PA$mean.norm, specificity.norm.df.purity.PA$factor.name, control = "dmso", conf.level = 0.95);
# Dunnett's test for comparing several treatments with a control :  
#     95% family-wise confidence level
# $dmso
#                           diff   lwr.ci   upr.ci   pval    
# propyl.acetate.A-dmso 8.206566 1.660719 14.75241 0.0147 *  
# propyl.acetate.B-dmso 8.207194 1.661347 14.75304 0.0147 *  
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

specificity.odor.plot.purity.MP <- ggplot(specificity.norm.df.purity.MP, aes(x = factor.name, y = mean.norm)) + # Feeds the data and the variables
  geom_jitter(data = specificity.norm.df.purity.MP, mapping = aes(x = factor.name, y = mean.norm), color = "808080", size = 1, inherit.aes = FALSE) +
  geom_boxplot(aes(fill = factor(factor.name)), outlier.shape = NA, alpha = 0.4) +
  scale_fill_manual(values = odor.palette) +
  scale_x_discrete(labels=odor.labels.MP) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  coord_cartesian(y = c(-15, 20)) +
  theme(legend.title = element_blank(), legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
  labs(y = expression(paste(rho," (%)"))) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1));
specificity.odor.plot.purity.MP;
ggsave("Output/specificity.odor.MP.pdf", plot = specificity.odor.plot.purity.MP, device = "pdf",  width = 25, height = 40, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

specificity.odor.plot.purity.PA <- ggplot(specificity.norm.df.purity.PA, aes(x = factor.name, y = mean.norm)) + # Feeds the data and the variables
  geom_jitter(data = specificity.norm.df.purity.PA, mapping = aes(x = factor.name, y = mean.norm), color = "808080", size = 1, inherit.aes = FALSE) +
  geom_boxplot(aes(fill = factor(factor.name)), outlier.shape = NA, alpha = 0.4) +
  scale_fill_manual(values = odor.palette) +
  scale_x_discrete(labels=odor.labels.PA) +
  scale_y_continuous(breaks=seq(-15,20,5)) +
  coord_cartesian(y = c(-15, 20)) +
  theme(legend.title = element_blank(), legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
  labs(y = expression(paste(rho," (%)"))) +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme_Publication() +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1));
specificity.odor.plot.purity.PA;
ggsave("Output/specificity.odor.PA.pdf", plot = specificity.odor.plot.purity.PA, device = "pdf",  width = 25, height = 40, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)

###------------------------------------
###--- Untransfected cells analysis ---
###------------------------------------

list.of.data.untr <- import_data("Input/ResponseProfile/Untransfected"); # Load data files as single list
list.of.metadata.untr <- import_metadata("Metadata/ResponseProfile/Untransfected"); # Load metadata as single list

untr.levels <- c("dmso", "3-octanol", "2-heptanone", "isobutyl.acetate", "hexyl.acetate", "3-methylthio-1-propanol",
                  "3-methylthio-1-propanol.B", "3-methylthio-1-propanol.A", "propyl.acetate", "propyl.acetate.A",
                  "propyl.acetate.B", "methyl.hexanoate", "butyl.acetate", "pentyl.acetate", "geosmin");


untr.labels <- c(
  "dmso" = "DMSO",
  "3-octanol" = "3-Octanol",
  "2-heptanone" = "2-Heptanone",
  "isobutyl.acetate" =  "Isobutyl Acetate",
  "hexyl.acetate" = "Hexyl Acetate",
  "3-methylthio-1-propanol" = "3-Methylthio-1-propanol",
  "3-methylthio-1-propanol.B" = "3-Methylthio-1-propanol",
  "3-methylthio-1-propanol.A" = "3-Methylthio-1-propanol",
  "propyl.acetate" = "Propyl Acetate",
  "propyl.acetate.B" = "Propyl Acetate",
  "propyl.acetate.A" = "Propyl Acetate",
  "methyl.hexanoate" = "Methyl Hexanoate",
  "butyl.acetate" = "Butyl Acetate",
  "pentyl.acetate" = "Pentyl Acetate",
  "geosmin" = "Geosmin"
);

untr.palette <- c(
  "dmso" = "#000000", # black
  "3-octanol" = "#994F00", # brown
  "2-heptanone" = "#E1BE6A",# light brown
  "isobutyl.acetate" =  "#56B4E9", # blue
  "hexyl.acetate" = "#56B4E9", # cyan
  "3-methylthio-1-propanol" = "#009E73", # green
  "3-methylthio-1-propanol.B" = "#009E73", # green
  "3-methylthio-1-propanol.A" = "#009E73", # green
  "propyl.acetate" = "#FFD700", # ocra
  "propyl.acetate.B" = "#FFD700", # ocra
  "propyl.acetate.A" = "#FFD700", # ocra
  "methyl.hexanoate" = "#ffa500", # orange
  "butyl.acetate" = "#D55E00", # vermillion
  "pentyl.acetate" = "#CC79A7", # reddish purple
  "geosmin" = "#FF00FF" # magenta
);

# Filter out non responding and spurious ROIs
specificity.filtered.untr <- lapply(list.of.data.untr, function(x) lapply(x, data_filtering.untr));

specificity.ts.untr <- mapply(function(x,y) ts_resp.calc(x, y, "odor", untr.levels), specificity.filtered.untr, list.of.metadata.untr);
specificity.df.untr <- specificity.df_mkr(specificity.ts.untr, untr.levels);

specificity.ts.plot.untr <- ggplot(data = specificity.df.untr, aes(x = time.sec, y = odor.mean)) + # Feeds the data and the variables
  geom_line(aes(colour= odor.name)) + # The line displays the mean values
  scale_colour_manual(values = untr.palette) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax=(odor.mean + odor.sd), ymin=(odor.mean - odor.sd), fill = odor.name), alpha = 0.4) +
  scale_fill_manual(values = untr.palette) +
  facet_wrap( ~ odor.name, ncol=2, labeller = labeller(odor.name = untr.labels)) + # Grouping according to the stimulus concentrations (groups)
  labs(x = "Time (s)", y =  "R (340/380 nm)") +
  #  scale_fill_gradient(name = "Concentration", trans = "log") + # Sets a log scale for the coloring of the different treatments
  scale_x_continuous(breaks=seq(0,200,50)) +
  #  scale_y_continuous(breaks=seq(0,1.75,0.25)) +
  coord_cartesian(x = c(0, 200), y = c(-10,160)) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
specificity.ts.plot.untr;
ggsave("Output/specificity.df.untr.pdf", plot = specificity.ts.plot.untr, device = "pdf",  width = 80, height = 100, units = "mm", useDingbats=FALSE); # Saves the plot

untr.resp <- mapply(function(x, y) norm_resp.untr(x, y, "odor", untr.levels), specificity.filtered.untr, list.of.metadata.untr);
untr.resp.df <- rbindlist(untr.resp, fill=TRUE);

untr.resp.df.plot <- which(untr.resp.df[,"factor.name"] == "3-methylthio-1-propanol" |
                           untr.resp.df[,"factor.name"] == "propyl.acetate");

untr.labels <- c(
  "dmso" = "DMSO",
  "3-octanol" = "3-Octanol",
  "2-heptanone" = "2-Heptanone",
  "isobutyl.acetate" =  "Isobutyl Acetate",
  "hexyl.acetate" = "Hexyl Acetate",
  "3-methylthio-1-propanol" = "3-Methylthio-1-propanol",
  "3-methylthio-1-propanol.B" = "3-Methylthio-1-propanol",
  "3-methylthio-1-propanol.A" = "3-Methylthio-1-propanol",
  "propyl.acetate" = "Propyl Acetate",
  "propyl.acetate.B" = "Propyl Acetate",
  "propyl.acetate.A" = "Propyl Acetate",
  "methyl.hexanoate" = "Methyl Hexanoate",
  "butyl.acetate" = "Butyl Acetate",
  "pentyl.acetate" = "Pentyl Acetate",
  "geosmin" = "Geosmin");

untr.resp.df.b <- untr.resp.df[-untr.resp.df.plot,];
untr.resp.df.dunn <- untr.resp.df.b[-76,]; # removing an outlier for geosmin

kruskal.test(mean.norm ~ factor.name, untr.resp.df.dunn);
# Kruskal-Wallis rank sum test
# data:  mean.norm by factor.name
# Kruskal-Wallis chi-squared = 27.115, df = 12, p-value = 0.007438

DunnettTest(untr.resp.df.dunn$mean.norm, untr.resp.df.dunn$factor.name, control = "dmso", conf.level = 0.95);
# Dunnett's test for comparing several treatments with a control :  
# 95% family-wise confidence level
# $dmso
#                                        diff       lwr.ci       upr.ci    pval    
# 3-octanol-dmso                 -0.020961250 -0.034814515 -0.007107986 0.00064 ***
# 2-heptanone-dmso               -0.005060096 -0.015532180  0.005411987 0.73434    
# isobutyl.acetate-dmso          -0.002414234 -0.012440491  0.007612022 0.99648    
# hexyl.acetate-dmso             -0.003466338 -0.013492595  0.006559918 0.94898    
# 3-methylthio-1-propanol.B-dmso  0.004914167 -0.005557916  0.015386251 0.76246    
# 3-methylthio-1-propanol.A-dmso  0.005862474 -0.004609609  0.016334557 0.57241    
# propyl.acetate.A-dmso           0.003889076 -0.006583008  0.014361159 0.92221    
# propyl.acetate.B-dmso           0.001451739 -0.009020344  0.011923823 0.99998    
# methyl.hexanoate-dmso          -0.001290304 -0.011316561  0.008735952 0.99999    
# butyl.acetate-dmso             -0.001716040 -0.011742296  0.008310217 0.99985    
# pentyl.acetate-dmso            -0.002963825 -0.011894447  0.005966798 0.96058    
# geosmin-dmso                    0.005285769 -0.005821553  0.016393091 0.74924    

specificity.untr.plot <- ggplot(untr.resp.df.b, aes(x = factor.name, y = mean.norm)) + # Feeds the data and the variables
  geom_jitter(data = untr.resp.df.b, mapping = aes(x = factor.name, y = mean.norm), color = "808080", size = 1, inherit.aes = FALSE) +
  geom_boxplot(aes(fill = factor(factor.name)), outlier.shape = NA, alpha = 0.4) +
  scale_fill_manual(values = untr.palette) +
  scale_x_discrete(labels=untr.labels) +
  scale_y_continuous(breaks=seq(-1,1,0.1)) +
  #  scale_y_continuous(breaks=seq(0,1.75,0.25)) +
  coord_cartesian(y = c(-0.1,0.2)) +
  theme(legend.title = element_blank(), legend.position="bottom", plot.title = element_text(hjust = 0.5)) +
  labs(y = expression(paste(rho," (%)"))) +
  theme_Publication() +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none") +
  theme(axis.title.x=element_blank()) +
  theme(axis.text.x = element_text(angle=45, vjust=1, hjust = 1));
specificity.untr.plot;
ggsave("Output/specificity.untr.pdf", plot = specificity.untr.plot, device = "pdf",  width = 60, height = 50, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off(); # To close the graphic device (plot)