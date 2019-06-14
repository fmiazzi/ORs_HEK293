### R and OS version used --------------
# >sessionInfo()
# R version 3.5.2 (2018-12-20)
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

### Set Working directory --------------
my_wd <- "/Users/basilio8/Documents/Postdoc/InsectORexpression/TagsTests/";
setwd(my_wd);

### Package loading --------------
library(ggplot2); # for graphs
library(gridExtra); # For custom grids
library(ggthemes); # Extra graph themes for ggplot
library(grid);
library(scales);
library(multcomp); # For Dunnet's post-hoc
library(extrafont);
library(stats); # for ecdf computing

### Functions --------------
# Loads as an element of a list each .csv file within the folder 'path' in the working directory
load_data <- function(path) { 
  files <- dir(path, pattern = '*.csv', full.names = TRUE); # Gets the file paths
  names <- dir(path, pattern = '*.csv', full.names = FALSE); # Gets the file names
  names.noext <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(names)); # Removes the extension from the file names
  output <- lapply(files, read.csv); # Imports all the .csv files as elements of a list
  names(output) <- names.noext; # Names each element using the file name w/out extension
  output; # Gives the renamed list as output
}

# Excludes ROIs with high Ca2+ levels and large st.dev on baselevels
polish_rois <- function(DataFrame) {
  # Eliminates ROIs with a basal [Ca2+]i > 200
  highCa.polish <- apply(DataFrame,2, function(x) mean(x[1:9])) > 200;
  highCa.polish[which(is.na(highCa.polish))] <- FALSE;
  output <- DataFrame[!highCa.polish];
  # Eliminates ROIs with a high st.dev before odor application (at 10th frame)
  highSDodor.polish <- apply(output,2, function(x) sd(x[1:9])) > 10;
  highSDodor.polish[which(is.na(highSDodor.polish))] <- FALSE;
  output <- output[!highSDodor.polish];
  # Eliminates ROIs with a high st.dev before VUAA1 application (at 70th frame)
  highSDvuaa1.polish <- apply(output,2, function(x) sd(x[61:69])) > 10;
  highSDvuaa1.polish[which(is.na(highSDvuaa1.polish))] <- FALSE;
  output <- output[!highSDvuaa1.polish];
}

# Calculates the mean for each row of a data.frame (the list element)
mean_calc <- function(DataFrame) {
  output <- apply(DataFrame, 1, function(x) mean(x, na.rm = TRUE));
}

# Calculates the standard deviation for each row of a data.frame (the list element)
sd_calc <- function(DataFrame) {
  output <- apply(DataFrame, 1, function(x) sd(x, na.rm = TRUE));
}

# Subtracts the background before the VUAA1 stimulation
VUAA1.bkg_calc <- function(DataFrame) {
  # Calculates the mean between the 61nd and the 68th row for each line of the data.frame
  bkg.vals <- apply(DataFrame[61:68,], 2, function(x) mean(x, na.rm = TRUE));
  # creates a data.frame with the same size of the input one filled with the bkg value for each column
  bkg.df <- data.frame(matrix(rep(bkg.vals, each = nrow(DataFrame)), nrow = nrow(DataFrame)));
  # subtracts the background from the input data.frame
  output <- DataFrame - bkg.df;
}

# Calculates and subtracts the background value for the data.frame
bkg_sub <- function(DataFrame) {
  # Calculates the mean between the 2nd and the 8th row for each line of the data.frame
  bkg.vals <- apply(DataFrame[2:8,], 2, function(x) mean(x, na.rm = TRUE));
  # creates a data.frame with the same size of the input one filled with the bkg value for each column
  bkg.df <- data.frame(matrix(rep(bkg.vals, each = nrow(DataFrame)), nrow = nrow(DataFrame)));
  # subtracts the background from the input data.frame
  output <- DataFrame - bkg.df;
}

# Creates the data.frame for the time series plot
time.series_mkr <- function(DataList) {
  # subtracts the background
  bkg.sub <- lapply(DataList, bkg_sub);
  # calculates the mean for each cuvette experiment and arranges all the values in a vector
  ts.mean <- unlist(lapply(bkg.sub, mean_calc), use.names = FALSE);
  # same as above but for the standard deviation
  ts.sd <- unlist(lapply(bkg.sub, sd_calc), use.names = FALSE);
  # calculates the time in seconds for each entry 
  ts.times <- rep(seq(0, (nrow(DataList[[1]])-1)*5, 5), times = length(names(DataList)));
  # gets the treatment name for each entry
  ts.treat <- rep(names(DataList), each = nrow(DataList[[1]]));
  ts.lab <- rep("generic.label", times = length(ts.treat));
  # creates the output data.frame
  output <- data.frame(ts.mean, ts.sd, ts.times, ts.treat, ts.lab);
  # arranges the output data.frame
  output$ts.treat <- factor(output$ts.treat, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
  output;
}

# Creates the data.frame for the statistical analysis of response intensities
resp.int_eval <- function(DataList) {
  # subtracts the background
  bkg.sub <- lapply(DataList, bkg_sub);
  # takes the values at frame=21 (100 seconds) for pentyl acetate plateau levels
  pentac.vals <- unlist(lapply(bkg.sub, function(x) x[21,]), use.names = FALSE);
  # length of the time.series
  pentac.num <- as.vector(unlist(lapply(DataList, ncol), use.names = FALSE));
  # assign the agonist for each measurement
  pentac.stim <- rep("pent.ac", length(pentac.vals));
  # assigns the treatment of each measuremnt
  pentac.treat <- rep(names(DataList), times = pentac.num);
  # creates the data.frame for the pentyl acetate values
  pentac.output <- data.frame(pentac.vals, pentac.stim, pentac.treat);
  # subtracts the baselevel before the VUAA1 stimulation
  VUAA1.bkg <- lapply(DataList, VUAA1.bkg_calc);
  # takes the values at frame=76 (375 seconds) for VUAA1 maximal levels
  VUAA1.vals <- unlist(lapply(VUAA1.bkg, function(x) x[76,]), use.names = FALSE);
  # length of the time.series
  VUAA1.num <- as.vector(unlist(lapply(DataList, ncol), use.names = FALSE));
  # assign the agonist for each measurement
  VUAA1.stim <- rep("VUAA1", length(VUAA1.vals));
  # assigns the treatment of each measuremnt
  VUAA1.treat <- rep(names(DataList), times = VUAA1.num);
  # creates the data.frame for the VUAA1 values
  VUAA1.output <- data.frame(VUAA1.vals, VUAA1.stim, VUAA1.treat);
  # renames the column names
  colnames(pentac.output) <- c("resp.vals", "resp.stim", "resp.treat");
  colnames(VUAA1.output) <- c("resp.vals", "resp.stim", "resp.treat");
  # merges the two data.frames
  output <- rbind(pentac.output, VUAA1.output);
  # re-orders the treatment level
  output$resp.treat <- factor(output$resp.treat, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
  output;
}

# Creates the data.frame for the barplot of the response intensities
resp.plot_eval <- function(DataList) {
  # Same procedure as 'resp.int_eval' but to obtain the mean and the st.dev of values for each treatment
  bkg.sub <- lapply(DataList, bkg_sub);
  pentac.means <- unlist(lapply(bkg.sub, function(x) mean(as.double(x[21,]))), use.names = FALSE);
  pentac.sd <- unlist(lapply(bkg.sub, function(x) sd(as.double(x[21,]))), use.names = FALSE);
  pentac.stim <- rep("pent.ac", length(pentac.means));
  pentac.treat <- names(DataList);
  pentac.output <- data.frame(pentac.means, pentac.sd, pentac.stim, pentac.treat);
  # Same procedure as 'resp.int_eval' but to obtain the mean and the st.dev of values for each treatment
  VUAA1.bkg <- lapply(DataList, VUAA1.bkg_calc);
  VUAA1.means <- unlist(lapply(VUAA1.bkg, function(x) mean(as.double(x[76,]))), use.names = FALSE);
  VUAA1.sd <- unlist(lapply(VUAA1.bkg, function(x) sd(as.double(x[76,]))), use.names = FALSE);
  VUAA1.stim <- rep("VUAA1", length(VUAA1.means));
  VUAA1.treat <- names(DataList);
  VUAA1.output <- data.frame(VUAA1.means, VUAA1.sd, VUAA1.stim, VUAA1.treat);
  # renames the column names
  colnames(pentac.output) <- c("resp.mean", "resp.sd", "resp.stim", "resp.treat");
  colnames(VUAA1.output) <- c("resp.mean", "resp.sd", "resp.stim", "resp.treat");
  # merges the two data.frames
  output <- rbind(pentac.output, VUAA1.output);
  # re-orders the treatment level (determines the order of appearance in the barplot)
  output$resp.treat <- factor(output$resp.treat, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
  output;
}

# Publication-ready theme for ggplot (Modified from Koundinya Desiraju)
# https://rpubs.com/Koundy/71792
theme_Publication <- function(base_size=7, base_family="Helvetica") {
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

# Assigns labels to each facet of a plot
fac_labeller <- function(variable){
  return(fac_treat[value])
}

# Subtracts the background from each ROI measurement after pentyl acetate stimulation
pentac_roi.resp <- function(DataFrame) {
  # Same procedure as 'resp.int_eval' but on single ROIs
  pentac.bkg <- unlist(apply(DataFrame, 2, function(x) mean(x[2:8])), use.names = FALSE);
  pentac.max <- unlist(DataFrame[21,], use.names = FALSE);
  pentac.delta <- pentac.max - pentac.bkg
}

# Creates the cumulative distribution of responding ROIs after stimulation with pentyl acetate
pentac_distr.mkr <- function(DataList, MetaList) {
  # ROI pentyl acetate plateau responses
  pentac.vals <- unlist(DataList, use.names = FALSE);
  # length of each response vector
  pentac.num <- unlist(lapply(DataList, length), use.names = FALSE);
  # treatment value for each measuremnt
  pentac.treat <- rep(MetaList$exp.treat, times = pentac.num);
  # transfection cuvette number for each measurement (same cuvette -> n = 1)
  pentac.cuv <- rep(MetaList$exp.cuvette, times = pentac.num);
  # merge as columns of a data.frame
  pentac.df <- data.frame(pentac.vals, pentac.treat, pentac.cuv);
  # remove NA values
  pentac.df[!is.na(pentac.df$pentac.vals),];
  # calculates the histogram for each cuvette; breaks determine the number and spacing of bins
  pentac.hist <- by(pentac.df$pentac.vals, pentac.df$pentac.cuv, function(x) hist(as.double(x), breaks = seq(-100, 450, by = 0.5)));
  # calculates the cumulative probability for ach bin; divided by 10 so that -> pmax = 1 (due to the chosen bin interval!)
  pentac.cum.prob <- lapply(pentac.hist, function(x) cumsum((x$density)/2));
  # binds column-wise in one data.frame all cuvettes belonging to the same treatment
  levs <- by(MetaList$exp.cuvette, MetaList$exp.treat, function(x) as.character(unique(x)));
  output <- lapply(levs, function(x) cbind.data.frame(pentac.cum.prob[x]));
  # adds the mid values of each bin
  xvals <- unlist(pentac.hist[[1]]["mids"], use.names = FALSE);
  output[["xvals"]] <- xvals;
  output;
}

# Prepares list for data analysis and plotting of cumulative distributions and percentage of responding cells
pentac_csplot.mkr <- function(CSdataList) {
  # Calculates the percentage of responding cells
  CSdata <- CSdataList[-length(CSdataList)]; # removes the bin mid-values column
  # Determination of the threshold value to consider a ROI as 'responding' ->
  CS.length <- unlist(lapply(CSdata, ncol), use.names = FALSE);
  CS.means <- lapply(CSdata, mean_calc);
  CS.sd <- lapply(CSdata, sd_calc);
  ctrl.length <- length(CS.means[["ctrl"]]);
  # Threshold value: the smallest bin where the probability of the (MEAN -2*ST.DEV) CTRL response value in higher that 99.5%.
  # this is a quite conservative definition
  resp.xthr <- CSdataList[["xvals"]][min(which((CS.means[["ctrl"]] - 2*CS.sd[["ctrl"]]) > 0.995))];
  # (A) data.frame for statistics on percentace of cells responding to pentyl acetate
  # gets the cum.prob as % of responding cells at the threshold level
  pr.vals <- round((1- unlist(lapply(CSdata, function(x) x[which(CSdataList[["xvals"]]==resp.xthr),]), use.names = FALSE))*100, 2);
  # assigns the treatment value to each measurement
  pr.treat <- rep(names(CSdata), times = CS.length);
  # assigns the agonist
  pr.agonist <- rep("pent.ac", length = length(pr.treat));
  # output data.frame for statistical analysis
  pentac.perc.resp <- data.frame(pr.vals, pr.agonist, pr.treat);
  # (B) data.frame for the barplot on the percentace of cells responding to pentyl acetate
  # gets the mean of the % of responding cells for each treatment
  pr.means <- unlist(by(pentac.perc.resp$pr.vals, pentac.perc.resp$pr.treat, function(x) round(mean(x), 2), simplify = FALSE), use.names = FALSE);
  # gets the st.dev of the % of responding cells for each treatment
  pr.sd <- unlist(by(pentac.perc.resp$pr.vals, pentac.perc.resp$pr.treat, function(x) round(sd(x), 2), simplify = FALSE), use.names = FALSE);
  # name of the treatment for each measurement
  pr.levels <- names(CSdata);
  # agonist for each measurement
  pr.agonist <- rep("pent.ac", times = length(pr.levels));
  pentac.pr.barplot <- data.frame(pr.means, pr.sd, pr.levels, pr.agonist);
  pentac.pr.barplot$pr.levels <- factor(pentac.pr.barplot$pr.levels, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
  # output data.frame for barplot
  pentac.pr.barplot
  # (C) data.frame for the distribution of cells responding to pentyl acetate
  # Calculates the mean response for each cuvette (as percentage of responding cells)
  cd.means <- (1- unlist(lapply(CSdata, mean_calc), use.names = FALSE))*100;
  # Calculates the s.d for each cuvette
  cd.sd <- (unlist(lapply(CSdata, sd_calc), use.names = FALSE))*100;
  cd.treat <- rep(names(CSdata), each = ctrl.length);
  cd.xthr <- rep(resp.xthr, times = length(cd.means));
  cd.xvals <- CSdataList[length(CSdataList)];
  names(cd.xvals) <- "cd.xvals";
  # assigns the mean percentage of responding cells for that treatment
  cd.mean.label <- rep(pr.means, each = ctrl.length);
  # assigns the st.dev. percentage of responding cells for that treatment
  cd.sd.label <- rep(pr.sd, each = ctrl.length);
  # assigns the agonist
  cd.agonist <- rep("Pentyl Acetate", times = length(cd.means));
  # output data.frame for cumulative plots
  pentac.cs.distr <- data.frame(cd.xvals, cd.xthr, cd.means, cd.sd, cd.mean.label, cd.sd.label, cd.treat, cd.agonist);
  pentac.cs.distr$cd.treat <- factor(pentac.cs.distr$cd.treat, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
  # Creates the output list of data.frames
  output <- list();
  output[["pentac.perc.resp"]] <- pentac.perc.resp;
  output[["pentac.pr.barplot"]] <- pentac.pr.barplot;
  output[["pentac.cs.distr"]] <- pentac.cs.distr;
  output;
  # Here is the code to get the p=0.5 values:
  #ec50.length <- nrow(data.subsets[[1]]);
  #ec50.means <- rep(unlist(lapply(data.subsets, function(x) x[which.min(abs(x$cd.means - 0.5)), 1]), use.names = FALSE), each = ec50.length);
  #ec50.highs <- rep(unlist(lapply(data.subsets, function(x)
  #  x[which.min(abs((x$cd.means - x$cd.sd) - 0.5)), 1]) , use.names = FALSE), each = ec50.length);
  #ec50.lows <- rep(unlist(lapply(data.subsets, function(x)
  #  x[which.min(abs((x$cd.means + x$cd.sd) - 0.5)), 1]) , use.names = FALSE), each = ec50.length);
  #ec50.treats <- rep(unique(levels(CSdata.frame$cd.treat)), each = ec50.length);
  #ec50df <- data.frame(ec50.means, ec50.lows, ec50.highs, ec50.treats);
  #ec50df$ec50.treats <- factor(ec50df$ec50.treats, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
}

# Subtracts the background from each ROI measurement after VUAA1 stimulation
# equivalent for VUAA1 to 'pentac_roi.resp'
VUAA1_roi.resp <- function(DataFrame) {
  VUAA1.bkg <- unlist(apply(DataFrame, 2, function(x) mean(x[61:68])), use.names = FALSE);
  VUAA1.max <- unlist(DataFrame[76,], use.names = FALSE);
  VUAA1.delta <- VUAA1.max - VUAA1.bkg
}

# Creates the cumulative distribution of responding ROIs after stimulation with VUAA1
# equivalent for VUAA1 to 'pentac_distr.mkr'
VUAA1_distr.mkr <- function(DataList, MetaList) {
  VUAA1.vals <- unlist(DataList, use.names = FALSE);
  VUAA1.num <- unlist(lapply(DataList, length), use.names = FALSE);
  VUAA1.treat <- rep(MetaList$exp.treat, times = VUAA1.num);
  VUAA1.cuv <- rep(MetaList$exp.cuvette, times = VUAA1.num);
  VUAA1.df <- data.frame(VUAA1.vals, VUAA1.treat, VUAA1.cuv);
  VUAA1.df[!is.na(VUAA1.df$VUAA1.vals),];
  VUAA1.hist <- by(VUAA1.df$VUAA1.vals, VUAA1.df$VUAA1.cuv, function(x) hist(as.double(x), breaks = seq(-50, 650, by = 0.5)));
  VUAA1.cum.prob <- lapply(VUAA1.hist, function(x) cumsum((x$density)/2));
  levs <- by(MetaList$exp.cuvette, MetaList$exp.treat, function(x) as.character(unique(x)));
  output <- lapply(levs, function(x) cbind.data.frame(VUAA1.cum.prob[x]));
  xvals <- unlist(VUAA1.hist[[1]]["mids"], use.names = FALSE);
  output[["xvals"]] <- xvals;
  output;
}

VUAA1.cs.dls <- VUAA1_csplot.mkr(VUAA1.cs);

# Prepares list for data analysis and plotting of cumulative distributions and percentage of responding cells for VUAA1
# equivalent for VUAA1 to 'pentac_csplot.mkr'
VUAA1_csplot.mkr <- function(CSdataList) {
  # Gives the percentage of responding cells
  CSdata <- CSdataList[-length(CSdataList)];
  CS.length <- unlist(lapply(CSdata, ncol), use.names = FALSE);
  CS.means <- lapply(CSdata, mean_calc);
  CS.sd <- lapply(CSdata, sd_calc);
  ctrl.length <- length(CS.means[["ctrl"]]);
  resp.xthr <- CSdataList[["xvals"]][min(which((CS.means[["ctrl"]] - 2*CS.sd[["ctrl"]]) > 0.995))];
  # data.frame for statistics on percentace of cells responding to VUAA1
  pr.vals <- round((1- unlist(lapply(CSdata, function(x) x[which(CSdataList[["xvals"]]==resp.xthr),]), use.names = FALSE))*100, 2);
  pr.treat <- rep(names(CSdata), times = CS.length);
  pr.agonist <- rep("VUAA1", length = length(pr.treat));
  VUAA1.perc.resp <- data.frame(pr.vals, pr.agonist, pr.treat);
  # data.frame for the barplot on the percentace of cells responding to VUAA1
  pr.means <- unlist(by(VUAA1.perc.resp$pr.vals, VUAA1.perc.resp$pr.treat, function(x) round(mean(x), 2), simplify = FALSE), use.names = FALSE);
  pr.sd <- unlist(by(VUAA1.perc.resp$pr.vals, VUAA1.perc.resp$pr.treat, function(x) round(sd(x), 2), simplify = FALSE), use.names = FALSE);
  pr.levels <- names(CSdata);
  pr.agonist <- rep("VUAA1", times = length(pr.levels));
  VUAA1.pr.barplot <- data.frame(pr.means, pr.sd, pr.levels, pr.agonist);
  VUAA1.pr.barplot$pr.levels <- factor(VUAA1.pr.barplot$pr.levels, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
  VUAA1.pr.barplot
  #data.frame for the distribution of cells responding to VUAA1 (as percentage of responding cells)
  cd.means <- (1 - unlist(lapply(CSdata, mean_calc), use.names = FALSE))*100;
  cd.sd <- unlist(lapply(CSdata, sd_calc), use.names = FALSE)*100;
  cd.treat <- rep(names(CSdata), each = ctrl.length);
  cd.xthr <- rep(resp.xthr, times = length(cd.means));
  cd.xvals <- CSdataList[length(CSdataList)];
  names(cd.xvals) <- "cd.xvals";
  cd.agonist <- rep("VUAA1", times = length(cd.means));
  cd.mean.label <- rep(pr.means, each = ctrl.length);
  cd.sd.label <- rep(pr.sd, each = ctrl.length);
  VUAA1.cs.distr <- data.frame(cd.xvals, cd.xthr, cd.means, cd.sd, cd.mean.label, cd.sd.label, cd.treat, cd.agonist);
  VUAA1.cs.distr$cd.treat <- factor(VUAA1.cs.distr$cd.treat, levels = c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a", "pDmel.or47a"));
  output <- list();
  output[["VUAA1.perc.resp"]] <- VUAA1.perc.resp;
  output[["VUAA1.pr.barplot"]] <- VUAA1.pr.barplot;
  output[["VUAA1.cs.distr"]] <- VUAA1.cs.distr;
  output;
}

### Main Analysis -------
# Input data loading
Data.values <- load_data("Input-CaConc"); # Fiji ouput .csv files loaded as a single list
Data.meta <- read.csv("HEKcellsMetaData.csv"); # Metadata for the input files

input.vals <- lapply(Data.values, function(x) x[-1]); # Remove the first column from each list (index column)
polished.vals <- lapply(input.vals, polish_rois); # Eliminates ROIs with high Ca2+ basal levels or high st.dev
list.of.exp <- list();
list.of.means <- list();
# Organizes the files in a single 3-level list where all means from the same trasfection cuvette (exp.cuvette)
# are grouped in a list within the list of the treatment (exp.treat) they belong to.
for (i in unique(Data.meta$exp.treat)) {
  sub.treat <- polished.vals[Data.meta[Data.meta$exp.treat == i,1]];
  sub.meta <- Data.meta[Data.meta$exp.treat == i,];
  sub.exp <- list();
  for (j in unique(sub.meta$exp.cuvette)) {
    df.exp <- do.call(cbind, sub.treat[as.character(sub.meta[sub.meta$exp.cuvette == j,1])]);
    sub.exp <- c(sub.exp,list(df.exp));
  }
  rm(df.exp, j);
  sub.means <-  as.data.frame(lapply(sub.exp, mean_calc));
  names(sub.exp) <- as.character(unique(sub.meta$exp.cuvette));
  names(sub.means) <- as.character(unique(sub.meta$exp.cuvette));
  list.of.exp <- c(list.of.exp, list(sub.exp));
  list.of.means <- c(list.of.means, list(sub.means))
}
rm(i, sub.treat, sub.meta, sub.exp, sub.means);
names(list.of.exp) <- as.character(unique(Data.meta$exp.treat));
names(list.of.means) <- as.character(unique(Data.meta$exp.treat));

# Color blind palette
cbPalette <- c("#000000", # black
               "#E69F00", # ocra
               "#56B4E9", # cyan
               "#009E73", # green
               "#F0E442", # yellow
               "#0072B2", # blue
               "#D55E00", # vermillion
               "#CC79A7");# reddish purple

# Arrows positions for 'plot.ts'
arrow_pos <- data.frame(x.vals = c(45, 345),
                        x.end = c(45, 345),
                        y.vals = c(90, 90),
                        y.end = c(80, 80));

arrow_pos.2 <- data.frame(x.vals = c(45, 345),
                        x.end = c(45, 345),
                        y.vals = c(230, 230),
                        y.end = c(200, 200));

# Labels of applied stimuli for 'plot.ts'
stim_labels <- data.frame(lab.x = c(45, 345),
                          lab.y = c(90, 90),
                          lab.names = c("PA", "VUAA1"));

stim_labels.2 <- data.frame(lab.x = c(45, 345),
                          lab.y = c(235, 235),
                          lab.names = c("PA", "VUAA1"));

# Facets names for treatment faceting
fac_treat <- c(ctrl = "EV",
               orco = "Orco",
               or47a = "hOr47a + Orco",
               ers.or47a = "E.hOr47a + Orco", # ERS = endoplasmic reticulum release signal
               dyn.or47a = "M.E.hOr47a + Orco") # MTS = Microtubule-mediated transport signal

# Facets names for treatment faceting
fac_treat.2 <- c(dyn.or47a = "M.E.hOr47a + Orco", # MTS = Microtubule-mediated transport signal
               pDmel.or47a = "pDmel-hOr47a"); #pDmel = optimized vector including DmelOrco

fac_treat.tot <- c(ctrl = "EV",
               orco = "Orco",
               or47a = "hOr47a + Orco",
               ers.or47a = "E.hOr47a + Orco", # ERS = endoplasmic reticulum release signal
               dyn.or47a = "M.E.hOr47a + Orco", # MTS = Microtubule-mediated transport signal
               pDmel.or47a = "pDmel-hOr47a"); #pDmel = optimized vector including DmelOrco

# Facets names for stimulus faceting
fac_stim <- c(pent.ac = "Pentyl Acetate",
              VUAA1 = "VUAA1");

fac_stim.tiny <- c(pent.ac = "PA",
              VUAA1 = "VUAA1");

# x-axis labels for barplots
stim_xaxis <- c("EV", "Orco", "hOr47a+Orco", "E.hOr47a+Orco", "M.E.hOr47a+Orco", "pDmel-Or47a");
stim_xaxis.2 <- c("M.E.hOr47a+Orco", "pDmel-hOr47a");


# Get the data.frame for the means plot
plot.ts <- time.series_mkr(list.of.means[c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a")]);
plot.ts.2 <- time.series_mkr(list.of.means[c("dyn.or47a", "pDmel.or47a")]);

# Palette for time series experiments
ts.pal <- c(cbPalette[1], cbPalette[7], cbPalette[2], cbPalette[4], cbPalette[3]);
names(ts.pal) <- as.character(unique(plot.ts$ts.treat));

ts.pal.2 <- c(cbPalette[3], cbPalette[8]);
names(ts.pal.2) <- as.character(unique(plot.ts.2$ts.treat));

# Time series plot
means.plot.1 <- ggplot(plot.ts, aes(x = ts.times, y = ts.mean)) + # Feeds the data and the variables
  geom_line(aes(colour= ts.treat)) + # The line displays the mean values
  scale_colour_manual(values = ts.pal) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax= ts.mean + ts.sd, ymin = ts.mean - ts.sd, fill = ts.treat), alpha = 0.4) +
  scale_fill_manual(values = ts.pal) +
  geom_segment(data = arrow_pos, aes(x = x.vals, xend = x.end, y = y.vals, yend = y.end), size = 0.5, 
               arrow=arrow(length=unit(0.7, "mm"), type="closed")) +
  geom_text(data = stim_labels, aes(x = lab.x, y = lab.y, label = lab.names), size = 2.45, vjust = 0, nudge_y = 3) + # size = 2.45 because geom_text measures in mm and not in pts! pts = mm*0.35
  theme_Publication() +
  facet_grid(ts.lab ~ ts.treat, labeller = labeller(ts.treat = fac_treat)) + # Grouping according to the stimulus concentrations (groups)
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  labs(title = "",x = "Time (s)", y = expression(paste(Delta,"[Ca"^"2+"*"]"[i]*" (nM)"))) +
  scale_y_continuous(breaks=seq(-100,200,20)) +
  coord_cartesian(y = c(0,100)) +
  theme(legend.position="none") +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank());
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
means.plot.1;
ggsave("Output/HEKStimRes_curves.pdf", plot = means.plot.1, device = "pdf", width = 180, height = 45, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

means.plot.2 <- ggplot(plot.ts.2, aes(x = ts.times, y = ts.mean)) + # Feeds the data and the variables
  geom_line(aes(colour= ts.treat)) + # The line displays the mean values
  scale_colour_manual(values = ts.pal.2) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax= ts.mean + ts.sd, ymin = ts.mean - ts.sd, fill = ts.treat), alpha = 0.4) +
  scale_fill_manual(values = ts.pal.2) +
  geom_segment(data = arrow_pos.2, aes(x = x.vals, xend = x.end, y = y.vals, yend = y.end), size = 0.5, 
               arrow=arrow(length=unit(0.7, "mm"), type="closed")) +
  geom_text(data = stim_labels.2, aes(x = lab.x, y = lab.y, label = lab.names), size = 2.45, vjust = 0, nudge_y = 3) + # size = 2.45 because geom_text measures in mm and not in pts! pts = mm*0.35
  theme_Publication() +
  facet_grid(ts.lab ~ ts.treat, labeller = labeller(ts.treat = fac_treat.2)) + # Grouping according to the stimulus concentrations (groups)
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  labs(title = "",x = "Time (s)", y = expression(paste(Delta,"[Ca"^"2+"*"]"[i]*" (nM)"))) +
  scale_y_continuous(breaks=seq(-50,400,50)) +
  coord_cartesian(y = c(0,250)) +
  theme(legend.position="none") +
  theme(strip.background.y = element_blank(), strip.text.y = element_blank());
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
means.plot.2;
ggsave("Output/HEKStimRes_curves_pDmel.pdf", plot = means.plot.2, device = "pdf", width = 70, height = 45, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

# Statistical analysis of max intensities
agonist.resp <- resp.int_eval(list.of.means);
agonist.resp.1 <- resp.int_eval(list.of.means[c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a")]);
agonist.resp.2 <- resp.int_eval(list.of.means[c("dyn.or47a", "pDmel.or47a")]);
# For pentyl acetate responses
pentac.resp <- agonist.resp[agonist.resp$resp.stim == "pent.ac", ]; # gets the data
pentac.resp.1 <- agonist.resp.1[agonist.resp.1$resp.stim == "pent.ac", ]; # gets the data
pentac.resp.2 <- agonist.resp.2[agonist.resp.2$resp.stim == "pent.ac", ]; # gets the data

# Shapiro normality tests for each treatment
pentac.resp.shapiro <- by(pentac.resp$resp.vals, pentac.resp$resp.treat, function(x) shapiro.test(x));
pentac.resp.shapiro[["orco"]]; # W = 0.94874, p-value = 0.7282
pentac.resp.shapiro[["ctrl"]]; # W = 0.94154, p-value = 0.6768
pentac.resp.shapiro[["or47a"]]; # W = 0.81772, p-value = 0.1121
pentac.resp.shapiro[["ers.or47a"]]; # W = 0.91468, p-value = 0.4962
pentac.resp.shapiro[["dyn.or47a"]]; # W = 0.77918, p-value = 0.05423
pentac.resp.shapiro[["pDmel.or47a"]]; # W = 0.99963, p-value = 0.963

# Barlett's test for homogeneity of variances
bartlett.test(resp.vals ~ resp.treat, pentac.resp.1);
#  Bartlett test of homogeneity of variances
#  data:  resp.vals by resp.treat
#  Bartlett's K-squared = 4.7001, df = 4, p-value = 0.3195

# Anova with Dunnett's post-hoc test to test response intensities vs the 'dyn.or47a' treatment
pentac.resp.1$resp.treat <- relevel(pentac.resp.1$resp.treat, ref = "dyn.or47a");
pentac.anova <- aov(resp.vals ~ resp.treat, data = pentac.resp.1);
pentac.dun <- glht(pentac.anova, linfct = mcp(resp.treat = "Dunnett"));
summary(pentac.dun);
#  Fit: aov(formula = resp.vals ~ resp.treat, data = pentac.resp)
#  Linear Hypotheses:
#  Estimate Std. Error t value Pr(>|t|)    
#  ctrl - dyn.or47a == 0        -6.955      0.537 -12.952  <0.001 ***
#  orco - dyn.or47a == 0        -6.773      0.537 -12.613  <0.001 ***
#  or47a - dyn.or47a == 0       -6.538      0.537 -12.175  <0.001 ***
#  ers.or47a - dyn.or47a == 0   -2.451      0.537  -4.564  <0.001 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#  (Adjusted p values reported -- single-step method)

# Repeating analysis with Welch's t-test and multiple correction
# test A
t.test(pentac.resp[pentac.resp$resp.treat == "dyn.or47a",1],
       pentac.resp[pentac.resp$resp.treat == "ers.or47a",1],
       var.equal = FALSE);
# data:  pentac.resp[pentac.resp$resp.treat == "dyn.or47a", 1] and pentac.resp[pentac.resp$resp.treat == "ers.or47a", 1]
# t = 3.5212, df = 5.7012, p-value = 0.0136

# test B
t.test(pentac.resp[pentac.resp$resp.treat == "dyn.or47a",1],
       pentac.resp[pentac.resp$resp.treat == "or47a",1],
       var.equal = FALSE);
# data:  pentac.resp[pentac.resp$resp.treat == "dyn.or47a", 1] and pentac.resp[pentac.resp$resp.treat == "or47a", 1]
# t = 9.5242, df = 5.4651, p-value = 0.0001316

# test C
t.test(pentac.resp[pentac.resp$resp.treat == "dyn.or47a",1],
       pentac.resp[pentac.resp$resp.treat == "orco",1],
       var.equal = FALSE);
# data:  pentac.resp[pentac.resp$resp.treat == "dyn.or47a", 1] and pentac.resp[pentac.resp$resp.treat == "orco", 1]
# t = 9.8063, df = 5.5701, p-value = 0.0001011

pentac.resp.pvals <- c(0.0136, 0.0001316, 0.0001011);
p.adjust(pentac.resp.pvals, method = "holm"); # A = 0.0136000; B = 0.0003033; C = 0.0003033;

# test D
t.test(pentac.resp[pentac.resp$resp.treat == "pDmel.or47a",1],
       pentac.resp[pentac.resp$resp.treat == "dyn.or47a",1],
       var.equal = FALSE);
# data:  pentac.resp[pentac.resp$resp.treat == "pDmel.or47a", 1] and pentac.resp[pentac.resp$resp.treat == "dyn.or47a", 1]
# t = 8.9005, df = 2.1312, p-value = 0.01013

# For VUAA1 responses
VUAA1.resp <- agonist.resp[agonist.resp$resp.stim == "VUAA1",];
VUAA1.resp.1 <- agonist.resp.1[agonist.resp.1$resp.stim == "VUAA1",];
VUAA1.resp.2 <- agonist.resp.2[agonist.resp.2$resp.stim == "VUAA1",];
# Shapiro normality tests for each treatment
VUAA1.resp.shapiro <- by(VUAA1.resp$resp.vals, VUAA1.resp$resp.treat, function(x) shapiro.test(x));
VUAA1.resp.shapiro[["ctrl"]]; # W = 0.88746, p-value = 0.3445
VUAA1.resp.shapiro[["orco"]]; # W = 0.81653, p-value = 0.1097
VUAA1.resp.shapiro[["or47a"]]; # W = 0.89685, p-value = 0.3927
VUAA1.resp.shapiro[["ers.or47a"]]; # W = 0.94031, p-value = 0.6681
VUAA1.resp.shapiro[["dyn.or47a"]]; # W = 0.84416, p-value = 0.1767
VUAA1.resp.shapiro[["pDmel.or47a"]]; # W = 0.7609, p-value = 0.02421

# Barlett's test for homogeneity of variances
bartlett.test(resp.vals ~ resp.treat, VUAA1.resp.1);
#  Bartlett test of homogeneity of variances
#  data:  resp.vals by resp.treat
#  Bartlett's K-squared = 27.546, df = 4, p-value = 1.541e-05

# Welch's t-tests with p-value correction
# test A
t.test(VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a",1],
       VUAA1.resp[VUAA1.resp$resp.treat == "ers.or47a",1],
       var.equal = FALSE);
# data:  VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a", 1] and VUAA1.resp[VUAA1.resp$resp.treat == "ers.or47a", 1]
# t = 2.4238, df = 6.0703, p-value = 0.05113
# test B
t.test(VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a",1],
       VUAA1.resp[VUAA1.resp$resp.treat == "or47a",1],
       var.equal = FALSE);
# data:  VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a", 1] and VUAA1.resp[VUAA1.resp$resp.treat == "or47a", 1]
# t = 7.1022, df = 5.2379, p-value = 0.0007058
# test C
t.test(VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a",1],
       VUAA1.resp[VUAA1.resp$resp.treat == "orco",1],
       var.equal = FALSE);
# data:  VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a", 1] and VUAA1.resp[VUAA1.resp$resp.treat == "orco", 1]
# t = 8.3258, df = 4.7581, p-value = 0.0005179
VUAA1.resp.pvals <- c(0.05113, 0.0007058, 0.0005179);
p.adjust(VUAA1.resp.pvals, method = "holm"); # A = 0.0511300; B = 0.0015537; C = 0.0015537;

wilcox.test(VUAA1.resp[VUAA1.resp$resp.treat == "pDmel.or47a",1],
            VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a",1],
            alternative = "two.sided");
# Wilcoxon rank sum test
# data:  VUAA1.resp[VUAA1.resp$resp.treat == "pDmel.or47a", 1] and VUAA1.resp[VUAA1.resp$resp.treat == "dyn.or47a", 1]
# W = 15, p-value = 0.03571
# alternative hypothesis: true location shift is not equal to 0

# Dummy panel to set the axis limits with facet_wrap scales = "free"
bkgr.panel.intstim <- data.frame(resp.mean = c(10, 80), resp.treat = c("ctrl", "ctrl"), resp.stim = c("pent.ac", "VUAA1"));

# Response intensities barplot
intstim.plot <- resp.plot_eval(list.of.means[c("ctrl", "orco", "or47a", "ers.or47a", "dyn.or47a")]);
bp.pal <- c(cbPalette[1], cbPalette[7], cbPalette[2], cbPalette[4], cbPalette[3]);
names(bp.pal) <- as.character(unique(intstim.plot$resp.treat));

int.plot <- ggplot(intstim.plot, aes(x = resp.treat, y = resp.mean, fill = resp.treat)) +
  scale_fill_manual(values = bp.pal) +
  geom_bar(stat = "identity", colour="black", width = 0.85, alpha = 0.4) +
  geom_jitter(data = agonist.resp.1, mapping = aes(x = resp.treat, y = resp.vals), color = "grey") +
  geom_errorbar(aes(ymin = resp.mean - resp.sd, ymax = resp.mean + resp.sd), width=.2) +
  geom_blank(data = bkgr.panel.intstim) +
  scale_x_discrete(labels = stim_xaxis) +
  theme_Publication() +
  labs(x = "", y = expression(paste(Delta,"[Ca"^"2+"*"]"[i]*" (nM)"))) +
  #scale_y_continuous(breaks=seq(0,100,10)) +
  # coord_cartesian(y = c(0,80)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(family = "Helvetica")) +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  facet_wrap(~ resp.stim, nrow = 1, scales = "free", labeller = labeller(resp.stim = fac_stim));
int.plot;
ggsave("Output/HEKStimRes_Max.pdf", plot = int.plot, device = "pdf", width = 70, height = 60, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

bkgr.panel.intstim.2 <- data.frame(resp.mean = c(50, 180), resp.treat = c("dyn.or47a", "dyn.or47a"), resp.stim = c("pent.ac", "VUAA1"));

intstim.plot.2 <- resp.plot_eval(list.of.means[c("dyn.or47a", "pDmel.or47a")]);
bp.pal.2 <- c(cbPalette[3], cbPalette[8]);
names(bp.pal.2) <- as.character(unique(plot.ts.2$ts.treat));

int.plot.2 <- ggplot(intstim.plot.2, aes(x = resp.treat, y = resp.mean, fill = resp.treat)) +
  scale_fill_manual(values = bp.pal.2) +
  geom_bar(stat = "identity", colour="black", width = 0.85, alpha = 0.4) +
  geom_jitter(data = agonist.resp.2, mapping = aes(x = resp.treat, y = resp.vals), color = "grey") +
  geom_errorbar(aes(ymin = resp.mean - resp.sd, ymax = resp.mean + resp.sd), width=.2) +
  geom_blank(data = bkgr.panel.intstim.2) +
  scale_x_discrete(labels = stim_xaxis.2) +
  theme_Publication() +
  labs(x = "", y = expression(paste(Delta,"[Ca"^"2+"*"]"[i]*" (nM)"))) +
  #scale_y_continuous(breaks=seq(0,100,10)) +
  # coord_cartesian(y = c(0,80)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(family = "Helvetica")) +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  facet_wrap(~ resp.stim, nrow = 1, scales = "free", labeller = labeller(resp.stim = fac_stim.tiny));
int.plot.2;
ggsave("Output/HEKStimRes_Max_pDmelOR.pdf", plot = int.plot.2, device = "pdf", width = 40, height = 55, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

# Distribution of response intensities for pentyl acetate
pentac.rois <- lapply(polished.vals, pentac_roi.resp);
pentac.cs <- pentac_distr.mkr(pentac.rois, Data.meta);

pentac.cs.dls <- pentac_csplot.mkr(pentac.cs);
pentac.pr.stat <- pentac.cs.dls[[1]];

pentac.pr.stat.1 <- pentac.pr.stat[-which(pentac.pr.stat[,"pr.treat"] == "pDmel.or47a"),]
pentac.pr.stat.2 <- pentac.pr.stat[c(which(pentac.pr.stat[,"pr.treat"] == "dyn.or47a"),
                                           which(pentac.pr.stat[,"pr.treat"] == "pDmel.or47a")),];

# Shapiro normality test for each treatment
pentac.pr.shapiro <- by(pentac.pr.stat$pr.vals, pentac.pr.stat$pr.treat, function(x) shapiro.test(x));
pentac.pr.shapiro[["ctrl"]]; # W = 0.73461, p-value = 0.0213
# These are the important treatments to compare!
pentac.pr.shapiro[["orco"]]; # W = 0.9756, p-value = 0.9098
pentac.pr.shapiro[["or47a"]]; # W = 0.92595, p-value = 0.569
pentac.pr.shapiro[["ers.or47a"]]; # W = 0.90624, p-value = 0.4453
pentac.pr.shapiro[["dyn.or47a"]]; # W = 0.87054, p-value = 0.2686
pentac.pr.shapiro[["pDmel.or47a"]]; # W = 0.95315, p-value = 0.5833

# Barlett's test for homogeneity of variances
pentac.pr.stat2 <- pentac.pr.stat[pentac.pr.stat$pr.treat != "ctrl",];
bartlett.test(pr.vals ~ pr.treat, pentac.pr.stat2);
#  Bartlett test of homogeneity of variances
#  data:  pr.vals by pr.treat
#  Bartlett's K-squared = 25.076, df = 4, p-value = 4.856e-05

# Welch's t-tests with p-value correction
# test A
t.test(pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a",1],
       pentac.pr.stat2[pentac.pr.stat2$pr.treat == "ers.or47a",1],
       var.equal = FALSE);
# data:  pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a", 1] and pentac.pr.stat2[pentac.pr.stat2$pr.treat == "ers.or47a", 1]
# t = 2.561, df = 6.7679, p-value = 0.03858
# test B
t.test(pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a",1],
       pentac.pr.stat2[pentac.pr.stat2$pr.treat == "or47a",1],
       var.equal = FALSE);
# data:  pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a", 1] and pentac.pr.stat2[pentac.pr.stat2$pr.treat == "or47a", 1]
# t = 9.2944, df = 4.2045, p-value = 0.0005866
# test C
t.test(pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a",1],
       pentac.pr.stat2[pentac.pr.stat2$pr.treat == "orco",1],
       var.equal = FALSE);
# data:  pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a", 1] and pentac.pr.stat2[pentac.pr.stat2$pr.treat == "orco", 1]
# t = 10.155, df = 4.0465, p-value = 0.0004992
pentac.pr.pvals <- c(0.03858, 0.0005866, 0.0004992);
p.adjust(pentac.pr.pvals, method = "holm"); # A = 0.0385800; B = 0.0014976; C = 0.0014976;

# test D
t.test(pentac.pr.stat2[pentac.pr.stat2$pr.treat == "pDmel.or47a",1],
       pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a",1],
       var.equal = FALSE);
# data:  pentac.pr.stat2[pentac.pr.stat2$pr.treat == "pDmel.or47a", 1] and pentac.pr.stat2[pentac.pr.stat2$pr.treat == "dyn.or47a", 1]
# t = 6.6831, df = 2.9391, p-value = 0.007285

pentac.cs.df <- pentac.cs.dls[[3]];

arrow_pos2 <- data.frame(x.vals = c(unique(pentac.cs.df$cd.xthr)),
                        x.end = c(unique(pentac.cs.df$cd.xthr)),
                        y.vals = 90,
                        y.end = 80);

cbPalette <- c("#000000", # black
               "#E69F00", # ocra
               "#56B4E9", # cyan
               "#009E73", # green
               "#F0E442", # yellow
               "#0072B2", # blue
               "#D55E00", # vermillion
               "#CC79A7");# reddish purple

ts.pal.2 <- c(cbPalette[1], cbPalette[3], cbPalette[4], cbPalette[2], cbPalette[7], cbPalette[8]);
names(ts.pal.2) <- as.character(unique(pentac.cs.df$cd.treat));

pentac.cs.plot <- ggplot(pentac.cs.df, aes(x = cd.xvals, y = cd.means)) + # Feeds the data and the variables
  geom_rect(aes(xmin = -Inf, xmax = unique(cd.xthr), ymin = -Inf, ymax = +Inf), fill="#D3D3D3", alpha=0.5) +
  geom_text(aes(x = 15, y = 95, label = paste(cd.mean.label,"\u00B1",cd.sd.label,"%", sep = " ")), vjust = 0, size = 2.45) +
  geom_segment(data = arrow_pos2, aes(x = x.vals, xend = x.end, y = y.vals, yend = y.end), size = 0.5, 
               arrow=arrow(length=unit(0.7, "mm"), type="closed")) +
  geom_line(aes(colour= cd.treat)) + # The line displays the mean values
  scale_colour_manual(values = ts.pal.2) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax= cd.means + cd.sd, ymin = cd.means - cd.sd, fill = cd.treat), alpha = 0.4) +
  scale_fill_manual(values = ts.pal.2) +
  theme_Publication() +
  labs(x = expression(paste(Delta,"[Ca"^"2+"*"]"[i]*" (\u2265 nM)")), y = "Analyzed cells (%)") + 
  scale_x_continuous(breaks=seq(-100,200,5)) +
  coord_cartesian(x = c(-5,25)) +
  facet_grid(cd.agonist ~ cd.treat, labeller = labeller(cd.treat = fac_treat.tot)) + # Grouping according to the stimulus concentrations (groups)
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
pentac.cs.plot;
ggsave("Output/HEKSCumsum_pentac.pdf", plot = pentac.cs.plot, device = "pdf", width = 180, height = 50, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

# Distribution of response intensities for VUAA1
VUAA1.rois <- lapply(polished.vals, VUAA1_roi.resp);
VUAA1.cs <- VUAA1_distr.mkr(VUAA1.rois, Data.meta);
VUAA1.cs.dls <- VUAA1_csplot.mkr(VUAA1.cs);
VUAA1.pr.stat <- VUAA1.cs.dls[[1]];

VUAA1.pr.stat.1 <- VUAA1.pr.stat[-which(VUAA1.pr.stat[,"pr.treat"] == "pDmel.or47a"),]
VUAA1.pr.stat.2 <- VUAA1.pr.stat[c(which(VUAA1.pr.stat[,"pr.treat"] == "dyn.or47a"),
                                     which(VUAA1.pr.stat[,"pr.treat"] == "pDmel.or47a")),];

# Shapiro normality test for each treatment
VUAA1.pr.shapiro <- by(VUAA1.pr.stat$pr.vals, VUAA1.pr.stat$pr.treat, function(x) shapiro.test(x));
VUAA1.pr.shapiro[["ctrl"]]; # W = 0.55218, p-value = 0.000131
# These are the important treatments to compare!
VUAA1.pr.shapiro[["orco"]]; # W = 0.72386, p-value = 0.01674 (!!!)
VUAA1.pr.shapiro[["or47a"]]; # W = 0.89694, p-value = 0.3932
VUAA1.pr.shapiro[["ers.or47a"]]; # W = 0.96476, p-value = 0.8407
VUAA1.pr.shapiro[["dyn.or47a"]]; # W = 0.95692, p-value = 0.7864
VUAA1.pr.shapiro[["pDmel.or47a"]]; # W = 0.90935, p-value = 0.4159
# Welch's t-tests with p-value correction
# test A
t.test(VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a",1],
       VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "ers.or47a",1],
       var.equal = FALSE);
# data:  VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a", 1] and VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "ers.or47a", 1]
# t = 3.376, df = 6.2405, p-value = 0.01406
# test B
t.test(VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a",1],
       VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "or47a",1],
       var.equal = FALSE);
# data:  VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a", 1] and VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "or47a", 1]
# t = 11.706, df = 5.5986, p-value = 3.791e-05
# test C
wilcox.test(VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a",1],
            VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "orco",1],
            var.equal = FALSE);
# data:  VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a", 1] and VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "orco", 1]
# W = 25, p-value = 0.007937
VUAA1.pr.pvals <- c(0.01406, 3.791e-05, 0.007937);
p.adjust(VUAA1.pr.pvals, method = "holm"); # A = 0.01587400; B = 0.00011373; C = 0.01587400;

# test D
t.test(VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "pDmel.or47a",1],
       VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a",1],
       var.equal = FALSE);
# data:  VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "pDmel.or47a", 1] and VUAA1.pr.stat[VUAA1.pr.stat$pr.treat == "dyn.or47a", 1]
# t = 7.5245, df = 3.6354, p-value = 0.002422

VUAA1.cs.df <- VUAA1.cs.dls[[3]];

arrow_pos3 <- data.frame(x.vals = c(unique(VUAA1.cs.df$cd.xthr)),
                         x.end = c(unique(VUAA1.cs.df$cd.xthr)),
                         y.vals = 110,
                         y.end = 100);

VUAA1.cs.plot <- ggplot(VUAA1.cs.df, aes(x = cd.xvals, y = cd.means)) + # Feeds the data and the variables
  geom_rect(aes(xmin = -Inf, xmax = unique(cd.xthr), ymin = -Inf, ymax = +Inf), fill="#d3d3d3", alpha=0.5) +
  geom_text(aes(x = 130, y = 90, label = paste(cd.mean.label,"\u00B1",cd.sd.label,"%", sep = " ")), vjust = 0, size = 2.45) +
  scale_colour_manual(values = ts.pal.2) +
  # The ribbon displays the mean ± st.dev values
  geom_ribbon(aes(ymax= cd.means + cd.sd, ymin = cd.means - cd.sd, fill = cd.treat), alpha = 0.4) +
  geom_line(aes(colour= cd.treat)) + # The line displays the mean values
  geom_segment(data = arrow_pos3, aes(x = x.vals, xend = x.end, y = y.vals, yend = y.end), size = 0.5, 
               arrow=arrow(length=unit(0.7, "mm"), type="closed")) +
  scale_fill_manual(values = ts.pal.2) +
  theme_Publication() +
  labs(x = expression(paste(Delta,"[Ca"^"2+"*"]"[i]*" (\u2265 nM)")), y = "Analyzed cells (%)") +
  scale_x_continuous(breaks=seq(-50,500,50)) +
  coord_cartesian(x = c(-5,250)) +
  facet_grid(cd.agonist ~ cd.treat, labeller = labeller(cd.treat = fac_treat.tot)) + # Grouping according to the stimulus concentrations (groups)
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0")) +
  theme(legend.position="none");
# Attention when saving as .pdf or other vectorial formats: it may not recognize the greek letters!
VUAA1.cs.plot;
ggsave("Output/HEKSCumsum_VUAA1.pdf", plot = VUAA1.cs.plot, device = "pdf", width = 180, height = 50, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

# Percentage of responding cells plot
# barplot data
pentac.pr.barplot <- pentac.cs.dls[[2]];
pentac.pr.barplot.1 <- pentac.pr.barplot[-which(pentac.pr.barplot[,"pr.levels"] == "pDmel.or47a"),]
pentac.pr.barplot.2 <- pentac.pr.barplot[c(which(pentac.pr.barplot[,"pr.levels"] == "pDmel.or47a"),
                                           which(pentac.pr.barplot[,"pr.levels"] == "dyn.or47a")),];

VUAA1.pr.barplot <- VUAA1.cs.dls[[2]];
VUAA1.pr.barplot.1 <- VUAA1.pr.barplot[-which(pentac.pr.barplot[,"pr.levels"] == "pDmel.or47a"),];
VUAA1.pr.barplot.2 <- VUAA1.pr.barplot[c(which(pentac.pr.barplot[,"pr.levels"] == "pDmel.or47a"),
                                         which(pentac.pr.barplot[,"pr.levels"] == "dyn.or47a")),];

pr.barplot.1 <- rbind(pentac.pr.barplot.1, VUAA1.pr.barplot.1);
# point jitter data
pr.jitterplot.1 <- rbind(pentac.pr.stat.1, VUAA1.pr.stat.1);
colnames(pr.jitterplot.1) <- c("pr.vals","pr.agonist","pr.levels");

bkgr.panel.pr <- data.frame(pr.means = c(80, 80), pr.levels = c("ctrl", "ctrl"), pr.agonist = c("pent.ac", "VUAA1"));

pr.plot <- ggplot(pr.barplot.1, aes(x = pr.levels, y = pr.means, fill = pr.levels)) +
  scale_fill_manual(values = bp.pal) +
  geom_bar(stat = "identity", colour="black", width = 0.85, alpha = 0.4) +
  geom_jitter(data = pr.jitterplot.1, mapping = aes(x = pr.levels, y = pr.vals), color = "grey") +
  geom_errorbar(aes(ymin = pr.means - pr.sd, ymax = pr.means + pr.sd), width=.2) +
  geom_blank(data = bkgr.panel.pr) +
  scale_x_discrete(labels = stim_xaxis) +
  facet_wrap(~ pr.agonist, nrow = 1, scales = "free", labeller = labeller(pr.agonist = fac_stim)) +
  theme_Publication() +
  labs(x = "", y = "Responding cells (%)") +
#  scale_y_continuous(breaks=seq(0,100,10)) +
#  coord_cartesian(y = c(0,80)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(family = "Helvetica")) +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0"));
pr.plot;
ggsave("Output/HEKPercRespCells.pdf", plot = pr.plot, device = "pdf", width = 70, height = 60, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();

pr.barplot.2 <- rbind(pentac.pr.barplot.2, VUAA1.pr.barplot.2);
# point jitter data
pr.jitterplot.2 <- rbind(pentac.pr.stat.2, VUAA1.pr.stat.2);
colnames(pr.jitterplot.2) <- c("pr.vals","pr.agonist","pr.levels");

bkgr.panel.pr.2 <- data.frame(pr.means = c(100, 100), pr.levels = c("dyn.or47a", "dyn.or47a"), pr.agonist = c("pent.ac", "VUAA1"));

bp.pal.2 <- c(cbPalette[8], cbPalette[3]);
names(bp.pal.2) <- as.character(unique(pr.barplot.2$pr.levels));

pr.plot.2 <- ggplot(pr.barplot.2, aes(x = pr.levels, y = pr.means, fill = pr.levels)) +
  scale_fill_manual(values = bp.pal.2) +
  geom_bar(stat = "identity", colour="black", width = 0.85, alpha = 0.4) +
  geom_jitter(data = pr.jitterplot.2, mapping = aes(x = pr.levels, y = pr.vals), color = "grey") +
  geom_errorbar(aes(ymin = pr.means - pr.sd, ymax = pr.means + pr.sd), width=.2) +
  geom_blank(data = bkgr.panel.pr.2) +
  scale_x_discrete(labels = stim_xaxis.2) +
  facet_wrap(~pr.agonist, ncol = 2, scales = "free", labeller = labeller(pr.agonist = fac_stim.tiny)) +
  theme_Publication() +
  labs(x = "", y = "Responding cells (%)") +
  #  scale_y_continuous(breaks=seq(0,100,10)) +
  #  coord_cartesian(y = c(0,80)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(axis.title.y = element_text(family = "Helvetica")) +
  theme(strip.background = element_rect(color = "black", fill="#e0e0e0"));
pr.plot.2;
ggsave("Output/HEKPercRespCells_pDmelOR.pdf", plot = pr.plot.2, device = "pdf", width = 40, height = 55, units = "mm", useDingbats=FALSE); # Saves the plot
dev.off();