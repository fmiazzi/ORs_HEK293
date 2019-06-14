/*  *********************
 *  *  Fura-2 Analyzer  *
 *  *********************
 *  
 *  DESCRIPTION
 *  This script performs:
 *  - background and flat-field correction
 *  - movement artifact correction
 *  - [Ca2+]i calculation
 *  for functional imaging stacks acquired using the Fura-2 calcium indicator
 *  and allows the user to selct regions of interest from which 
 *  extract intensity values for subsequent statistical analysis
 *  
 *  INPUT
 *  - the path to the folder where the:
 *  		> .tiff stack with the intensity values acquired with a 340nm stimulation
 *  		> .tiff stack with the intensity values acquired with a 380nm stimulation
 *  		> .tiff image of the cells under transmission light
 *  are located;
 *  - a .tiff stack with > 340nm background fluorescence image 
 *  					 > 380nm background fluorescence image 
 *  					 > 340nm stimulation for flat-field correction
 *  					 > 380nm stimulation for flat-field correction;
 *  - the path to the folder where to store the output files;
 *  
 *  OUTPUT
 *  - .tif stacks for the background and flat-field corrected 340nm and 380nm measurements
 *  - .tif stacks for the stabilized 340nm and 380nm measurements
 *  - a .txt file containing the coefficients used by the Image Stabilizer plugin
 *  - a .zip file containing the ROI used for background subtraction (for the [Ca2+]i values)
 *  - .csv files containing the mean and standard deviation values used for background subtraction
 *    for the 340nm and 380nm measurements (for the [Ca2+]i values)
 *  - a .tif hyperstack with the calibrated [Ca2+]i values
 *  - a .zip file containg the ROI used for mean gray values extraction
 *  - .csv files containing the mean gray values throughout the stabilized 340nm 
 *   and 380nm stacks of the selected ROIs
 *  - a .txt file documenting the image elaboration process performed by the script
 *  
 *  REQUIREMENTS 
 *  This script requires the "Image Stabilizer" Plugin by Kang Li (2008)
 *  http://www.cs.cmu.edu/~kangli/code/Image_Stabilizer.html
 *  And requires that the [Ca2+]i, frame size and frame rate calibration constants
 *  are entered in the 'Constants' section of this script!
 */

/* Import of Classes */
importClass(Packages.java.io.File);
importClass(Packages.java.util.Calendar);
importClass(Packages.ij.IJ);
importClass(Packages.ij.Prefs);
importClass(Packages.ij.WindowManager)
importClass(Packages.ij.ImagePlus);
importClass(Packages.ij.ImageStack);
importClass(Packages.ij.gui.Roi);
importClass(Packages.ij.gui.GenericDialog);
importClass(Packages.ij.gui.NonBlockingGenericDialog);
importClass(Packages.ij.gui.WaitForUserDialog);
importClass(Packages.ij.gui.Overlay);
importClass(Packages.ij.gui.Plot);
importClass(Packages.ij.measure.ResultsTable);
importClass(Packages.ij.plugin.ImageCalculator);
importClass(Packages.ij.plugin.Duplicator);
importClass(Packages.ij.plugin.frame.RoiManager);
importClass(Packages.ij.plugin.HyperStackConverter);
importClass(Packages.ij.process.ImageStatistics);
importClass(Packages.ij.process.ImageConverter);
importClass(Packages.ij.plugin.filter.BackgroundSubtracter);

/* Constants */
/*
 * These values are used to calibrate the [Ca2+]i time series and the ROI plot.
 * Please check whether they are valid for your experiments; in case you changed
 * setup or protocol conditions, please update these constants here.
 */
// Constants for intracellular calcium concentration calculation
var Rmin = 0.1; //var Rmin = 0.1734732;
var Rmax = 5.1; // var Rmax = 5.8363955;
var Fratio = 13.4415395;
var Keff = 130; //nM [Ca2+]i
var Ktot = Math.round(Fratio*Keff);
//var Ktot = 1300;

// Time interval between frames
var deltaT = 5; //seconds
var tScale = "sec";

// Calibration (40x Water objective, Dieter's optical imaging lab)
var pixDist = 37; // distance in pixels
var umDist = 10; // known distance in micrometers
var Scale = "µm"; // measure scale

/* Check requirements */

// Check Software Information
IJ.log("****************************");
IJ.log("** Fura-2 Analyzer (v0.1) **");
IJ.log("****************************");
IJ.log("\n")

var ijVersion = IJ.getVersion(); // Checks the current version of Fiji ImageJ
var lastIJtestedVersion = "2.0.0-rc-69/1.52n"; // This is the last ImageJ version this script was checked for
var javaRE = java.lang.System.getProperty("java.version"); // Checks Java Realtime Environment version running this script
var lastJREtestedVersion = "1.8.0_66"; // This is the last JavaRE version this script was checked for

// Message error in case of different ImageJ version than the one tested
if (ijVersion !== lastIJtestedVersion) {
IJ.error("This script was tested using Fiji ImageJ " + lastIJtestedVersion + "\n \
your current version is the " + ijVersion + "\n \
The script may work incorrectly!");

	IJ.log("********** Attention! ***********");
	IJ.log("This script was tested using:");
	IJ.log("Fiji ImageJ "+ lastIJtestedVersion +",");
	IJ.log("your current version is the " + ijVersion);
	IJ.log("This script may work incorrectly!");
	IJ.log("*********************************");
	IJ.log("\n");
}

// Message error in case of different JavaRE version than the one tested
if (javaRE !== lastJREtestedVersion) {
IJ.error("This script was tested using JavaRE " + lastJREtestedVersion + "\n \
your current version is the " + javaRE + "\n \
The script may work incorrectly!");

	IJ.log("********** Attention! ***********");
	IJ.log("This script was tested using:");
	IJ.log("Java Runtime Environment " + lastJREtestedVersion + ",");
	IJ.log("your current version is the "+ javaRE);
	IJ.log("This script may work incorrectly!");
	IJ.log("*********************************");
	IJ.log("\n");
}

// Print Software and Time Information
IJ.log("This script is running on:");
var osName = java.lang.System.getProperty("os.name"); // Operative system name
var osArch = java.lang.System.getProperty("os.arch"); // Operative system architecture
var osVersion = java.lang.System.getProperty("os.version"); // Operative system version
IJ.log(osName +" "+ osVersion +" "+ osArch);
IJ.log("Fiji ImageJ "+ ijVersion); // ImageJ version
var javaREv = java.lang.System.getProperty("java.vendor"); // Java Realtime Environment source
var javaVM = java.lang.System.getProperty("java.vm.specification.version"); // Java Virtual Machine specification version
var javaVMi = java.lang.System.getProperty("java.vm.version");  // Java Virtual Machine version
var javaVMv = java.lang.System.getProperty("java.vm.specification.vendor"); // Java Virtual Machine source
IJ.log("Java Runtime Environment "+ javaRE +", by "+ javaREv);
IJ.log("Java Virtual Machine "+ javaVM +" ("+ javaVMi +"), by "+ javaVMv);
var cal = Calendar.getInstance(); // Get the current date and time
IJ.log("@ "+ cal.getTime().toString());
IJ.log("\n");

/* Get input files */
batchMode = true; // prevents the refresh of the open images to speed up the computation

// Get the Fura-2 input file paths
// useJFileChooser is used to display title bars from file open and directory chooser dialogs on Mac OS X 10.11 and later
Prefs.useJFileChooser = true; //
var inputDir = IJ.getDirectory("Select the folder containing the input files for this measurement");
// Keep here for debugging 
//IJ.log("Input directory: ");
//IJ.log(inputDir);
//IJ.log("\n");
Prefs.useJFileChooser = false; //

var inputList = getTiffs(inputDir); //Gets the list of tiff files in Input Directory
IJ.log("* List of input files *");
if (inputList.length !== 3) {
	IJ.error("The input folder must contain exactly 3 .tiff files. \n" +
	"Please check your input folder and the script for specifications.");
	throw new Error("The input folder must contain only 3 .tif or .tiff files");
}
var fura340 = "";
var fura380 = "";
var trImage = "";
for (var i = 0; i<inputList.length; i++) {
	var fileStr = String(inputList[i]);
	if ((fileStr.endsWith(".tif") || fileStr.endsWith(".tiff")) == true) { // Analyzes only .tif(f) files!
		fileStrNoExt = fileStr.replace(/\.[^/.]+$/, "");
		var fileType = fileStrNoExt.slice(-5);
		switch (fileType) {
		case "340nm":
			fura340 = inputList[i];
			IJ.log("340nm stack: " + fura340);
		break;
		case "380nm":
			fura380 = inputList[i];
			IJ.log("380nm stack: " + fura380);
		break;
		case "trans":
			trImage = inputList[i];
			IJ.log("Transmission image: " + trImage);
		break;
		default:
		IJ.error( fileStr + " is an unidentified input file in the " + String(inputDir).match(/([^\/]*)\/*$/)[1] + " input directory");
		throw new Error("Unidentified input file in the " + String(inputDir).match(/([^\/]*)\/*$/)[1] + " input directory");
		}
	}	else {
		IJ.log(fileStr +" cannot be processed, only .tif and .tiff files are accepted!");
	}
}

// Get the file with the correction stack
Prefs.useJFileChooser = true;
var ffCorr = IJ.getFilePath("Select the .tiff stack for background and flat-field correction");
var ffCorrName = ffCorr.match(/([^\/]*)\/*$/)[1];
IJ.log("Stack for background and flat-field correction: " + ffCorrName);
IJ.log("\n");
Prefs.useJFileChooser = false;

// Get the path of the output directory
Prefs.useJFileChooser = true; //
var outputDir = IJ.getDirectory("Save output files in...");
// Keep here for debugging 
//IJ.log("Output directory: ");
//IJ.log(outputDir);
//IJ.log("\n");
Prefs.useJFileChooser = false; //

/* Main analysis */

// Creates a new directory where to store the output for the processed file
var fileName = inputDir.match(/([^\/]*)\/*$/)[1];
var stackOutDir = outputDir + fileName + "/";
var bln = new File(stackOutDir).mkdir();

var indexCorr = ffCorr_fun(inputDir, fura340, fura380, ffCorr, stackOutDir); // Background and flat-field correction
if (indexCorr !== true) {
	IJ.error("The ffCorr_fun() function did not work correctly on " + fileName);
}

var indexStab = stabIm(stackOutDir, fura340, fura380); // Image Stabilization
if (indexStab !== true) {
	IJ.error("The stabIm() function did not work correctly on " + fileName);
}

var indexCaCalc = caCalc(inputDir, stackOutDir, fura340, fura380, trImage, Rmin, Rmax, Ktot, deltaT, tScale, pixDist, umDist, Scale); // [Ca2+]i calculation
if (indexCaCalc !== true) {
	IJ.error("The indexCaCalc() function did not work correctly on " + fileName);
}

var indexSegm = segmIm(inputDir, stackOutDir, fura340, trImage); // Image segmentation
if (indexSegm !== true) {
	IJ.error("The indexSegm() function did not work correctly on " + fileName);
}

var indexRes = resultsMaker(stackOutDir, fura340, fura380); // Writes the CSV file and saves all result files
if (indexRes !== true) {
	IJ.error("The indexRes() function did not work correctly on " + fileName);
}

// Log the parameters used during the analysis
IJ.log("** Constants used **");
IJ.log("\n");
IJ.log("Calibration bar: " + pixDist/umDist + " pix/" + Scale);
IJ.log("\n");
IJ.log("Time series: " + deltaT + " " + tScale + "/frame");
IJ.log("\n");
IJ.log("[Ca2+]i calculation formula: [Ca2+]i = Keff(R - Rmin)/(Rmax - R)");
IJ.log("with R = f340nm/f380nm");
IJ.log("Keff = " + Ktot + " nM");
IJ.log("Rmin = " + Rmin);
IJ.log("Rmax = " + Rmax);
IJ.log("\n");
IJ.log("*** " + fileName + " analysis successful! ***");
// Save and close log file
IJ.selectWindow("Log");
IJ.saveAs("Text", stackOutDir + fileName + "_log.txt");
//IJ.run("Close");


/* Functions */

// Gets a list of .tif files from the Input Directory
function getTiffs(dirst) { 
        var dir = new java.io.File(dirst); 
        var names = dir.list();
        var fileNameList = [];
//       print("nr. files in the folder: " + names.length);
        for (var i = 0 ; i < names.length; i++) {
        var fileName = names[i].match(/([^\/]*)\/*$/)[1];
		var isHidden = fileName.substring(0,2);
        	if (isHidden == "._") {
//        		print(names[i]);
        		names[i] = [];
        	}
        }
        for (var i = 0 ; i < names.length; i++) {
   			if ((names[i].endsWith(".tif") || names[i].endsWith(".tiff")) == true) {
            	fileNameList.push(names[i]);
            }  
        }
//        print(fileNameList);
        return fileNameList;
}

// Background and flat-field correction
function ffCorr_fun(inDir, f340, f380, furaCorr, outDir) {
	batchMode = true;
	// Get file names without extension
	var file340NoExt = f340.replace(/\.[^/.]+$/, "");
	var file380NoExt = f380.replace(/\.[^/.]+$/, "");
	var fileName = file380NoExt.slice(0, -6);
	// Open input files
	var image340 = new ImagePlus();
	var image380 = new ImagePlus();
	var imageCorr = new ImagePlus();
	image340 = IJ.openImage(inDir + f340);
//	image340.show();
	image380 = IJ.openImage(inDir + f380);
//	image380.show();
	imageCorr = IJ.openImage(furaCorr);
//	imageCorr.show();
	// Get the individual background and flat-field correction images
	var bk340 = new ImagePlus();
	var bk380 = new ImagePlus();
	var ff340 = new ImagePlus();
	var ff380 = new ImagePlus();
	imageCorr.setSlice(1);
	bk340 = imageCorr.crop();
//	bk340.show();
	imageCorr.setSlice(2);
	bk380 = imageCorr.crop();
//	bk380.show();
	imageCorr.setSlice(3);
	ff340 = imageCorr.crop();
//	ff340.show();
	imageCorr.setSlice(4);
	ff380 = imageCorr.crop();
//	ff380.show();
	imageCorr.changes = false;
	imageCorr.close();
	// Background and flat-field correction for the 340nm file
	var stack340bc = new ImageCalculator().run("Subtract create stack", image340, bk340);
	stack340bc.setTitle("Stack340BkgCorr");
	IJ.run(stack340bc, "16 colors", "");
	IJ.run(stack340bc, "Enhance Contrast", "saturated=0.1");
	stack340bc.show();
	var ff340bc = new ImageCalculator().run("Subtract create", ff340, bk340);
	ff340bc.setTitle("Ff340BkgCorr");
	// Mean(S - B) is equivalent to Mean(S) - Mean(B) as: Sum(Sij - Bij) = Sum(Sij) - Sum(Bij);
	var meanFF340 = ff340bc.getStatistics().mean;
	IJ.run("Conversions...", " "); // Prevents scaling when converting between 32- and 16-bit images and viceversa
	new ImageConverter(stack340bc).convertToGray32();
	new ImageConverter(ff340bc).convertToGray32();
	ff340bc.show();
	IJ.run("Calculator Plus", "i1=" + stack340bc.getTitle() + " i2=" + ff340bc.getTitle() + 
	" operation=[Divide: i2 = (i1/i2) x k1 + k2] k1=" + meanFF340 + " k2=0 create");
	stack340ff = WindowManager.getCurrentImage();
	IJ.run(stack340ff, "16-bit", "");
	IJ.saveAs(stack340ff, "Tiff", outDir + file340NoExt + "_FlatField.tif");
	stack340bc.changes = false;
	stack340bc.close();
	ff340bc.changes = false;
	ff340bc.close();
	stack340ff.close();
	// Background and flat-field correction for the 380nm file
	var stack380bc = new ImageCalculator().run("Subtract create stack", image380, bk380);
	stack380bc.setTitle("Stack380BkgCorr");
	// Subtract background ROI from the 380nm measurements
	IJ.run(stack380bc, "16 colors", "");
	IJ.run(stack380bc, "Enhance Contrast", "saturated=0.1");
	stack380bc.show();
	var ff380bc = new ImageCalculator().run("Subtract create", ff380, bk380);
	ff380bc.setTitle("Ff380BkgCorr");
	// Mean(S - B) is equivalent to Mean(S) - Mean(B) as: Sum(Sij - Bij) = Sum(Sij) - Sum(Bij);
	var meanFF380 = ff380bc.getStatistics().mean;
	new ImageConverter(stack380bc).convertToGray32();
	new ImageConverter(ff380bc).convertToGray32();
	stack380bc.show();
	ff380bc.show();
	IJ.run("Calculator Plus", "i1=" + stack380bc.getTitle() + " i2=" + ff380bc.getTitle() + 
	" operation=[Divide: i2 = (i1/i2) x k1 + k2] k1=" + meanFF380 + " k2=0 create");
	stack380ff = WindowManager.getCurrentImage();
	IJ.run(stack380ff, "16-bit", "");
	IJ.run("Conversions...", "scale weighted"); // Sets back to scaling options necessary for later...
	IJ.saveAs(stack380ff, "Tiff", outDir + file380NoExt + "_FlatField.tif");
	stack380bc.changes = false;
	stack380bc.close();
	ff380bc.changes = false;
	ff380bc.close();
	stack380ff.close();
	batchMode = false;
	IJ.run("Collect Garbage");
	return true;
}

// Image stabilization
function stabIm(fileDir, stack340name, stack380name) {
	var file340NoExt = stack340name.replace(/\.[^/.]+$/, "");
	var file380NoExt = stack380name.replace(/\.[^/.]+$/, "");
	var fileName = file380NoExt.slice(0, -6);
	var stack380 = IJ.openImage(fileDir + file380NoExt + "_FlatField.tif");
	// Apply the Image Stabilizer plugin to the 380nm stack and save stabilized stack
	stack380.show();
	IJ.run(stack380, "Image Stabilizer", "transformation=Translation maximum_pyramid_levels=1 " + 
	"template_update_coefficient=0.90 maximum_iterations=200 error_tolerance=0.0000001 " +
	"log_transformation_coefficients");
	// For the stacks with strong movement ->
	/*
	 * IJ.run(imp, "Image Stabilizer", "transformation=Translation maximum_pyramid_levels=2 " +
	 * "template_update_coefficient=0.90 maximum_iterations=500 error_tolerance=0.0000001 " +
	 * "log_transformation_coefficients");
	 */
	IJ.saveAs(stack380, "Tiff", fileDir + file380NoExt + "_FlatFieldStabilized.tif"); 
	stack380.close();
	// Apply stabilization coefficients to 340nm stack and save stabilized stack
	var stack340 = IJ.openImage(fileDir + file340NoExt + "_FlatField.tif");
	stack340.show();
	IJ.run(stack340, "Image Stabilizer Log Applier", " ");
	IJ.saveAs(stack340, "Tiff", fileDir + file340NoExt + "_FlatFieldStabilized.tif");
	stack340.close();
	// Saves and then closes the coefficient log file
	IJ.selectWindow(file380NoExt + "_FlatField.log");
	IJ.saveAs("Text", fileDir + fileName + "_StabCoeff.txt"); 
	IJ.run("Close");
	IJ.run("Collect Garbage");
	return true
}

// [Ca2+]i visualization
function caCalc(inputDir, outputDir, stack340name, stack380name, transmissionImage, Rmin, Rmax, Keff, deltaT, tScale, pixDist, umDist, Scale) {
	// Open 380nm stack
	batchMode = true;
	var file340NoExt = stack340name.replace(/\.[^/.]+$/, "");
	var file380NoExt = stack380name.replace(/\.[^/.]+$/, "");
	var fileName = file340NoExt.slice(0, -6);
	var stack380 = IJ.openImage(outputDir + file380NoExt + "_FlatFieldStabilized.tif");
	IJ.run(stack380, "HiLo", "");
	IJ.run(stack380, "Enhance Contrast", "saturated=0.1");
	stack380.setSlice(1);
	batchMode = false;
	//stack380.show();
	// Create the mask for background elimination
	var satCa = new Boolean (false); // Is set to false until the user is happy with the overall result
	while (satCa == false) {
	var satMask = new Boolean(false); // Is set to false until the user is happy with the masking result
	while (satMask == false) {
		batchMode = true;
		IJ.setTool("line");
		var mask380 = stack380.crop();
		IJ.run("Conversions...", "scale"); // Scale the values in the 16- to 8-bit conversion
		IJ.run(mask380, "8-bit", "");
		IJ.run("Conversions...", " ");
		mask380.show();
		new WaitForUserDialog("Feature selection", "Use the 'lines' tool to determine the dimension of the biggest feature to be extracted").show();
		IJ.run(mask380, "Subtract Background...", ""); // Invokes the background subtraction filter
		IJ.run(mask380, "Subtract...", "")// Subtracts a value from all the pixels
		IJ.run(mask380, "Median...", ""); // Invokes a median filter (to smoothen areas belonging to the same feature (e.g. the same ROI)
		Prefs.blackBackground = false;
		IJ.setThreshold(0,0);
		IJ.run(mask380, "Convert to Mask", "");
		IJ.run(mask380, "Open", ""); // Removes small particles
		batchMode = false;
		var trIm = IJ.openImage(inputDir + transmissionImage);
		trIm.show();
		// Asks the user for confirmation of the ROI thresholding
		var sf = new GenericDialog("Thresholding confirmation");
		var choice = ["Yes, proceed", "No, repeat"];
		sf.addMessage("Are you satisfied with the thresholding results?");
		sf.addChoice("", choice, "Yes, proceed");
		sf.showDialog();
		var usChoice = sf.getNextChoice();
		switch (usChoice) {
			case "Yes, proceed":
//			mask380 = WindowManager.getCurrentImage();
			IJ.run("Conversions...", "scale"); // Scale the values in the 16- to 8-bit conversion
			IJ.run(mask380, "16-bit", "");
			IJ.run(mask380, "Multiply...", "value=10000");
			IJ.run("Conversions...", " ");
			IJ.saveAs(mask380, "Tiff", outputDir + fileName + "_BkgMask.tif");
			satMask = true;
			break;
			case "No, repeat":
//			mask380 = WindowManager.getCurrentImage(); // If not, closes the image and starts again
			trIm.changes = false;
			trIm.close();
			mask380.changes = false;
			mask380.close();
			break;
			default:
			IJ.error("Unexpected input in 'Thresholding confirmation' dialog");
		}
		IJ.run("Collect Garbage");
		}
	batchMode = true;
	var stack340 = IJ.openImage(outputDir + file340NoExt + "_FlatFieldStabilized.tif");
	IJ.run(stack340, "16 colors", ""); // Specifies the LUT
	IJ.run(stack340, "Enhance Contrast", "saturated=0.1"); // Sets the color scale for better visualization
	stack340.show();
	stack380.show();
	// [Ca2+]i calculation
	/* Calculates the intracellular concentration of Calcium using the formula
	 * [Ca2+]i = Keff(R - Rmin)/(Rmax - R)")
	 */
	var size340 = stack340.getNSlices();
	var ic = new ImageCalculator();
	var ratio = ic.run("Divide create 32-bit stack", stack340, stack380); // Calculates the ratio "R" between the 340 and 380nm measurements
	stack340.changes = false;
	stack340.close();
	stack380.hide();
	var num = ratio.duplicate();
	num.setTitle("Numerator");
	var denom = ratio.duplicate();
	denom.setTitle("Denominator");
	IJ.run(num, "Subtract...", "value=" + Rmin + " stack"); // Calculates the numerator of the [Ca2+]i formula (R - Rmin)
	IJ.run(denom, "Subtract...", "value=" + Rmax + " stack"); // Together they calculate the denominator
	IJ.run(denom, "Multiply...", "value=-1 stack"); // of the [Ca2+]i formula (Rmax - R)
	num.show();
	denom.show();
// Creates the [Ca2+]i stack
	var ic2 = new ImageCalculator();
	var freeCa = ic2.run("Divide create 32-bit stack", num, denom);
	IJ.run(freeCa, "Multiply...", "value=" + Keff + " stack");
	IJ.run("Conversions...", " "); 
	IJ.run(freeCa, "16-bit", "");
	num.changes = false;
	num.close();
	denom.changes = false;
	denom.close();
// Converts it into a hyperstack (with time frames) and calibrates pixel size and time interval
	IJ.run(freeCa, "Set Scale...", "distance="+ pixDist + " known=" + umDist + " pixel=1 unit=" + Scale);
	freeCa.show();
	var ic3 = new ImageCalculator();
	var freeCa2 = ic3.run("Subtract create stack", freeCa, mask380);
	freeCa2 = HyperStackConverter.toHyperStack(freeCa2, 1, 1, size340, "default", "Grayscale");
	IJ.run(freeCa2, "Properties...", "channels=1 slices=1 frames="+ size340 + 
	" unit=" + Scale + " pixel_width=" + umDist/pixDist + " pixel_height=" + umDist/pixDist + 
	" voxel_depth=0 frame=[" + deltaT + " " + tScale + "]");
	freeCa.changes = false;
	freeCa.close();
	freeCa2.setDisplayRange(0, 400);
	IJ.run(freeCa2, "royal", "");
	freeCa2.setTitle("Intracellular free Calcium (nM)");
	freeCa2.show();
	var sf2 = new NonBlockingGenericDialog("Background masking");
		var choice = ["Yes, proceed", "No, repeat"];
		sf2.addMessage("Are you satisfied with the background mask result?");
		sf2.addChoice("", choice, "Yes, proceed");
		sf2.showDialog();
		var usChoice2 = sf2.getNextChoice();
		switch (usChoice2) {
			case "Yes, proceed":
			satCa = true;
			break;
			case "No, repeat":
			freeCa2 = WindowManager.getCurrentImage(); // If not, closes the image and starts again
			freeCa2.changes = false;
			freeCa2.close();
			mask380.changes = false;
			mask380.close();
			trIm.changes = false;
			trIm.close();
			break;
			default:
			IJ.error("Unexpected input in 'Thresholding confirmation' dialog");
		}
		IJ.run("Collect Garbage");
		}
	IJ.saveAs(freeCa2, "Tiff", outputDir + fileName + "_FreeCa.tif"); 
	freeCa2.close();
	stack380.changes = false;
	stack380.close();
	mask380.changes = false;
	mask380.close();
	trIm.changes = false;
	trIm.close();
	batchMode = false;
	IJ.run("Collect Garbage");
	return true;
}

// ROI selection
function segmIm(inputDir, outputDir, stack340name, transmissionImage) {
	batchMode = true;
	var qualityCheck = new Boolean(false); // This value determines when the analysis is finished (see later)
	// Open the image files
	var file340NoExt = stack340name.replace(/\.[^/.]+$/, "");
	var fileName = file340NoExt.slice(0, -6);
	var trIm = IJ.openImage(inputDir + transmissionImage);
	var caIm = IJ.openImage(outputDir + fileName + "_FreeCa.tif");
	trIm.show();
	caIm.show();
// Starts the dialog for the selection of the ROIs
IJ.setTool("freehand");
//new RoiManager(true);
var manager = RoiManager.getRoiManager();
manager.runCommand(caIm,"Show All with labels");
new WaitForUserDialog("Calcium imaging ROI selection", "Select the ROIs for data analysis:\n" +
"- use the freehand selection tool to select a ROI\n" +
"- press 't' to add the ROI to the ROI manager\n" +
"- repeat for all ROIs\n" +
"- press OK when ready").show();
// These settings are required to plot correctly the Z-axis profile

manager.runCommand("Associate", "false");
manager.runCommand("Centered", "false");
manager.runCommand("UseNames", "false");
manager.runCommand("Show All");
// This loop gives the possibility to the user make modifications to the ROI selection
while (qualityCheck == false) {
	manager.runCommand("Deselect");
	var nROIs = manager.getCount();
	IJ.run("Set Measurements...", "mean redirect=None decimal=3");
	var caVal = manager.multiMeasure(caIm); // Measures the mean gray value of each ROI through the hyperstack
	// Plots the fluorescence profile throughout the hyperstack
	var plot = new Plot("Intracellular Free Calcium", "Time (s)", "[Ca2+]i (nM)");
	var N = caIm.getNFrames();
	var frNr = Array.apply(null, Array(N)).map(function (_, i) {return i*deltaT;}); // Sets the x-axis numeric array 
	for (var i = 0; i < nROIs; i++) {
		var roiF = caVal.getColumnAsDoubles(i);
		plot.add("line", frNr, roiF);
	}
	plot.setLimits(NaN, NaN, NaN, NaN);
	batchMode = false;
	var plotWin = plot.show();
	// Asks for confirmation about the quality of the ROI selection process
	var qc = new NonBlockingGenericDialog("ROI measurement confirmation");
	var choice = ["Yes, confirm", "No, make manual changes"];
	qc.addMessage("Are you satisfied with the overall result?");
	qc.addChoice("", choice, "Yes, confirm");
	qc.showDialog();
	var qcChoice = qc.getNextChoice();
	switch (qcChoice) {
		case "Yes, confirm": // To save the ROI values without further changes
		qualityCheck = true;
		break;
		case "No, make manual changes": // If manual changes are needed
		caVal.reset();
		IJ.selectWindow("Intracellular Free Calcium");
		IJ.run("Close");
		// Attention: ROIs must be manually added to the ROI manager, an overlay is not enough!
		new WaitForUserDialog("Make manual changes", "- delete/add ROIs to the ROI manager using the manager commands\n" + 
		"- select a ROI and go to Stacks > 'Plot Z-axis profile' to plot a single ROI\n" +
		"- press OK when ready").show();
		qualityCheck = false;
		break;
		default: // Default error message in case of weird behavior
		IJ.error("Unexpected input in 'Make manual changes' dialog");
	}
}
var imPlot = new ImagePlus();
imPlot = plot.getImagePlus();
IJ.saveAs(imPlot, "Tiff", outputDir + fileName + "_RoiPlot.tif");
imPlot.close();
manager.runCommand("Deselect");
manager.runCommand("Save", outputDir + fileName + "_RoiSet.zip");
manager.runCommand(caIm,"Delete");
manager.close();
caIm.changes = false;
caIm.close();
trIm.chnages = false;
trIm.close();
batchMode = false;
IJ.run("Collect Garbage");
return true;
}

// Extract and save results
function resultsMaker(fileDir, stack340name, stack380name) {
	batchMode = true;
	IJ.run("Set Measurements...", "mean redirect=None decimal=3");
	var file340NoExt = stack340name.replace(/\.[^/.]+$/, "");
	var file380NoExt = stack380name.replace(/\.[^/.]+$/, "");
	var fileName = file340NoExt.slice(0, -6);
	IJ.open(fileDir + fileName + "_RoiSet.zip");
	var stack340 = IJ.openImage(fileDir + file340NoExt + "_FlatFieldStabilized.tif");
	stack340.show();
	IJ.run("From ROI Manager", "");
	var manager = RoiManager.getInstance();
	manager.runCommand("Associate", "false");
	manager.runCommand("Centered", "false");
	manager.runCommand("UseNames", "false");
	manager.runCommand("Show All");
	var nROIs = manager.getCount();
	var rt = manager.multiMeasure(stack340);
	rt.save(fileDir + file340NoExt + "_Results.csv");
	manager.runCommand(stack340,"Deselect");
	IJ.run("Clear Results");
	stack340.changes = false;
	stack340.close();
	var stack380 = IJ.openImage(fileDir + file380NoExt + "_FlatFieldStabilized.tif");
	stack380.show();
	IJ.run("From ROI Manager", "");
	var manager = RoiManager.getInstance();
	manager.runCommand("Associate", "false");
	manager.runCommand("Centered", "false");
	manager.runCommand("UseNames", "false");
	manager.runCommand("Show All");
	var nROIs = manager.getCount();
	var rt = manager.multiMeasure(stack380);
	rt.save(fileDir + file380NoExt + "_Results.csv");
	manager.runCommand(stack380,"Deselect");
	IJ.run("Clear Results");
	stack380.changes = false;
	stack380.close();
	var stackCa = IJ.openImage(fileDir + fileName + "_FreeCa.tif");
	stackCa.show();
	IJ.run("From ROI Manager", "");
	var manager = RoiManager.getInstance();
	manager.runCommand("Associate", "false");
	manager.runCommand("Centered", "false");
	manager.runCommand("UseNames", "false");
	manager.runCommand("Show All");
	var nROIs = manager.getCount();
	var rt = manager.multiMeasure(stackCa);
	rt.save(fileDir + fileName + "_CaValResults.csv");
	manager.runCommand(stackCa,"Deselect");
	IJ.run("Clear Results");
	manager.runCommand(stackCa,"Delete");
	manager.close();
	stackCa.changes = false;
	stackCa.close();
	IJ.run("Collect Garbage");
	batchMode = false;
	return true;
}