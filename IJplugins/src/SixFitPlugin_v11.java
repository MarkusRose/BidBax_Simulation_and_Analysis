/*-------------------------------------------------------------------\
Version 3: Aug. 2, 2011
- Amalgamation of versions 1 and 2.
- Order of code reflects the chronological order of analysis steps.
- Variables listed by order of use.

Version 4: Aug. 12, 201110

- shows copy image with detected spots regions erased during analysis loop

Version 5: Aug. 28, 2011
- moved radius (xradius, yradius) for ROI fitting into analysis loop
- Draws boxes around spot, measures total intensity of boxed area until brightest area found
- Sept. 24: outputs data to txt file

Version 6: Oct. 18, 2011
- separates output data into streaks and spots, sorts spots based on wellness of fit
- no cross drawn for streaks
- brightest pixel of spot indicated in blue 

Version7: July 9, 2012
- spots/streaks separated by w0x, w0y values from fit
- streaks sorted good/bad based on gaussian fit 
- no ROIs detected in outer row of pixels

Version8: Sept 10, 2012 - new author Marty Kurylowicz
-automated for command-line batch processing (in Linux)
-all pop-ups removed, parameter settings fixed in code: either blue or red channel
-shell script loops through all tif stacks in a directions, calling a macro which compiles and runs the plugin

SixFitGreenv2: June, 2014 - Modified by Cecile Fradin and Kelly Cathcart
- Now starts at image 1.
- Use pixel size, pixel dwell time and CPP
- Return the real w0
- Classify as spot or streak
- Give real particle position
- Calculate apparent stoichiometry
- Calculate correlation coefficient
To do next:
- Count spots, count streaks, count correlated spots, correlated streaks, and generate a result file (with everEst, noisest).
- Adapt for red to green images

SixFitGreenv3: May-August, 2015 - Modified by Sheldon Winkel
- Corrected error in calculation for noiseEst
	- changed "if (i+k>-1 && i+k<width && j+l>-1 && j+l<height) {noiseEst = noiseEst + mirror[i+k][j+l]/Math.pow(2*(rad+1),2);}"
		   to "if (i+k>-1 && i+k<width && j+l>-1 && j+l<height) {noiseEst = noiseEst + mirror[i+k][j+l]/Math.pow(2*rad+1,2);}"
- Program runs twice, first with the Bax (green) image then the Bid (red) image as the initial image, with the position selected by the green(Bax) image
	- check performed to determine whether bound protein classifications agree in both cases
- Removed unused "draw" functions
- Added if statement to remove unrealistically large ROIs and added test for ellipticity to differentiate between streaks and spots
- Compare results from both runs to see if bound protein classifications agree with each other
- Ripley's K analysis for Bid and Bax proteins with respect to same protein with edge effects correction
- Ripley's K analysis with Bid proteins as point of interest and Bax proteins as surrounding particles, with edge effects correction

SixFitPlugin Version 11: September, 2020 - Modified by Markus Rose
- Set new detection intensity and size thresholds, optimized for Bid and Bax detection. 
- Added functionality to run automatically for all images in directory with file "SixFitFolderAnalysis_.java"
- 
	
\-------------------------------------------------------------------*/ 

import java.awt.Font;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import ij.*;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.*;
import ij.measure.ResultsTable;
import ij.plugin.filter.*;
import ij.process.ImageProcessor;

public class SixFitPluginMarkus_v1 implements PlugInFilter 
	{
	ImagePlus imp;
	ImageStack stack;
	public static String DIR;
	public static int MaxCounter;

/*-------------------------------------------------------------------\ 
Setup method: Getting info about the PlugIn
\-------------------------------------------------------------------*/

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		return DOES_16;
	}

		
	/*-------------------------------------------------------------------\ 
	Run method: Main program
	\-------------------------------------------------------------------*/
	
	@SuppressWarnings("deprecation")
	public void run(ImageProcessor ip) {
	
		for(int run = 1; run <= 2; run++) {		// Program runs twice, first analyzing the (green) Bax image, then the (red) Bid image
		
			/* Important quantities */
			double pixelsize =100;					// Pixel size in nm
			double pixeltime = 1; 					// Pixel dwell time in ms
			
			double CPP = 6/pixeltime; 				// Specific molecular brightness in photon/pixel - CPP for Bax fluorophore
			if(run == 2) {
				CPP = 11.5/pixeltime; 				//CPP for Bid fluorophore
			}
			double wnotuser = 300/pixelsize;        // Expected size of Point Spread Function (PSF) in pixel

			/* Declaration of variables */
	
			int loopcount=0; 	//MK
			int maxcounter=0; 	//MK
			int xint=0;
			int yint=0;
			int xint0=0;
			int yint0=0;
			int w0pix=0;

			int size=imp.getStackSize();			//getting info about the image stack
			int width = ip.getWidth();
			int height = ip.getHeight();
			String title1 = imp.getTitle();
			String title2;
	
			if(run==1){
				title2 = "shifted-0003.tif";	//
			}
			else {
				title2 = "shifted-0001.tif";	//Added for case where Bid/red image is initially analyzed - SW
			}
			
			String dir = IJ.getDirectory("image");
			setDir(dir);
			//IJ.write("size: " + size + " width: " + width + " title: " + title1);

			double SumIm[] = new double[size+1]; 	//variables in loop finding brightest slice
			double MaxSum = 0;	
			
			int position = 0; 				       // slice number of brightest slice
		     		
			double mirror1[][]= new double [width][height];        // array for saving original image data
			double mirror2[][]= new double [width][height];      // array for saving the corresponding red image
			
			double maxint = 0;				// variables for creating copy image with 10X resolution (1000x1000 pixels) 
			double gardpixel = 0;
			int c = 0;

			double minint1 = 4095;				 // variables for finding estimated noise
			double minint2 = 4095;
			double noiseEst1 = 0; 
			double noiseEst2 = 0;
			int noiseCounter1 = 0;
			int noiseCounter2 = 0;

			double pThreshold = 0;				// variable for threshold used in analysis				
			
			double zSum = 0; //MK

			double a  = 0;				// variables for results of fit			
			double z  = 0;
			double x0  = 0;
			double y0  = 0;
			double wx  = 0;
			double wy  = 0;
			double chi  = 0;

			int endloop = 0;				// variable for looping analysis loop 

			int xmaxint = 0;			      	// variables for finding pixels of highest intensity in copy image
			int ymaxint = 0;	
			
			int xmin = 0;					// variables defining ROI for Gaussian fit
			int xmax = 0;
			int ymin = 0;
			int ymax = 0;
			int xmin4 = 0;					// variables defining region to erase from mirror, on spurious maxima
			int xmax4 = 0;
			int ymin4 = 0;
			int ymax4 = 0;

			int MaxROIx=0;				// variables for finding index of intensity maximum in ROI
			int MaxROIy=0;
			int MaxROIint=0;
			int pixvalue=0;
			
			int countcall = 0;
			
			double boundTable1[][] = new double [6][8];		//contains information for bound proteins - Bax protein for run 1, Bid protein for run 2
			double boundTable2[][] = new double [6][8];		//contains information for bound proteins - Bax protein for run 2
			int numberBound = 0;				// number of bound particles found during analysis
			double unboundTable1[][] = new double [101][30];	//number of "spots" bound or unbound - Bax protein for run 1, Bid protein for run 2
			double unboundTable1b[][] = new double [101][30];	//Used to hold Ripley's K data for same proteins (Bax in run 1, Bid in run 2)
			double unboundTable2[][] = new double [101][30];	//number of "spots" bound or unbound - Bax protein for run 2
			int numberUnbound = 0;								//number of unbound "spots"

			int stoploop = 300;                    // Maximum number of objects that can be found in the image (should be 300)

			int wnot = (int) wnotuser;   // Calculation of correlation coefficient
			double iFirst = 0;
			double iSecond = 0;
			double corrX = 0;
			double varFirst = 0;
			double varSecond = 0;
			double chiX = 0;

			int rad = 1;
			
			String conclusion;

			ResultsTable rt = Analyzer.getResultsTable();
			rt = null;

			/* Find the brightest image in the stack  */

			for (int k=2; k<size+1; k++) { // For a fake start to skip bright blobs start at k = 2 
				IJ.write("k = " + k);
				imp.setSlice(k);	    
				double Sum = 0;
				for (int i=0; i<width; i++) {	
					for (int j=0; j<height; j++) {
						Sum = Sum + ip.getPixelValue(i,j);
					}
				}
				SumIm[k]=Sum;
				IJ.write("Slice: "+ k +", Intensity: " + Sum);
				if (Sum>MaxSum) {
					MaxSum = Sum; 
					position = k;		//Change to specific number if position is known but not selected by plugin
				}
			}
			
			imp.setSlice(position);

			String filename = title1 + "slice#" + position;

			/* Copies image to array (the "mirror image") */

			for (int i=0; i<width; i++) { 
				for (int j=0; j<height; j++)  
				{mirror1[i][j]=ip.getPixelValue(i,j);} 
			}


			/* Create a copy image with a 10 times higher resolution */

			ImagePlus Result1_image = NewImage.createRGBImage("Result1", width*10, height*10, 1, NewImage.FILL_RAMP);

			ImageProcessor Result1_ip = Result1_image.getProcessor();

			maxint = 0;
			for (int i=0; i<width; i++) {                           // Find the highest intensity value in the image
				for (int j=0; j<height; j++) {
					if (mirror1[i][j]>maxint) 
					{maxint=mirror1[i][j];}
				} 
			}

			for (int i=0; i<width; i++) {                      // Fill the high resolution colored image with intensity values normalized by the highest intensity 
				for (int j=0; j<height; j++) {
					gardpixel=(int) (ip.getPixelValue(i,j)/maxint*255);
			
					if(run==1){
						c=((0 & 0xff)<<16) | (((int) gardpixel & 0xff)<<8) | 0 & 0xff;      //Green for Bax
					}
					else if(run==2){
						c=(((int) gardpixel & 0xff)<<16) | ((0 & 0xff)<<8) | 0 & 0xff;           //Red for Bid
					}
			
//			c = (int) ip.getPixelValue(i,j);  		                                          // Grayscale
					for (int k=0; k<10; k++) {  		
						for (int l=0; l<10; l++) {
							Result1_ip.putPixel(10*i+k,10*j+l,c);
						}
					}
				}
			}

			Result1_image.show();                                   // Shows the enhanced first image

			double maxmax=maxint;

			/* Open the corresponding slice in the red channel */
	
			IJ.write("Directory: "+ dir +", redtitle: " + title2);	//Troubleshooting null pointer exception in line 281

			ImagePlus Second_image = IJ.openImage(dir + title2);
			ImageProcessor Second_ip;
			Second_ip = Second_image.getProcessor();
			Second_image.setSlice(position);                                  
			for (int i=0; i<width; i++) {                                           // We probably don't need that we can work on the ip directly
				for (int j=0; j<height; j++) {
					mirror2[i][j]=Second_ip.getPixelValue(i,j);
				} 
			}
			Second_image.show();														// Shows the untouched red image

			// Generate a red result image at 10 times the resolution and plot it.

			ImagePlus Result2_image = NewImage.createRGBImage("Result2", width*10, height*10, 1, NewImage.FILL_RAMP);
			ImageProcessor Result2_ip = Result2_image.getProcessor();

			maxint = 0;
			for (int i=0; i<width; i++) {                       // Find the highest intensity value in the image
				for (int j=0; j<height; j++) {
					if (mirror2[i][j]>maxint) {
						maxint=mirror2[i][j];
					}
				} 
			}

			for (int i=0; i<width; i++) {                          // Fill the high resolution colored image with intensity values normalized by the highest intensity 
				for (int j=0; j<height; j++) {
					gardpixel=(int) (Second_ip.getPixelValue(i,j)/maxint*255);
			
					if(run==1){
						c=(((int) gardpixel & 0xff)<<16) | ((0 & 0xff)<<8) | 0 & 0xff;		//Red for Bid
					}
					else if(run==2){
						c=((0 & 0xff)<<16) | (((int) gardpixel & 0xff)<<8) | 0 & 0xff;      //Green for Bax
					}
			
					for (int k=0; k<10; k++) {		
						for (int l=0; l<10; l++) {
							Result2_ip.putPixel(10*i+k,10*j+l,c);
						}
					}
				}
			}

			Result2_image.show();                                   // Shows the enhanced red image


			/* Calculate noise floor and number of pixels with lowest intensity for both images (draw crosses) */
	
			rad = 2;	

			for (int i=0; i<width; i++)
				{ for (int j=0; j<height; j++) { if (mirror1[i][j]<minint1) {minint1=mirror1[i][j]; } } }     // Find the lowest intensity value in the image

			for (int i=0; i<width; i++) {                             // removed -> Draw a cross around those points with low intensity value and calculate an average background value from the area around them
				for (int j=0; j<height; j++) {
					if (mirror1[i][j] == minint1) {					// if there is another pixel with a minimum value
						noiseCounter1++;
						//drawCross(i*10+5, j*10+5, 1, 1, Result_ip);
						for (int k=-rad; k<rad+1; k++){
							for (int l=-rad; l<rad+1; l++) {
								//if the area falls within the boundaries of the image
								if (i+k>-1 && i+k<width && j+l>-1 && j+l<height) {noiseEst1 = noiseEst1 + mirror1[i+k][j+l]/Math.pow(2*rad + 1,2);} 
							}
						}
					}
				}
			}


			noiseEst1=noiseEst1/noiseCounter1;			    // estimated noise from points of low intensity

			pThreshold = noiseEst1 + 0.2*CPP; //everEst;      						// threshold value for analysis  (used to be noiseEst+(0.5*CPP); 
			IJ.saveAs(Result1_image, ".png", dir + filename);


			/* Same for red image */

			for (int i=0; i<width; i++)                              
			{ for (int j=0; j<height; j++) { if (mirror2[i][j]<minint2) {minint2=mirror2[i][j];} } }  // Find the lowest intensity value in the image

			for (int i=0; i<width; i++) {                             // removed -> Draw a cross around those points with low intensity value and calculate an average background value from the area around them
				for (int j=0; j<height; j++) { 
					if (mirror2[i][j] == minint2) {
						noiseCounter2++;
						//drawCross(i*10+5, j*10+5, 1, 1, Result_ip);
						for (int k=-rad; k<rad+1; k++) {
							for (int l=-rad; l<rad+1; l++) {
								if (i+k>-1 && i+k<width && j+l>-1 && j+l<height) {noiseEst2 = noiseEst2 + mirror2[i+k][j+l]/Math.pow(2*rad + 1,2);} 
							}
						}
					}
				}
			}

			noiseEst2=noiseEst2/noiseCounter2;			           // estimated noise from points of low intensity


			/* Save the results window */

			IJ.selectWindow("Results");
			IJ.save(dir + filename + "_intensities.txt");
				
			/* Set values for ROI1 and ROI2 */
	
			int radius1 = 8;  // size of exclusion area from border
			int radius2 = 5;  // size of both ROIs, mirror and export 

			//	w0pix = (int) Math.rint(wnotuser/pixelsize)/2; // minimum number of pixels to remove around maxint region
			w0pix = 1;

			/* Analysis loop=================================================================== */

			while (endloop == 0) {
				loopcount=loopcount+1;
				if (loopcount==stoploop) {
					endloop=1;
					//		IJ.write("Warning! The search as stopped because the maximum allowed number of particles, " + stoploop + ", has been reached!");
				}	
	
				/* Find the ROI(s) - point(s) of highest intensity in the image */
				maxint = 0;
				for (int i=radius1; i<width-radius1; i++) {// skips data within radius1 of edge (larger than radius2 to give buffer for second pass)
					for (int j=radius1; j<height-radius1; j++) {// skips data within radius1 of edge
						if (mirror1[i][j]>maxint) {
							maxint=mirror1[i][j];	// record highest intensity value
							xmaxint=i;ymaxint=j;	// record its coordinates
						}						// compare with next pixel etc.
					}
				}

				if (maxint>pThreshold) { 	// If local max is larger than noise floor (pThreshold = avg pixel intensity + CPP/5, or average image intensity, or other)
							// else endloop = 1 -> no more ROIs in the image

					maxcounter = maxcounter+1;
					setMaxCounter(maxcounter);
					Result1_image.draw();	
					
					/* Determine the boundary of the ROI(s)  */
					xmin = xmaxint - radius2; 
					xmax = xmaxint + radius2;
					ymin = ymaxint - radius2; 
					ymax = ymaxint + radius2;

					/* Export first ROI sub-image  */ 

					ImagePlus MaxROI0_image = NewImage.createByteImage("MaxROI0", (xmax-xmin), (ymax-ymin), 1, NewImage.FILL_RAMP);
					ImageProcessor MaxROI0_ip = MaxROI0_image.getProcessor();	
					MaxROIx=0;
					MaxROIy=0;
					MaxROIint=0;
					for (int i=xmin; i<xmax+1; i++) {
						//for (int i=xmin; i<xmax; i++)
						for (int j=ymin; j<ymax+1; j++) {
							//for (int j=ymin; j<ymax; j++)
							pixvalue = (int) ip.getPixelValue(i,j);
							MaxROI0_ip.putPixel(i-xmin,j-ymin, pixvalue);
							if (pixvalue > MaxROIint) {
								MaxROIint=pixvalue; //Max pixel value in ROI (not mirror) -> due to preivously erased ROIs
								MaxROIx=i;
								MaxROIy=j;
							}
						}
					}
		
					//MaxROI0_image.show();

					if (xmaxint==MaxROIx && ymaxint==MaxROIy) {	//if max point is highest within ROI
         														//else spurious max -> ROI previously analyzed
						IJ.saveAs(MaxROI0_image, ".fits", dir + "MaxROI0/" + filename + ".MAXROI" + (maxcounter));       //BLUE
						//drawmaxpointBlue(xmaxint, ymaxint, Result_ip);
						//drawBoxBlue(xmin, xmax, ymin, ymax, Result_ip);
						//drawspotnumberBlue(maxcounter, xdraw, ydraw, Result_ip);
				
						//call FIRST PASS best fitting parameters from sixFittingMethod		
						countcall = countcall+1;

						sixFittingMethod1 sixFitting0 = new sixFittingMethod1();
						double[] results0=sixFitting0.theFitter(MaxROI0_ip,maxcounter);
						double z0=results0[0];
						double a0=results0[1];
						double x00=results0[2];
						double y00=results0[3];
						double wx0=results0[4];
						double wy0=results0[5];
						//create a guassian matrix from best fitting parameters and save model image
						double[][] gaussianMatrix0=gaussianFunction(xmax-xmin,z0,a0,x00,y00,wx0,wy0);			
						ImagePlus gaussianFunction0_image = NewImage.createShortImage("Model0", (xmax-xmin), (ymax-ymin), 1, NewImage.FILL_RAMP);
						ImageProcessor gaussianFunction0_ip = gaussianFunction0_image.getProcessor();
						for (int ii=0; ii<xmax-xmin; ii++) {   		
							for (int jj=0; jj<ymax-ymin; jj++) {
								gaussianFunction0_ip.putPixelValue(ii,jj,gaussianMatrix0[ii][jj]);
							}
						}			
						//IJ.saveAs(gaussianFunction0_image, ".fits", dir + "Model0/" + filename + ".Model" + (maxcounter));

						// Run SECOND PASS through fitting, using fit x0 and y0 to reset ROI position

						/* Determine the boundary of the SECOND ROI(s)  */

						xint = (int) Math.rint(x00);
						yint = (int) Math.rint(y00);

						xmin = xmaxint + xint - radius2; 
						xmax = xmaxint + xint + radius2;

						ymin = ymaxint + yint - radius2; 
						ymax = ymaxint + yint + radius2;

						ImagePlus MaxROI_image = NewImage.createByteImage("MaxROI", (xmax-xmin), (ymax-ymin), 1, NewImage.FILL_RAMP);
						ImageProcessor MaxROI_ip = MaxROI_image.getProcessor();	
			    
						/* Export SECOND PASS ROI, smaller than first  */

						MaxROIx=0;
						MaxROIy=0;
						MaxROIint=0;
						for (int i=xmin; i<xmax+1; i++) {
							//for (int i=xmin; i<xmax; i++)
							for (int j=ymin; j<ymax+1; j++) {
								//for (int j=ymin; j<ymax; j++)
								pixvalue = (int) ip.getPixelValue(i,j);
								MaxROI_ip.putPixel(i-xmin,j-ymin, pixvalue);
								if (pixvalue > MaxROIint)
								{
									MaxROIint=pixvalue; //Max pixel value in ROI (not mirror)
									MaxROIx=i;
									MaxROIy=j;
								}
							}
						}
						IJ.saveAs(MaxROI_image, ".fits", dir + "MaxROI/" + filename + ".MAXROI" + (maxcounter));
						//RED
						//drawmaxpointRed(xmaxint + xint, ymaxint+yint, Result_ip);
						drawBoxBlue(xmin, xmax, ymin, ymax, Result1_ip);                                         //    Draw blue box to show the last considered ROI
						//drawspotnumberRed(maxcounter, 10*(xmaxint+xint), 10*(ymaxint+yint), Result_ip);
				
						//call SECOND PASS  sixFittingMethod 		
						sixFittingMethod2 sixFitting= new sixFittingMethod2();
						// x0 and y0 are zero again since ROI has been re-centred at {x0,y0}, other parameters come from first fit
						double[] results2=sixFitting.theFitter(MaxROI_ip,maxcounter,a0,z0,0,0,wx0,wy0); 
						z=results2[0];
						a=results2[1];
						x0=results2[2];
						y0=results2[3];
						wx=results2[4];
						wy=results2[5];
						chi=results2[6];
			    
						//Check if dimensions of ROI are reasonable

						//Update pThreshold from running average of fit noise level z	
						zSum = zSum + z;
						pThreshold = zSum/countcall + 0.2*CPP;       // + 0.5 * CPP;

						// draw fit w0x and w0y borders
						xint0 =  (int) Math.rint(x00);                   //SixFit delivers wx=w0x/2
						yint0 = (int) Math.rint(y00);
						xint =  (int) Math.rint(x0);                     //SixFit delivers wx=w0x/2
						yint = (int) Math.rint(y0);
						//if the ROI is greater than the initial box, do not count - SW
						if(radius2 > 0.5*((xmaxint+x0+xint0+2*wx)-(xmaxint+x0+xint0-2*wx)) && radius2 > 0.5*((ymaxint+y0+yint0+2*wy)-(ymaxint+y0+yint0-2*wy))) {

							drawBoxYellow2(xmaxint+x0+xint0-2*wx, xmaxint+x0+xint0+2*wx, ymaxint+y0+yint0-2*wy, ymaxint+y0+yint0+2*wy, Result1_ip);
							drawBoxYellow2(xmaxint+x0+xint0-2*wx, xmaxint+x0+xint0+2*wx, ymaxint+y0+yint0-2*wy, ymaxint+y0+yint0+2*wy, Result2_ip);                                                         // draw fit w0x and w0y borders
							drawmaxpointYellow2(xmaxint+x0+xint0, ymaxint+y0+yint0, Result1_ip);                                // The real center
							drawspotnumberYellow(maxcounter, 10*(xmaxint+xint+xint0), 10*(ymaxint+yint+yint0), Result1_ip);

							// Calculate apparent stoichiometry

							double stoi = 0;
							stoi = a/CPP;

							// Is it a good fit?
							int goodfit = 0;
							if (chi<2) {goodfit = 1;}

							// Is the particle above threshold?

							double realthreshold = 0;
							String thr = "no";

							realthreshold = pThreshold + Math.sqrt(pThreshold);
							if (a > realthreshold) {thr = "yes";}

							// Is the particle a spot? A streak?

							String spot = "Undefined";

							if (wy*2 < wnotuser/2) {spot = "streak";}
							//Check for PSF and ellipticity as well
							if ( Math.sqrt(Math.pow((wx*2-wnotuser),2) + Math.pow((wy*2-wnotuser),2)) < 0.75*wnotuser && (((Math.abs(wx-wy))/(wx+wy)) < 0.3)) {spot = "spot";}

							// Calculate the correlation coefficient for the data around the particle in a box of radius wnotuser. 
							// This is done as in Friaa et al., BRL 2014, except that estimated noise is subtracted.
							// Also, obviously, particle detection is different.
							// As a test, when the same image is given for correlation, chi ~ 1 for all particles.

							wnot = (int) wnotuser;
							iFirst = 0;
							iSecond = 0;
							corrX = 0;
							varFirst = 0;
							varSecond = 0;
							chiX = 0;

							for (int i=xmaxint + xint + xint0-wnot; i<xmaxint + xint + xint0+wnot+1; i++) {   		
								for (int j=ymaxint + yint + yint0 -wnot; j<ymaxint + yint + yint0+wnot+1; j++) {
									iFirst = iFirst + ip.getPixelValue(i,j)/Math.pow((2*wnot+1),2);
									iSecond = iSecond + Second_ip.getPixelValue(i,j)/Math.pow((2*wnot+1),2);
								}
							}

							for (int i=xmaxint + xint + xint0-wnot; i<xmaxint + xint + xint0+wnot+1; i++) {   		
								for (int j=ymaxint + yint + yint0-wnot; j<ymaxint + yint + yint0+wnot+1; j++) {
									corrX = corrX + (ip.getPixelValue(i,j)-iFirst)*(Second_ip.getPixelValue(i,j)-iSecond)/Math.pow((2*wnot+1),2);
									varFirst = varFirst +Math.pow((ip.getPixelValue(i,j)-(iFirst-noiseEst1)),2)/Math.pow((2*wnot+1),2);
									varSecond = varSecond +Math.pow((Second_ip.getPixelValue(i,j)-(iSecond-noiseEst2)),2)/Math.pow((2*wnot+1),2);
								}
							}
			    
							chiX = corrX/Math.sqrt(varFirst-iFirst)/Math.sqrt(varSecond-iSecond);
							String bound = "free";
							if (chiX > 0.5) {
								bound = "bound";
			    	
								numberBound = numberBound + 1;
								boundTable1[numberBound][1] = run;
								boundTable1[numberBound][2] = maxcounter;
								boundTable1[numberBound][3] = xmaxint+x0+xint0;
								boundTable1[numberBound][4] = ymaxint+y0+yint0;
								if(spot == "streak") {
									boundTable1[numberBound][5] = 1;
									boundTable1[numberBound][6] = 0;
								}
								if(spot == "spot") {
									boundTable1[numberBound][5] = 0;
									boundTable1[numberBound][6] = 1;
								}					    
							}
			    
							if (spot == "spot") {// && bound == "free") {	//save "free spots" in a different table
								numberUnbound = numberUnbound + 1;
								unboundTable1[numberUnbound][1] = run;
								unboundTable1[numberUnbound][2] = maxcounter;
								unboundTable1[numberUnbound][3] = xmaxint+x0+xint0;
								unboundTable1[numberUnbound][4] = ymaxint+y0+yint0;
								if(spot == "streak") {//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
									unboundTable1[numberUnbound][5] = 1;
									unboundTable1[numberUnbound][6] = 0;
								}
								if(spot == "spot") {
									unboundTable1[numberUnbound][5] = 0;
									unboundTable1[numberUnbound][6] = 1;
								}
							}
			    				    
							conclusion = spot + "-" + bound;
			    
							//save fitting parameters into a results table
							if (rt == null) {
								rt = new ResultsTable();
								Analyzer.setResultsTable(rt);
			    				}
							rt.incrementCounter();
							rt.addValue("ROI", maxcounter);
							rt.addValue("z", z);
							rt.addValue("a", a);
							rt.addValue("stoichiometry", stoi);
							rt.addValue("xpos", xmaxint+x0+xint0);     // Absolute position in original image
							rt.addValue("ypos", ymaxint+y0+yint0);  
							rt.addValue("x0", x0);            		
							rt.addValue("y0", y0);	
							rt.addValue("Chi^2", chi);
							rt.addValue("Good?", goodfit);	
							rt.addValue("wx", wx*2*pixelsize);   //  Now give the real w0 in nm
							rt.addValue("wy", wy*2*pixelsize);
							rt.addValue("streak or spot?", spot);
							rt.addValue("Correlation Coefficient", chiX);
							rt.addValue("free or bound?", bound);
							rt.addValue("iFirst", iFirst);
							rt.addValue("iSecond", iSecond);
							rt.addValue("pthresh", pThreshold);
							rt.addValue("realthreshold", realthreshold);
							rt.addValue("Above threshold?", thr);
							rt.addValue("Shape-Binding", conclusion);
							rt.show("Results");
			    
							//create a guassian matrix from best fitting parameters and save model image
							double[][] gaussianMatrix=gaussianFunction(xmax-xmin,z,a,x0,y0,wx,wy);			
							ImagePlus gaussianFunction_image = NewImage.createShortImage("Model", (xmax-xmin), (ymax-ymin), 1, NewImage.FILL_RAMP);
							ImageProcessor gaussianFunction_ip = gaussianFunction_image.getProcessor();
							for (int ii=0; ii<xmax-xmin; ii++) {
								for (int jj=0; jj<ymax-ymin; jj++) {
									gaussianFunction_ip.putPixelValue(ii,jj,gaussianMatrix[ii][jj]);
								}
							}		
							IJ.saveAs(gaussianFunction_image, ".fits", dir + "Model/" + filename + ".Model" + (maxcounter));
			    
						}

						// Erase the area corresponding to the particle in the mirror image

						for (int i=xmin; i<xmax+1; i++) {	       // xmin = xmaxint + x00 - radius2
							for (int j=ymin; j<ymax+1; j++) {
                                if (mirror1.length <= i || mirror1[i].length <= j){ continue; }
                                if (i < 0 || j < 0){continue;}
								mirror1[i][j]=0;
							}
						}

					}             // end if: intensity is max within ROI

					else    {   //else there is a spurious maximum elswhere in ROI

						//		   IJ.write("Spurious maximum!");
						xmin4 = xmaxint - w0pix;
						xmax4 = xmaxint + w0pix;
						ymin4 = ymaxint - w0pix;
						ymax4 = ymaxint + w0pix;
		    
						//drawmaxpointBlue(xmaxint, ymaxint, Result_ip);	
						//drawBoxBlack(xmin4, xmax4, ymin4, ymax4, Result_ip);    // draw box for blackout
						//			IJ.write("xmin4: " + xmin4 + ", xmax4: " + xmax4);
						//			IJ.write("ymin4: " + ymin4 + ", ymax4: " + ymax4);
						for (int i=xmin4; i<xmax4; i++) {
							for (int j=ymin4; j<ymax4; j++)	{
								mirror1[i][j]=0;
							}
						}
					} //end if else: int is not maximum within ROI

				} //end if maxint      

				else {
					endloop = 1;
					//		IJ.write("Maxint small: " + maxint);
				}
			}

			//show mirror (erased areas)

			ImagePlus Erase_image = NewImage.createRGBImage("Erased spots", width*10, height*10, 1, NewImage.FILL_RAMP);
			ImageProcessor Erase_ip = Erase_image.getProcessor();

			for (int i=0; i<width; i++)	{   		
				for (int j=0; j<height; j++) {
					gardpixel=(int) (mirror1[i][j]/maxmax*255);
					c=((0 & 0xff)<<16) | (((int) gardpixel & 0xff)<<8) | 0 & 0xff;	
					for (int k=0; k<10; k++) {  		
						for (int l=0; l<10; l++) {
							Erase_ip.putPixel(10*i+k,10*j+l,c);
						}
					}
				}
			}

	
			Erase_image.show();
			IJ.saveAs(Erase_image, ".png", dir + filename + ".erasedareas");
  
			IJ.selectWindow("Results");
			IJ.save(dir + filename + ".txt");
			rt = new ResultsTable();
			Analyzer.setResultsTable(rt);	//clear columns
			//	rt.saveAs(dir + filename + "2.csv");
			//int li = rt.getCounter();;
			//for (int i = li; i>-2; i--) {rt.deleteRow(i);}
			rt.show("Results");
			
			IJ.write("Directory: " + dir);
			IJ.write("First stack name: " + title1);
			IJ.write("Second stack name: " + title2);
			IJ.write("Stack size: " + size);
			IJ.write("Brightest Stack : " + position);
			IJ.write("Image width: " + width);
			IJ.write("Image height: "+ height);
			IJ.write("Estimated noise: " + noiseEst1);
			IJ.write("Estimated noise in red image: " + noiseEst2);
			IJ.write("Pixel dwell time (ms): " + pixeltime);
			IJ.write("Pixel size (nm): " + pixelsize);
			IJ.write("CPP (kHz): " + CPP);
			
			String filename2 = (title2 + "slice#" + position);
		
			IJ.selectWindow("Results");	
			IJ.save(dir + filename + "_info.txt");
			//li = rt.getCounter();;
			//for (int i = li; i>-2; i--) {rt.deleteRow(i);}
			rt.show("Results");
	
			/*double boundTable1[][] = new double [numberBound][6];				//Table to compare if "bound" classifications agree, assumes less than 10 bound particles
			for(int i = 0; i < numberBound; i++) {
				for(int j = 0; j <= 5; j++) {
				boundTable1[i][j] = boundTable[i][j];
				}
			}*/
	
			//save information for bound proteins in file
			rt = new ResultsTable();		//clears the columns, separate file for bound particles
			Analyzer.setResultsTable(rt);
			for(int i = 1; i <= numberBound; i++) {
				rt.incrementCounter();
				rt.addValue("run", boundTable1[i][1]);			  
				rt.addValue("ROI", boundTable1[i][2]);
				rt.addValue("xpos", boundTable1[i][3]);
				rt.addValue("ypos", boundTable1[i][4]);
				rt.addValue("streak", boundTable1[i][5]);     // Absolute position in original image
				rt.addValue("spot", boundTable1[i][6]);
				rt.show("Results");
			}
			IJ.selectWindow("Results");	
			IJ.save(dir + filename + "_bound_info.txt");
	
			//save information for unbound "spots" in file
			rt = new ResultsTable();		//clears the columns, separate file for bound particles
			Analyzer.setResultsTable(rt);
			for(int i = 1; i <= numberUnbound; i++) {
				rt.incrementCounter();
				rt.addValue("run", unboundTable1[i][1]);			  
				rt.addValue("ROI", unboundTable1[i][2]);
				rt.addValue("xpos", unboundTable1[i][3]);
				rt.addValue("ypos", unboundTable1[i][4]);
				rt.addValue("streak", unboundTable1[i][5]);     // Absolute position in original image
				rt.addValue("spot", unboundTable1[i][6]);
				rt.show("Results");
			}
			IJ.selectWindow("Results");	
			IJ.save(dir + filename + "_unbound_info.txt");
			//} //end if
	
			//read bound particles from previous run, save in boundTable2 and compare
			if(run == 2) {
		
				BufferedReader br = null;
		 
				try {	//table for bound values
 
					String sCurrentLine;
					int lineNumber = 0;

					br = new BufferedReader(new FileReader(dir + filename2 + "_bound_info.txt"));
					while ((sCurrentLine = br.readLine()) != null) {
						lineNumber++;
						String[] numbers = sCurrentLine.split("\\s+");
						double[] answer = new double[numbers.length];
						//System.out.println(sCurrentLine);
				
						if(lineNumber > 1) {	//exclude first line

							for (int i = 0; i < numbers.length; i++)
								answer[i] = Double.parseDouble(numbers[i]);
							//System.out.println(Arrays.toString(answer));
							//System.out.println(answer.length);
							for(int j = 1; j <= 6; j++) {
								boundTable2[lineNumber-1][j] = answer[j];
							}
						}
					}
 
				}
		
				catch (IOException e) {
					e.printStackTrace();
				} finally {
					try {
						if (br != null)br.close();
					} catch (IOException ex) {
						ex.printStackTrace();
					}
				}	//end buffered reader
		
				try {	//table for unbound values
			 
					String sCurrentLine;
					int lineNumber = 0;
 
					/*br = new BufferedReader(new FileReader(dir + filename2 + "_bound_info.txt"));
					while ((sCurrentLine = br.readLine()) != null) {	//count number of lines
					numberOfLines++;
					}*/
			
					br = new BufferedReader(new FileReader(dir + filename2 + "_unbound_info.txt"));
					while ((sCurrentLine = br.readLine()) != null) {
						lineNumber++;
						String[] numbers = sCurrentLine.split("\\s+");
						double[] answer = new double[numbers.length];
						//System.out.println(sCurrentLine);
				
						if(lineNumber > 1) {	//exclude first line

							for (int i = 0; i < numbers.length; i++)
								answer[i] = Double.parseDouble(numbers[i]);
							//System.out.println(Arrays.toString(answer));
							//System.out.println(answer.length);
							for(int j = 1; j <= 6; j++) {
								unboundTable2[lineNumber-1][j] = answer[j];
							}
						}
					}
 
				}
				catch (IOException e) {
					e.printStackTrace();
				} finally {
					try {
						if (br != null)br.close();
					} catch (IOException ex) {
						ex.printStackTrace();
					}
				}	//end buffered reader
		
				//Compare values of boundTable1 and boundTable2 to determine whether the classified bound particles from one image are the same as the bound particles from the other image
				//Array of bound particles stored in boundTable1 & 2
				//Check if position of bound particle in bT1 matches the position of a bound particle in bT2
		
				for(int i = 1; i <= numberBound; i++) {
					for(int j = 1; j <= numberBound; j++) {
						if((boundTable2[i][3] < (boundTable1[j][3] + 1)) && (boundTable2[i][3] > (boundTable1[j][3] - 1)) && (boundTable2[i][4] < (boundTable1[j][4] + 1)) && (boundTable2[i][4] > (boundTable1[j][4] - 1))) {
							boundTable2[i][7] = boundTable1[j][2];
							boundTable1[j][7] = boundTable2[i][2];
						}
					}
				}
		
			} //end if(run==2) - Now have arrays "boundTable1" and "boundTable2" containing information for bound particles
	

			IJ.saveAs(Result1_image, ".png", dir + filename + ".objectsfound");
			IJ.saveAs(Result2_image, ".png", dir + filename + ".redcorrelation");
	
			//   	IJ.run("Quit"); // necessary for looping through files from shell script (macro won't wait for completion of plugin if IJ is open
	
			//close all images without prompt to save	
			ImagePlus img;									
			while (null != WindowManager.getCurrentImage()) {
				img = WindowManager.getCurrentImage();
				img.changes = false;
				img.close();
			}

			//open second image (red Bid image) for analysis for next loop
			if(run==1) {
				Second_image = IJ.openImage(dir + title2);
				Second_image.show();														// Shows the untouched red image
				imp = Second_image;// Sets as active image
				ip = imp.getProcessor();
			}
	
	
			/*********************************************************************************************/
			//Ripley's K for proteins of the same type
			unboundTable1b = unboundTable1;
			calculateRipleysKb(unboundTable1b);
	
	
			rt = new ResultsTable();		//clears the columns, separate file for bound particles
			Analyzer.setResultsTable(rt);
			for(int i = 1; i <= numberUnbound; i++) {
			    rt.incrementCounter();
			    rt.addValue("run", unboundTable1b[i][1]);			  
			    rt.addValue("ROI", unboundTable1b[i][2]);
			    rt.addValue("xpos", unboundTable1b[i][3]);
			    rt.addValue("ypos", unboundTable1b[i][4]);
			    rt.addValue("streak", unboundTable1b[i][5]);     // Absolute position in original image
			    rt.addValue("spot", unboundTable1b[i][6]);
			    rt.addValue("N1", unboundTable1b[i][8]);
			    rt.addValue("N2", unboundTable1b[i][9]);
			    rt.addValue("N3", unboundTable1b[i][10]);
			    rt.addValue("N4", unboundTable1b[i][11]);
			    rt.addValue("N5", unboundTable1b[i][12]);
			    rt.addValue("N6", unboundTable1b[i][13]);
			    rt.addValue("N7", unboundTable1b[i][14]);
			    rt.addValue("N8", unboundTable1b[i][15]);
			    rt.addValue("N9", unboundTable1b[i][16]);
			    rt.addValue("N10", unboundTable1b[i][17]);
			    rt.addValue("RK 1", unboundTable1b[i][18]);
			    rt.addValue("RK 2", unboundTable1b[i][19]);
			    rt.addValue("RK 3", unboundTable1b[i][20]);
			    rt.addValue("RK 4", unboundTable1b[i][21]);
			    rt.addValue("RK 5", unboundTable1b[i][22]);
			    rt.addValue("RK 6", unboundTable1b[i][23]);
			    rt.addValue("RK 7", unboundTable1b[i][24]);
			    rt.addValue("RK 8", unboundTable1b[i][25]);
			    rt.addValue("RK 9", unboundTable1b[i][26]);
			    rt.addValue("RK 10", unboundTable1b[i][27]);

			    rt.show("Results");	
			    
			    IJ.selectWindow("Results");	
			    IJ.save(dir + filename + "Ripley's_K.txt");
			}

			/***************************************************************************************************/
	
			//Ripley's K for Bid proteins with respect to Bax
			if (run==2) {
				//drawBoxYellow2(xmaxint+x0+xint0-2*wx, xmaxint+x0+xint0+2*wx, ymaxint+y0+yint0-2*wy, ymaxint+y0+yint0+2*wy, Result1_ip);
				calculateRipleysK(unboundTable1, unboundTable2);
			}
	
			if(run==2){
				rt = new ResultsTable();		//clears the columns, separate file for bound particles
				Analyzer.setResultsTable(rt);
				for(int i = 1; i <= numberBound; i++) {
				    rt.incrementCounter();
				    rt.addValue("run1", boundTable2[i][1]);			  
				    rt.addValue("ROI1", boundTable2[i][2]);
				    rt.addValue("xpos1", boundTable2[i][3]);
				    rt.addValue("ypos1", boundTable2[i][4]);
				    rt.addValue("streak1", boundTable2[i][5]);     // Absolute position in original image
				    rt.addValue("spot1", boundTable2[i][6]);
				    rt.addValue("Corresponding ROI1", boundTable2[i][7]);
				    rt.addValue("run2", boundTable1[i][1]);			  
				    rt.addValue("ROI2", boundTable1[i][2]);
				    rt.addValue("xpos2", boundTable1[i][3]);
				    rt.addValue("ypos2", boundTable1[i][4]);
				    rt.addValue("streak2", boundTable1[i][5]);     // Absolute position in original image
				    rt.addValue("spot2", boundTable1[i][6]);
				    rt.addValue("Corresponding ROI2", boundTable1[i][7]);

				    rt.show("Results");	
				    
				    IJ.selectWindow("Results");	
				    IJ.save(dir + filename + "_matched_bound_proteins.txt");
				}
				
				rt = new ResultsTable();		//clears the columns, separate file for bound particles
				Analyzer.setResultsTable(rt);
				for(int i = 1; i <= numberUnbound; i++) {
				    rt.incrementCounter();
				    rt.addValue("run1", unboundTable1[i][1]);			  
				    rt.addValue("ROI1", unboundTable1[i][2]);
				    rt.addValue("xpos1", unboundTable1[i][3]);
				    rt.addValue("ypos1", unboundTable1[i][4]);
				    rt.addValue("streak1", unboundTable1[i][5]);     // Absolute position in original image
				    rt.addValue("spot1", unboundTable1[i][6]);
				    rt.addValue("N1", unboundTable1[i][8]);
				    rt.addValue("N2", unboundTable1[i][9]);
				    rt.addValue("N3", unboundTable1[i][10]);
				    rt.addValue("N4", unboundTable1[i][11]);
				    rt.addValue("N5", unboundTable1[i][12]);
				    rt.addValue("N6", unboundTable1[i][13]);
				    rt.addValue("N7", unboundTable1[i][14]);
				    rt.addValue("N8", unboundTable1[i][15]);
				    rt.addValue("N9", unboundTable1[i][16]);
				    rt.addValue("N10", unboundTable1[i][17]);
				    rt.addValue("RK 1", unboundTable1[i][18]);
				    rt.addValue("RK 2", unboundTable1[i][19]);
				    rt.addValue("RK 3", unboundTable1[i][20]);
				    rt.addValue("RK 4", unboundTable1[i][21]);
				    rt.addValue("RK 5", unboundTable1[i][22]);
				    rt.addValue("RK 6", unboundTable1[i][23]);
				    rt.addValue("RK 7", unboundTable1[i][24]);
				    rt.addValue("RK 8", unboundTable1[i][25]);
				    rt.addValue("RK 9", unboundTable1[i][26]);
				    rt.addValue("RK 10", unboundTable1[i][27]);
				    
				    rt.show("Results");	
				    
				    IJ.selectWindow("Results");	
				    IJ.save(dir + filename + "_Ripley's_K_for_Bid_wrt_Bax.txt");

				}   // show values of boundTable2
			}
	
			/*************************************************************************/
	

	
			if(run == 1) {
				rt = new ResultsTable(); //Clear the results table
				Analyzer.setResultsTable(rt);
				rt.show("Results");
			}
	
		}	//end for(run = 1,2)
	}
	
/*   ------------End of main program----------------- */





public  double[][] gaussianFunction(int matrixLength, double noise, double amplitude, double centreX, double centreY, double majorAxis, double minorAxis){
	
	int imgwidth=matrixLength;
	double z = noise;
	double a = amplitude;
	double x0 = centreX;
	double y0 = centreY;
	double w1 = majorAxis;
	double w2 = minorAxis;
	
	//<<< createing the two independant varibles x and y >>>//
	double[] x = new double[imgwidth];
        double[] y= new double[imgwidth];

        for(int ii= 0; ii< imgwidth; ii++){
       		x[ii]=(-imgwidth/2)+ii;
		y[ii]=(-imgwidth/2)+ii;
        }//end for ii
	

	//<<< creating the guassian Matrix>>>//
	double[][] gaussian= new double[imgwidth][imgwidth];
	
	for(int i =0; i<imgwidth; i++){
	 	for(int j=0; j<imgwidth;j++){
			gaussian[i][j] = z + a * Math.exp(-0.5*Math.pow(((x[i]-x0)/w1),2) -0.5*Math.pow(((y[j]-y0)/w2), 2));
		}//end for j
				
	}//end for i
	return gaussian;
	
 	   }

public static void drawmaxpointYellow2(double xmaxint, double ymaxint, ImageProcessor ip)
	{
	int paintvalue=  ((255 & 0xff)<<16) | ((255 & 0xff)<<8) | 0 & 0xff;	int xx = (int) (xmaxint*10+5);
	int yy = (int) (ymaxint*10+5);
	ip.setValue(paintvalue);
	ip.setLineWidth(5);
	ip.drawDot(xx, yy);
	}

/* Method for drawspotnumberYellow - draws numbers above spots that are found*/
public static void drawspotnumberYellow(int counter, int xpos, int ypos, ImageProcessor ip)
	{
	int paintvalue= ((255 & 0xff)<<16) | ((255 & 0xff)<<8) | 0 & 0xff;
	ip.setValue(paintvalue);
	Font f;
	f = new Font ("SansSerif", Font.BOLD, 22);
	ip.setFont(f);
	ip.drawString(Integer.toString(counter), xpos, ypos+2);
	}

// Method for drawBoxYellow - draws boxes around areas of highest total intensity after gaussian blur
public static void drawBoxYellow(double roixmingard, double roixmaxgard, double roiymingard, double roiymaxgard, ImageProcessor ip)
	{
	int paintvalue= ((255 & 0xff)<<16) | ((255 & 0xff)<<8) | 0 & 0xff;
	ip.setValue(paintvalue);
	ip.setLineWidth(3);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	}

/* Method for drawBoxYellow2 - draws boxes around areas of highest total intensity after gaussian blur AT THE REAL POSITION*/
public static void drawBoxYellow2(double roixmingard, double roixmaxgard, double roiymingard, double roiymaxgard, ImageProcessor ip)
	{
	int paintvalue= ((255 & 0xff)<<16) | ((255 & 0xff)<<8) | 0 & 0xff;
	ip.setValue(paintvalue);
	ip.setLineWidth(3);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	}

// Method for drawBoxBlue - draws boxes around areas of highest total intensity after gaussian blur
public static void drawBoxBlue(double roixmingard, double roixmaxgard, double roiymingard, double roiymaxgard, ImageProcessor ip)
	{
	int paintvalue= ((0 & 0xff)<<16) | ((0 & 0xff)<<8) | 255 & 0xff;
	ip.setValue(paintvalue);
	ip.setLineWidth(1);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmingard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	ip.drawLine(((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymaxgard*10)))+5, ((int)(Math.rint(roixmaxgard*10)))+5, ((int)(Math.rint(roiymingard*10)))+5);
	}			

//set global variable DIR to obtain directory
public static void setDir(String dir){
		DIR=dir;
	}
  
public static String getDir(){                       // This one is unused as far as I can tell. CF.
		return DIR;
	}                         
////////////////////////////////////////////////////////

//set global variable maxCount 
	public static void setMaxCounter(int counter){
	MaxCounter=counter;
	}
	public static int getMaxCounter(){
		return MaxCounter;
	}

//method to save image of the matrix created from the gaussian function
public static void saveImageMatrix(double[][] matrix){
		int maxcounter=getMaxCounter();
		String counterString=Integer.toString(maxcounter);
		double[][] copiedMatrix=matrix;
		int matrixDimension = copiedMatrix.length;
		ImagePlus model_Image=NewImage.createShortImage("Model", matrixDimension, matrixDimension, 1, NewImage.FILL_RAMP);
		ImageProcessor model_ImageProcessor=model_Image.getProcessor();
 		for (int z=0; z<matrixDimension; z++){   		
			for (int x=0; x<matrixDimension; x++){
				model_ImageProcessor.putPixelValue(z,x,copiedMatrix[z][x]);
				}
			}
		 IJ.saveAs(model_Image, ".tif", "/home/nehad/Desktop/test" + "/Model"+counterString);
	}
/*Method for write - writes data to text file*/

protected void write(String dir, String filename, int width, int height, String channel, double MaxSum, int positionuser, double wnotuser, double[][] spotArray, double[][] streakArray, int position, double noiseEst, int counter, int badcounter, int kmax, int linter, int mmax, int ninter, int noisecounter, double everEst, double pThreshold, double CPP, String location)
	{
	try
		{
		BufferedWriter bw = new BufferedWriter(new FileWriter(dir + filename));
		bw.write(filename + "\n");
		bw.write("Image dimensions are " + width + "x" + height + " pixels" + "\n");
		bw.write("The brightest image is slice #" + position + ", Average brightness: " + MaxSum/(width*height) + " ph/pxl" + "\n");
		bw.write("Analysis conditions - channel: " + channel + "     theoretical w0: " + wnotuser + "     CPP: " + CPP + "     threshold set for: " + location + "\n");
		bw.write("The slice analyzed was # " + positionuser + "\n");
		bw.write("average pixel brightness: " + everEst + "\n"); 
		bw.write("noise: " + noiseEst + " ,calculated from " + noisecounter + " points" + "\n");
		bw.write("threshold used: " + pThreshold + "\n"); 
		bw.write("spots found: " + (counter-1) + " streaks founds: " + (badcounter-1) + "\n");
		
		bw.write("fit range intensity: " + -kmax + "-" + kmax + " noise: 0-" + linter + " wnotfitx: " + -mmax + "-" + mmax + " wnotfity: 0-" + ninter + "\n");
		bw.write("N  MaxIntensity(ph/pxl)   xint(pxl)   yint(pxl)   xcenter(pxl)   ycenter(pxl)  fitintensity(ph/pxl)  fitnoise(ph/pxl)  fitw0x(nm)   fitw0y(nm)   totalIntensity(ph) intensity  noise  w0x   w0y   chisquare \n");
					
		bw.write("streaks" + "\n");
		for (int ii = 0; ii<badcounter-1; ii++)
			{
			for (int jj = 0; jj<16; jj++)
				{
				bw.write(Double.toString(streakArray[jj][ii]));
				if (jj<16){ bw.write("  ");}
				}
			bw.write("\n");
			}
				
		bw.write("good streaks" + "\n");
		for (int ii = 0; ii<badcounter-1; ii++)
			{
			if (streakArray[11][ii]>-kmax && streakArray[11][ii]<kmax && streakArray[12][ii]>0 && streakArray[12][ii]<linter && streakArray[9][ii]<=wnotuser)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(streakArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");	
				}
			} 
		bw.write("intensity out of bounds" + "\n");
		for (int ii = 0; ii<badcounter-1; ii++)
			{
			if (streakArray[11][ii]==-kmax || streakArray[11][ii]==kmax)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(streakArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");	
				}
			} 
		bw.write("background out of bounds" + "\n");
		for (int ii = 0; ii<badcounter-1; ii++)
			{
			if (streakArray[12][ii]==0 || streakArray[12][ii]==linter)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(streakArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");
				}
			} 
		bw.write("w0y too large" + "\n");
		for (int ii = 0; ii<badcounter-1; ii++)
			{
			if (streakArray[9][ii]>wnotuser)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(streakArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");
				}
			} 
		bw.write("spots" + "\n");
		for (int ii = 0; ii<counter-1; ii++)
			{
			for (int jj = 0; jj<16; jj++)
				{
				bw.write(Double.toString(spotArray[jj][ii]));
				if (jj<16){ bw.write("  ");}
				}
			bw.write("\n");
			}

		bw.write("good guassian fits" + "\n");
		for (int ii = 0; ii<counter-1; ii++)
			{
			if (spotArray[11][ii]>-kmax && spotArray[11][ii]<kmax && spotArray[12][ii]>0 && spotArray[12][ii]<linter)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(spotArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");	
				}
			} 

		bw.write("intensity out of bounds" + "\n");
		for (int ii = 0; ii<counter-1; ii++)
			{
			if (spotArray[11][ii]==-kmax || spotArray[11][ii]==kmax)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(spotArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");	
				}
			} 

		bw.write("background out of bounds" + "\n");
		for (int ii = 0; ii<counter-1; ii++)
			{
			if (spotArray[12][ii]==0 || spotArray[12][ii]==linter)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(spotArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");
				}
			} 

		bw.write("wnot x out of bounds" + "\n");
		for (int ii = 0; ii<counter-1; ii++)
			{
			if (spotArray[13][ii]==-mmax || spotArray[13][ii]==mmax)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(spotArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");
				}
			} 

		bw.write("wnot y out of bounds" + "\n");
		for (int ii = 0; ii<counter-1; ii++)
			{
			if (spotArray[14][ii]==0 || spotArray[14][ii]==ninter)
				{
				for (int jj = 0; jj<16; jj++)
					{
					bw.write(Double.toString(spotArray[jj][ii]));
					if (jj<16){ bw.write("  ");}
					}
				bw.write("\n");
				}
			} 
		bw.close();
		}
	catch (Exception e) 
		{
		IJ.error("Simple ASCII Writer", e.getMessage());
		return;
		}
	}	

public static void calculateRipleysK(double unboundTable1[][], double unboundTable2[][]) {
	//Calculates Ripley's K for each bound Bid protein, with respect to bound Bax proteins

	int N = 0;
	for(int i = 0; i <= 100; i++) {
		if (unboundTable2[i][6] == 1) {
			N = N + 1;
		}
	}
		double area = 84*84;
		double lambda = N/area;
	
		double count = 0;
		double xpos1 = 0;
		double ypos1 = 0;
		double xpos2 = 0;
		double ypos2 = 0;
		double deltaX = 0;
		double deltaY = 0;
		double distance = 0;
		
		double radius = 0;
		double xBase = 0;
		double yHeight = 0;
		double yBase = 0;
		double xHeight = 0;
		double areaTriangles = 0;
		double areaSquare = 0;
		double angleCircle = 0;
		double totalArea = 0;
		double pi = 3.14159;
		double omega = 0;
		
		
		for(int t = 1; t <= 10; t++) {

				
			for(int i = 1; i <= N; i++) {
				//if blah blah = 0, quit
				count = 0;
				xBase = 0;
				yHeight = 0;
				yBase = 0;
				xHeight = 0;
				areaTriangles = 0;
				areaSquare = 0;
				angleCircle = 0;
				totalArea = 0;
				omega = 0;
				
				xpos1 = unboundTable1[i][3];
				ypos1 = unboundTable1[i][4];
				radius = t*5;
				
				//Determine weight function value for edge correction
				
	/*************Begin if statement to determine weight function value for edge correction***************/
				
				//check if particle is too close to any edge
				if( ((xpos1 - radius) < 8) || ((xpos1 + radius) > 92) || ((ypos1 - radius) < 8) || ((ypos1 + radius) > 92) ) {
					
					//check if particle is too close to left or right side
					if(((xpos1 - radius) < 8) || ((xpos1 + radius) > 92)) {
						
						//check if particle is too close to left side
						if((xpos1 - radius) < 8) {
							xBase = xpos1 - 8;
						}
						
						//check if particle is too close to right side
						if((xpos1 + radius) > 92) {
							xBase = 92 - xpos1;
						}
						
						yHeight = Math.sqrt(Math.pow(radius, 2) - Math.pow(xBase, 2));
					}
					
					//check if particle is too close to the top or bottom
					if(((ypos1 - radius) < 8) || ((ypos1 + radius) > 92)) {
						
						//check if particle is too close to top
						if((ypos1 - radius) < 8) {
							yBase = ypos1 - 8;
						}
					
						//check if particle is too close to bottom
						if((ypos1 + radius) > 92) {
							yBase = 92 - ypos1;
					}
					
					xHeight = Math.sqrt(Math.pow(radius, 2) - Math.pow(yBase, 2));
					}
					
					//check if the portion of the circle outside the boundaries is one piece and extends over two sides
					if((Math.sqrt(Math.pow(xBase,  2) + Math.pow(yBase,  2))) > radius) {
						areaTriangles = 0.5*xBase*yHeight + 0.5*yBase*xHeight;
						areaSquare = xBase*yBase;
						angleCircle = 360 - Math.abs((180/pi)*Math.asin(yHeight/radius)) - Math.abs((180/pi)*Math.asin(xHeight/radius)) - 90;
					}
					
					//otherwise, portion of circle outside the boundaries crosses two sides and are unattached or only one boundary is crossed
					if((Math.sqrt(Math.pow(xBase,  2) + Math.pow(yBase,  2))) < radius) {
						areaTriangles = 2*(0.5*xBase*yHeight) + 2*(0.5*yBase*xHeight);
						areaSquare = 0;
						angleCircle = 360 - Math.abs((180/pi)*2*Math.asin(yHeight/radius)) - Math.abs((180/pi)*2*Math.asin(xHeight/radius));
					}
					
					
					totalArea = areaTriangles + areaSquare + (angleCircle/360)*pi*Math.pow(radius, 2);
					omega = (pi*Math.pow(radius, 2))/totalArea;
				} else {
					omega = 1;
				}

				
	/****************** End if statement for edge correction **************************************************/
				
				
				count = 0;
				
					for(int j = 1; j <= 100; j++) {
						xpos2 = unboundTable2[j][3];
						ypos2 = unboundTable2[j][4];
						deltaX = Math.abs(xpos2 - xpos1);
						deltaY = Math.abs(ypos2 - ypos1);
						distance = Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
			
						if(distance < (radius)) {
							count = count + 1;
						}
						
						if(unboundTable2[j+1][4] == 0) {
							j = 101;
						}
					}
				
						unboundTable1[i][t+7] = omega*count;
						unboundTable1[i][t+17] = (count*omega)/(lambda*pi*Math.pow(radius, 2));
			}
		}
	}
	


public static void calculateRipleysKb(double unboundTable1b[][]) {
	//Calculates Ripley's K for each bound protein with respect to bound proteins of the same type

	int N = 0;
	for(int i = 0; i <= 100; i++) {
		if (unboundTable1b[i][6] == 1) {
			N = N + 1;
		}
	}
		double area = 84*84;
		double lambda = N/area;
		
		double count = 0;
		double xpos1 = 0;
		double ypos1 = 0;
		double xpos2 = 0;
		double ypos2 = 0;
		double deltaX = 0;
		double deltaY = 0;
		double distance = 0;
		double radius = 0;
		
		double xBase = 0;
		double yHeight = 0;
		double yBase = 0;
		double xHeight = 0;
		double areaTriangles = 0;
		double areaSquare = 0;
		double angleCircle = 0;
		double totalArea = 0;
		double pi = 3.14159;
		double omega = 0;
		
		for(int t = 1; t <= 10; t++) {

			for(int i = 1; i <= N; i++) {
				count = 0;
				xBase = 0;
				yHeight = 0;
				yBase = 0;
				xHeight = 0;
				areaTriangles = 0;
				areaSquare = 0;
				angleCircle = 0;
				totalArea = 0;
				omega = 0;
				xpos1 = unboundTable1b[i][3];
				ypos1 = unboundTable1b[i][4];
				radius = t*5;
				
	/*************Begin if statement to determine weight function value for edge correction***************/

				//check if particle is too close to any edge
				if( ((xpos1 - radius) < 8) || ((xpos1 + radius) > 92) || ((ypos1 - radius) < 8) || ((ypos1 + radius) > 92) ) {
					
					//check if particle is too close to left or right side
					if( ((xpos1 - radius) < 8) || ((xpos1 + radius) > 92) ) {
						
						//check if particle is too close to left side
						if((xpos1 - radius) < 8) {
							xBase = xpos1 - 8;
						}
						
						//check if particle is too close to right side
						if((xpos1 + radius) > 92) {
							xBase = 92 - xpos1;
						}
						
						yHeight = Math.sqrt(Math.pow(radius, 2) - Math.pow(xBase, 2));
					}
					
					//check if particle is too close to the top or bottom
					if( ((ypos1 - radius) < 8) || ((ypos1 + radius) > 92) ) {
						
						//check if particle is too close to top
						if((ypos1 - radius) < 8) {
							yBase = ypos1 - 8;
						}
					
						//check if particle is too close to bottom
						if((ypos1 + radius) > 92) {
							yBase = 92 - ypos1;
					}
					
					xHeight = Math.sqrt(Math.pow(radius, 2) - Math.pow(yBase, 2));
					}
					
					//check if the portion of the circle outside the boundaries is one piece and extends over two sides
					if((Math.sqrt(Math.pow(xBase,  2) + Math.pow(yBase,  2))) > radius) {
						areaTriangles = 0.5*xBase*yHeight + 0.5*yBase*xHeight;
						areaSquare = xBase*yBase;
						angleCircle = 360 - Math.abs((180/pi)*Math.asin(yHeight/radius)) - Math.abs((180/pi)*Math.asin(xHeight/radius)) - 90;
					}
					
					//otherwise, portion of circle outside the boundaries crosses two sides and are unattached or only one boundary is crossed
					else {
						areaTriangles = 2*(0.5*xBase*yHeight) + 2*(0.5*yBase*xHeight);
						angleCircle = 360 - Math.abs((180/pi)*2*Math.asin(yHeight/radius)) - Math.abs((180/pi)*2*Math.asin(xHeight/radius));
					}
					
					totalArea = areaTriangles + areaSquare + (angleCircle/360)*pi*Math.pow(radius, 2);
					omega = (pi*Math.pow(radius, 2))/totalArea;
				}
				
				else {
					omega = 1;
				}

	/****************** End if statement for edge correction **************************************************/
				
				count = 0;
				
					for(int j = 1; j <= 100; j++) {
						xpos2 = unboundTable1b[j][3];
						ypos2 = unboundTable1b[j][4];
						if(j != i) {
							deltaX = Math.abs(xpos2 - xpos1);
							deltaY = Math.abs(ypos2 - ypos1);
							distance = Math.sqrt(Math.pow(deltaX, 2) + Math.pow(deltaY, 2));
			
							if(distance < (t*5)) {
								count = count + 1;
							}
						}
						
						if(unboundTable1b[j+1][4] == 0) {
							j = 101;
						}
					}
				
						unboundTable1b[i][t+7] = omega*count;
						unboundTable1b[i][t+17] = (count*omega)/(lambda*pi*Math.pow(radius, 2));
			}
		}
	}
}
	
