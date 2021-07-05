/* modified version of SixFittingMethod Nov,26th 2012 
Modified by Marty*/
import ij.process.ImageProcessor;
import ij.ImagePlus;
import ij.*;

public class sixFittingMethod1 {

	ImagePlus imp;

	public static int[][] imageMatrix;

	public double[] theFitter(ImageProcessor ip, int iROI) {

		imageMatrix = ip.getIntArray();
		setEqualToImageMatrix(imageMatrix);
		int imgWidth = imageMatrix.length;
		double abest = obtainMax(imageMatrix);
		double zbest = estimateNoise(imageMatrix);
		double xbest = 0.0;
		double ybest = 0.0;
		double wxbest = 1.0;
		double wybest = 1.0;

		// div = 4 corresponds to adding bestvalue + (-50%, -25%, 0, 25%,
		// +50%)*bestvalue for i=(-2,-1,0,1,2)
		abest = abest - zbest; // abest is relative not absolute height
		double div = 4.00;
		double da = 1.00 / div;
		double dz = 1.00 / div;
		double dx = 4.00 / div; // xbest-2 to xbest+2 (-2 to 2 pixels)
		double dy = 4.00 / div;
		double dwx = 2.00 / div; // wxbest=2: {1 to 5}, wxbest=1: {0.5 to 2.5},
		double dwy = 2.00 / div;

		// Use best estimate parameters (zbest, abest, xbest...) to start fits

		int dimension = 5;
		double[][][][][][] fit1 = new double[dimension][dimension][dimension][dimension][dimension][dimension];

		for (int i = 0; i < dimension; i++) {
			for (int j = 0; j < dimension; j++) {
				for (int k = 0; k < dimension; k++) {
					for (int l = 0; l < dimension; l++) {
						for (int m = 0; m < dimension; m++) {
							for (int n = 0; n < dimension; n++) {

								fit1[i][j][k][l][m][n] = fitGaussian(zbest
										+ (i - 2) * dz * zbest, abest + (j - 2)
										* da * abest, xbest + (k - 2) * dx,
										ybest + (l - 2) * dy, wxbest + (m - 1)
												* dwx * wxbest, wybest
												+ (n - 1) * dwy * wybest);

							}// end n
						}// end m
					}// end for l
				}// end for k
			}// end for j
		}// end for i

		int[] valuesFit1 = locateSmallest6D(fit1);

		double z1 = (valuesFit1[0] - 2) * dz * zbest + zbest;
		double a1 = (valuesFit1[1] - 2) * da * abest + abest;
		double x1 = (valuesFit1[2] - 2) * dx + xbest;
		double y1 = (valuesFit1[3] - 2) * dy + ybest;
		double wx1 = (valuesFit1[4] - 1) * dwx * wxbest + wxbest;
		double wy1 = (valuesFit1[5] - 1) * dwy * wybest + wybest;

		String ztext = Double.toString(z1);
		String atext = Double.toString(a1);
		String xtext = Double.toString(x1);
		String ytext = Double.toString(y1);
		String wxtext = Double.toString(wx1);
		String wytext = Double.toString(wy1);
		String iROItext = Integer.toString(iROI);
		IJ.log(iROItext + ", " + ztext + ", " + atext + ", " + xtext + ", "
				+ ytext + ", " + wxtext + ", " + wytext);

		div = 2 * div; // double denominator: divide twice as fine each time
		da = 1 / div;
		dz = 1 / div;
		dx = 4 / div;
		dy = 4 / div;
		dwx = 2 / div;
		dwy = 2 / div;

		double[][][][][][] fit2 = new double[dimension][dimension][dimension][dimension][dimension][dimension];
		for (int ii = 0; ii < dimension; ii++) {
			for (int jj = 0; jj < dimension; jj++) {
				for (int kk = 0; kk < dimension; kk++) {
					for (int ll = 0; ll < dimension; ll++) {
						for (int mm = 0; mm < dimension; mm++) {
							for (int nn = 0; nn < dimension; nn++) {

								// if(z1==0.0){z1=z1-2*dz*z1;}

								fit2[ii][jj][kk][ll][mm][nn] = fitGaussian(z1
										+ (ii - 2) * dz * z1, a1 + (jj - 2)
										* da * a1, x1 + (kk - 2) * dx, y1
										+ (ll - 2) * dy, wx1 + (mm - 1) * dwx
										* wx1, wy1 + (nn - 1) * dwy * wy1);
							}// end nn
						}// end mm
					}// end for ll
				}// end for kk
			}// end for jj
		}// end for ii

		int[] valuesFit2 = locateSmallest6D(fit2);

		double z2 = (valuesFit2[0] - 2) * dz * z1 + z1;
		double a2 = (valuesFit2[1] - 2) * da * a1 + a1;
		double x2 = (valuesFit2[2] - 2) * dx + x1;
		double y2 = (valuesFit2[3] - 2) * dy + y1;
		double wx2 = (valuesFit2[4] - 1) * dwx * wx1 + wx1;
		double wy2 = (valuesFit2[5] - 1) * dwy * wy1 + wy1;

		ztext = Double.toString(z2);
		atext = Double.toString(a2);
		xtext = Double.toString(x2);
		ytext = Double.toString(y2);
		wxtext = Double.toString(wx2);
		wytext = Double.toString(wy2);
		IJ.log(iROItext + ", " + ztext + ", " + atext + ", " + xtext + ", "
				+ ytext + ", " + wxtext + ", " + wytext);

		div = 2 * div; // double denominator: divide twice as fine each time
		da = 1 / div;
		dz = 1 / div;
		dx = 4 / div;
		dy = 4 / div;
		dwx = 2 / div;
		dwy = 2 / div;

		double[][][][][][] fit3 = new double[dimension][dimension][dimension][dimension][dimension][dimension];

		for (int iii = 0; iii < dimension; iii++) {
			for (int jjj = 0; jjj < dimension; jjj++) {
				for (int kkk = 0; kkk < dimension; kkk++) {
					for (int lll = 0; lll < dimension; lll++) {
						for (int mmm = 0; mmm < dimension; mmm++) {
							for (int nnn = 0; nnn < dimension; nnn++) {

								fit3[iii][jjj][kkk][lll][mmm][nnn] = fitGaussian(
										z2 + (iii - 2) * dz * z2, a2
												+ (jjj - 2) * da * a2, x2
												+ (kkk - 2) * dx, y2
												+ (lll - 2) * dy, wx2
												+ (mmm - 1) * dwx * wx2, wy2
												+ (nnn - 1) * dwy * wy2);

							}// end nn
						}// end mm
					}// end for ll
				}// end for kk
			}// end for jj
		}// end for ii

		int[] valuesFit3 = locateSmallest6D(fit3);

		double z3 = (valuesFit3[0] - 2) * dz * z2 + z2;
		double a3 = (valuesFit3[1] - 2) * da * a2 + a2;
		double x3 = (valuesFit3[2] - 2) * dx + x2;
		double y3 = (valuesFit3[3] - 2) * dy + y2;
		double wx3 = (valuesFit3[4] - 1) * dwx * wx2 + wx2;
		double wy3 = (valuesFit3[5] - 1) * dwy * wy2 + wy2;

		ztext = Double.toString(z3);
		atext = Double.toString(a3);
		xtext = Double.toString(x3);
		ytext = Double.toString(y3);
		wxtext = Double.toString(wx3);
		wytext = Double.toString(wy3);
		IJ.log(iROItext + ", " + ztext + ", " + atext + ", " + xtext + ", "
				+ ytext + ", " + wxtext + ", " + wytext);

		div = 2 * div; // double denominator: divide twice as fine each time
		da = 1 / div;
		dz = 1 / div;
		dx = 4 / div;
		dy = 4 / div;
		dwx = 2 / div;
		dwy = 2 / div;

		double[][][][][][] fit4 = new double[dimension][dimension][dimension][dimension][dimension][dimension];

		for (int iiii = 0; iiii < dimension; iiii++) {
			for (int jjjj = 0; jjjj < dimension; jjjj++) {
				for (int kkkk = 0; kkkk < dimension; kkkk++) {
					for (int llll = 0; llll < dimension; llll++) {
						for (int mmmm = 0; mmmm < dimension; mmmm++) {
							for (int nnnn = 0; nnnn < dimension; nnnn++) {

								fit4[iiii][jjjj][kkkk][llll][mmmm][nnnn] = fitGaussian(
										z3 + (iiii - 2) * dz * z3, a3
												+ (jjjj - 2) * da * a3, x3
												+ (kkkk - 2) * dx, y3
												+ (llll - 2) * dy, wx3
												+ (mmmm - 1) * dwx * wx3, wy3
												+ (nnnn - 1) * dwy * wy3);
							}// end nn
						}// end mm
					}// end for ll
				}// end for kk
			}// end for jj
		}// end for ii

		int[] valuesFit4 = locateSmallest6D(fit4);

		double z4 = (valuesFit4[0] - 2) * dz * z3 + z3;
		double a4 = (valuesFit4[1] - 2) * da * a3 + a3;
		double x4 = (valuesFit4[2] - 2) * dx + x3;
		double y4 = (valuesFit4[3] - 2) * dy + y3;
		double wx4 = (valuesFit4[4] - 1) * dwx * wx3 + wx3;
		double wy4 = (valuesFit4[5] - 1) * dwy * wy3 + wy3;

		ztext = Double.toString(z4);
		atext = Double.toString(a4);
		xtext = Double.toString(x4);
		ytext = Double.toString(y4);
		wxtext = Double.toString(wx4);
		wytext = Double.toString(wy4);
		IJ.log(iROItext + ", " + ztext + ", " + atext + ", " + xtext + ", "
				+ ytext + ", " + wxtext + ", " + wytext);

		div = 2 * div; // double denominator: divide twice as fine each time
		da = 1 / div;
		dz = 1 / div;
		dx = 4 / div;
		dy = 4 / div;
		dwx = 2 / div;
		dwy = 2 / div;

		double[][][][][][] fit5 = new double[dimension][dimension][dimension][dimension][dimension][dimension];

		for (int iiiii = 0; iiiii < dimension; iiiii++) {
			for (int jjjjj = 0; jjjjj < dimension; jjjjj++) {
				for (int kkkkk = 0; kkkkk < dimension; kkkkk++) {
					for (int lllll = 0; lllll < dimension; lllll++) {
						for (int mmmmm = 0; mmmmm < dimension; mmmmm++) {
							for (int nnnnn = 0; nnnnn < dimension; nnnnn++) {

								fit5[iiiii][jjjjj][kkkkk][lllll][mmmmm][nnnnn] = fitGaussian(
										z4 + (iiiii - 2) * dz * z4, a4
												+ (jjjjj - 2) * da * a4, x4
												+ (kkkkk - 2) * dx, y4
												+ (lllll - 2) * dy, wx4
												+ (mmmmm - 1) * dwx * wx4, wy4
												+ (nnnnn - 1) * dwy * wy4);

							}// end nn
						}// end mm
					}// end for ll
				}// end for kk
			}// end for jj
		}// end for ii

		int[] valuesFit5 = locateSmallest6D(fit5);

		double z5 = (valuesFit5[0] - 2) * dz * z4 + z4;
		double a5 = (valuesFit5[1] - 2) * da * a4 + a4;
		double x5 = (valuesFit5[2] - 2) * dx + x4;
		double y5 = (valuesFit5[3] - 2) * dy + y4;
		double wx5 = (valuesFit5[4] - 1) * dwx * wx4 + wx4;
		double wy5 = (valuesFit5[5] - 1) * dwy * wy4 + wy4;

		ztext = Double.toString(z5);
		atext = Double.toString(a5);
		xtext = Double.toString(x5);
		ytext = Double.toString(y5);
		wxtext = Double.toString(wx5);
		wytext = Double.toString(wy5);
		IJ.log(iROItext + ", " + ztext + ", " + atext + ", " + xtext + ", "
				+ ytext + ", " + wxtext + ", " + wytext);

		div = 2 * div; // double denominator: divide twice as fine each time
		da = 1 / div;
		dz = 1 / div;
		dx = 4 / div;
		dy = 4 / div;
		dwx = 2 / div;
		dwy = 2 / div;

		double[][][][][][] fit6 = new double[dimension][dimension][dimension][dimension][dimension][dimension];

		for (int iiiiii = 0; iiiiii < dimension; iiiiii++) {
			for (int jjjjjj = 0; jjjjjj < dimension; jjjjjj++) {
				for (int kkkkkk = 0; kkkkkk < dimension; kkkkkk++) {
					for (int llllll = 0; llllll < dimension; llllll++) {
						for (int mmmmmm = 0; mmmmmm < dimension; mmmmmm++) {
							for (int nnnnnn = 0; nnnnnn < dimension; nnnnnn++) {

								fit6[iiiiii][jjjjjj][kkkkkk][llllll][mmmmmm][nnnnnn] = fitGaussian(
										z5 + (iiiiii - 2) * dz * z5, a5
												+ (jjjjjj - 2) * da * a5, x5
												+ (kkkkkk - 2) * dx, y5
												+ (llllll - 2) * dy, wx5
												+ (mmmmmm - 1) * dwx * wx5, wy5
												+ (nnnnnn - 1) * dwy * wy5);

							}// end nn
						}// end mm
					}// end for ll
				}// end for kk
			}// end for jj
		}// end for ii

		int[] valuesFit6 = locateSmallest6D(fit6);

		double z6 = (valuesFit6[0] - 2) * dz * z5 + z5;
		double a6 = (valuesFit6[1] - 2) * da * a5 + a5;
		double x6 = (valuesFit6[2] - 2) * dx + x5;
		double y6 = (valuesFit6[3] - 2) * dy + y5;
		double wx6 = (valuesFit6[4] - 1) * dwx * wx5 + wx5;
		double wy6 = (valuesFit6[5] - 1) * dwy * wy5 + wy5;

		ztext = Double.toString(z6);
		atext = Double.toString(a6);
		xtext = Double.toString(x6);
		ytext = Double.toString(y6);
		wxtext = Double.toString(wx6);
		wytext = Double.toString(wy6);

		IJ.log(iROItext + ", " + ztext + ", " + atext + ", " + xtext + ", "
				+ ytext + ", " + wxtext + ", " + wytext);
		IJ.log("");

		double chiSquared = fitGaussian(z6, a6, x6, y6, wx6, wy6);
		double[] fitParameters = new double[7];
		fitParameters[0] = z6;
		fitParameters[1] = a6;
		fitParameters[2] = x6;
		fitParameters[3] = y6;
		fitParameters[4] = wx6;
		fitParameters[5] = wy6;
		fitParameters[6] = chiSquared;

		return fitParameters;

	}// end run method

	// method to estimate the noise by taking the average of the pixel at the
	// edge of the image
	public double estimateNoise(int[][] matrix) {
		int width = matrix.length;
		System.out.println("");

		int noiseCount = 0;
		for (int i = 0; i < width; i++) {

			noiseCount = noiseCount + matrix[i][0] + matrix[i][width - 1];
		}

		for (int j = 1; j < width - 1; j++) {
			noiseCount = noiseCount + matrix[0][j] + matrix[width - 1][j];
		}

		double averageNoise = (double) noiseCount / ((width - 1) * 4);
		for (int i = 0; i < width; i++) {
			System.out.println("matrix " + Integer.toString(matrix[i][0]));
		}

		System.out.println("AverageNoise " + Double.toString(averageNoise));
		return averageNoise;
	}// end

	public static void setEqualToImageMatrix(int[][] matrix) {
		imageMatrix = matrix;
	}

	public static int[][] getImageMatrix() {
		return imageMatrix;
	}

	public double obtainMax(int[][] matrix) {
		int maxValue = 0;
		for (int i = 0; i < matrix.length; i++) {
			for (int j = 0; j < matrix.length; j++) {
				if (matrix[i][j] > maxValue) {
					maxValue = matrix[i][j];
				}
			}
		}
		double maximum = maxValue;
		return maximum;

	}

	public double fitGaussian(double noise, double amplitude, double centreX,
			double centreY, double majorAxis, double minorAxis) {
		int[][] imgMatrix = getImageMatrix();// imgToMatrix(imp);
		int imgwidth = imgMatrix.length;
		double z = noise;
		double a = amplitude;
		double x0 = centreX;
		double y0 = centreY;
		double w1 = majorAxis;
		double w2 = minorAxis;

		// <<< creating the two independant varibles x and y >>>//
		double[] x = new double[imgwidth];
		double[] y = new double[imgwidth];

		for (int ii = 0; ii < imgwidth; ii++) {
			x[ii] = (-imgwidth / 2) + ii;
			y[ii] = (-imgwidth / 2) + ii;
		}// end for ii

		// <<< creating the guassian Matrix>>>//
		double[][] gaussian = new double[imgwidth][imgwidth];
		double[][] difference = new double[imgwidth][imgwidth];
		double sum = 0; // sum of the difference squared
		double chi;
		for (int i = 0; i < imgwidth; i++) {
			for (int j = 0; j < imgwidth; j++) {
				gaussian[i][j] = z
						+ a
						* Math.exp(-0.5 * Math.pow(((x[i] - x0) / w1), 2) - 0.5
								* Math.pow(((y[j] - y0) / w2), 2));
			}// end for j

		}// end for i

		// <<< subtracting the guassian matrix from the image matrix >>>//
		for (int q = 0; q < imgwidth; q++) {
			for (int r = 0; r < imgwidth; r++) {
				difference[q][r] = Math.pow((imgMatrix[q][r] - gaussian[q][r]),
						2) / gaussian[q][r];
				sum = sum + difference[q][r];
			}// end for r
		}// end for q
		chi = sum / (Math.pow(imgwidth, 2) - 6);
		return chi;// return the sum of differences squared.
	}

	public static int[] locateSmallest6D(double[][][][][][] a) {
		int[] location = { 0, 0, 0, 0, 0, 0 };
		double min = a[0][0][0][0][0][0];

		// what comes after my for statements
		for (int row = 0; row < a.length; row++) {
			for (int column = 0; column < a[row].length; column++) {
				for (int height = 0; height < a[row][column].length; height++) {
					for (int step = 0; step < a[row][column][height].length; step++) {
						for (int deep = 0; deep < a[row][column][height][step].length; deep++) {
							for (int vall = 0; vall < a[row][column][height][step][deep].length; vall++) {

								if (a[row][column][height][step][deep][vall] < min) {
									location[0] = row;
									location[1] = column;
									location[2] = height;
									location[3] = step;
									location[4] = deep;
									location[5] = vall;

									min = a[row][column][height][step][deep][vall];
								} // end if

							}// end vall
						}// end deep
					} // end step
				} // end height
			} // end column
		} // end row
		return location;
	}// end locateSmallest6D Mehod

}// end finalGuassianFitter class

