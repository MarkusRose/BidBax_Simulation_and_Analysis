/*-------------------------------------------------------------------------------/
 * Plugin for ImageJ to run "SixFitPlugin_v11.class" for all images in directory.
 * Markus Rose
 * Sep 2020
 * -----------------------------------------------------------------------------*/

import java.io.File;

import ij.IJ;
import ij.ImagePlus;
import ij.io.Opener;
import ij.plugin.PlugIn;
import ij.process.ImageProcessor;

public class SixFitFolderAnalysis_ implements PlugIn {


	public void run(String arg0) {
		// TODO Auto-generated method stub
		String outDir = IJ.getDir("myDir");
		IJ.log(IJ.getDir("current"));

		File dir = new File(outDir);
		File[] filelist = dir.listFiles();
		for (int i = 0; i < filelist.length; i++) {
			if (filelist[i].isDirectory()) {
				String imageFilename = filelist[i].getAbsolutePath() +"/"+ filelist[i].getName() + ".bin.counts_ch0.stack.tif";
				File imageFile = new File(imageFilename);
				if (imageFile.isFile()) {
					System.out.println(imageFilename);
					Opener opener = new Opener();  
					ImagePlus imp = opener.openImage(imageFilename);
					ImageProcessor ip = imp.getProcessor(); // ImageProcessor from ImagePlus 
					imp.show();
					try {
					SixFitPlugin_v3 myAnalyzer = new SixFitPlugin_v3();
					myAnalyzer.setup(arg0, imp);
					myAnalyzer.run(ip);
					} catch (ArrayIndexOutOfBoundsException e) {
						continue;
					} catch (NullPointerException e) {
						continue;
					}
				}
			}
		}
		
	}

}
