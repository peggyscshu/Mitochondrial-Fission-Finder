/* This is designed to measure the fission signature of mitochondria based on the study of this paper (Nature, 593, p435, 2021).
 * Edited by Shao-Chun, Peggy, Hsu.
 * E-mail peggyschsu@gmail.com
 * 
 * The size of cropped image is defined as xResult and yResult with 200*200 pxs.
 * 
 * Date: 2025/5/6-5/9
 * 1. Demo image structure: Mito matrix, Drp1, biosensor from LSM980 with LSM plus processing.
 * 2. Compatible file formats: .czi, .lsm, .tif.
 * 3. The input file should be saved on the local computer first.
 * 4. Prepare a csv file for computing mito.
 * 5. Analyze fission pairs.
 * 
 * Date: 2025/5/11
 * 1. Paste fission pairs image from different channels.
 * 
 * Date:2025/5/12
 * 1. Confirm the sequential images.
 * 
 * Date: 2025/5/13
 * 1. Mosaic the images.
 * 2. Use "synchronize window to check Ch123_seq.tif"
 * 3. It took 783.871 sec to analyze 621 pairs from the dataset with 8 frames with pixel size 2024*2024. 
 * 
 * Date 2025/5/13A
 * 1. Remove the mito pair with intensity lower than 500 in the mother mito.
 * 
 * Date 2025/6/3
 * 1. Automatically generates a folder to store the intermediate table for computing.
 * 
 * Date 2025/6/27
 * 1. Mute the steps for checking codes.
 * 2. Clear the intermediate data in the hard disc.
 */


setBatchMode("hide");
//Preparation
    //Import image data	
		tifName = getTitle();
		if (endsWith(tifName, ".czi")) {
   			Name = replace(tifName, ".czi", "");
		}
		if (endsWith(tifName, ".lsm")) {
   			Name = replace(tifName, ".lsm", "");
		}
		if (endsWith(tifName, ".tif")) {
   			Name = replace(tifName, ".tif", "");
		}
		dir1 = getDirectory("image");
		getDimensions(width, height, channels, slices, frames);
		BD = bitDepth();
		rename("Raw");
	//Generate a folder to store results
	    ResultPath = dir1 + File.separator + "Results";
		File.makeDirectory(ResultPath);
    //Prepare an intermediate table for computing
    	dir2 = dir1 + File.separator + "MitoIDcsv";
    	File.makeDirectory(dir2);
    	Table.create("MitoID");
    	path2 = dir2 + File.separator + "MitoID.csv";
    	saveAs("Results", dir2 + File.separator + "MitoID.csv");
	//Miscelleneous
		startTime = getTime();
		run("Options...", "iterations=1 count=1 black do=Nothing");
		run("Set Measurements...", "mean redirect=None decimal=1");
	//Create Arrays for data collection
	    DaughterNoArray = newArray();
		ParentLabelArray = newArray();
		DaughterALabelArray = newArray();
		DaughterBLabelArray = newArray();
        DaughterCLabelArray = newArray();
        DaughterDLabelArray = newArray();
//Create a summmary table
		Table.create("Fission analysis");
		Table.setColumn("Daughter No.", DaughterNoArray);
		Table.setColumn("Parent label", ParentLabelArray);
		Table.setColumn("DaughterA label", DaughterALabelArray);
		Table.setColumn("DaughterB label", DaughterBLabelArray);
		Table.setColumn("DaughterC label", DaughterCLabelArray);
		Table.setColumn("DaughterD label", DaughterDLabelArray);

//Generate Mito ID img
	//Generate mito mask with bleach corrected brightness but the same contour
   	 selectWindow("Raw");
   	 run("Duplicate...", "title=Mito duplicate channels=1");
   	 run("Bleach Correction", "correction=[Histogram Matching]");
   	 selectWindow("DUP_Mito");
   	 run("Grays");
   	 rename("BleachCorrectedMito");
 	 run("8-bit");
     run("Median...", "radius=2 stack");
   	 run("Subtract Background...", "rolling=50 stack");
     run("Auto Local Threshold", "method=Phansalkar radius=15 parameter_1=0 parameter_2=0 white stack");
     selectWindow("Mito");
     run("Duplicate...", "title=[Mito_Simple mask] duplicate");
     selectWindow("Mito_Simple mask");
     setAutoThreshold("Default dark");
	 setOption("BlackBackground", false);
	 run("Convert to Mask", "background=Dark black");
	 run("Divide...", "value=255 stack");
   	 imageCalculator("Multiply stack", "BleachCorrectedMito","Mito_Simple mask");
    //Clear
      selectWindow("Mito");
      close();
      selectWindow("Mito_Simple mask");
      close();
   //Identify individual mito
     selectWindow("BleachCorrectedMito");
	 setThreshold(129, 255);
	 run("Analyze Particles...", "size=0.009-Infinity show=[Count Masks] add stack");
	 selectWindow("Count Masks of BleachCorrectedMito");
	 rename("MitoID");
	 run("glasbey on dark");
//Get mito parental and daughter ID
    //Get the mito number at t+1 from the roi at t
	 r = roiManager("count");
     for (i = 0; i < r; i++) {
		roiManager("select", i);
		ParentalLabel = i+1;
		a = getSliceNumber();
		nextT = a+1;
		if(nextT<=frames){
			selectWindow("MitoID");
			roiManager("select", i);
			setSlice(nextT);
            run("Save XY Coordinates...", "save=[path2]");
            open(path2);
            //Get unique Daughter Label
            	selectWindow("MitoID.csv");
            	Array_DaughterLabel = Table.getColumn("Value");
				Array_DaughterLabel = Array.sort(Array_DaughterLabel);
            	//Array.print(Array_DaughterLabel);
     			unique = newArray();
     			unique = Array.concat(unique,Array_DaughterLabel[0]);
				for (k = 0; k < Array_DaughterLabel.length-1; k++) {
					if (Array_DaughterLabel[k+1] != Array_DaughterLabel[k]) {
        				unique = Array.concat(unique,Array_DaughterLabel[k+1]);
    				}
				}
				unique = Array.deleteValue(unique, 0);
				//Array.print(unique);
				uqCount = unique.length;
				
				if (uqCount==2) {
					DaughterNoArray = Array.concat(DaughterNoArray,uqCount);
					ParentLabelArray = Array.concat(ParentLabelArray,ParentalLabel);
					DaughterALabelArray = Array.concat(DaughterALabelArray,unique[0]);
					DaughterBLabelArray = Array.concat(DaughterBLabelArray,unique[1]);
					DaughterCLabelArray = Array.concat(DaughterCLabelArray,0);
					DaughterDLabelArray = Array.concat(DaughterDLabelArray,0);
					selectWindow("Fission analysis");
					Table.setColumn("Daughter No.", DaughterNoArray);
					Table.setColumn("Parent label", ParentLabelArray);
					Table.setColumn("DaughterA label", DaughterALabelArray);	
					Table.setColumn("DaughterB label", DaughterBLabelArray);
					Table.setColumn("DaughterC label", DaughterCLabelArray);
					Table.setColumn("DaughterD label", DaughterDLabelArray);
				}
				if (uqCount==3) {
					DaughterNoArray = Array.concat(DaughterNoArray,uqCount);
					ParentLabelArray = Array.concat(ParentLabelArray,ParentalLabel);
					DaughterALabelArray = Array.concat(DaughterALabelArray,unique[0]);
					DaughterBLabelArray = Array.concat(DaughterBLabelArray,unique[1]);
					DaughterCLabelArray = Array.concat(DaughterCLabelArray,unique[2]);
					DaughterDLabelArray = Array.concat(DaughterDLabelArray,0);
					selectWindow("Fission analysis");
					Table.setColumn("Daughter No.", DaughterNoArray);
					Table.setColumn("Parent label", ParentLabelArray);
					Table.setColumn("DaughterA label", DaughterALabelArray);	
					Table.setColumn("DaughterB label", DaughterBLabelArray);
					Table.setColumn("DaughterC label", DaughterCLabelArray);
					Table.setColumn("DaughterD label", DaughterDLabelArray);
				}
				if (uqCount==4) {
					DaughterNoArray = Array.concat(DaughterNoArray,uqCount);
					ParentLabelArray = Array.concat(ParentLabelArray,ParentalLabel);
					DaughterALabelArray = Array.concat(DaughterALabelArray,unique[0]);
					DaughterBLabelArray = Array.concat(DaughterBLabelArray,unique[1]);
					DaughterCLabelArray = Array.concat(DaughterCLabelArray,unique[2]);
					DaughterDLabelArray = Array.concat(DaughterDLabelArray,unique[3]);
					selectWindow("Fission analysis");
					Table.setColumn("Daughter No.", DaughterNoArray);
					Table.setColumn("Parent label", ParentLabelArray);
					Table.setColumn("DaughterA label", DaughterALabelArray);	
					Table.setColumn("DaughterB label", DaughterBLabelArray);
					Table.setColumn("DaughterC label", DaughterCLabelArray);
					Table.setColumn("DaughterD label", DaughterDLabelArray);
				}
			}
		}

	saveAs("tif", ResultPath + File.separator + Name + "_Mito Label");
	roiManager("Deselect");
	roiManager("Save", ResultPath + File.separator + Name + "_MitoIndex.zip");
	selectWindow("MitoID.csv");
	run("Close");

//Generate result image
	//Prepare the result image stack with the slice number of pair count 
		selectWindow("Fission analysis");
		FissionPairNo = Table.size;
		xResult = 200;
		yResult = 200;
		newImage("Ch1_t-1", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch1_Parent", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch1_Daughters", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch2_t-1", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch2_Parent", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch2_Daughters", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch3_t-1", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch3_Parent", BD + "-bit black", xResult, yResult, FissionPairNo);
		newImage("Ch3_Daughters", BD + "-bit black", xResult, yResult, FissionPairNo);
	//Prepare images for duplication
		selectWindow("Raw");
	    run("Split Channels");
	//Paste pairs to result images
		for (i = 0; i < FissionPairNo; i++) {
			ParentROIindex = parseInt(Table.get("Parent label", i))-1;
			DaughterROIArray = newArray();
			Aindex = parseInt(Table.get("DaughterA label", i))-1;
			Bindex = parseInt(Table.get("DaughterB label", i))-1;
			DaughterROIArray = Array.concat(DaughterROIArray, Aindex);
			DaughterROIArray = Array.concat(DaughterROIArray, Bindex);
			//print(ParentROIindex);
			//Array.print(DaughterROIArray);
		//Locate Parental Mito------------------------------------------------------------------------
			selectWindow("C1-Raw");
			run("Duplicate...", "title=C1-Raw_Draw duplicate");
			selectWindow("C1-Raw_Draw");
			roiManager("select", ParentROIindex);
			sl = getSliceNumber();
			getSelectionBounds(x, y, width, height);
			cenX = x + width/2 ;
			cenY = y + height/2 ;
			setForegroundColor(255, 255, 255);
            run("Draw", "slice");
			//Paste Ch1 before fission
				selectWindow("C1-Raw_Draw");				
				run("Specify...", "width=xResult height=yResult x=cenX y=cenY slice=sl centered");
				run("Copy");
				selectWindow("Ch1_Parent");
				setSlice(i+1);
				run("Paste");
				selectWindow("C1-Raw_Draw");
				close();
			//Paste Ch2 before fission
				selectWindow("C2-Raw");
				run("Duplicate...", "title=C2-Raw_Draw duplicate");
				selectWindow("C2-Raw_Draw");
				roiManager("select", ParentROIindex);
				setForegroundColor(255, 255, 255);
                run("Draw", "slice");
				run("Specify...", "width=xResult height=yResult x=cenX y=cenY slice=sl centered");
				run("Copy");
				selectWindow("Ch2_Parent");
				setSlice(i+1);
				run("Paste");
				selectWindow("C2-Raw_Draw");
				close();
			//Paste Ch3 before fission
				selectWindow("C3-Raw");
				run("Duplicate...", "title=C3-Raw_Draw duplicate");
				selectWindow("C3-Raw_Draw");
				roiManager("select", ParentROIindex);
				setForegroundColor(255, 255, 255);
                run("Draw", "slice");
				run("Specify...", "width=xResult height=yResult x=cenX y=cenY slice=sl centered");
				run("Copy");
				selectWindow("Ch3_Parent");
				setSlice(i+1);
				run("Paste");
				selectWindow("C3-Raw_Draw");
				close();
				//print(sl);
		//Locate Daughter Mito------------------------------------------------------------------------
			selectWindow("C1-Raw");
			run("Duplicate...", "title=C1-Raw_Draw duplicate");
			selectWindow("C1-Raw_Draw");
			roiManager("show none");
			roiManager("select", DaughterROIArray);
			sldaughter = getSliceNumber();
			if (sldaughter != frames) {
					sldaughter = sl +1;
				}
				if (sl == frames ) {
					sldaughter = sldaughter;
				}
			setSlice(sldaughter);
			roiManager("Combine");
			setForegroundColor(255, 255, 255);
            run("Draw", "slice");
			getSelectionBounds(x, y, width, height);
			cenXdaughter = x + width/2 ;
			cenYdaughter = y + height/2 ;
			//Paste Ch1 after fission
				selectWindow("C1-Raw_Draw");	
				run("Specify...", "width=xResult height=yResult x=cenXdaughter y=cenYdaughter slice=sldaughter centered");
				run("Copy");
				selectWindow("Ch1_Daughters");
				setSlice(i+1);
				run("Paste");
				selectWindow("C1-Raw_Draw");
				close();
			//Paste Ch2 after fission
				selectWindow("C2-Raw");
				run("Duplicate...", "title=C2-Raw_Draw duplicate");
				selectWindow("C2-Raw_Draw");
				roiManager("show none");
				roiManager("select", DaughterROIArray);
				setSlice(sldaughter);
				roiManager("Combine");
				setForegroundColor(255, 255, 255);
                run("Draw", "slice");
				run("Specify...", "width=xResult height=yResult x=cenXdaughter y=cenYdaughter slice=sldaughter centered");
				run("Copy");
				selectWindow("Ch2_Daughters");
				setSlice(i+1);
				run("Paste");
				selectWindow("C2-Raw_Draw");
				close();
			//Paste Ch3 after fission
				selectWindow("C3-Raw");
				run("Duplicate...", "title=C3-Raw_Draw duplicate");
				selectWindow("C3-Raw_Draw");
				roiManager("show none");
				roiManager("select", DaughterROIArray);
				setSlice(sldaughter);
				roiManager("Combine");
				setForegroundColor(255, 255, 255);
                run("Draw", "slice");
				run("Specify...", "width=xResult height=yResult x=cenXdaughter y=cenYdaughter slice=sldaughter centered");
				run("Copy");
				selectWindow("Ch3_Daughters");
				setSlice(i+1);
				run("Paste");
				selectWindow("C3-Raw_Draw");
				close();
				//print(sldaughter);
		//Locate Parental Mito at t-1------------------------------------------------------------------------
			//Paste Ch1 before fission t-1 --------------------------------------------------------------------------------------
				selectWindow("C1-Raw");
				run("Duplicate...", "title=C1-Raw_Draw duplicate");
				selectWindow("C1-Raw_Draw");
				roiManager("select", ParentROIindex);
				sl = getSliceNumber();
				if (sl != 1) {
					slbefore = sl -1;
				}
				if (sl == 1 ) {
					slbefore = sl;
				}
				setSlice(slbefore);
				setForegroundColor(255, 255, 255);
                run("Draw", "slice");
				run("Specify...", "width=xResult height=yResult x=cenX y=cenY slice=slbefore centered");
				run("Copy");
				selectWindow("Ch1_t-1");
				setSlice(i+1);
				run("Paste");
				selectWindow("C1-Raw_Draw");
				close();
			//Paste Ch2 before fission t-1 --------------------------------------------------------------------------------------
				selectWindow("C2-Raw");
				run("Duplicate...", "title=C2-Raw_Draw duplicate");
				selectWindow("C2-Raw_Draw");
				roiManager("select", ParentROIindex);
				sl = getSliceNumber();
				if (sl != 1) {
					slbefore = sl -1;
				}
				if (sl == 1 ) {
					slbefore = sl;
				}
				setSlice(slbefore);
				setForegroundColor(255, 255, 255);
                run("Draw", "slice");
				run("Specify...", "width=xResult height=yResult x=cenX y=cenY slice=slbefore centered");
				run("Copy");
				selectWindow("Ch2_t-1");
				setSlice(i+1);
				run("Paste");
				selectWindow("C2-Raw_Draw");
				close();
			//Paste Ch3 before fission t-1 --------------------------------------------------------------------------------------
				selectWindow("C3-Raw");
				run("Duplicate...", "title=C3-Raw_Draw duplicate");
				selectWindow("C3-Raw_Draw");
				roiManager("select", ParentROIindex);
				sl = getSliceNumber();
				if (sl != 1) {
					slbefore = sl -1;
				}
				if (sl == 1 ) {
					slbefore = sl;
				}
				setSlice(slbefore);
				setForegroundColor(255, 255, 255);
                run("Draw", "slice");
				run("Specify...", "width=xResult height=yResult x=cenX y=cenY slice=slbefore centered");
				run("Copy");
				selectWindow("Ch3_t-1");
				setSlice(i+1);
				run("Paste");
				selectWindow("C3-Raw_Draw");
				close();
				//print(slbefore);
		}
		//Adjust B&C
			selectWindow("C1-Raw");//--------------------------
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(0, 19207);
			selectWindow("Ch1_t-1");
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(0, 19207);
			selectWindow("Ch1_Parent");
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(0, 19207);
			selectWindow("Ch1_Daughters");
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(0, 19207);
			selectWindow("Ch2_t-1");//------------------------
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(34, 24364);
			selectWindow("Ch2_Parent");
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(34, 24364);
			selectWindow("Ch2_Daughters");
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(34, 24364);
			selectWindow("Ch3_t-1");//-------------------------
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(116, 1603);
			selectWindow("Ch3_Parent");
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(116, 1603);
			selectWindow("Ch3_Daughters");
			run("Select None");
			run("Rainbow RGB");
			setMinAndMax(116, 1603);
		//Save individual sequence images
			selectWindow("Ch1_t-1");
			saveAs("tif", ResultPath + File.separator + Name + "_Ch1_t-1");
			rename("Ch1_t-1");
			selectWindow("Ch1_Parent");
			saveAs("tif", ResultPath + File.separator + Name + "Ch1_Parent");
			rename("Ch1_Parent");
			selectWindow("Ch1_Daughters");
			saveAs("tif", ResultPath + File.separator + Name + "Ch1_Daughters");
			rename("Ch1_Daughters");
			selectWindow("Ch2_t-1");
			saveAs("tif", ResultPath + File.separator + Name + "_Ch2_t-1");
			rename("Ch2_t-1");
			selectWindow("Ch2_Parent");
			saveAs("tif", ResultPath + File.separator + Name + "Ch2_Parent");
			rename("Ch2_Parent");
			selectWindow("Ch2_Daughters");
			saveAs("tif", ResultPath + File.separator + Name + "Ch2_Daughters");
			rename("Ch2_Daughters");
			selectWindow("Ch3_t-1");
			saveAs("tif", ResultPath + File.separator + Name + "_Ch3_t-1");
			rename("Ch3_t-1");
			selectWindow("Ch3_Parent");
			saveAs("tif", ResultPath + File.separator + Name + "Ch3_Parent");
			rename("Ch3_Parent");
			selectWindow("Ch3_Daughters");
			saveAs("tif", ResultPath + File.separator + Name + "Ch3_Daughters");
			rename("Ch3_Daughters");
		//Mosaic the images
			run("Combine...", "stack1=[Ch1_t-1] stack2=[Ch1_Parent] combine");
            run("Combine...", "stack1=[Combined Stacks] stack2=[Ch1_Daughters] combine");
			selectWindow("Combined Stacks");
			rename("Ch1_seq");
			run("Combine...", "stack1=[Ch2_t-1] stack2=[Ch2_Parent] combine");
            run("Combine...", "stack1=[Combined Stacks] stack2=[Ch2_Daughters] combine");
			selectWindow("Combined Stacks");
			rename("Ch2_seq");
			run("Combine...", "stack1=[Ch3_t-1] stack2=[Ch3_Parent] combine");
            run("Combine...", "stack1=[Combined Stacks] stack2=[Ch3_Daughters] combine");
			selectWindow("Combined Stacks");
			rename("Ch3_seq");
       //Save individual sequence images
			selectWindow("Ch1_seq");
			saveAs("tif", ResultPath + File.separator + "_Ch1_seq_" + Name);
			selectWindow("Ch2_seq");
			saveAs("tif", ResultPath + File.separator + "_Ch2_seq_" + Name);
  			selectWindow("Ch3_seq");
			saveAs("tif", ResultPath + File.separator + "_Ch3_seq_" + Name);
//Save Fission anaylysis result table
	selectWindow("Fission analysis");
	saveAs("Results", ResultPath + File.separator + Name + "_Fission analysis.csv");
//Clear show up
	roiManager("reset");
	run("Close All");
//Clear the hard disc
	File.delete(path2);
	File.delete(dir2);
	File.delete(ResultPath + File.separator + Name + "_Ch1_t-1.tif");
	File.delete(ResultPath + File.separator + Name + "_Ch2_t-1.tif");
	File.delete(ResultPath + File.separator + Name + "_Ch3_t-1.tif");
	File.delete(ResultPath + File.separator + Name + "Ch1_Daughters.tif");
	File.delete(ResultPath + File.separator + Name + "Ch2_Daughters.tif");
	File.delete(ResultPath + File.separator + Name + "Ch3_Daughters.tif");
	File.delete(ResultPath + File.separator + Name + "Ch1_Parent.tif");
	File.delete(ResultPath + File.separator + Name + "Ch2_Parent.tif");
	File.delete(ResultPath + File.separator + Name + "Ch3_Parent.tif");
	
//Report
    endTime = getTime();
    processTime = (endTime-startTime)/1000/60;
    print("It took " + processTime + " min to analyze " + Name + " with " + frames + " frames.");
setBatchMode("show");	
 
	 





		
		
		



