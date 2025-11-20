// This macro will skeletonize the MAP2 image to get the branch length


// WARNING: Install "hiPNAT" plugin for the "Summarize Skeleton" function

imageList = getList("image.titles");
	
	for (i = 0; i < imageList.length; i++) {
    	title = imageList[i];
    	
    	
    // apply on C2 image, assuming  MAP2 is C2. !!! Change here if the MAP2 is in another channel.
    if (indexOf(title, "C2") != -1) {
        selectWindow(title);
        setOption("BlackBackground", true); // white is foreground
        run("Skeletonize");   
		run("Analyze Particles...", "size=0.05-Infinity show=Masks display clear summarize overlay");
		run("Summarize Skeleton");
    }
}
