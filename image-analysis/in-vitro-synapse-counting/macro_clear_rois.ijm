// This macro will clear all the information inside all the ROIs in ROI manager of all open images


// !! First draw the ROIs around the somas and save them

// get the open images
titles = getList("image.titles");

// get the number of ROIs in the ROI Manager
roiManager("Select", 0); // this initializes the manager
roiCount = roiManager("count");

// loop through each image
for (i = 0; i < titles.length; i++) {
    if (startsWith(titles[i], "C")) { // restrict to split channels
        selectImage(titles[i]);

        // apply each ROI and clear inside it
        for (j = 0; j < roiCount; j++) {
            roiManager("Select", j);
            run("Clear");
        }
    }
}