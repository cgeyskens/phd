// This macro will perform a Z-stack max projection, split the channels and apply LUT grey


// get the title of the opened image
imgTitle = getTitle();

// selecting that image
selectImage(imgTitle);

// run Z-projection
run("Z Project...", "projection=[Max Intensity]");

// close the original image
selectImage(imgTitle);
close();

// run split channelss
run("Split Channels");

// get the open images
titles = getList("image.titles");

// change the open images color to grey
for (i = 0; i < titles.length; i++) {
    if (startsWith(titles[i], "C")) {
        selectImage(titles[i]);
        run("Grays"); // Just change LUT, not bit depth
    }
}