// Macro to close everything related to the current image off


// closing all images
while (nImages > 0) {
    selectImage(nImages);
    close();
}

// Closing all summary tables
if (isOpen("Results")) close("Results");
if (isOpen("Summary")) close("Summary");
if (isOpen("Skeleton Stats")) close("Skeleton Stats");
roiManager("reset");