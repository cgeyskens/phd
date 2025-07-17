// This macro count synapses in an entire FOV


imageList = getList("image.titles");
vgat_image = "";
geph_image = "";

// check the vgat and geph image

for (i = 0; i < imageList.length; i++) {
    title = imageList[i];
    if (indexOf(title, "C1") != -1)
        geph_image = title;
    else if (indexOf(title, "C3") != -1)
        vgat_image = title;
}

if (geph_image == "" || vgat_image == "") {
    showMessage("Could not find both C1 and C3 images.");
    exit();
}

// first create the AND image between pre and post
imageCalculator("AND create", geph_image, vgat_image);
andResult = "Result of " + geph_image;

// Pre-synaptic (GEPH)
selectImage(geph_image);
run("Analyze Particles...", "size=0.05-Infinity show=Masks display clear summarize overlay");

// Post-synaptic (VGAT) 
selectImage(vgat_image);
run("Analyze Particles...", "size=0.05-Infinity show=Masks display clear summarize overlay");

// Synapses (overlapping areas)
selectImage(andResult);
run("Analyze Particles...", "size=0.05-Infinity show=Masks display clear summarize overlay");
