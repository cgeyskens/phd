// This macro count synapses in an entire FOV


imageList = getList("image.titles");
pre_image = "";
post_image = "";

// check the vgat and geph image

for (i = 0; i < imageList.length; i++) {
    title = imageList[i];
    if (indexOf(title, "C1") != -1)
        post_image = title;
    else if (indexOf(title, "C3") != -1)
        pre_image = title;
}

if (post_image == "" || pre_image == "") {
    showMessage("Could not find both C1 and C3 images.");
    exit();
}

// first create the AND image between pre and post
imageCalculator("AND create", post_image, pre_image);
andResult = "Result of " + post_image;

// Pre-synaptic
selectImage(post_image);
run("Analyze Particles...", "size=0.05-Infinity show=Masks display clear summarize overlay");

// Post-synaptic
selectImage(pre_image);
run("Analyze Particles...", "size=0.05-Infinity show=Masks display clear summarize overlay");

// Synapses (pre + post)
selectImage(andResult);
run("Analyze Particles...", "size=0.05-Infinity show=Masks display clear summarize overlay");
