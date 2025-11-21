// This Macro will set a threshold to each of the three open images.


// !!! Adjust this threshold to the mean threshold that you calculated
C3_threshold = 19373; // assuming presynapse is in C3
C2_threshold = 21308; // assuming MAP2 is in C2
C1_threshold = 23221; // assuming postsynapse is in C1

list = getList("image.titles");

// create static list for looping
copyList = newArray(list.length);
for (i = 0; i < list.length; i++) {
    copyList[i] = list[i];
}

for (i = 0; i < copyList.length; i++) {
    selectWindow(copyList[i]);
    title = getTitle();

    // determine threshold based on filename
    if (indexOf(title, "C1") != -1) {
        setThreshold(C1_threshold, 65535); 
    } else if (indexOf(title, "C2") != -1) {
        setThreshold(C2_threshold, 65535);
    } else if (indexOf(title, "C3") != -1) {
        setThreshold(C3_threshold, 65535);
    } 
    setOption("BlackBackground", true);
    selectWindow(copyList[i]);
    run("Convert to Mask");
}  
