
run("Stack to Hyperstack...", "order=xytzc channels=1 slices=300 frames=7 display=Color");
run("Enhance Contrast", "saturated=0.35");
run("Grouped Z Project...", "projection=[Sum Slices] all");
run("Enhance Contrast", "saturated=0.35");
