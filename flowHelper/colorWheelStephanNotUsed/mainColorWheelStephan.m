% parameters for color coded flow visualization
ccopt=struct();
ccopt.ordered=1;
ccopt.overlap=0;
ccopt.stepsize=5;
ccopt.arrowtype=0;
ccopt.usecolor=1;
ccopt.maxspeed=5;


ccode = computeColor2(u,u,ccopt);
imshow(ccode);