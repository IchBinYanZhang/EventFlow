close all;
clear;
clc;

n = 100;
xs = -n:n;
ys = xs;

[X, Y] = meshgrid(xs, ys);
anglesPlot = getColorWheelCoded(Y, X, n);

subplot 211;
imshow(anglesPlot);

subplot 212;
imshow(reshape(hsv(256), [1, 256, 3]));
xlabel('-\pi to \pi');



