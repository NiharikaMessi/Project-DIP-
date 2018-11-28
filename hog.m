%% HOG Computation 

function m=hog(im)
im=im2double(im);
F = extractHOGFeatures(im,'CellSize',[16 16],'BlockSize',[4 4]);
m=1000*mean(F(:));
end