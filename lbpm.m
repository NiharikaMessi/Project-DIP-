%% LBPM 
function m=lbpm(im)
im=im2double(im);
F=extractLBPFeatures(im,'CellSize',[32 32],'Normalization','None');
m=var(F(:));
end