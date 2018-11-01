% Fourier Analyis 
% Algortim :
% 1. Patch 2.Grad of the Patch 3. Fourier Transform 4.Low pass filtering 
%5.Removal of Low amplitude signals 6.Non Maximal supression 7. Alignment 
% 6. Count the Peaks 7.Reconstruct the image
function fa(im)
[M,N,K]=size(im);
if(K>0)
    im=rgb2gray(im);
end

% Convert to proper uint8 class type 
im=im2uint8(im);

% Compute the gradient of the image 
[Gmag,Gdir]=imgradient(im,'central');

% Fourier transform of the image  & Applying Low pass filters 
h= fspecial('gaussian',3,2);

PQ = paddedsize(size(Gmag));
F = fft2(double(Gmag),PQ(1),PQ(2));
H = fft2(double(h),PQ(1),PQ(2));
Flowpass = H.*F;
Flpmag=abs(Flowpass);
Flph=angle(Flowpass);

% Removing Low Amplitude signals 
thresh=multithresh(Flpmag);
Flpmag(Flpmag<=thresh*0.1)=0; 
Fimage=Flpmag;

% Computing Inverse Image

Frecon=abs(Fimage).*exp(1i.*Flph);
ffi = ifft2(Frecon);
ffi= ffi(2:size(im,1)+1, 2:size(im,2)+1);

% Applying Non Maximal Supression 
Him=nonmaxsup(ffi);

figure
imshow(abs(ffi),[])
title('Inverse Reconstruction');
figure
imshow(fftshift(log(Fimage+1)),[])
title('Low pass frouerir');
figure
imshow(Him,[])
title('Non Max supresssion');

end

function PQ = paddedsize(AB, CD)
%PADDEDSIZE Computes padded sizes useful for FFT-based filtering.
%   PQ = PADDEDSIZE(AB), where AB is a two-element size vector,
%   computes the two-element size vector PQ = 2*AB.
if nargin == 1
   PQ = 2*AB;
elseif nargin == 2 & ~ischar(CD)
   PQ = AB + CD - 1;
   PQ = 2 * ceil(PQ / 2);
elseif nargin == 2
   m = max(AB); % Maximum dimension.

   % Find power-of-2 at least twice m.
   P = 2^nextpow2(2*m);
   PQ = [P, P];
elseif nargin == 3
   m = max([AB CD]); %Maximum dimension.
   P = 2^nextpow2(2*m);
   PQ = [P, P];
else
   error('Wrong number of inputs.')
end
end

function edgeMap = nonmaxsup(im)
        [mag,dir]=mag_dir(im);
        [x, y] = size(mag);
        edgeMap = mag;
        % For-loop to iterate through size of image, with 1 pixel padding
        for i=2:x-1
            for j=2:y-1
                % Find co-ordinates of the next and prev pixels depening
                % on the direction of the gradient
                switch(dir(i,j))
                    case 0
                        prevPixel = mag(i,j-1);
                        nextPixel = mag(i,j+1);
                    case 45
                        prevPixel = mag(i-1,j-1);
                        nextPixel = mag(i+1,j+1);
                    case 90
                        prevPixel = mag(i-1,j);
                        nextPixel = mag(i+1,j);
                    case 135
                        prevPixel = mag(i+1,j-1);
                        nextPixel = mag(i-1,j+1);
                    otherwise
                        error('Invalid Direction');
                end
                % If either neighbours are bigger, turn off candidate pixel
                if(mag(i,j) <= prevPixel || mag(i,j) < nextPixel)
                    edgeMap(i,j) = 0;
                end
            end
        end
        imshow(im2uint8(edgeMap));
end
    
 function [mag, dir] = mag_dir(fimg)
        fimg=double(fimg);
        fimg = padarray(fimg,[1,1], 'replicate','post');

        % Applies matrix values to get dy/dx
        dy = fimg(:,2:end) - fimg(:,1:end-1);
        dy = padarray(dy,[0,1], 'replicate','post');
        dx = fimg(2:end,:) - fimg(1:end-1,:);
        dx = padarray(dx,[1,0], 'replicate','post');

        dir = atan2d(dx,dy);
        mag = sqrt((dy.^2) + (dx.^2));

        % Applies roundUp function to the direction matrix
        dir = roundUp(dir);
    end

% Rounds up the angle to either 0, 45, 90, or 135
    function dir = roundUp(dir)
        dir(dir < 22.5 & dir >= -22.5) = 0;
        dir(dir < 67.5 & dir >= 22.5 | dir >= -67.5 & dir < -22.5) = 45;
        dir(dir < 112.5 & dir >= 67.5 | dir >= -112.5 & dir < -67.5) = 90;
        dir(dir < 157.5 & dir >= 112.5 | dir >= -157.5 & dir < -112.5) = 135;
        dir(dir >= 157.5 | dir < -112.5) = 0;
    end
