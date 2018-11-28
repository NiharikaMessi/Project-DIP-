
%% Support vector Regressor 
function J=SVR(im)
%im=rgb2gray(im);
im=im2double(im);
load CC2.txt
load FAD.txt
load meanhog.tx1t
load varlbp.txt

%Each data 
hogdata=meanhog(:,1);
lbpm1=varlbp(:);
fa1=FAD(:,1);
fa2=FAD(:,2);
fa3=FAD(:,3);
fa4=FAD(:,4);
ac=CC2(:,1);

% MOdel For each case 
% Md1=fitrsvm(HogFeatures,CC2);
% Md2=fitrsvm(FAD,CC2);
% Md3=fitrsvm(lbp,CC2);

% All 4 data fitting 

Table=[fa1 fa2 fa3 fa4  lbpm1];
Model=fitrsvm(Table,ac);

% Calculation for Image Parameter 
T1=lbpm(im); % 1 number 
T2=hog(im);    % 1 number 
T3=finaltokyo(im); % 4 cols nos 

Tdata=[T3(1) T3(2) T3(3) T3(4) T2 T1 ];
yfit=predict(Model,Tdata);
J=yfit;

end


%% LBPM 
function M1=lbpm(im)
im=im2double(im);
F=extractLBPFeatures(im);
M1=var(F(:));
end

%% HOG Computation 

function M2=hog(im)
im=im2double(im);
F = extractHOGFeatures(im,'CellSize',[16 16],'BlockSize',[4 4]);
M2=1000*mean(F(:));
end
%% Fourier Analysis 
function Z=finaltokyo(im)
% Input Image patch and find its size 
[M,N]=size(im);
% Convert to proper uint8 class type 
im=im2uint8(im);
im=im2double(im);

% Removing Noise in the image 
im=wiener(im);

% Compute the gradient of the image 
imgrad=grad1(im);

%Fourier Transform of the image gradient  
PQ=paddedsize(size(imgrad));
Ift1= fft2(double(imgrad),PQ(1),PQ(2));
Fmag=abs(Ift1);
Fphase=angle(Ift1);

% Low pass fitering with Gaussian Filter of sigma=3 
Ik=gaussianfil(imgrad);

% Reconstruction of the image 
Frecon=abs(Ik).*exp(1i.*Fphase);
ffi = ifft2(Frecon);
ffi= ffi(2:size(im,1)+1, 2:size(im,2)+1);
Ifinal=abs(ffi);

% Non-Maxima Supression 
Ifnm=nonmaximasup(Ifinal);

%Computing the Variance,Kurtosis,Skewness & Energy of the image patch 
Z=parameters(Ifnm);
end

%% ALL FUNCTIONS USED ABOVE 

function P=grad1(im)
[P,Q]=imgradient(im,'Sobel');
end

%% GAUSSIAN LOW PASS CORRECTED
function Z=gaussianfil(I)
pq=paddedsize(size(I));
h= fspecial('gaussian',19,3);
H = fft2(double(h),pq(1),pq(2));
I = fft2(double(I),pq(1),pq(2));
Z = I.*H;
end

%% Removal of NOISE 
function J=wiener(I)
    m=5;
    n=5;
    J=wiener2(I,[m n]);
end

%% %% Non Maximal Supression
%%Compute the gradient magnitufde based on derivatives in x and y:

    function I=nonmaximasup(im)
    [H,W]=size(im);   
     size_of_kernel = 6*3+1;   
    GRADIENT=zeros(H,W);
    %Derivatives in x and y   
    [derivative_x,derivative_y] = imgradientxy(im,'sobel');   
 
    %%Compute the gradient magnitufde based on derivatives in x and y:
        for r=1+ceil(size_of_kernel/2):H-ceil(size_of_kernel/2)  
            for c=1+ceil(size_of_kernel/2):W-ceil(size_of_kernel/2)  
                GRADIENT(r,c) = sqrt((derivative_x(r,c)^2) + (derivative_y(r,c)^2));
            end
        end   
    non_max =GRADIENT ;
    for r=1+ceil(size_of_kernel/2):H-ceil(size_of_kernel/2) 
        for c=1+ceil(size_of_kernel/2):W-ceil(size_of_kernel/2) 
            %%quantize:
            if (derivative_x(r,c) == 0) 
                tangent = 5;       
            else
                tangent = (derivative_y(r,c)/derivative_x(r,c));
            end
            if (-0.4142<tangent && tangent<=0.4142)
                if(GRADIENT(r,c)<GRADIENT(r,c+1) || GRADIENT(r,c)<GRADIENT(r,c-1))
                    non_max(r,c)=100;
                end
            end
        if (0.4142<tangent & tangent<=2.4142)
            if(GRADIENT(r,c)<GRADIENT(r-1,c+1) || GRADIENT(r,c)<GRADIENT(r+1,c-1))
                non_max(r,c)=100;
            end
        end
        if ( abs(tangent) >2.4142)
            if(GRADIENT(r,c)<GRADIENT(r-1,c) || GRADIENT(r,c)<GRADIENT(r+1,c))
                non_max(r,c)=100;
            end
        end
        if (-2.4142<tangent && tangent<= -0.4142) 
            if(GRADIENT(r,c)<GRADIENT(r-1,c-1) || GRADIENT(r,c)<GRADIENT(r+1,c+1))
                non_max(r,c)=100;
            end
        end
        end
    end
    I=non_max;
    end
    
    %% %% Calculation of Parameters 
    function R=parameters(im)
        [M,N]=size(im);
         Isqr =im.^2;
         im=abs(im2double(im));
         V=sum(var(im(:)));
         E=sum(Isqr(:));
         S =sum(skewness(im(:)));
         K= sum(kurtosis(im(:)));
         R=[V E S K ];
         R=round(R);

    end
    %% %% Padded SIZE 
 function PQ = paddedsize(AB, CD)
%PADDEDSIZE Computes padded sizes useful for FFT-based filtering.
%   PQ = PADDEDSIZE(AB), where AB is a two-element size vector,
%   computes the two-element size vector PQ = 2*AB.
if nargin == 1
   PQ = 2*AB;
elseif nargin == 2 && ~ischar(CD)
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