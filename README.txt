UCF DATA : https://drive.google.com/open?id=1LOefxfxSRke1DhNYZGLG-C5DrCCkRbTZ


Function & its usage procedure : 

1.finaltokyo(im) 

* Takes input as an image , and gives out a 4*1 vector - which consits of 
Variance,Energy,Skeness & Kurtosis respectively . 

2.Hog(im)
Retrurms the final HOG feature value after detecting the features. 

3.lbpm(im)
Returns the LBP patch's feature vector in decimal form  (instead of structure).


4.SVR(im)

Takes the input and also loads data : Ground data of the UCF DATA ,
along with each approch's output feature vector .

FAD.txt - Fourier Analysis 
CC2.txt - Consists all the text output 
varHog.txt -Consists all the variance of each image(HOG).
meanLbpm.txt - Mean values of LBPM vector .

5. getpatch.m
Gives  series of patches for the input image .
----END-----

5.hog_feature_vector.m 
Consists of whole code of HOG .

6.fa.m
Older version of the fouirer analysis code .
