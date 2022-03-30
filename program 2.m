%%
   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\1.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I1), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
k=BW;
k1=mat2gray(double(k))
figure,imshow(k1),title('double');
 
% calculating centroid  and pixel value

s = regionprops(BW, k1, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k1, {'Centroid','PixelValues','BoundingBox'});
imshow(k1),title('Standard Deviation of Regions');
hold on
for k = 1 : numObj
s(k).StandardDeviation = std(double(s(k).PixelValues));
text(s(k).Centroid(1),s(k).Centroid(2), ...
sprintf('%2.1f', s(k).StandardDeviation), ...
'EdgeColor','b','Color','r');
end
hold off
th1 = finalpolcor(k1)

%% target 1

   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\a.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
k=BW;
k1=mat2gray(double(k))
figure,imshow(k1),title('double');

imshow(k1)
BW = k1 > 0;
s = regionprops(BW, k1, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k1, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th1 = finalpolcor(k1)


%% target 2

   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\b.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 

BW = k3 > 0;
s = regionprops(BW, k3, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k3, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th2=finalpolcor(k3)


%% target 3


   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\c.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
BW = k5 > 0;
s = regionprops(BW, k5, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off


s = regionprops(BW, k5, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th3=finalpolcor(k5)


%% target 4


   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\d.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
k=BW;
k1=mat2gray(double(k))
figure,imshow(k1),title('double');
 
% calculating centroid  and pixel value

BW = k7 > 0;
s = regionprops(BW, k7, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k7, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th4=finalpolcor(k7)

%% target 5


   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\e.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
BW = k9 > 0;
s = regionprops(BW, k9, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k9, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th5=finalpolcor(k9)

%% target 6



   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\f.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
k=BW;
k1=mat2gray(double(k))
figure,imshow(k1),title('double');
 
% calculating centroid  and pixel value

BW = k11 > 0;
s = regionprops(BW, k11, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k11, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th6=finalpolcor(k11)


%% target 7



   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\g.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
BW = k13 > 0;
s = regionprops(BW, k13, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k13, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th7=finalpolcor(k13)

%% target 8

   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\admin\Desktop\h.jpg');
   img=img_orig; %copy of original image
   hsv=rgb2hsv(img);
   h=hsv(:,:,1);
   s=hsv(:,:,2);
   
   [r c v]=find(h>0.25 | s<=0.3 | s>0.9); %non skin
   numid=size(r,1);
   
   for i=1:numid
       img(r(i),c(i),:)=0;
   end
   
   
   imshow(img);
   
   
   %%
   %ycbcr segmentation
   img_ycbcr=img;  %image from the previous segmentation
   ycbcr=rgb2ycbcr(img_ycbcr);
   cb=ycbcr(:,:,2);
   cr=ycbcr(:,:,3);
   
    
    %Detect Skin
    %[r,c,v] = find(cb>=77 & cb<=127 & cr>=133 & cr<=173);
     [r c v] = find(cb<=77 | cb >=127 | cr<=133 | cr>=173);
    numid = size(r,1);
    
    %Mark Skin Pixels
    for i=1:numid
        img_ycbcr(r(i),c(i),:) = 0;
       % bin(r(i),c(i)) = 1;
    end
    
    figure
    title('ycbcr segmentation');
   imshow(img_ycbcr);
   
   
   
  %%
  %rgb segmentation

img_rgb=img_ycbcr;
r=img_rgb(:,:,1);
g=img_rgb(:,:,2);
b=img_rgb(:,:,3);


[row col v]= find(b>0.79*g-67 & b<0.78*g+42 & b>0.836*g-14 & b<0.836*g+44 ); %non skin pixels
numid=size(row,1);

for i=1:numid
    img_rgb(row(i),col(i),:)=0;
end


imshow(img_ycbcr);

imshow(img_rgb);

 %%
%segmentation 

E = entropyfilt(img_ycbcr);
Eim = mat2gray(E);
imshow(Eim);
BW1 = im2bw(Eim, .8);
imshow(BW1);
 imshow(img_ycbcr);
BWao = bwareaopen(BW1,2000);
imshow(BWao);
nhood = true(9);
closeBWao = imclose(BWao,nhood);
imshow(closeBWao)
roughMask = imfill(closeBWao,'holes');
imshow(roughMask);
imshow(img_ycbcr);
I2 = img_ycbcr;
I2(roughMask) = 0;
imshow(I2);
E2 = entropyfilt(I2);
E2im = mat2gray(E2);
imshow(E2im);
BW2 = im2bw(E2im,graythresh(E2im));
imshow(BW2)
imshow(img_ycbcr);
mask2 = bwareaopen(BW2,1000);
imshow(mask2);
texture1 = img_ycbcr;
texture1(~mask2) = 0;
texture2 = img_ycbcr;
texture2(mask2) = 0;

figure, imshow(texture2),title ('segmented');
boundary = bwperim(mask2);
segmentResults = img_ycbcr;
segmentResults(boundary) = 255;
imshow(segmentResults),title ('segmented result')



% filteration for noise removal
h = fspecial('unsharp');
I2 = imfilter(mask2,h);
imshow(mask2), title('Original Image')
figure, imshow(I2), title('Filtered Image')

% binarization
BW = I2 > 0;
figure,imshow(BW)
title('Binary Image');
% converting to double 
k=BW;
k1=mat2gray(double(k))
figure,imshow(k1),title('double');
 
% calculating centroid  and pixel value


BW = k17 > 0;
s = regionprops(BW, k17, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k17, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th9=finalpolcor(k17)

th1
th2
th3
th4
th5
th6
th7
th8
th9





x1=th1-th2
x2=th1-th3
x3=th1-th4
x4=th1-th5
x5=th2-th3
x6=th2-th4
x7=th2-th5
x8=th3-th4
x9=th3-th5
x10=th4-th5

a=x1(2)
b=x2(2)
c=x3(2)
d=x4(2)
e=x5(2)
f=x6(2)
g=x7(2)
h=x8(2)
i=x9(2)
j=x10(2)

A=[a b c d e f g h i j]

B=min(A)