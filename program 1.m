%%
   clear all;
  close all;
  clc;
  
   %% 
   %hsv segmentation
   img_orig=imread('C:\Documents and Settings\Admin\Desktop\1.jpg');
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


[row col v]= find(b>0.79 & b<0.78 & b>0.836 & b<0.840 ); %non skin pixels
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

J1=imread('C:\Users\admin\Desktop\a.jpg');
k2=(rgb2gray(J1));
k3=mat2gray(double(k2))
imshow(k3)
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

J2=imread('C:\Users\admin\Desktop\b.jpg');
k4=(rgb2gray(J2));
k5=mat2gray(double(k4))
imshow(k5)
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

J3=imread('C:\Users\admin\Desktop\c.jpg');
k6=(rgb2gray(J3));
k7=mat2gray(double(k6))
imshow(k7)
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

J4=imread('C:\Users\admin\Desktop\d.jpg');
k8=(rgb2gray(J4));
k9=mat2gray(double(k8))
imshow(k9)
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

J5=imread('C:\Users\admin\Desktop\e.jpg');
k10=(rgb2gray(J5));
k11=mat2gray(double(k10))
imshow(k11)
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

J6=imread('C:\Users\admin\Desktop\f.jpg');
k12=(rgb2gray(J6));
k13=mat2gray(double(k12))
imshow(k13)
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


J7=imread('C:\Users\admin\Desktop\g.jpg');
k14=(rgb2gray(J7));
k15=mat2gray(double(k14))
imshow(k15)
BW = k15 > 0;
s = regionprops(BW, k15, {'Centroid','WeightedCentroid'});
hold on
numObj = numel(s);
for k = 1 : numObj
    plot(s(k).WeightedCentroid(1), s(k).WeightedCentroid(2), 'r*');
    plot(s(k).Centroid(1), s(k).Centroid(2), 'bo');
end
hold off
s = regionprops(BW, k15, {'Centroid','PixelValues','BoundingBox'});
hold on
for k = 1 : numObj
    s(k).StandardDeviation = std(double(s(k).PixelValues));
    text(s(k).Centroid(1),s(k).Centroid(2), ...
        sprintf('%2.1f', s(k).StandardDeviation), ...
        'EdgeColor','b','Color','r');
end
hold off
th8=finalpolcor(k15)



J8=imread('C:\Users\admin\Desktop\h.jpg');
k16=(rgb2gray(J8));
k17=mat2gray(double(k16))
imshow(k17)
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