# hand_gesture_recognition_using_nntool_

Prerequsities:- MATLAB, Dataset

In this project, we have discussed how hand gestures recognition can be done using 
MATLAB and optimization of results is done using NNTool in MATLAB. NNTool is 
based on Artificial Neural Network concept. We are using camera as a detecting device as 
well as input device for Reality System. One of the most effective of software computing 
techniques is Artificial Neural Networks that has many applications on hand gesture 
recognition problem. The eight target images of different orientations (0 degree, 45 
degree, 90 degree, 135 degree, 180 degree, 225 degree, 270 degree, 315 degree) and input 
samples are captured using image acquisition toolbox and went under skin segmentation 
where skin detection algorithm is applied on them. The image is converted in skin pixel 
and non skin pixels and further the preprocessing is done on the skin detected images like 
edge detection, background removal, noise removal, and binarization. Then feature 
extraction must be done, different methods can be used like geometric features or 
nongeometric features, geometric features that use angles and orientations, palm center.
Non geometric such as color, silhouette and textures, but they are in adequate in recognition.
Here we have used geometric feature extraction using various geometric algorithms. Then 
the geometric feature of input sample images were compared with the eight target images 
of different orientations (0 degree, 45 degree, 90 degree, 135 degree, 180 degree, 225 
degree, 270 degree, 315 degree). Centroid method of geometrical feature extraction is 
used and thus theta value is calculated for each target and input samples. And these theta 
values are given as input and target values for creating 8 different networks giving using 
us the network errors such that we can decide the orientation of the input sample with 
respect to different targets images. The network giving the minimum error is decided as 
49
the best network for the input sample and the target values. Here the second network 
created gives us the minimum error and thus we can conclude that the input image falls in 
the category of target image 2 of orientation of 45 degrees. Thus this approach of using 
NNTool can be used for different image orientations analysis
