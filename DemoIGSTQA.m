% This code is written by: S. Alireza Golestaneh;
% Email: sgolest1@asu.edu;
% March 1 2017; for CVPRW 2018; Paper: Synthesized Texture Quality Assessment 
% via Multi-scale Spatial and Statistical Texture Attributes of Image and Gradient Magnitude Coefficients
% S. Alireza Golestaneh and Lina Karam

% @InProceedings{Golestaneh_2018_CVPR,
% author = {Alireza Golestaneh, S. and Karam, Lina J.},
% title = {Spatially-Varying Blur Detection Based on Multiscale Fused and Sorted Transform Coefficients of Gradient Magnitudes},
% booktitle = Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition Workshops},
% month = {July},
% year = {2018}
% }
% This code was tested on MATLAB2017a
%=====================================================================

%%
clc; clear all;
close all;
path(pathdef);

Ref_Img=imread('bananas.png');
Synth_Img=imread('bananas_Alg1.png');
QualityScore = IGSTQA(Ref_Img,Synth_Img)

% Important Note
% For inputs:
% Ref_Img=imread('bananas.png');
% Synth_Img=imread('bananas_Alg1.png');
% you must get the QualityScore =    17.6479; otherwise, there is something
% different in your setting than the setting that we
% used for our experiment/analysis
