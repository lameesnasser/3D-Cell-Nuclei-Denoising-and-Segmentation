
%################# Example for Image Denoising  #########################################################################
% This script aims at demonstrating how to produce the denoised image and simultaneously the potential locations
% of cell nuclei in 2D/3D images as described in "A novel generic dictionary-based
% denoising method for improving noisy and
% densely packed nuclei segmentation in 3D time-lapse fluorescence microscopy
% images" Scientific Reports9 (2019): 5654 (2019)
% doi:10.1038/s41598-019-41683-3
% The script is implemented in MATLAB release R2017b by
% Lamees Nasser
% April 2019
%###################################################################################################################
clear all  % Close all figures
close all  % Clear  workspace
clc        % Clear command window.

%############################ Set the folder path contains the images ##############################################
Path='E:\GITHUB\KSVD_Matlab_ToolBox -3D\Images\3D\'
f= dir (Path)
% ############################  Select the image to be processed ##################################################
fname1='Noisy(SNR-7).tif'
info1 = imfinfo([Path,fname1]);
%Check if image 2D or 3D
z=length(info1)

% ############################ Read Tiff images ##################################################################
if z<=1  %2D images
    IMin=im2double((imread([Path,fname1],1,'Info', info1)));
else     %3D images
    for count=1:length(info1)
        IMin0=im2double((imread([Path,fname1],count,'Info', info1)));
        IMin(:,:,count)=(IMin0);
    end
end
%############################ Initialization parameters ##########################################################
param.K=64                          % - The number of dictionary atoms to train.
param.numIteration =15;             % - The number of K-SVD training iterations
param.L=3;                          % - The number of nonzero elements (used atoms from the dictionary)
%   for the sparse representation coefficient
param.blocksize=[15 15 5];          % - The size of the blocks the algorithm
%                                       works. Example in 2D imgages p=[N N]. 3D images p= [N N M].
param.trainnum=10^4;                % - The number of blocks to train on.
%############################# Find the denoising image and Detection Map ##########################################
[filter_image,map_image,Max_image,Dictionary] = ImageswWithKSVD(IMin,param);
%############################# Display the denoised mage and detection map ##########################################
close all
if z<=1
    figure;
    subplot(2,2,1); imshow(IMin,[]); title('Noisy image');
    subplot(2,2,2); imshow(filter_image,[]); title('Denoised image');
    subplot(2,2,3); imshow(map_image,[]); title('Detection map');
    subplot(2,2,4); imshow(Max_image,[]); title('Maximum response image');
else
    
    for i=1:length(info1)
% figure('name',['Z=' num2str(i)]);
figure;
set(gcf, 'name', ['Frame(Z)=' num2str(i)])
        subplot(2,2,1); imshow(IMin(:,:,i),[]); title('Noisy image');
        subplot(2,2,2); imshow(filter_image(:,:,i),[]); title('Denoised image');
        subplot(2,2,3); imshow(map_image(:,:,i),[]); title('Detection map');
        subplot(2,2,4); imshow(Max_image(:,:,i),[]); title('Maximum response image');
    end
end
%######################################################################################################################



