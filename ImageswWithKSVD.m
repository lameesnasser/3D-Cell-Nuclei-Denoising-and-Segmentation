
function [Denoised_Image,Detection_Map,Max_image,Dictionary] = ImageswWithKSVD(Image,param);
%####################################################################################################################################################
%   Image Denoising Using Sparse Representation Model and Dictionary
%                Learning on a Noisy Image
%####################################################################################################################################################
% ImageswWithKSVD used for denoising and detecting of 2D/3D cell nuclei images
% based on a sparse representation model. A more detailed description can
% be found in:
% "A novel generic dictionary-based denoising method for improving noisy and
% densely packed nuclei segmentation in 3D time-lapse fluorescence microscopy
% images" Scientific Reports9 (2019): 5654 (2019)
% doi:10.1038/s41598-019-41683-3
% The script is implemented in MATLAB release R2017b by
%  Lamees Nasser
%  April 2019

% %##################################################################################################################################################

% INPUT ARGUMENTS :  Image               -  The noisy image in gray-level scale.
%                    param.K             -  The number of dictionary atoms to train.
%                    param.numIteration  -  The number of K-SVD training iterations
%                    param.blockSize     -  The size of the blocks the algorithm
%                                        works. Example in 2D imgages p=[N N]. 3D images p= [N N M].
%                    param.L             - The number of nonzero elements (used atoms from the dictionary)
%                                         for the sparse representation coefficient
%                    param.trainnum      -  The number of blocks to train on.


% OUTPUT ARGUMENTS :  Denoised_Image - The cleaned image.
%
%                     Detection_Map  - Image indicates the potential
%                                      locations of cell nuclei.

%                     Dictionary     -  Updated dictionary.

%%############################# Check If 2D or 3D image##########################################################################
[x,y,z]=size(Image);

if z>1
    p=3; %3D Image
    param.stepsize=[1 1 1];  % Step size of sliding window
    
else
    p=2; %2D Image
    param.stepsize=[1 1];    % Step size of sliding window
    
end

%% ################################## An initial dictionary is constructed by selecting random patches from the patches extracted from the noisy image##########################################################################
%  among those having intensities greater than the obtained average intensity.
ids = cell(p,1);
[ids{:}] = reggrid(size(Image)-param.blocksize+1, param.trainnum, 'eqdist');
blkMatrix = sampgrid(Image,param.blocksize,ids{:});
summation=sum(blkMatrix);
[val ind]=find(summation>mean(summation));
param.initialDictionary =  blkMatrix(:,ind(1:param.K ));
%% Waiting for Training the dictionary
counterForWaitBar = param.numIteration+1;
h = waitbar(0,'Trainning Dictionary In Process ...');
param.waitBarHandle = h;
param.counterForWaitBar = counterForWaitBar;
%% ################################# Update and obtain the final dictionary by using K-SVD  algorithm from the paper "The K-SVD: An Algorithm for
% Designing of Overcomplete Dictionaries for Sparse Representation", written  by M. Aharon, M. Elad, and A.M. Bruckstein
%and appeared in the IEEE Trans. On Signal Processing, Vol. 54, no. 11, pp. 4311-4322, November 2006.
% The used K-SVD toolbox can be find in https://elad.cs.technion.ac.il/software/
[Dictionary,output] = KSVD(blkMatrix,param);
output.D = Dictionary;
Dictionary(Dictionary>1)=1;
Dictionary(Dictionary<0)=0;
disp('Finished Trainning Dictionary');

%%###############################  Reconstruction of the denoised image and detection map##########################################################################
counterForWaitBar = 1
h1 = waitbar(0,'Images Reconstrucion In Process ...')

stepsize=param.stepsize;
blocksize=param.blocksize;
xount=0;
nz=0;
y_clean = zeros(size(Image));
y_Map = zeros(size(Image));
if p==3  %3D Images
    for k = 1:stepsize(3):size(y_clean,3)-blocksize(3)+1
        for j = 1:stepsize(2):size(y_clean,2)-blocksize(2)+1
            xount=xount+1;
            blocks = im2colstep(Image(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1),blocksize,stepsize);
            gamma  = full(mexOMP(blocks, Dictionary, param));
            Mapblocks= ones(prod(param.blocksize),param.K)*gamma;
            cleanblocks= Dictionary*gamma;
            nz = nz + nnz(gamma);
            cleanvol = col2imstep(cleanblocks,[size(y_clean,1) blocksize(2:3)],blocksize,stepsize);
            y_clean(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) = y_clean(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) + cleanvol;
            Mapvol = col2imstep(Mapblocks,[size(y_Map,1) blocksize(2:3)],blocksize,stepsize);
            y_Map(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) =  y_Map(:,j:j+blocksize(2)-1,k:k+blocksize(3)-1) + Mapvol;
        end
        
        param.waitBarHandle = h1;
        waitbar(k/param.counterForWaitBar);
        param.counterForWaitBar = counterForWaitBar;
        
    end
    
elseif p==2 %2D Images
    
    blocks = im2colstep(Image,blocksize,stepsize);
    gamma  = full(mexOMP(blocks, Dictionary, param));
    blocks1 =Dictionary*gamma;
    blocks_map = ones( prod(param.blocksize),param.K)*gamma;
    y_clean = col2imstep(blocks1,[size(Image)],blocksize,stepsize);
    y_Map = col2imstep(blocks_map,[size(Image)],blocksize,stepsize);
      param.waitBarHandle = h1;
        waitbar(param.counterForWaitBar);
        param.counterForWaitBar = counterForWaitBar;
end
lambda=0;
cnt = countcover(size(Image),blocksize,stepsize);
Denoised_Image = (y_clean+lambda*Image)./(1*(cnt + lambda));
Detection_Map= (y_Map+lambda*Image)./(1*(cnt + lambda));
Max_image=Denoised_Image.*Detection_Map;
disp('Finished Images Reconstrucion')


