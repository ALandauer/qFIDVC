function [large_bead] = generate_single_speckle_3d(transform_matrix,bead_size)
%Generate bead image by densely sampling the bead representation, the psf, and
%convolving the two
%
% input: transformation matrix [3x3] with which to warp the bead
% output: the image of the tranformed bead
%
%AL, 5/10/18

%make sure the input transform matrix is [3x3]
if size(transform_matrix,1) == 2
    transform_matrix(1,1) = transform_matrix(1,1)+1;
    transform_matrix(2,2) = transform_matrix(2,2)+1;
    transform_matrix(1,3) = 0;
    transform_matrix(2,3) = 0;
    transform_matrix(3,3) = 1;
else
    transform_matrix(1,1) = transform_matrix(1,1)+1;
    transform_matrix(2,2) = transform_matrix(2,2)+1;
    transform_matrix(1,2) = transform_matrix(1,2);
    transform_matrix(2,1) = transform_matrix(1,2);
    transform_matrix(1,3) = transform_matrix(1,3);
    transform_matrix(2,3) = transform_matrix(2,3);
    transform_matrix(3,3) = transform_matrix(3,3)+1;
    transform_matrix(3,1) = transform_matrix(1,3);
    transform_matrix(3,2) = transform_matrix(2,3);
end
if size(transform_matrix,2) == 2
    transform_matrix(3,1) = 0;
    transform_matrix(3,2) = 0;
    transform_matrix(3,3) = 1;
end

for ii = 1:4
    transform_matrix(4,ii) = 0;
    transform_matrix(ii,4) = 0;
    transform_matrix(4,4) = 1;
end

T = affine3d(transform_matrix);

%create a bead image and psf
r = bead_size;
a1 = r;
a2 = r;
a3 = r;

num_px = 20*bead_size+1;
k = 1000;
v = linspace(-3*bead_size,3*bead_size,num_px);

%compute the correct st dev of the Gaussian approx to the PSF
L = 532*10^-9*(bead_size/(0.75*10^-6)); %~actual bead size for 40x imaging
NA = 0.75;% ~for 40x
s1 = 0.42*L*1/(2*NA);
s2 = 0.42*L*1/(2*NA);
s3 = 0.84*L*1/(2*NA);

%compute each gaussian
gaussian_psf_x = exp(-((v.^2)/(2*s1^2)))';
gaussian_psf_y = exp(-((v.^2)/(2*s2^2)))';
gaussian_psf_z = exp(-((v.^2)/(2*s3^2)))';

%generate the large bead on a meshgrid of v
[xg,yg,zg] = meshgrid(v,v,v);
large_bead_sharp = 1./(1+exp(-2*k*(1-((xg).^2/a1^2)-...
                ((yg).^2/a2^2)-((zg).^2/a3^2))));

%since the Gaussianis seperable, do the convolution of the PSF and the
%bead as a series of nx1 convolutions
large_bead_ = convn(gaussian_psf_x,large_bead_sharp);
large_bead_ = convn(gaussian_psf_y,rot90_3D(large_bead_(floor(num_px/2):num_px+floor(num_px/2)-1,:,:),2,1));
large_bead_ = convn(gaussian_psf_z,rot90_3D(large_bead_(floor(num_px/2):num_px+floor(num_px/2)-1,:,:),3,1));
large_bead_ = rot90_3D(rot90_3D(large_bead_(floor(num_px/2):num_px+floor(num_px/2)-1,:,:),3,-1),2,1);
% fig,imagesc3d(large_bead_),axis image

large_bead_warp = imwarp(large_bead_,T);

large_bead = large_bead_warp;%(num_px-floor(num_px/2):num_px+floor(num_px/2),...
%     num_px-floor(num_px/2):num_px+floor(num_px/2),...
%     num_px-floor(num_px/2):num_px+floor(num_px/2));

end

function M = rot90_3D(M, Dim, NumRot)
% Extends the rot90 command to 3-Dimensions.
% The input matrix A_Mat is assumed to be a 3D matrix. This matrix is %
% rotated by 90degrees (or 180 or 270) about any of the 3 axes.
% Dim is the dimension to rotate normal to.
% NumRot is the number of 90-degree rotations.
% A_Mat must be specified as input matrix.
% Dim should be an integer in the range 1..3 for 3 dimensions.
% if Dim is not specified, then Dim=1. 
% If NumRot not specified then NumRot=1
 
%
% EXAMPLE
% here, the input matrix is smoki and we rotate it 90 degrees about the
% z-axis.
% smoki(1:3,1,1)=1;
% smoki(1:3,2,2)=1;
% smoki(1:3,3,3)=1;
% yoki=rot90_3D(smoki,3,1)
% EXAMPLE
% smoki(1:3,1,1)=1;
% smoki(1:3,3,3)=2;
% yoki=rot90_3D(smoki,2,1);
% 
% Tested: Matlab 7.7
% Author: Matt Fetterman, mattinjersey@yahoo.com . Thanks to Jan Simon for
% modifications.
aSize = size(M);

switch Dim
  case 1
      switch NumRot
        case 1
            X = permute(M, [1, 3, 2]);
            M = X(:, aSize(3):-1:1, :);
        case 2
            M = M(:, aSize(2):-1:1, aSize(3):-1:1);
        case 3
            M = permute(M(:, aSize(2):-1:1, :), [1, 3, 2]);
      end
      
  case 2
      switch NumRot
        case 1
            X = permute(M, [3, 2, 1]);
            M = X(aSize(3):-1:1, :, :);
        case 2
            M = M(aSize(1):-1:1, :, aSize(3):-1:1);
        case 3
            M = permute(M(aSize(1):-1:1, :, :), [3, 2, 1]);
      end
      
  case 3
      switch NumRot
        case 1
            X = permute(M, [2, 1, 3]);
            M = X(aSize(2):-1:1, :, :);
        case 2
            M = M(aSize(1):-1:1, aSize(2):-1:1, :);
        case 3
            M = permute(M(aSize(1):-1:1, :, :), [2, 1, 3]);
      end
      
  otherwise
      error('Dim must be 1, 2 or 3');
end

return;
end
