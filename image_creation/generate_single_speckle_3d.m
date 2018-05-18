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
s3 = 0.8*L*1/(2*NA);

large_bead_sharp = zeros(num_px,num_px,num_px);
gaussian_psf_x = zeros(num_px,1);
gaussian_psf_y = zeros(1,num_px);
gaussian_psf_z = zeros(1,1,num_px);

for ii = 1:num_px
    x = v(ii);
    for jj = 1:num_px
        y = v(jj);
        for kk = 1:num_px
            z = v(kk);
            large_bead_sharp(ii,jj,kk) = 1./(1+exp(-2*k*(1-((x).^2/a1^2)-...
                ((y).^2/a2^2)-((z).^2/a3^2))));
            if jj == 1
                gaussian_psf_z(kk) = exp(-((z.^2)/(2*s3^2)));
            end
        end
        if ii == 1
            gaussian_psf_y(jj) = exp(-((y.^2)/(2*s2^2)));
        end
    end
    gaussian_psf_x(ii) = exp(-((x.^2)/(2*s1^2)));
    
end

%since the Gaussianis seperable, do the convolution of the PSF and the
%bead as a series of nx1 convolutions
large_bead_ = convn(gaussian_psf_x,large_bead_sharp);
large_bead_ = convn(gaussian_psf_y,large_bead_);
large_bead_ = convn(gaussian_psf_z,large_bead_);

large_bead_warp = imwarp(large_bead_,T);

large_bead = large_bead_warp(num_px-floor(num_px/2):num_px+floor(num_px/2),...
    num_px-floor(num_px/2):num_px+floor(num_px/2),...
    num_px-floor(num_px/2):num_px+floor(num_px/2));

end

