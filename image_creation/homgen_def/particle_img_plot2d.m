%This function plots the deformed and undeformed configuration of the whole
%seeded image as a centerline 2d slice
close all
clearvars -except def_img undef_img

if ~exist('undef_img','var')
load 100000beads_noisy_0p05step_lambda_1_static
undef_img = I0;
end

if ~exist('def_img','var')
load('100000beads_noisy_0p05step_lambda_0.25_static.mat')
def_img = I0;%/max(I0(:));
end

clear I0 I1

X(1) = 460;
X(2) = 2050;
Y(1) = 1580;
Y(2) = 2120;

example_img_undef_ = imread('c095p1_1_25ms_000135','tif');
size_img = size(example_img_undef_);
example_img_undef = example_img_undef_(Y(1):Y(2),X(1):X(2));

example_img_def_ = imread('c095p1_1_25ms_000356','tif');
size_img = size(example_img_def_);
example_img_def = example_img_def_(Y(1):Y(2),X(1):X(2));

sizeI = size(def_img);

undef_img_slice = undef_img(:,:,floor(sizeI(3)/2));
def_img_slice = def_img(:,:,floor(sizeI(3)/2));

% figure
% imshow3d(undef_img)
% figure
% imshow3d(def_img)

undef_slice_norm = uint8(255*undef_img_slice/(max(undef_img_slice(:))));
def_slice_norm = uint8(255*def_img_slice/(max(def_img_slice(:))));

figure
subplot(2,4,1)
imshow(undef_slice_norm(1:130,:))
% imshow(undef_slice_norm(:,:))
title('Synthetic, undeformed')
set(gca,'fontsize',16)
subplot(2,4,2)
imshow(def_slice_norm(1:130,:))
set(gca,'fontsize',16)
title('Synthetic, 75% compressed')
subplot(2,4,3)
imshow(example_img_undef)
set(gca,'fontsize',16)
title('Experimental, undeformed')
subplot(2,4,4)
imshow(example_img_def)
title('Experimental, 75% compressed')
set(gca,'fontsize',16)

subplot(2,4,5)
histogram(undef_slice_norm(1:130,:),20)
axis([0 300 0 9250])
title('Histograms')
set(gca,'fontsize',22)
subplot(2,4,6)
histogram(def_slice_norm(1:130,:),20)
axis([0 300 0 9250])
set(gca,'fontsize',22)
subplot(2,4,7)
histogram(example_img_undef(:,:),20)
axis([0 300 0 125000])
set(gca,'fontsize',22)
subplot(2,4,8)
histogram(example_img_def(:,:),20)
axis([0 300 0 125000])
set(gca,'fontsize',22)
% 
% figure
% subplot(1,2,1)
% imagesc(undef_img_slice)
% subplot(1,2,2)
% imshow(def_img_slice)



