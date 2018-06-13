%compose image with many beads
clearvars -except bead_img;

disp('locations')
%get the original location for each bead
sizeI = [2047,511,191]; %final image size less one
spacing = 4; %partical min. spacing
nPts = round(2.75*185000); %number of particals
try
    load(['seed_points_',num2str(nPts)])
catch
    x0 = poissonDisc(sizeI,spacing,nPts,0); %generate locations
    save(['seed_points_',num2str(nPts)],'x0')
end


% x0 = [100,sizeI(2)/2,sizeI(3)/2
%       200,sizeI(2)/2,sizeI(3)/2
%       300,sizeI(2)/2,sizeI(3)/2
%       400,sizeI(2)/2,sizeI(3)/2
%       500,sizeI(2)/2,sizeI(3)/2];


sn = 25; %signal to noise ratio, 10 compares reasonably to exp images.
bl = 5; %black level, i.e. background brightness onto which,
% a small amount of noise is superimosed on this
% gaussian noise is applied

nu = 0.0; %Possion's ratio

% x0 = [8,8,8];

%define stretch levels to use
% stretch_lvls = [1:-0.05:0.24];
stretch_lvls = 1:-0.02:0.95;
% translation = [0:0.2:2];
int_decrement = linspace(1.0,0.99,length(stretch_lvls));

if ~exist('bead_img','var')
    load large_bead_img_series
end

%set up parallel pool
% pool = gcp('nocreate');
% if isempty(pool)
%     curCluster = parcluster('local');
%     curCluster.NumWorkers = 16;
%     saveProfile(curCluster);
%     pool = parpool(16);
% end

% stretch_lvls = 0.25;

% s1 = 2;
% s2 = 2;
% s3 = 2;

%%

%final bead size (with 1x bead size margin on each side, undeformed)
bead_size = 17;
disp('images generation')
%generate small bead images
small_bead_img = small_bead_synth(bead_img,stretch_lvls,bead_size);


t0 = tic;
% I = cell(1,length(translation));
% I = cell(1,length(stretch_lvls));

%loop through to compose the images
for stretch = 1:length(stretch_lvls)
    
    t1 = tic;
    
    cur_trans = 0;
    cur_stretch = stretch_lvls(stretch);
    
    %decrease the bead intensity slightly with increasing strain
%     small_bead_img{stretch} = int_decrement(stretch).*small_bead_img{stretch};
    
    %grab bead image
    small_bead = small_bead_img{stretch}(1:end,1:end,1:end);
    
    %compute stretch in each dir
    %     l(1) = cur_stretch;
    %     l(2) = 1/sqrt(cur_stretch);
    %     l(3) = 1/sqrt(cur_stretch);
    
    %compute the strain field
    e_mag = cur_stretch-1;
    e = [e_mag      0            0
        0     -nu*(e_mag)       0
        0          0       -nu*(e_mag)];
    
    %impose motion
    %     x1 = x0 + (x0-sizeI(1)/2)*e;
    x1 = zeros(size(x0));
    %     x1 = x0 + (x0(:,1)-sizeI(1)/2)*e(1,:);
    x1 = x0 + x0(:,1)*e(1,:);
    x1 = x1 + (x0(:,2)-sizeI(2)/2)*e(2,:);
    x1 = x1 + (x0(:,3)-sizeI(3)/2)*e(3,:);
    
    %     x1(:,1) = cur_trans+x0(:,1);
    %     x1(:,2) = x0(:,2);
    %     x1(:,3) = x0(:,3);
    
    disp('bead seeding')
    %seed in beads
    I{1} = seedDefBeads(small_bead,sizeI,x1,bead_size);
    
    I0_ = I{1}(1:512,1:512,1:192);
    
    %     clk = clock;
    %     save(['bead_image_noPSF_output_lambda_',num2str(cur_stretch),'_',...
    %         num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    %         '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'cur_img')
    %
    %     save(['bead_image_noPSF_output_trans_',num2str(cur_trans),'_',...
    %         num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    %         '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'I0')
    %
    I0 = I0_;
%     save(['./test_imgs',num2str(nPts),'_beads_0p05step_lambda_',num2str(cur_stretch),'.mat'],'I0')
    %     save(['5beads_lambda_',num2str(cur_stretch),'.mat'],'I0')
    
    disp('noise generation')
    %rescale to match noise snr and set bl (using a smallish Gaussian)
    bl_mat = zeros(size(I0_));
    for n = 1:size(I0_,3)
        seed_mat = awgn(bl*ones(size(I0_(:,:,n))),0.1);
        seed_mat(seed_mat < 0.01) = 0;
        bl_mat(:,:,n) = seed_mat;
    end
    sig = sn^2;
    I0_ = I0_*sig/max(I0_(:)) + bl.^2;
    
    %resample each pt from the Poisson distro centered on the pt, to mimick
    %photon capture
    [rng1,rng2] = RandStream.create('mrg32k3a','Seed','shuffle','NumStreams',2);
    RandStream.setGlobalStream(rng1);
    I0 = poissrnd(I0_);
    RandStream.setGlobalStream(rng2);
%     tic
    I1 = poissrnd(I0_);
%     toc
    
    %renormalize to 256
    I0(I0>300) = 300; %flatten anything over 118%
    I1(I1>300) = 300;
    I0 = uint8(255*I0/max(I0(:)));
    I1 = uint8(255*I1/max(I1(:)));
    
    save(['./test_imgs',num2str(nPts),'beads_noisy_0p05step_lambda_nu0',...
        num2str(cur_stretch),'_static.mat'],'I0','I1')
    
    clear vol
    vol{1} = I0;
    
    save(['vol_series_e1_',num2str(stretch+1000)],'vol')
    
    toc(t1)
end
toc(t0)

% % for stretch = 1:19
% %     cur_stretch = stretch_lvls(stretch);
% %     load(['./test_imgs',num2str(nPts),'beads_noisy_0p05step_lambda_',...
% %         num2str(cur_stretch),'_static.mat'])
% %         clear vol
% %     vol{1} = I0;
% %     save(['vol_series_',num2str(stretch+1000)],'vol')
% % end

disp('done')
clear *_

% clk = clock;
% save(['bead_image_65000p_output_',...
%     num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
%     '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'I','-v7.3')

% figure
% imshow3D(I0)


%old
%load in the bead images for this deformation case
%     load(['bead_PSF_stretch_',num2str(cur_stretch),'.mat'])

%     cur_trans = translation(stretch);
%
%     %     cur_bead_idx = find(bead_imgs.stretch == cur_stretch);
%     %     cur_bead_img = bead_imgs.bead(cur_bead_idx);
%
%
%     %create a bead image with no psf
%     r = 5;
%     a(1) = r*cur_stretch;
%     a(2) = r/sqrt(cur_stretch);
%     a(3) = r/sqrt(cur_stretch);
%     num_px = 10*bead_size+1;
%     k = 3;
%     v = linspace(-2*r,2*r,num_px);
%
% %     large_bead = zeros(num_px,num_px,num_px);
%     large_bead_sharp = zeros(num_px,num_px,num_px);
%
%     for ii = 1:num_px
%         x = v(ii);
%         for jj = 1:num_px
%             y = v(jj);
%             for kk = 1:num_px
%                 z = v(kk);
%                 large_bead_sharp(ii,jj,kk) = 1./(1+exp(-2*k*(1-((x).^2/a(1)^2)-...
%                     ((y).^2/a(2)^2)-((z).^2/a(3)^2))));
%                 gaussian_psf = exp(-((x.^2)/(2*s1^2) + (y.^2)/(2*s2^2) + (z.^2)/(2*s1^2)));
%             end
%         end
%     end


