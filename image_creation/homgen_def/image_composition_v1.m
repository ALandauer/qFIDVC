%compose image with many beads
clear;

%get the original location for each bead
sizeI = [511,511,191];
spacing = 6;
nPts = 65000;
x0 = poissonDisc(sizeI,spacing,nPts,0);
x1 = zeros(size(x0));

% x0 = [8,8,8];

%define stretch levels to use
stretch_lvls = [1:-0.01:0.25];
translation = [0:0.2:2];


load large_bead_img_series

% s1 = 2;
% s2 = 2;
% s3 = 2;

%%

%final bead size
bead_size = 17;
t0 = tic;
% I = cell(1,length(translation));
I = cell(1,length(stretch_lvls));

%loop through to compose the images
for stretch = 1:length(stretch_lvls)
%     load(['bead_PSF_stretch_',num2str(cur_stretch),'.mat'])
    
    t1 = tic;
    defC = sizeI;
    
        cur_trans = 0;
        cur_stretch = stretch_lvls(stretch);
    
%     cur_stretch = 1;
    
    large_bead = bead_img{stretch}(1:end,1:end,1:end);
    
    %compute stretch in each dir
    e(1) = cur_stretch;
    e(2) = 1/sqrt(cur_stretch);
    e(3) = 1/sqrt(cur_stretch);
    
    %impose motion
        x1(:,1) = e(1).*x0(:,1);
        x1(:,2) = e(2).*x0(:,2);
        x1(:,3) = e(3).*x0(:,3);
    
%     x1(:,1) = cur_trans+x0(:,1);
%     x1(:,2) = x0(:,2);
%     x1(:,3) = x0(:,3);
    
    %seed in beads
    I{stretch} = seedDefBeads(large_bead,sizeI,x1,bead_size);
    
    I0 = I{stretch};
    
    clk = clock;
    %     save(['bead_image_noPSF_output_lambda_',num2str(cur_stretch),'_',...
    %         num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    %         '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'cur_img')
%     
%     save(['bead_image_noPSF_output_trans_',num2str(cur_trans),'_',...
%         num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
%         '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'I0')
%     
    save(['bead_image_65000p_lambda_',num2str(cur_stretch),'.mat'],'I0')
    
    toc(t1)
end
toc(t0)

clk = clock;
save(['bead_image_65000p_output_',...
    num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'I','-v7.3')

figure
imshow3D(I0)


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


