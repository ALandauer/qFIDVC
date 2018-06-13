%Generate bead image by densely sampling the bead representation, the psf, and
%convolving the two

%A Landauer, 3/11

clear,clc

stretch_lvls = [1:-0.01:0.25];

%set the size of the bead
bead_size = 9;
t0 = tic;
I = cell(1,length(stretch_lvls));

% %set up parallel pool
% pool = gcp('nocreate');
% if isempty(pool)
%     curCluster = parcluster('local');
%     curCluster.NumWorkers = 6;
%     saveProfile(curCluster);
%     pool = parpool(6);
% end


%loop through to compose the images
parfor nn = 1:length(stretch_lvls)
    t1 = tic;
    
    cur_stretch = stretch_lvls(nn);
    
    %create a bead image and psf
    r = bead_size;
    a1 = r*cur_stretch; %Vol. preserving uniaxial def.
    a2 = r/sqrt(cur_stretch);
    a3 = r/sqrt(cur_stretch);
    num_px = 20*bead_size+1;
    k = 1000;
    v = linspace(-(2*r/sqrt(min(stretch_lvls))),2*r/sqrt(min(stretch_lvls)),num_px);
    
    %compute the correct st dev of the Gaussian approx to the PSF
    L = 532*10^-9*(bead_size/(0.5*10^-6)); %~actual bead size for 40x imaging
    NA = 0.7;% ~for 40x
    s1 = 0.42*L*1/(2*NA);
    s2 = 0.42*L*1/(2*NA);
    s3 = 2*(0.42*L*1/(2*NA)); %blurring factor of 2x along the axis of the lens
    
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
                %large bead shape bead on volume-preserving deformation of
                %a sphere
                large_bead_sharp(ii,jj,kk) = 1./(1+exp(-2*k*(1-((x).^2/a1^2)-...
                    ((y).^2/a2^2)-((z).^2/a3^2))));
                %                 gaussian_psf(ii,jj,kk) = exp(-((x.^2)/(2*s1^2) + (y.^2)/(2*s2^2)...
                %                     + (z.^2)/(2*s1^2)));
                if ii == 1
                    %Gaussian approx. to the Bessel function
                    gaussian_psf_z(1,1,kk) = exp(-((z.^2)/(2*s3^2)));
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
    
%     tic
%     gx = gpuArray(gaussian_psf_x);
%     gy = gpuArray(gaussian_psf_y);
%     gz = gpuArray(gaussian_psf_z);
%     lbs = gpuArray(large_bead_sharp);
%     
%     lb_ = convn(gx,lbs);
%     lb_ = convn(gy,lb_);
%     lb_ = convn(gz,lb_);
%     
%     large_bead_ = gather(lb_);
%     toc
    
    large_bead_ = convn(gaussian_psf_x,large_bead_sharp);
    large_bead_ = convn(gaussian_psf_y,large_bead_);
    large_bead_ = convn(gaussian_psf_z,large_bead_);
    
    large_bead = large_bead_(num_px-floor(num_px/2):num_px+floor(num_px/2),...
        num_px-floor(num_px/2):num_px+floor(num_px/2),...,...
        num_px-floor(num_px/2):num_px+floor(num_px/2))/max(large_bead_(:));%trim back down
    
    I{nn} = large_bead;
    I_2d{nn} = large_bead(:,:,round(size(large_bead,3)/2)); 
    %     hold on
    %     imshow3d(large_bead)
    %     drawnow
    %     clk = clock;
    
    %     save(['bead_PSF_stretch_',num2str(cur_stretch),'_',...
    %         num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    %         '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'large_bead')
    
    % save(['bead_PSF_stretch_',num2str(cur_stretch),'.mat'],'large_bead')
    
    toc(t1)
end
toc(t0)

clear *_

clk = clock;
save(['bead_images_output_',...
    num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'])

bead_img = I;
save('large_bead_img_series','I','-v7.3');
save('large_bead_img_series_2d','I_2d','-v7.3');


% % %% PLOTTING
% % % figure
% % vid_out = VideoWriter('./large_bead_compression.mp4');%,'MPEG-4');
% % 
% % % Set the frame rate
% % set(vid_out,'FrameRate',8);
% % % Open the object for writing
% % open(vid_out);
% % framepar.resolution = [1920,1080];
% % % Create a new figure
% % % f1 = figure();
% % hFigRotated = figure;
% % scrsz = [1 1 1920 1080];% get(0,'ScreenSize');
% % set(hFigRotated,'position',[50,50,scrsz(3)-150,scrsz(4)-150],'visible','on');
% % 
% % I_extrema = [min(min(min(cell2mat(I)))),max(max(max(cell2mat(I))))];
% % 
% % colors = [0.4039    0.0000    0.1216
% %     0.6801    0.0627    0.1607
% %     0.8169    0.3268    0.2720
% %     0.9302    0.5762    0.4449
% %     0.9876    0.7904    0.6813
% %     0.9907    0.9392    0.9073
% %     0.9156    0.9520    0.9741
% %     0.7412    0.8585    0.9179
% %     0.4869    0.7230    0.8424
% %     0.2210    0.5427    0.7503
% %     0.1247    0.3834    0.6603
% %     0.0196    0.1882    0.3804];
% % 
% % for ii = 1:length(I)
% %     ii
% %         
% %     sizeOut = size(I{ii});
% %     I_xy = squeeze(I{ii}(:,:,ceil(sizeOut(3)/2)));
% %     I_yz = squeeze(I{ii}(ceil(sizeOut(1)/2),:,:));
% %     I_xz = squeeze(I{ii}(:,ceil(sizeOut(2)/2),:));
% %     
% %     I_xy(I_xy < I_extrema(2)*0.0001) = nan; %define an intensity threashold of 0.01%
% %     I_yz(I_yz < I_extrema(2)*0.0001) = nan;
% %     I_xz(I_xz < I_extrema(2)*0.0001) = nan;
% %     
% %     %do the contour plots
% %     subplot(1,3,1)
% %     %         hold on
% %     contourf(I_xy,'LineStyle','none')
% %     title(['Intensity on x-y, \lambda = ',num2str(stretch_lvls(ii))])
% %     colormap(colors)
% %     caxis(I_extrema)
% %     colorbar
% %     axis image
% %     
% %     subplot(1,3,2)
% %     %        hold on
% %     contourf(I_yz,'LineStyle','none')
% %     title(['Intensity on y-z'])
% %     colormap(colors)
% %     caxis(I_extrema)
% %     colorbar
% %     axis image
% %     
% %     
% %     subplot(1,3,3)
% %     hold on
% %     contourf(I_xz,'LineStyle','none')
% %     title(['Intensity on x-z'])
% %     colormap(colors)
% %     caxis(I_extrema)
% %     colorbar
% %     axis image
% %    
% %     %     hAxRotated  = axes;
% %     %     slice(double(I{ii}),sizeOut(2)/2,sizeOut(1)/2,sizeOut(3)/2,'method','nearest');
% %     %     grid on, colormap gray, axis image
% %     %     shading interp
% %     drawnow
% %     
% %     F1 = fig2frame(hFigRotated,framepar); %construct a frame, from fig2frame
% %     % https://github.com/cmacminn/fig2frame
% %     writeVideo(vid_out,F1)
% %     clear F1
% %     
% % end
% % close(vid_out);

