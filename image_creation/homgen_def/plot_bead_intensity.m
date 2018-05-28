%this function generates line plots of intensity and intensity change
%versus the spatial coordinate for the series of bead images
%A Landauer, 1/10/17
close all
% clearvars -except bead_img
clear

% if ~exist('bead_img','var')
%     load large_bead_img_series
% end

if ~exist('bead_img','var')
    load small_bead_img_series
    bead_img = small_bead;
    clear small_bead
end

str_lvls = 1:-0.01:0.25;

%colors from the 'RdBu' color map of "BrewerMap" from the File Exchange
% https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer--
%   attractive-and-distinctive-colormaps
[colors1,~,~] = brewermap(length(str_lvls),'Blues');
[colors2,~,~] = brewermap(length(str_lvls),'GnBu');
% colors = colors1+colors2(end:-1:1,:);
% colors = colors/(1.3*max(colors(:)));
colors = colors2.*repmat(linspace(0.5,1.0,size(colors2,1))',1,3);
colors = colors/(1.08*max(colors(:)));
% colors = jet(length(str_lvls));

f1 = figure;
cnt = 0;
for ii = 1:76
    cnt = cnt+1;
    
    sizeOut = size(bead_img{ii});
    I_xy{cnt} = squeeze(bead_img{ii}(:,:,floor(sizeOut(3)/2)));
    I_yz{cnt} = squeeze(bead_img{ii}(floor(sizeOut(1)/2),:,:));
    I_xz{cnt} = squeeze(bead_img{ii}(:,floor(sizeOut(2)/2),:));
    
    I_xy_x{cnt} = I_xy{cnt}(:,floor(sizeOut(1)/2));
    gradI_xy_x{cnt} = gradient(I_xy_x{cnt});
    
    I_xy_y{cnt} = I_xy{cnt}(floor(sizeOut(1)/2),:);
    gradI_xy_y{cnt} = gradient(I_xy_y{cnt});
    
    I_xz_x{cnt} = I_xz{cnt}(:,floor(sizeOut(1)/2));
    gradI_xz_x{cnt} = gradient(I_xz_x{cnt});
    
    I_xz_z{cnt} = I_xz{cnt}(floor(sizeOut(1)/2),:);
    gradI_xz_z{cnt} = gradient(I_xz_z{cnt});
    
    I_yz_r1{cnt} = I_yz{cnt}(:,floor(sizeOut(1)/2));
    I_yz_r2{cnt} = I_yz{cnt}(floor(sizeOut(1)/2),:);
    
    I_yz_r{cnt} = (I_yz_r1{cnt}(:)+I_yz_r2{cnt}(:))/2;
    gradI_yz_r{cnt} = gradient(I_yz_r{cnt});
    
    subplot(2,3,1)
    hold on
    plot(I_xy_x{cnt}(:),'color',colors(cnt,:))
    xlabel('distance, px','Fontsize',18)
    ylabel('Normalized intensity','Fontsize',18)
    title('xy-slice, x-dir')
    set(gca,'Fontsize',16)
    %     axis image
    
    subplot(2,3,2)
    hold on
    plot(I_xy_y{cnt}(:),'color',colors(cnt,:))
    xlabel('distance, px','Fontsize',18)
    ylabel('Normalized intensity','Fontsize',18)
    title('xy-slice, y-dir')
    set(gca,'Fontsize',16)
    %     axis image
    
    subplot(2,3,3)
    hold on
    plot(I_yz_r{cnt}(:),'color',colors(cnt,:))
    xlabel('distance, px','Fontsize',18)
    ylabel('Normalized intensity','Fontsize',18)
    title('yz-slice, radial-dir (radially symmetric)')
    set(gca,'Fontsize',16)
    
    subplot(2,3,4)
    hold on
    plot(gradI_xy_x{cnt}(:),'color',colors(cnt,:))
    xlabel('distance, px','Fontsize',18)
    ylabel('$\frac{dI}{dr}$','interpreter','latex','Fontsize',18)
    title('gradient xy_x')
    set(gca,'Fontsize',16)
    %     axis image
    
    subplot(2,3,5)
    hold on
    plot(gradI_xy_y{cnt}(:),'color',colors(cnt,:))
    xlabel('distance, px','Fontsize',18)
    ylabel('$\frac{dI}{dr}$','interpreter','latex','Fontsize',18)
    title('gradient xy_y')
    set(gca,'Fontsize',16)
    %     axis image
    
    subplot(2,3,6)
    hold on
    plot(gradI_yz_r{cnt}(:),'color',colors(cnt,:))
    xlabel('distance, px','Fontsize',18)
    ylabel('$\frac{dI}{dr}$','interpreter','latex','Fontsize',18)
    title('gradient yz_r')
    set(gca,'Fontsize',16)
    %     axis image
    
end
% legend('\lambda = 1.0','\lambda = 0.85','\lambda = 0.70','\lambda = 0.55',...
%     '\lambda = 0.40','\lambda = 0.25')
set(gca,'Fontsize',16)

%%

figure;

for ii = 1:76
%     if ii == 1
%         cnt = 1;
%         f2 = figure;
%     else
%         cnt = 6;
%         f3 = figure;
%     end
    cnt = ii;
% sizeOut = size(I{ii});
%     I_xy = squeeze(I{ii}(:,:,ceil(sizeOut(3)/2)));
%     I_yz = squeeze(I{ii}(ceil(sizeOut(1)/2),:,:));
%     I_xz = squeeze(I{ii}(:,ceil(sizeOut(2)/2),:));
    
%     I_xy(I_xy < I_extrema(2)*0.0001) = nan; %define an intensity threashold of 0.01%
%     I_yz(I_yz < I_extrema(2)*0.0001) = nan;
%     I_xz(I_xz < I_extrema(2)*0.0001) = nan;
    
    %do the contour plots
%     subplot(1,3,1)
    %         hold on
    imagesc(I_xy{cnt})%,'LineStyle','none')
%     title(['Intensity on x-y'])%, \lambda = ',num2str(stretch_lvls(ii))])
%     colormap(colors)
%     caxis(I_extrema)
%     colorbar
%     axis image
    
% %     subplot(1,3,2)
% %     %        hold on
% %     imshow(I_yz{cnt});
% % %     contourf(I_yz,'LineStyle','none')
% % %     title(['Intensity on y-z'])
% % %     colormap(colors)
% % %     caxis(I_extrema)
% % %     colorbar
% % %     axis image
% %     
% %     subplot(1,3,3)
% %     imshow(I_xz{cnt})
% % %     contourf(I_xz,'LineStyle','none')
% % %     title(['Intensity on x-z'])
% % %     colormap(colors)
% % %     caxis(I_extrema)
% % %     colorbar
% % %     axis image
drawnow
pause(0.1)
end




