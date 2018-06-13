function [u_masked,u_imp_masked,u_err,l2err,std_err,mean_err,imp_strain,nan_mask] = ...
    error_mapping_homogen_3d(u,dm,m,bndry_wd)

% try
%     load(['.././results_hybrid_',num2str(numImages)],'nan_mask_hybrid')
%     load(['.././results_inc_',num2str(numImages)],'nan_mask_inc')
%     load(['.././results_cum_',num2str(numImages)],'nan_mask_cum')
%     nan_mask = nan_mask_hybrid.*nan_mask_inc.*nan_mask_cum;
% catch
% end

nu = 0.0; %from image_composition.m

stretch_lvls = 0.98:-0.02:0.25;

% seedPoints = load(['.',filesep,'seed_points_508750.mat']);
% x0 = seedPoints.x0;

[xq,yq,zq] = ndgrid(1:512,1:512,1:192);
    
x0_1 = reshape(xq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
    numel(xq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end])),1);
x0_2 = reshape(yq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
    numel(yq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end])),1);
x0_3 = reshape(zq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
    numel(zq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end])),1);

x0 = [x0_1,x0_2,x0_3];

sizeI = [2047,511,191];
final_size = [512,512,192];

% sizeI = [512,512];

% x0 = [];
% for ii = 1:(sizeI(1)-0)
%     x_ = [1:(sizeI(2)-0);ii*ones(1,sizeI(2))]';
%     x0 = cat(1,x0,x_);
% end

% if ~exist('nan_mask','var')
nan_mask = ones(final_size/8+1);
nan_mask(1:10,:,:) = nan;
nan_mask(:,1:10,:) = nan;
nan_mask(:,:,1:4) = nan;
nan_mask(end-9:end,:,:) = nan;
nan_mask(:,end-9:end,:) = nan;
nan_mask(:,:,end-3:end) = nan;
% end

for cnt = 1:length(u)
    
    
    disp(['disp num ',num2str(cnt)])
    
    
    cur_stretch = stretch_lvls(cnt);
    e_mag = cur_stretch-1;
    e = [-nu*(e_mag)      0            0
        0     -nu*(e_mag)      0
        0          0       e_mag];
    
    %imposed motion
    % %
    % %     x1 = x0 + (x0(:,1)-sizeI(1)/2)*e(1,:);
    % %     x1 = x1 + (x0(:,2)-sizeI(2)/2)*e(2,:);
    % %     x1 = x1 + (x0(:,3)-sizeI(3)/2)*e(3,:);
    
    x1 = x0 + x0(:,1)*e(1,:);
    x1 = x1 + (x0(:,2)-sizeI(2)/2)*e(2,:);
    x1 = x1 + (x0(:,3)-sizeI(3)/2)*e(3,:);
    
    u1 = x1(:,1)-x0(:,1);
    u2 = x1(:,2)-x0(:,2);
    u3 = x1(:,3)-x0(:,3);
    
    %     [xq,yq,zq] = ndgrid(round(final_size(2)/2+1):round(3/2*final_size(2)),...
    %         (floor(sizeI(1)/2)+1):(final_size(1)+floor(sizeI(1)/2)),...
    %         round(final_size(3)/2+1):round(3/2*final_size(3)));
    
    [xq,yq,zq] = ndgrid(1:512,1:512,1:192);
    
    U1 = scatteredInterpolant(x0(:,1),x0(:,2),x0(:,3),u1,'natural');
    Uq1 = U1(xq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
        yq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
        zq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]));
    
    U2 = scatteredInterpolant(x0(:,1),x0(:,2),x0(:,3),u2,'natural');
    Uq2 = U2(xq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
        yq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
        zq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]));
    
    U3 = scatteredInterpolant(x0(:,1),x0(:,2),x0(:,3),u3,'natural');
    Uq3 = U3(xq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
        yq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]),...
        zq([1:dm:end,end],[1:dm:end,end],[1:dm:end,end]));
    
    u_imp{cnt}{2} = Uq1;
    u_imp{cnt}{1} = Uq2;
    u_imp{cnt}{3} = Uq3;
    
    %     u_imp{cnt}{4} = sqrt(u_([1:dm:end,end],1).^2 ...
    %         + u_([1:dm:end,end],2).^2 + u_([1:dm:end,end],3).^2);
    
    
    %     close all
    %     figure,imagesc(u_new{1}{1}),colorbar
    %     figure,imagesc(u_new{1}{2}),colorbar
    %     figure,imagesc(u_new{1}{3}),colorbar
    %
    %     figure,imagesc(u_imp{1}{1}),colorbar
    %     figure,imagesc(u_imp{1}{2}),colorbar
    %     figure,imagesc(u_imp{1}{3}),colorbar
    
    for ii = 1:3
        % %         figure
        % %
        % % %         u_err_std{cnt}{ii} = u{cnt}{ii} - u_imp{cnt}{ii};
        
        u_imp_masked{cnt}{ii} = nan_mask.*u_imp{cnt}{ii};
        u_masked{cnt}{ii} = nan_mask.*u{cnt}{ii};
        
        u_err{cnt}{ii} = nan_mask.*u{cnt}{ii} - nan_mask.*u_imp{cnt}{ii};
        
        % % %         figure,imagesc(u_err_std{1}{ii}(dm:end-dm,dm:end-dm)),colorbar
        % % %         title('old')
        % %         subplot(1,2,1)
        % %         imagesc(u_err_new{1}{ii}(dm:end-dm,dm:end-dm)),colorbar
        % %         title(['new, u',num2str(ii)])
        % %         axis image
        % %
        % % %         figure,histogram(u_err_std{1}{ii}(:),1000)
        % %         subplot(1,2,2)
        % %         histogram(u_err_new{1}{ii}(:),1000)
        % %         title(['Histogram of errors, u',num2str(ii)])
        % % %         std(u_err_std{1}{ii}(:))
        % %         disp(['std of error, u',num2str(ii)])
        % %         std(u_err_new{1}{ii}(:))
        
        %compute L2 norm...
        %         u_err_2norm{cnt}(ii) = sqrt(nansum(nansum(u_err{cnt}{ii}.*u_err{cnt}{ii})));
        %         imp_strain{cnt}(ii) = e(ii,ii);
        
        l2err(cnt,ii) = sqrt(nansum(nansum(u_err{cnt}{ii}...
            (1+bndry_wd:end-bndry_wd).*u_err{cnt}{ii}(1+bndry_wd:end-bndry_wd))));
        std_err(cnt,ii) = nanstd(u_err{cnt}{ii}(1+bndry_wd:end-bndry_wd));
        mean_err(cnt,ii) = sqrt(nansum(nansum((...
            u_err{cnt}{ii}(1+bndry_wd:end-bndry_wd)).^2)))/...
            nansum(nansum(isfinite(u_err{cnt}{ii}(1+bndry_wd:end-bndry_wd))));
        imp_strain(cnt,ii) = e(ii,ii);
        
    end
    
    
    
end



% figure,
% subplot(2,2,1)
% contourf(u_imp_masked{1}{1}),colorbar
% axis image
% title('u1_{imp}')
% subplot(2,2,2)
% contourf(u_imp_masked{end}{1}),colorbar
% axis image
% subplot(2,2,3)
% contourf(u_imp_masked{1}{2}),colorbar
% axis image
% title('u2_{imp}')
% subplot(2,2,4)
% contourf(u_imp_masked{end}{2}),colorbar
% axis image
% saveas(gcf,['u_imp',datestr(now,'mmddTHHMMSS')])
% % set(gcf,'Position',[3 3 15 6]);
%
% figure
% subplot(2,2,1)
% contourf(u_masked{1}{1}),colorbar
% axis image
% title('u1_{DIC}')
% subplot(2,2,2)
% contourf(u_masked{end}{1}),colorbar
% axis image
% subplot(2,2,3)
% contourf(u_masked{1}{2}),colorbar
% axis image
% title('u2_{DIC}')
% subplot(2,2,4)
% contourf(u_masked{end}{2}),colorbar
% axis image
% saveas(gcf,['u_DIC',datestr(now,'mmddTHHMMSS')])
% % set(gcf,'Position',[3 3 15 6]);
%
% figure
% subplot(2,2,1)
% contourf(u_err{1}{1}),colorbar
% axis image
% title('u1_{Err}')
% subplot(2,2,2)
% contourf(u_err{end}{1}),colorbar
% axis image
% subplot(2,2,3)
% contourf(u_err{1}{2}),colorbar
% axis image
% title('u2_{Err}')
% subplot(2,2,4)
% contourf(u_err{end}{2}),colorbar
% axis image
% saveas(gcf,['u_err',datestr(now,'mmddTHHMMSS')])
% % set(gcf,'Position',[3 3 15 6]);
%
% figure
% subplot(2,2,1)
% histogram(u_err{1}{1})
% title('u_{err} {1} {1}')
% subplot(2,2,2)
% histogram(u_err{1}{2})
% title('u_{err} {1} {2}')
% subplot(2,2,3)
% histogram(u_err{end}{1})
% title('u_{err} {2} {1}')
% subplot(2,2,4)
% histogram(u_err{end}{2})
% title('u_{err} {2} {2}')
% saveas(gcf,['u_err_hist',datestr(now,'mmddTHHMMSS')])

disp('error compute done')





