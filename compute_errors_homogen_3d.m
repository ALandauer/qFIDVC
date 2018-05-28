% clear,close all

load('results_qDVC_uniform_e3_inc')

if runMode(1) == 'i'
    u_c = inc2cum(u,dm,m,'cubic');
end

%%
bndry_wd = 0;
[u_masked_hyb,u_imp_masked_hyb,u_err_hyb,l2err_hyb,std_err_hyb,mean_err_hyb,imp_strain_hyb,nan_mask]...
    = error_mapping_homogen_3d(u_c,dm,m,bndry_wd);
1
%%
img = 1;
comp = 1;
fig
imagesc3d(nan_mask.*u_err_hyb{img}{comp})
%%
fig
subplot(1,2,2)
plot(std_err_hyb)
title('st dev')
legend('1','2','3')

subplot(1,2,1)
plot(mean_err_hyb)
title('mean')

% % %%
% % fig
% % plane = 5;
% % subplot(1,2,1)
% % imagesc(nan_mask(:,:,plane).*u_masked_hyb{img}{comp}(:,:,plane)),colorbar,axis image
% % subplot(1,2,2)
% % imagesc(nan_mask(:,:,plane).*u_imp_masked_hyb{img}{comp}(:,:,plane)),colorbar,axis image
% % %%
% % fig
% % plane = 22;
% % subplot(1,2,1)
% % imagesc(squeeze(nan_mask(:,plane,:)).*squeeze(u_masked_hyb{img}{comp}(:,plane,:))),colorbar,axis image
% % subplot(1,2,2)
% % imagesc(squeeze(nan_mask(:,plane,:)).*squeeze(u_imp_masked_hyb{img}{comp}(:,plane,:))),colorbar,axis image
% % %%
% % fig
% % plane = 22;
% % subplot(1,2,1)
% % imagesc(squeeze(nan_mask(plane,:,:)).*squeeze(u_masked_hyb{img}{comp}(plane,:,:))),colorbar,axis image
% % subplot(1,2,2)
% % imagesc(squeeze(nan_mask(plane,:,:)).*squeeze(u_imp_masked_hyb{img}{comp}(plane,:,:))),colorbar,axis image




