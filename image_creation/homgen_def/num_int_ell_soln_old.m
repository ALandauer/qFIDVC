%This function computes the convolution of the bead intensity and the point
%spread function of pinhole optics

clear all
close all
clc

addpath('/gpfs/data/cf5/alandaue/particle_deform')
cd /gpfs/data/cf5/alandaue/particle_deform

%% CONVOLUTION
% set run parameters
l = 532*10^-9; %incident wavelength
d = 1*10^-8;   %pinhole radius
% d^2/l
I_0 = 255; %reference intensity

r = 0.2; %bead radius
stretch = 1; %uniaxial stretch level

%compute semi-axis lengths
a(1) = r*stretch;
a(2) = r/sqrt(stretch);
a(3) = r/sqrt(stretch);

%compute the convolution
st = -1.0;
sp = 0.1;
ed = 1.0;
% I_field = zeros(4,length(st:sp:ed)*length(st:sp:ed)*length(st:sp:ed));
cnt = 0;
t0 = tic;
t1 = tic;
% for m = st:sp:ed
%     for n = st:sp:ed
%         for p = st:sp:ed
%             cnt  = cnt+1;
%
%             %             FofI = @(x,y,z) I_0*exp(-(x.^2 + y.^2)./...
%             %                 (0.3528*((l/d).*(sqrt(x.^2+y.^2+z.^2))).^2))./...
%             %                 (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));
%             FofI = @(x,y,z) bead_img_conv(x,y,z,m,n,p,a,d,l,I_0);
%
%
%             I_cur = integral3(FofI,st,ed,st,ed,st,ed);
%
%             I_field(:,cnt) = [m,n,p,I_cur];
%
%             if ~mod(cnt,100)
%                 fprintf('\n %i of %i and loop time %f',cnt,length(I_field),toc(t1))
%                 t1 = tic;
%             end
%
%         end
%     end
% end



curCluster = parcluster('local');
curCluster.NumWorkers = 16;
saveProfile(curCluster);
pool = parpool(16);

t0 = tic;
t1 = tic;
v = st:sp:ed;
I_field = zeros(4,length(st:sp:ed)*length(st:sp:ed)*length(st:sp:ed));
for ii = 1:length(v) %st:sp:ed
    m = v(ii);
    fprintf('\n %i of %i and loop time %f',cnt,length(I_field),toc(t1));
    t1 = tic;
    for jj = 1:length(v) %st:sp:ed
        n = v(jj);
        parfor kk = 1:length(v) %st:sp:ed

            p = v(kk);
            FofI = @(x,y,z) bead_img_conv(x,y,z,m,n,p,a,d,l,I_0);
            
            I_cur = integral3(FofI,st,ed,st,ed,st,ed);
            
            I_field(:,cnt+kk) = [m,n,p,I_cur];
            
        end
        cnt = cnt+length(v);
    end
end
fprintf('\nTotal integration time: %f\n',toc(t0))

save beam_image_output.mat

%% PLOTTING
% % 
% % %set up for plotting
% % xp = I_field(1,:)';
% % yp = I_field(2,:)';
% % zp = I_field(3,:)';
% % Ii = I_field(4,:)';
% % 
% % Ixy = Ii(zp==min(abs(zp))); %intensity on the xy plane
% % Ixz = Ii(yp==min(abs(yp))); %on xz
% % Iyz = Ii(xp==min(abs(xp))); %on yz
% % 
% % xc2 = yp(xp==min(abs(xp))); %the y-coords of points on the x=0 plane
% % xc3 = zp(xp==min(abs(xp))); %z-coords of x=0
% % 
% % yc1 = xp(yp==min(abs(yp))); %the x-coords of points on the y=0 plane
% % yc3 = zp(yp==min(abs(yp))); %z-coords
% % 
% % zc1 = xp(zp==min(abs(zp))); %the x-coords of points on the z=0 plane
% % zc2 = yp(zp==min(abs(zp))); %y-coords
% % 
% % % Ixy = Ii(zp==0.5); %intensity on the xy plane
% % % Ixz = Ii(yp==0.5); %on xz
% % % Iyz = Ii(xp==0.5); %on yz
% % % 
% % % xc2 = yp(xp==0.5); %the y-coords of points on the x=0 plane
% % % xc3 = zp(xp==0.5); %z-coords of x=0
% % % 
% % % yc1 = xp(yp==0.5); %the x-coords of points on the y=0 plane
% % % yc3 = zp(yp==0.5); %z-coords
% % % 
% % % zc1 = xp(zp==0.5); %the x-coords of points on the z=0 plane
% % % zc2 = yp(zp==0.5); %y-coords
% % 
% % [x_cord,y_cord] = meshgrid(min(zc1):sp:max(zc1),min(zc2):sp:max(zc2));
% % Ixy_grid_ = griddata(zc1,zc2,Ixy,x_cord,y_cord,'cubic');
% % 
% % [x_cord,z_cord] = meshgrid(min(yc1):sp:max(yc1),min(yc3):sp:max(yc3));
% % Ixz_grid_ = griddata(yc1,yc3,Ixz,x_cord,z_cord,'cubic');
% % 
% % [y_cord,z_cord] = meshgrid(min(xc2):sp:max(xc2),min(xc3):sp:max(xc3));
% % Iyz_grid_ = griddata(xc2,xc3,Iyz,y_cord,z_cord,'cubic');
% % 
% % %pad for better viewing
% % Ixy_grid = Ixy_grid_;%padarray(Ixy_grid_,[4,4],0);
% % Ixz_grid = Ixy_grid_;%padarray(Ixz_grid_,[4,4],0);
% % Iyz_grid = Ixy_grid_;%padarray(Iyz_grid_,[4,4],0);
% % 
% % %from the 'BuRd' color map of "BrewerMap" from the File Exchange
% % % https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer--
% % %   attractive-and-distinctive-colormaps
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
% % I_extrema = [min(Ii(:)),max(Ii(:))];
% % 
% % %do the contour plots
% % % figure
% % % contour3()
% % 
% % 
% % figure
% % subplot(3,1,1)
% % hold on
% % contourf(x_cord,z_cord,Ixy_grid)%,'LineStyle','none')
% % title('I on x-y')
% % colormap(colors)
% % caxis(I_extrema)
% % colorbar
% % axis equal
% % 
% % subplot(3,1,2)
% % hold on
% % contourf(x_cord,z_cord,Ixz_grid)%,'LineStyle','none')
% % title('I on x-z')
% % colormap(colors)
% % caxis(I_extrema)
% % colorbar
% % axis equal
% % 
% % subplot(3,1,3)
% % hold on
% % contourf(y_cord,z_cord,Iyz_grid)%,'LineStyle','none')
% % title('I on y-z')
% % colormap(colors)
% % caxis(I_extrema)
% % colorbar
% % axis equal
% % 

%% testing
% clear all
% tic
% l = 532*10^-9;
% d = 1*10^-8;
% % d^2/l
% k = 100; %large number
% I_0 = 1; %reference intensity
%
% r = 0.5;
% stretch = 1;
%
% a(1) = r*stretch;
% a(2) = r/sqrt(stretch);
% a(3) = r/sqrt(stretch);
%
% m = 0.1;
% n=0.1;
% p=0.1;
%
% FofI = @(x,y,z) I_0*exp(-(x.^2 + y.^2)./...
%                 (0.3528*((l/d).*(sqrt(x.^2+y.^2+z.^2))).^2))./...
%                 (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));
%
% % FofI = @(x,y,z) bead_img_conv(x,y,z,m,n,p,a,d,l,I_0);
%
% ic = 0;
% jc = 0;
% kc = 0;
% for ii = -1:0.01:1
%     ic = ic+1
%     jc = 0;
%     for jj = -1:0.01:1
%         jc = jc+1;
%         kc = 0;
%         for kk = -1:0.01:1
%             kc = kc+1;
%             I_(ic,jc,kc) = FofI(ii,jj,kk);
%         end
%     end
% end
%
% xv = -1:0.01:1;
% yv = -1:0.01:1;
% zv = -1:0.01:1;
% toc
%
% surf(xv,yv,I_(:,:,101)),axis image
%
%
%
