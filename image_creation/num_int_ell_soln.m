%This function computes the convolution of the bead intensity and the point
%spread function of pinhole optics

clear all
close all
clc

%% CONVOLUTION
stretch = 0.5; %uniaxial stretch level

% for stretch = 1:-0.15:0.25

% set run parameters
L = 532*10^-9; %incident wavelength

% d^2/l
I_0 = 255; %reference intensity

r = 0.3; %bead radius

dz = 1.2; %spacing of image plane from pinhole

d = sqrt(0.61*L*dz);   %pinhole radius, based on Young, 1971 "pinhole optics"


%compute semi-axis lengths
a(1) = r*stretch;
a(2) = r/sqrt(stretch);
a(3) = r/sqrt(stretch);

%compute the convolution
st = -3.0;
sp = 0.2;
ed = 3.0;

pool = gcp('nocreate');

if isempty(pool)
    curCluster = parcluster('local');
    curCluster.NumWorkers = 16;
    saveProfile(curCluster);
    pool = parpool(16);
end

sigma = [0.5,0.5,1]; %3d gaussian distros

k = 500;

t0 = tic;
t1 = tic;
v = st:sp:ed;
I_field = zeros(4,length(st:sp:ed)*length(st:sp:ed)*length(st:sp:ed));
I_mat = zeros(length(v),length(v),length(v));

cnt = 0;
for ii = 1:length(v) %st:sp:ed
    m = v(ii);
    fprintf('\n %i of %i and loop time %f',cnt,length(I_field),toc(t1));
    t1 = tic;
    for jj = 1:length(v) %st:sp:ed
        n = v(jj);
%         t2 = tic;
        parfor kk = 1:length(v) %st:sp:ed
            
            p = v(kk)+dz;
            FofI = @(x,y,z) I_0*exp(-(x.^2 + y.^2 + z.^2)./...
                (0.3528*((L/d).*(sqrt(x.^2+y.^2+z.^2))).^2))./...
                (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));
%             FofI = @(x,y,z) bead_img_conv(x,y,z,m,n,p,a,d,l,I_0,sigma);
            
            I_cur = integral3(FofI,-1.5*sigma(1),1.5*sigma(1),-1.5*sigma(2),1.5*sigma(2),...
                -(1.5*sigma(3)+dz),(1.5*sigma(3)+dz));%,'AbsTol',5e-4,'RelTol',5e-3);

%             I_cur = triplequad(FofI,-1.5*sigma(1),1.5*sigma(1),-1.5*sigma(2),1.5*sigma(2),...
%                 -(1.5*sigma(3)+dz),(1.5*sigma(3)+dz));
            
%             for qq = 1:8
%                 I_cur_(qq) = integral3(FofI,o{qq});
%             end
%             I_cur = sum(I_cur_);

%             I_cur = triplequad(FofI,st,ed,st,ed,st+5,ed+5);
            
            
            I_field(:,cnt+kk) = [m,n,p,I_cur];
            
            I_mat(ii,jj,kk) = I_cur;
            
        end
%         toc(t2)
        cnt = cnt+length(v);
    end
end
fprintf('\nTotal integration time: %f\n',toc(t0))

clk = clock;
% save(['bead_image_output_lambda_',num2str(stretch),'_',num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
%     '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'])

% end

%% PLOTTING

% m = -1:0.1:1;
% n=-1:0.1:1;
% p=-1:0.1:1;
% cnt = 0;
% for ii = 1:length(m)
%     for jj = 1:length(n)
%         for kk = 1:length(p)
%             cnt = cnt+1;
% bead_i = 1./(1+exp(-2*k*(1-((m(ii)).^2/a(1)^2)-((n(jj)).^2/a(2)^2)-((p(kk)).^2/a(3)^2))));
%         bead_field(:,cnt) = [ii,jj,kk,bead_i];
%         end
%     end
% end
%
% clear I_field

% I_field = bead_field;

figure
imshow3D(I_mat,[])

%set up for plotting
xp = I_field(1,:)';
yp = I_field(2,:)';
zp = I_field(3,:)'-dz;
Ii = I_field(4,:)'/max(I_field(4,:));

Ixy = Ii(zp.^2==min(zp.^2)); %intensity on the xy plane
Ixz = Ii(yp.^2==min(abs(yp).^2)); %on xz
Iyz = Ii(xp.^2==min(abs(xp).^2)); %on yz

xc2 = yp(xp.^2==min(abs(xp).^2)); %the y-coords of points on the x=0 plane
xc3 = zp(xp.^2==min(abs(xp).^2)); %z-coords of x=0

yc1 = xp(yp.^2==min(abs(yp).^2)); %the x-coords of points on the y=0 plane
yc3 = zp(yp.^2==min(abs(yp).^2));%z-coords

zc1 = xp(zp.^2==min(abs(zp).^2)); %the x-coords of points on the z=0 plane
zc2 = yp(zp.^2==min(abs(zp).^2)); %y-coords

% Ixy = Ii(zp==-3.05); %intensity on the xy plane
% Ixz = Ii(yp==-3.05); %on xz
% Iyz = Ii(xp==-3.05); %on yz
% 
% xc2 = yp(xp==-3.05); %the y-coords of points on the x=0 plane
% xc3 = zp(xp==-3.05); %z-coords of x=0
% 
% yc1 = xp(yp==-3.05); %the x-coords of points on the y=0 plane
% yc3 = zp(yp==-3.05); %z-coords
% 
% zc1 = xp(zp==-3.05); %the x-coords of points on the z=0 plane
% zc2 = yp(zp==-3.05); %y-coords

[x_cord,y_cord] = meshgrid(min(zc1):sp:max(zc1),min(zc2):sp:max(zc2));
Ixy_grid_ = griddata(zc1,zc2,Ixy,x_cord,y_cord,'cubic');

[x_cord,z_cord] = meshgrid(min(yc1):sp:max(yc1),min(yc3):sp:max(yc3));
Ixz_grid_ = griddata(yc1,yc3,Ixz,x_cord,z_cord,'cubic');

[y_cord,z_cord] = meshgrid(min(xc2):sp:max(xc2),min(xc3):sp:max(xc3));
Iyz_grid_ = griddata(xc2,xc3,Iyz,y_cord,z_cord,'cubic');

%pad for better viewing
Ixy_grid = Ixy_grid_;%padarray(Ixy_grid_,[4,4],0);
Ixz_grid = Ixy_grid_;%padarray(Ixz_grid_,[4,4],0);
Iyz_grid = Ixy_grid_;%padarray(Iyz_grid_,[4,4],0);

%from the 'BuRd' color map of "BrewerMap" from the File Exchange
% https://www.mathworks.com/matlabcentral/fileexchange/45208-colorbrewer--
%   attractive-and-distinctive-colormaps
colors = [0.4039    0.0000    0.1216
    0.6801    0.0627    0.1607
    0.8169    0.3268    0.2720
    0.9302    0.5762    0.4449
    0.9876    0.7904    0.6813
    0.9907    0.9392    0.9073
    0.9156    0.9520    0.9741
    0.7412    0.8585    0.9179
    0.4869    0.7230    0.8424
    0.2210    0.5427    0.7503
    0.1247    0.3834    0.6603
    0.0196    0.1882    0.3804];

I_extrema = [min(Ii(:)),max(Ii(:))];

%do the contour plots
% figure
% contour3()

figure
sizeOut = size(I_mat);
hFigRotated = figure;
hAxRotated  = axes;
slice(double(I_mat),sizeOut(2)/2,sizeOut(1)/2,sizeOut(3)/2);
grid on, shading interp, colormap gray, axis image

figure
subplot(3,1,1);
hold on
contourf(x_cord,z_cord,Ixy_grid);%,'LineStyle','none')
title('I on x-y');
colormap(colors);
caxis(I_extrema);
colorbar
axis equal

subplot(3,1,2);
hold on
contourf(x_cord,z_cord,Ixz_grid);%,'LineStyle','none')
title('I on x-z');
colormap(colors);
caxis(I_extrema);
colorbar
axis equal

subplot(3,1,3);
hold on
contourf(y_cord,z_cord,Iyz_grid);%,'LineStyle','none')
title('I on y-z');
colormap(colors);
caxis(I_extrema);
colorbar
axis equal


% % %% testing
% % clear all
% % % tic
% % l = 532*10^-9;
% % d = 1*10^-8;
% % % d^2/l
% % k = 10000; %large number
% % I_0 = 1; %reference intensity
% % 
% % r = 0.5;
% % stretch = 0.5;
% % 
% % a(1) = r*stretch;
% % a(2) = r/sqrt(stretch);
% % a(3) = r/sqrt(stretch);
% % 
% % % m = 0;
% % % n = 0;
% % % p = 0;
% % % cnt = 0;
% % % for m = -1:0.2:1
% % %     for n = -1:0.2:1
% % %         for p = -1:0.2:1
% % %             cnt = cnt+1
% % % FofI = @(x,y,z) I_0*exp(-(x.^2 + y.^2)./...
% % %                 (0.3528*((l/d).*(sqrt(x.^2+y.^2+z.^2))).^2))./...
% % %                 (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));
% % 
% % m = 0;
% % n = 0;
% % p = 0;
% % sigma = [0.5,0.5,1];
% % % FofI = @(x,y,z) exp(-((x.^2)/(2*sigma(1)^2) + y.^2/(2*sigma(2)^2) +z.^2./(2*sigma(3)^2)));
% % % FofI = @(x,y,z) exp(-(x.^2 + y.^2)./(0.3528*((l/d).*(sqrt(x.^2+y.^2+z.^2))).^2));
% % FofI = @(x,y,z) I_0*exp(-((x.^2)/(2*sigma(1)^2) + y.^2/(2*sigma(2)^2) +z.^2./(2*sigma(3)^2)))./...
% %     (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));
% % 
% % % FofI = @(x,y,z) bead_img_conv(x,y,z,m,n,p,a,d,l,I_0);
% % 
% % 
% % m = -2:0.1:2;
% % n = -2:0.1:2;
% % p = -2:0.1:2;
% % for ii = 1:length(m)
% %     for jj = 1:length(n)
% %         for kk = 1:length(p)
% %             bead(ii,jj,kk) = 1./(1+exp(-2*k*(1-((m(ii)).^2/a(1)^2)-((n(jj)).^2/a(2)^2)-((p(kk)).^2/a(3)^2))));
% %         end
% %     end
% % end
% % 
% % ic = 0;
% % jc = 0;
% % kc = 0;
% % for ii = -1:0.05:1
% %     ic = ic+1;
% %     jc = 0;
% %     for jj = -1:0.05:1
% %         jc = jc+1;
% %         kc = 0;
% %         for kk = -1:0.05:1
% %             kc = kc+1;
% %             I_(ic,jc,kc) = FofI(ii,jj,kk);
% %         end
% %     end
% % end
% % 
% % % f(cnt) = sum(sum(sum(isnan(I_))));
% % %         end
% % %     end
% % % end
% % % xv = -1:0.01:1;
% % % yv = -1:0.01:1;
% % % zv = -1:0.01:1;
% % % toc
% % 
% % figure
% % sizeOut = size(bead);
% % hFigRotated = figure;
% % hAxRotated  = axes;
% % slice(double(bead),sizeOut(2)/2,sizeOut(1)/2,sizeOut(3)/2);
% % grid on, shading interp, colormap gray, axis image
% % 
% % figure
% % imshow3D(I_,[])
% % % surf(xv,yv,I_(:,:,101),'linestyle','none'),axis image
% % %
% % %
% % %

% % % % % % % t0 = tic;
% % % % % % % t1 = tic;
% % % % % % % for m = st:sp:ed
% % % % % % %     for n = st:sp:ed
% % % % % % %         for p = st:sp:ed
% % % % % % %             cnt  = cnt+1;
% % % % % % %
% % % % % % %             %             FofI = @(x,y,z) I_0*exp(-(x.^2 + y.^2)./...
% % % % % % %             %                 (0.3528*((l/d).*(sqrt(x.^2+y.^2+z.^2))).^2))./...
% % % % % % %             %                 (1+exp(-2*k*(1-((m-x).^2/a(1)^2)-((n-y).^2/a(2)^2)-((p-z).^2/a(3)^2))));
% % % % % % %             FofI = @(x,y,z) bead_img_conv(x,y,z,m,n,p,a,d,l,I_0);
% % % % % % %
% % % % % % %
% % % % % % %             I_cur = integral3(FofI,st,ed,st,ed,st,ed);
% % % % % % %
% % % % % % %             I_field(:,cnt) = [m,n,p,I_cur];
% % % % % % %
% % % % % % %             if ~mod(cnt,100)
% % % % % % %                 fprintf('\n %i of %i and loop time %f',cnt,length(I_field),toc(t1))
% % % % % % %                 t1 = tic;
% % % % % % %             end
% % % % % % %
% % % % % % %         end
% % % % % % %     end
% % % % % % % end
% % % % % % % lower = [-1.0001,-0.0001];
% % % % % % % upper = [.0001,1.0001];
% % % % % % % 
% % % % % % % o{1} = [upper,upper,upper];
% % % % % % % o{2} = [lower,upper,upper];
% % % % % % % o{3} = [upper,lower,upper];
% % % % % % % o{4} = [upper,upper,lower];
% % % % % % % o{5} = [lower,lower,upper];
% % % % % % % o{6} = [lower,upper,lower];
% % % % % % % o{7} = [upper,lower,lower];
% % % % % % % o{8} = [lower,lower,lower];
