function [] = image_comp_generalDef(folder)

% folder = ['./generated_datasets'];
%compose image with many beads
format short

%specific strains or displacements
strORdisp = 's';

% folder = './2d_bead_images_testing/';
if exist(folder,'dir') ~= 7
    mkdir(folder);
end
%Make a new output folder if none exists
if exist(['./',folder(1:end-1),'_tifs/'],'dir') ~= 7
    mkdir(['./',folder(1:end-1),'_tifs/']);
end

disp('locations')
%get the original location for each bead
sizeI = [511,511,195]; %final image size less one
spacing = 4; %particle min. spacing
nPts = round(2500); %number of particles
try
    load([folder,'seed_points_3d_',num2str(nPts)])
catch
    x0 = poissonDisc(sizeI,spacing,nPts,0); %generate locations
    save([folder,'seed_points_3d_',num2str(nPts)],'x0')
end

% x0 = [100,sizeI(2)/2,sizeI(3)/2
%       200,sizeI(2)/2,sizeI(3)/2
%       300,sizeI(2)/2,sizeI(3)/2
%       400,sizeI(2)/2,sizeI(3)/2
%       500,sizeI(2)/2,sizeI(3)/2];

sn = 30; %signal to noise ratio, 10 compares reasonably to exp images.
bl = 5; %black level, i.e. background brightness onto which,
% a small amount of noise is superimposed on this
% gaussian noise is applied

%%

%final bead size (with 1x bead size margin on each side, undeformed)
bead_size = 15; %needs to be odd (?)
disp('image generation')

% %generate small bead images
% small_bead_img = small_bead_synth(defMat,bead_size);

%steps to progress the deformation through
dSteps = 0:0.3:1;

t0 = tic;
% I = cell(1,length(translation));
I = cell(1,length(dSteps));

%loop through to compose the images
for step = 1:length(dSteps)
    %     load(['bead_PSF_stretch_',num2str(cur_stretch),'.mat'])
    
    t1 = tic;
    
    cur_stretch = dSteps(step);
    
    I0_ = zeros(sizeI+1);
    x1 = zeros(size(x0));
    small_bead_img = cell(1,length(x0));
    
    %for each bead, compute the local deformation gradient (do in parallel)
    parfor ii = 1:length(x0)
        
        %compute the deformation matrix to impose on the bead
        if strORdisp == 's'
            
            defMat = strain_field(x0(ii,1),x0(ii,2),x0(ii,3),step);
            
        elseif strORdisp == 'd'
            
            [xl,yl] = meshgrid(linspace(x0(ii,1)-2,x0(ii,1)+2,10),...
                linspace(x0(ii,2)-2,x0(ii,2)+2,10));
            
            u_loc = disp_field(xl,yl,step,sizeI);
            
            %calculate the linearized local deformation gradient
            defMat = disp_grad_planefit(u_loc);
            
        end
        
        %generate small bead image
        small_bead_img{ii} = small_bead_synth_genDef_3d(defMat,bead_size);
        
        %decrease the bead intensity slightly with increasing volumetric strain
        small_bead_img{ii} = 1/(det(defMat))*small_bead_img{ii};
        
        x_ = x0(ii,1);
        y_ = x0(ii,2);
        z_ = x0(ii,3);
        %find the bead displacement
        u_inc = disp_field(x_,y_,z_,step,sizeI);
        
        if sum((abs(u_inc(:)))) == 0
            small_bead_img{ii} = 0*small_bead_img{ii};
        end
        
        %new bead locations
        x1(ii,:) = x0(ii,:) + u_inc;
        
    end
   
    %place bead into image
    I0_ = seedGenDefBeadSeries_3d(I0_,small_bead_img,sizeI,x1,bead_size);
    
    cT = toc(t1);
    disp(cT)
    
    I{step} = I0_;
    
    clk = clock;
    %     save(['bead_image_noPSF_output_lambda_',num2str(cur_stretch),'_',...
    %         num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    %         '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'cur_img')
    %
    %     save(['bead_image_noPSF_output_trans_',num2str(cur_trans),'_',...
    %         num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
    %         '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'I0')
    %
    I0 = I0_;
    save([folder,num2str(nPts),'noNoise_0p01step_lambda_',...
        num2str(cur_stretch,3),'.mat'],'I0')
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
    I1 = poissrnd(I0_);
    
    %renormalize to 256
    
    I0 = (255+bl)*I0/max(I0(:));
    I1 = (255+bl)*I1/max(I1(:));
    
    I0(I0 > (255)) = 255;
    I1(I1 > (255)) = 255;
    
    I0 = uint8(255*I0/max(I0(:)));
    I1 = uint8(255*I1/max(I1(:)));
    
    save([folder,num2str(nPts),'noisyImage_0p01step_lambda_',...
        num2str(cur_stretch,3),'.mat'],'I0')
    
%     imwrite(I0(1:512,1:512),[folder(1:end-1),'_tifs/',num2str(step+999),'_',...
%         'noisyImage_0p01step_lambda_',...
%         num2str(cur_stretch,3),'.tif'],...
%         'tif','Compression','none');
    
    save([folder,num2str(nPts),'beads_noisy_0p01step_lambda_',...
        num2str(cur_stretch,3),'_static.mat'],'I0','I1')
    
    toc(t1)
end
toc(t0)

disp('done')
clear *_
end


function [F_total] = disp_grad_planefit(u)

u_size = size(u{1}{1});
u_x = (0:(u_size(1)-1)) + 1;
u_y = (0:(u_size(2)-1)) + 1;

%Set up the data structures for surface fitting of u1a and u2
[fitting_params.xData, fitting_paramas.yData, ...
    fitting_params.zData_u2] = prepareSurfaceData( u_x(5:end-5),...
    u_y(5:end-5), rot90(u{ii}{2}(5:end-5,5:end-5),1));

[~, ~, fitting_params.zData_u1] = prepareSurfaceData( u_x(5:end-5),...
    u_y(5:end-5), rot90(u{ii}{1}(5:end-5,5:end-5),3));

% Set up fittype and options
ft = fittype( 'poly11' );
opts = fitoptions( 'Method', 'LinearLeastSquares' );
opts.Robust = 'Bisquare';

% Fit model to data
[fitting_params.fitresult_u2, fitting_params.gof_u2] = ...
    fit( [fitting_params.xData, fitting_params.yData], ...
    fitting_params.zData_u2, ft, opts );
[fitting_params.fitresult_u1, fitting_params.gof_u1] = ...
    fit( [fitting_params.xData, fitting_params.yData], ...
    fitting_params.zData_u1, ft, opts );

%extract fitted parameters of the planes
plane_coeffs{1}(1) = fitting_params.fitresult_u2.p10;
plane_coeffs{1}(2) = fitting_params.fitresult_u2.p01;
plane_coeffs{1}(3) = fitting_params.fitresult_u2.p00;

plane_coeffs{2}(1) = fitting_params.fitresult_u1.p10;
plane_coeffs{2}(2) = fitting_params.fitresult_u1.p01;
plane_coeffs{2}(3) = fitting_params.fitresult_u1.p00;

%Find incremental deformation gradient from the plane fit params
H_inc(1,1) = plane_coeffs{1}(1);
H_inc(1,2) = plane_coeffs{1}(2);
H_inc(2,1) = plane_coeffs{2}(1);
H_inc(2,2) = plane_coeffs{2}(2);

%Convert
F_total(1,1) = H_inc(1,1)+1;
F_total(1,2) = H_inc(1,2);
F_total(2,1) = H_inc(2,1);
F_total(2,2) = H_inc(2,2)+1;

end




function [u] = disp_field(x,y,z,step,sizeI)
%This function takes in the curent location, and returns the cartisean
%displacement at that location.

nu = 0.3;
E = 1e6;

%change the coord system
x = x - sizeI(1)/2;
y = y - sizeI(2)/2;
z = z - sizeI(2)/2;

a = 10;
r = sqrt(x.^2 + y.^2);
t = atan2(y,x);

P = step*10000;

Sinf = step*10000;
%
% %the line-load in a half-space.
% u_r = -(2*(1 - nu^2)/(pi*E))*P*cos(t).*log(r)-(1+nu)*(1-2*nu)/(pi*E)*P.*t.*sin(t);
% u_t = (2*(1 - nu^2)/(pi*E))*P.*sin(t).*log(r)+(1+nu)/(pi*E)*P*sin(t)-...
%     2*(1+nu)*(1-2*nu)/(pi*E)*P.*t.*cos(t);
%
% u_x = u_r*cos(u_t);
% u_y = u_r*sin(u_t);
%
% u = [u_x,u_y];

%the lame problem
if r >= a
    Srr_int = Sinf/2*(r*t-2*a*t+(r+4*a-(9*a^4/r^3)+(4*a^2/r))*sin(t)*cos(t));
    Stt_int = Sinf/2*(r*t-a^2/r-(r-(9*a^4/r^3)-10*a)*sin(t)*cos(t));
    Srt_int = Sinf/4*(r+(9*a^4/r^3)-a^2/r-8*a)*cos(2*t);
    
    u_r = 1/E*(Srr_int-nu*Stt_int);
    u_t = 1/E*(Stt_int-nu*Srr_int);
    u_z = 1/E*(-nu*(Srr_int+Stt_int));
%     urt = 1/E*(1+nu)*Srt_int;
else
    u_r = 0;
    u_t = 0;
    u_z = 0;
%     urt = 0;
end

u_x = u_r*cos(u_t);
u_y = u_r*sin(u_t);

u = [u_x,u_y,u_z];

end

function [e] = strain_field(x,y,z,step)
%this function take in the current location and step and give back the
%strain tensor.  Currently set up for the hole-plate (lame) problem.


nu = 0.3;
E = 1e6;
Sinf = step*1e3;
a = 10;

r = sqrt(x.^2 + y.^2);
t = atan2(y,x);
z = z;

if r >= a
    Srr = Sinf/2*(1-a^2/r^2) + Sinf/2*(1+3*a^4/r^4-4*a^2/r^2)*cos(2*t);
    Stt = Sinf/2*(1+a^2/r^2) - Sinf/2*(1+3*a^4/r^4)*cos(2*t);
    Srt = -Sinf/2*(1-3*a^4/r^4+2*a^2/r^2)*sin(2*t);
    
    err = 1/E*(Srr-nu*Stt);
    ett = 1/E*(Stt-nu*Srr);
    ert = 1/E*(1+nu)*Srt;
    ezz = 1/E*(-nu*(Srr+Stt));
    
    e = [err ert 0; ert ett 0; 0 0 ezz];
    
else
    
    for ii = 1:2
        for jj = 1:2
            e(ii,jj) = 0;
        end
    end
end



end

% clk = clock;
% save(['bead_image_65000p_output_',...
%     num2str(fix(clk(3))),'_',num2str(fix(clk(4))),...
%     '_',num2str(fix(clk(5))),'_',num2str(fix(clk(6))),'.mat'],'I','-v7.3')

% figure
% imshow3D(I0)





%================================old===================================

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


