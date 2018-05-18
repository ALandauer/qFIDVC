%% Generate test volume
clear all; close all;
%%
% desired standard deviation of bead
sigma       = [0.9 0.9 0.9];
% nSigma      = 3*sigma;

% desired output size.
sizeI           = [512,512,192];

% beads per pixel
seedingDensity  = 0.006;

% number of beads
nBeads          = round(seedingDensity*prod(sizeI));

% generate seed points from uniform distribution and scale to image size
x0{1} =  rand(nBeads,3);
x0{1} = [x0{1}(:,1)*(sizeI(1) - 1) + 1, ...
    x0{1}(:,2)*(sizeI(2) - 1) + 1,x0{1}(:,3)*(sizeI(3) - 1) + 1];
% [x0] = poissonDisc(sizeI,5,nBeads,1);
I0 = seedBeadsN(sigma,x0{1},sizeI);

%% DEFORM IMAGE
nDeformations = 8; % number of deformation increments

% stretch
L(1) = 1.05;
L(2) = 1/L(1);
L(3) = L(2);

e_mag = linspace(0,2.1,nDeformations);

I{1} = I0;
u = cell(1,length(2:nDeformations));% zeros(size(x0));
for i = 2:nDeformations
    fprintf('Deformation step: %i / %i \n',i,nDeformations)
    
    step = e_mag(i);
    
%     xC = (sizeI+1)/2; % center of deformation

%     u{i-1}(:,1) = step*(L(1) - 1)*(x0{1}(:,1) - xC(1));
%     u{i-1}(:,2) = step*(L(2) - 1)*(x0{1}(:,2) - xC(2));
%     u{i-1}(:,3) = step*(L(3) - 1)*(x0{1}(:,3) - xC(3));
    
%     x1{i-1} = x0{1} + u{i-1};
    
    x = x0{1}(:,1);
    y = x0{1}(:,2);
    z = x0{1}(:,3);
    
    I{i} = zeros(sizeI);
    for ii = 1:length(x)
        u{i-1}(ii,:) = disp_field(x(ii),y(ii),z(ii),step,sizeI);
    end
    x1{i-1} = x0{1} + u{i-1};
    I{i} = seedBeadsN(sigma,x1{i-1},sizeI);
end

for i = 1:length(I)
vol{1} = I{i};
save(['vol_series_',num2str(1000+i),'.mat'],'vol');
end
save(['imposed_disp_series','.mat'],'u','x0');

% vol{1} = I{2};
% save('vol01.mat','vol')

for ii = 1:length(I)
    fig
    imagesc3d(I{ii}),colorbar
end








