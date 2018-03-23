%% Generate test volume
clear all; close all;
%%
% desired standard deviation of bead
sigma       = [1 1 1];
% nSigma      = 3*sigma;

% desired output size.
sizeI           = [512,512,192];

% beads per pixel
seedingDensity  = 0.005;

% number of beads
nBeads          = round(seedingDensity*prod(sizeI));

% generate seed points from uniform distribution and scale to image size
x0 =  rand(nBeads,3);
x0 = [x0(:,1)*(sizeI(1) - 1) + 1, x0(:,2)*(sizeI(2) - 1) + 1,x0(:,3)*(sizeI(3) - 1) + 1];
% [x0] = poissonDisc(sizeI,5,nBeads,1);
I0 = seedBeadsN(sigma,x0,sizeI);

%% DEFORM IMAGE
nDeformations = 5; % number of deformation increments



% stretch
`L(1) = 1.01;
L(2) = 1/L(1);
L(3) = L(2);


I{1} = I0;
u = zeros(size(x0));
for i = 2:nDeformations
    fprintf('Deformation step: %i / %i \n',i,nDeformations)
    
    step = (i)/nDeformations;
    
    xC = (sizeI+1)/2; % center of deformation

    u(:,1) = step*(L(1) - 1)*(x0(:,1) - xC(1));
    u(:,2) = step*(L(2) - 1)*(x0(:,2) - xC(2));
    u(:,3) = step*(L(3) - 1)*(x0(:,3) - xC(3));
    
    x1 = x0 + u;
    I{i} = seedBeadsN(sigma,x1,sizeI);
end

for i = 1:length(I)
vol{i} = I{i};
save(['vol0',num2str(i),'.mat','vol']);
end

% vol{1} = I{2};
% save('vol01.mat','vol')









