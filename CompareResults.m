%% Compare displacement

sizeI0 = [512,512,101];
% % % sSize = [64 64 32];
% % % 
% % % %Predict the displacement grid points :P
% % % sSpacing = sSize/2; 
% % % sizeI = ceil(sizeI0./sSpacing).*sSpacing;
% % % prePad = ceil((sizeI - sizeI0)/2);
% % % postPad = floor((sizeI - sizeI0)/2);
% % % 
% % % for i = 1:3, m{i} = 1:dm:(sizeI(i)+1); end
% % % 
% % % [m{:}] = meshgrid(m{:});

load('resultsFIDVC.mat')
% Compute displacement at the grid points
xC = ([513,513,193])/2; % center of deformation
% stretch
L(1) = 1.01;
L(2) = 1/L(1);
L(3) = L(2);
step = 1;
v{1} = step*(L(1) - 1)*(m{1} - xC(1));
v{2} = step*(L(2) - 1)*(m{2} - xC(2));
v{3} = step*(L(3) - 1)*(m{3} - xC(3));


temp = u{1}{2} - v{1};

figure;
contourf(squeeze(temp(6:end-5,6:end-5,7)))
axis image
colorbar