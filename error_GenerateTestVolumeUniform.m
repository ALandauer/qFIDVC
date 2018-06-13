
%% load data
clear
imposed = load('imposed_disp_series.mat');
% measured = load('results_qDVC_ss16_inc.mat');
measured = load('results_qDVC_uniform_side_ss32_inc.mat');

a = 25;
dm = 8;
sizeI = [512,512,192];
%% map imposed and measured data onto same mesh

% grid = permute(measured.m,[213]);
% [gridPoints_m{1},gridPoints_m{2},gridPoints_m{3}] = meshgrid(measured.m{1},...
%     measured.m{2},measured.m{3});
[gridPoints{1},gridPoints{2},gridPoints{3}] = ndgrid(measured.m{1}(:,1,1),...
    measured.m{2}(1,:,1),measured.m{3}(1,1,:));

gridPoints_m{1} = measured.m{2};
gridPoints_m{2} = measured.m{1};
gridPoints_m{3} = measured.m{3};
e_mag = linspace(0,2.1,8);

for ii= 1:length(measured.u)
    
    
    % %     for xx = 1:size(gridPoints_m{1},1)
    % %         for yy = 1:size(gridPoints_m{2},2)
    % %             for zz = 1:size(gridPoints_m{3},3)
    % %                 [u] = disp_field(gridPoints_m{1}(xx,yy,zz),gridPoints_m{2}(xx,yy,zz),gridPoints_m{3}(xx,yy,zz),ii+1,sizeI);
    % %                 imposed.u_grid{ii}{1}(xx,yy,zz) = u(1);
    % %                 imposed.u_grid{ii}{2}(xx,yy,zz) = u(2);
    % %                 imposed.u_grid{ii}{3}(xx,yy,zz) = u(3);
    % %
    % %             end
    % %         end
    % %     end
    
    for jj = 1:3
        imposed.u_grid{ii}{jj} = griddata(imposed.x0{1}(:,1),imposed.x0{1}(:,2), ...
            imposed.x0{1}(:,3),-imposed.u{ii}(:,jj),...
            gridPoints_m{1},gridPoints_m{2},gridPoints_m{3});
        
        nan_mask{ii}{jj} = measured.u{ii}{jj}./measured.u{ii}{jj};
% %         nan_mask{ii}{jj}(round(size(measured.u{ii}{jj},1)/2-a/dm):round(size(measured.u{ii}{jj},1)/2+a/dm),...
% %             round(size(measured.u{ii}{jj},2)/2-a/dm):round(size(measured.u{ii}{jj},2)/2+a/dm),...
% %             round(size(measured.u{ii}{jj},3)/2-a/dm):round(size(measured.u{ii}{jj},3)/2+a/dm)) = nan;
        
    end
    
    imposed.u_grid{ii}{4} = sqrt(imposed.u_grid{ii}{1}.^2 + ...
        imposed.u_grid{ii}{2}.^2 + imposed.u_grid{ii}{3}.^2);
    
end

% imposed.u_grid = inc2cum(imposed.u_grid,measured.dm,gridPoints,'spline');
measured.u_c = inc2cum(measured.u,measured.dm,gridPoints,'linear');
1
%% plotting

comp = 2;
img = 7;

fig
% subplot(1,2,1)
imagesc3d(nan_mask{img}{1}.*imposed.u_grid{img}{comp})
% subplot(1,2,2)
fig
imagesc3d(nan_mask{img}{1}.*measured.u_c{img}{comp})

fig
imagesc3d(nan_mask{img}{1}.*measured.u_c{img}{comp} - nan_mask{img}{1}.*imposed.u_grid{img}{comp})

2

