function I_bead = seedDefBeads(cur_bead_img,sizeI,x1_,bead_size_small)
%This function places beads defined by the volume matrix cur_bead_img, into
%and image of size sizeI, at the location specified by x1_, the bead is
%downsampled to the size given in bead_size.

composed_img = zeros(sizeI(1)+2*bead_size_small,sizeI(2)+2*bead_size_small,sizeI(3)+2*bead_size_small);

sizeB = size(cur_bead_img);

%only keep beads that stay within the image
x1_exclude_(:,1) = (x1_(:,1))>sizeI(1) | (x1_(:,1))<1;
x1_exclude_(:,2) = (x1_(:,2))>sizeI(2) | (x1_(:,2))<1;
x1_exclude_(:,3) = (x1_(:,3))>sizeI(3) | (x1_(:,3))<1;

x1_exclude = x1_exclude_(:,1) | x1_exclude_(:,2) | x1_exclude_(:,3);

x1(:,1) = x1_(~x1_exclude,1);
x1(:,2) = x1_(~x1_exclude,2);
x1(:,3) = x1_(~x1_exclude,3);

%set up parallel pool
pool = gcp('nocreate');
if isempty(pool)
    curCluster = parcluster('local');
    curCluster.NumWorkers = 6;
    saveProfile(curCluster);
    pool = parpool(6);
end

[px_x,px_y,px_z] = ndgrid(1:sizeB(1),1:sizeB(2),1:sizeB(3));

small_bead = cell(1,size(x1,1));
parfor bead_num = 1:size(x1,1)
    
    
    
    %find the center loc
    center = x1(bead_num,:) + bead_size_small;
    %     center = [10.5,10.5,10.5];
    
    
    center_int = floor(center); %interger loc
    center_offset = (center - center_int)*(sizeB(1)/bead_size_small); %decimal part
    
    %interpolate the bead image onto the pixel grid
    cur_bead_shift = interpn(px_x+center_offset(1),px_y+center_offset(2),...
        px_z+center_offset(3),cur_bead_img,px_x,px_y,px_z,'linear');
    
    %     cur_bead_shift = mirt3D_mexinterp(cur_bead_img,px_x+center_offset(1),px_y+center_offset(2),...
    %         px_z+center_offset(3));
    
    % % % %     %offset in each dir
    % % % %     del_x = center_offset(1);%s/(num_px/final_res);
    % % % %     del_y = center_offset(2);%/(num_px/final_res);
    % % % %     del_z = center_offset(3);%/(num_px/final_res);
    % % % %
    %build the summation kernal
    kernel = ones(floor(sizeB(1)/bead_size_small)+1,floor(sizeB(2)/bead_size_small)+1,...
        floor(sizeB(3)/bead_size_small)+1);
    % % % %     kernel(1,:,:) = del_x;
    % % % %     kernel(end,:,:) = 1-del_x;
    % % % %
    % % % %     kernel(:,1,:) = del_y;
    % % % %     kernel(:,end,:) = 1-del_y;
    % % % %
    % % % %     kernel(:,:,1) = del_z;
    % % % %     kernel(:,:,end) = 1-del_z;
    % % % %
    % % % %     %the plan is to find the average intensity in each pixel, wieghted by
    % % % %     %the kernel function, to corrected account for the sub-pixel shift.
    % % % %
    %do the summation
    x=0;
    small_bead{bead_num} = zeros(bead_size_small,bead_size_small,bead_size_small);
    for ii = 1:size(kernel,1)-1:size(cur_bead_shift,1)-1
        x=x+1;
        y=0;
        for jj = 1:size(kernel,2)-1:size(cur_bead_shift,2)-1
            y=y+1;
            z=0;
            for kk = 1:size(kernel,3)-1:size(cur_bead_shift,3)-1
                z=z+1;
                %                 small_bead{bead_num}(x,y,z) = ...
                %                     sum(sum(sum(kernel.*cur_bead_shift(ii:ii+size(kernel,1)-1,...
                %                     jj:jj+size(kernel,2)-1,kk:kk+size(kernel,3)-1))))/numel(kernel);
                
                small_bead{bead_num}(x,y,z) = ...
                    sum(sum(sum(cur_bead_shift(ii:ii+size(kernel,1)-1,...
                    jj:jj+size(kernel,2)-1,kk:kk+size(kernel,3)-1))))/numel(kernel);
                
            end
        end
    end
    
end

for bead_num = 1:size(x1,1)
    
    %find the center loc
    center = x1(bead_num,:) + bead_size_small;
    %     center = [10.5,10.5,10.5];
    
    center_int = floor(center); %interger loc
    
    %     %put the bead into the image
    %     composed_img(center_int(1)-floor(bead_size/2):center_int(1)+floor(bead_size/2),...
    %         center_int(2)-floor(bead_size/2):center_int(2)+floor(bead_size/2),...
    %         center_int(3)-floor(bead_size/2):center_int(3)+floor(bead_size/2)) = ...
    %         composed_img(center_int(1)-floor(bead_size/2):center_int(1)+floor(bead_size/2),...
    %         center_int(2)-floor(bead_size/2):center_int(2)+floor(bead_size/2),...
    %         center_int(3)-floor(bead_size/2):center_int(3)+floor(bead_size/2)) + small_bead{bead_num};
    
    %add the bead into the image
    composed_img(center_int(1)-fix(bead_size_small/2):center_int(1)+floor(bead_size_small/2),...
        center_int(2)-fix(bead_size_small/2):center_int(2)+floor(bead_size_small/2),...
        center_int(3)-fix(bead_size_small/2):center_int(3)+floor(bead_size_small/2)) = ...
        composed_img(center_int(1)-fix(bead_size_small/2):center_int(1)+floor(bead_size_small/2),...
        center_int(2)-fix(bead_size_small/2):center_int(2)+floor(bead_size_small/2),...
        center_int(3)-fix(bead_size_small/2):center_int(3)+floor(bead_size_small/2)) + ...
        small_bead{bead_num}(2:end,2:end,2:end);
    
end

I_bead = composed_img(bead_size_small:(end-bead_size_small),bead_size_small:(end-bead_size_small),...
    bead_size_small:(end-bead_size_small));

%
% %display options
% figure
% imshow3D(composed_img,[])
%
% figure
% imshow3D(small_bead,[])
%
% figure
% sizeOut = size(small_bead);
% hFigRotated = figure;
% hAxRotated  = axes;
% slice(double(small_bead),sizeOut(2)/2,sizeOut(1)/2,sizeOut(3)/2);
% grid on, shading interp, colormap gray, axis image
%
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function I = seedDefBeads(cur_bead_img,sizeI,x1)

% % %set up variables
% % sizeB = size(cur_bead_img);
% % I = zeros(sizeI(1)+2*sizeB(1),sizeI(2)+2*sizeB(2),sizeI(3)+2*sizeB(3));
% % [px_x,px_y,px_z] = ndgrid(1:sizeB(1),1:sizeB(2),1:sizeB(3));
% %
% % %loop through to add each bead
% % for ii = 1:length(x1)
% %
% %     %find the center loc
% %     center = x1(ii,:) + sizeB;
% %     center_int = floor(center); %interger loc
% %
% %     center_offset = center - center_int; %decimal
% %
% %     %interpolate the bead image onto the pixel grid
% %     bead_shifted = interpn(px_x+center_offset(1),px_y+center_offset(2),px_z+center_offset(3),...
% %         cur_bead_img,px_x,px_y,px_z,'spline');
% %
% %     %add the bead into the image
% %     I((center_int(1)-ceil(sizeB(1)/2)):(center_int(1)+floor(sizeB(1)/2)),...
% %         (center_int(2)-ceil(sizeB(2)/2)):(center_int(2)+floor(sizeB(2)/2)),...
% %         (center_int(3)-ceil(sizeB(3)/2)):(center_int(3)+floor(sizeB(3)/2))) = ...
% %         bead_shifted;
% %
% % end

%
%bead image

%      x=0;
%     n = 0;
%     small_bead{bead_num} = zeros(final_res,final_res,final_res);
%     sub_vol = zeros(size(kernel,1),size(kernel,2),size(kernel,3),final_res^3);
%     for ii = 1:size(kernel,1)-1:size(cur_bead_img,1)-1
%         x=x+1;
%         y=0;
%         for jj = 1:size(kernel,2)-1:size(cur_bead_img,2)-1
%             y=y+1;
%             z=0;
%             for kk = 1:size(kernel,3)-1:size(cur_bead_img,3)-1
%                 z=z+1;
%                 n = n+1;
%                 sub_vol(:,:,:,n) = cur_bead_img(ii:ii+size(kernel,1)-1,...
%                     jj:jj+size(kernel,2)-1,kk:kk+size(kernel,3)-1);
%
%             end
%         end
%     end
%     toc(t1)
%     t2 = tic;
%     small_bead{bead_num}(:) = squeeze(sum(sum(sum(bsxfun(@times,sub_vol,kernel),1),2),3));
%     toc(t2)

