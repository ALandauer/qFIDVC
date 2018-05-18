function I_bead = seedGenDefBeadSeries_3d(cur_image,bead_imgs,sizeI,x1,bead_size_small)
%This function places beads defined by the bead_imgs, into
%and image of size sizeI, at locations specified by x1, the bead is
%downsampled to the size given in bead_size.
%
% A Landauer, Franck Lab, 5/10/18

%Params
margin_offset = bead_size_small;
composed_img = zeros(sizeI(1)+2*margin_offset,sizeI(2)+2*margin_offset,sizeI(3)+2*margin_offset);
composed_img(bead_size_small:(end-bead_size_small),...
    bead_size_small:(end-bead_size_small),bead_size_small:(end-bead_size_small)) = cur_image;
sizeB = size(bead_imgs{1}); %should be n x n, n=bead_size_small

[px_x,px_y,px_z] = ndgrid(1:sizeB(1),1:sizeB(2),1:sizeB(3));

cur_bead_shift = cell(1,size(x1,1));

%interpolate each bead to its subpixel location in parallel
parfor bead_num = 1:size(x1,1)
    %find the center loc
    center = x1(bead_num,:) + margin_offset;
    %     center = [10.5,10.5,10.5];
    
    center_int = floor(center); %interger center loc
    center_offset = (center - center_int).*(sizeB(1)./bead_size_small); %fractional part
    
    %interpolate the bead image onto the pixel grid
    cur_bead_shift{bead_num} = interpn(px_x+center_offset(1),px_y+center_offset(2),...
        px_z+center_offset(3),bead_imgs{bead_num},px_x,px_y,px_z,'spline');
end

for bead_num = 1:size(x1,1)
    %find the center loc
    center = x1(bead_num,:) + bead_size_small;
    
    center_int = floor(center); %interger loc
    %only keep beads that stay within the image
    if (center_int(1)+bead_size_small) > sizeI(1)
        I_bead = cur_image;
    elseif (center_int(1)-bead_size_small) < 1
        I_bead = cur_image;
    elseif (center_int(2)+bead_size_small) > sizeI(2)
        I_bead = cur_image;
    elseif (center_int(1)-bead_size_small) <1
        I_bead = cur_image;
    elseif (center_int(3)+bead_size_small) > sizeI(3)
        I_bead = cur_image;
    elseif (center_int(3)-bead_size_small) < 1
        I_bead = cur_image;
    else
        %add the bead into the image
        composed_img(center_int(1)-fix(bead_size_small/2):center_int(1)+floor(bead_size_small/2),...
            center_int(2)-fix(bead_size_small/2):center_int(2)+floor(bead_size_small/2),...
            center_int(3)-fix(bead_size_small/2):center_int(3)+floor(bead_size_small/2)) ...
            = ...
            composed_img(center_int(1)-fix(bead_size_small/2):center_int(1)+floor(bead_size_small/2),...
            center_int(2)-fix(bead_size_small/2):center_int(2)+floor(bead_size_small/2),...
            center_int(3)-fix(bead_size_small/2):center_int(3)+floor(bead_size_small/2))...
            + ...
            cur_bead_shift{bead_num};
    end
end

%set output
I_bead = composed_img(bead_size_small:(end-bead_size_small),...
    bead_size_small:(end-bead_size_small),bead_size_small:(end-bead_size_small));
end



