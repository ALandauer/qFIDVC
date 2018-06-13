function I_bead = seedGenDefBead_2d(cur_image,cur_bead_img,sizeI,x1,bead_size_small)
%This function places a bead defined by the cur_bead_img, into
%and image of size sizeI, at the location specified by x1_, the bead is
%downsampled to the size given in bead_size.
%
% A Landauer, Franck Lab, 4/5/17

%Params
margin_offset = bead_size_small;
composed_img = zeros(sizeI(1)+2*margin_offset,sizeI(2)+2*margin_offset);
composed_img(bead_size_small:(end-bead_size_small),bead_size_small:(end-bead_size_small)) = cur_image;
sizeB = size(cur_bead_img); %should be n x n, n=bead_size_small

%only keep beads that stay within the image
if (x1(1))>sizeI(1)
        I_bead = cur_image;
elseif (x1(1))<1 
        I_bead = cur_image;
elseif (x1(2))>sizeI(2)
        I_bead = cur_image;
elseif (x1(2))<1
    I_bead = cur_image;
else
        
    [px_x,px_y] = ndgrid(1:sizeB(1),1:sizeB(2));
    
    cur_bead_shift = cell(1,size(x1,1));
    
    %interpolate each bead to its subpixel location in parallel
    
        %find the center loc
        center = x1 + margin_offset;
        %     center = [10.5,10.5,10.5];
        
        center_int = floor(center); %interger center loc
        center_offset = (center - center_int)*(sizeB(1)/bead_size_small); %fractional part
        
        %interpolate the bead image onto the pixel grid
        cur_bead_shift = interpn(px_x+center_offset(1),px_y+center_offset(2),...
            cur_bead_img,px_x,px_y,'spline');
        
        %find the center loc
        center = x1 + bead_size_small;
        
        center_int = floor(center); %interger loc
        
        %add the bead into the image
        composed_img(center_int(1)-fix(bead_size_small/2):center_int(1)+floor(bead_size_small/2),...
            center_int(2)-fix(bead_size_small/2):center_int(2)+floor(bead_size_small/2)) ...
            = ...
            composed_img(center_int(1)-fix(bead_size_small/2):center_int(1)+floor(bead_size_small/2),...
            center_int(2)-fix(bead_size_small/2):center_int(2)+floor(bead_size_small/2))...
            + ...
            cur_bead_shift;
    
    %set output
    I_bead = composed_img(bead_size_small:(end-bead_size_small),bead_size_small:(end-bead_size_small));
end

end



