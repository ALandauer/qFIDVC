function small_bead = small_bead_synth_genDef(def_matrix,bead_size_small)
%This function downsamples a high-res bead image to a given size, using
%subset volume-intensity summation to simulate a sensor.
%
%A Landauer, 3/30/17 latest

%generate the large bead based on the given deformation matrix
cur_bead_img = generate_single_speckle(def_matrix);

%set up bead size and summation kernel
sizeB = size(cur_bead_img);
%     kernel = ones(floor(sizeB(1)/bead_size_small)+1,floor(sizeB(2)/bead_size_small)+1,...
%         floor(sizeB(3)/bead_size_small)+1);

sum_width(1) = floor(sizeB(1)/bead_size_small);
offset(1) = ceil(rem(sizeB(1),bead_size_small)/2);
sum_width(2) = floor(sizeB(2)/bead_size_small);
offset(2) = ceil(rem(sizeB(2),bead_size_small)/2);

if length(size(cur_bead_img)) > 2
    sum_width(3) = floor(sizeB(3)/bead_size_small);
    offset(3) = ceil(rem(sizeB(3),bead_size_small)/2);
end

%do the summation
x=0;
if length(size(cur_bead_img)) > 2
    small_bead = zeros(bead_size_small,bead_size_small,bead_size_small);
else
    small_bead = zeros(bead_size_small,bead_size_small);
end

for ii = offset(1):sum_width(1):size(cur_bead_img,1)-offset(1)-1
    x=x+1;
    y=0;
    for jj = offset(2):sum_width(2):size(cur_bead_img,2)-offset(1)-1
        y=y+1;
        
        if length(size(cur_bead_img)) > 2
            z=0;
            for kk = offset(3):sum_width(3):size(cur_bead_img,3)-offset(1)-1
                z=z+1;
                
                small_bead(x,y,z) = ...
                    sum(sum(sum(cur_bead_img(ii:ii+sum_width(1)-1,...
                    jj:jj+sum_width(2)-1,kk:kk+sum_width(3)-1))));
                
            end %kk
        else
            small_bead(x,y) = ...
                sum(sum(cur_bead_img(ii:ii+sum_width(1)-1,...
                jj:jj+sum_width(2)-1)));
        end
        
    end %jj
end %ii

small_bead = 255*small_bead/max(small_bead(:));

end %function

% save('small_bead_img_series','small_bead')