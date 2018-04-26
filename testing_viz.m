%%
[x,y,z] = ndgrid(linspace(1,65,512),linspace(1,65,512),linspace(1,25,192));
[x1,y1,z1] = ndgrid(1:65,1:65,1:25);

qf1 = interpn(x1,y1,z1,cc{1, 1}{1, 4}.alpha_mask.qf1,x,y,z,'nearest');
qf2 = interpn(x1,y1,z1,cc{1, 1}{1, 4}.alpha_mask.qf2,x,y,z,'nearest');

figure,imagesc3D(qf1)
figure,imagesc3D(qf2)