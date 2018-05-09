%%
load('vol00.mat')
I0 = vol{1};
load('vol01.mat')
I1 = vol{1};

sizeI = [512,512,101];
I0 = I0(1:sizeI(1),1:sizeI(2),1:sizeI(3));
I1 = I1(1:sizeI(1),1:sizeI(2),1:sizeI(3));

vol{1} = I0;
save('test00.mat','vol')
vol{1} = I1;
save('test01.mat','vol')