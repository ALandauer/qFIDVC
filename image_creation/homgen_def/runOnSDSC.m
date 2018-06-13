processors=32;
clusterHost='comet.sdsc.edu';
ppn=16;
username='alandaue';
% account='sds154';
account='brn110';
queue='compute';
time='00:10:00';
% DataLocation='C:\Users\FRUVR\Google Drive\ParticleDeformDIC\mfiles\image_synth';
DataLocation='C:\Users\FRUVR\Documents\syncfiles\ParticleDeformDIC\mfiles\image_synth';
RemoteDataLocation='/oasis/scratch/comet/alandaue/temp_project/parallel_testing';
% RemoteDataLocation='/home/alandaue/parallel_testing';
keyfile='C:\Users\FRUVR\Documents\syncfiles\Research\ssh\ssh_key_FRUVR';
matlabRoot='/opt/matlab/2016b';
cluster = getCluster(username,account,clusterHost,ppn,queue,time,DataLocation,...
RemoteDataLocation,keyfile,matlabRoot);
j = createCommunicatingJob(cluster);
j.AttachedFiles={'testparfor2.m'};
set(j,'NumWorkersRange',[1 processors]);
set(j,'Name','Test');
t = createTask(j,@testparfor2,1,{processors});
submit(j);
% wait(j);
% pause(30);
o=j.fetchOutputs;
o{:}
   