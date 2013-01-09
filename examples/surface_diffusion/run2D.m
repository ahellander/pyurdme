% Run2D. Runs the mean hitting time problem on the surface of a sphere. 
%
% A. Hellander, 05/10/2012. 
% 

clear all; close all;
load surfacemeshes.mat;

N = 10;
nmsh = 4;
mean_hitting = zeros(N,nmsh);
comptime     = zeros(N,nmsh);

fem02.xmesh = meshextend(fem02);
for i=1:N
    tic;
    fem = urdme(fem02,@hitting2D,{'Propensities','hitting2D'});
    mean_hitting(i,1) = mht(fem);
    comptime(i,1)=toc;
    disp(['Dataset 1, realization ' num2str(i) 'of ' num2str(N)  ' complete.']);
end

fem01.xmesh = meshextend(fem01);
for i=1:N
    tic;
    fem = urdme(fem01,@hitting2D,{'Propensities','hitting2D'});
    mean_hitting(i,2) = mht(fem);
    comptime(i,2)=toc;
    disp(['Dataset 2, realization ' num2str(i) 'of ' num2str(N)  ' complete.']);

end

fem008.xmesh = meshextend(fem008);
for i=1:N
    tic;
    fem = urdme(fem008,@hitting2D,{'Propensities','hitting2D'});
    mean_hitting(i,3) = mht(fem);
    comptime(i,3)=toc;
    disp(['Dataset 3, realization ' num2str(i) 'of ' num2str(N)  ' complete.']);

end

fem006.xmesh = meshextend(fem006);
for i=1:N
    tic;
    fem = urdme(fem006,@hitting2D,{'Propensities','hitting2D'});
    mean_hitting(i,4) = mht(fem);
    comptime(i,4) = toc;
    disp(['Dataset 4, realization ' num2str(i) 'of ' num2str(N)  ' complete.']);
end

voxels(1) = size(fem02.mesh.p,2);
voxels(2) = size(fem01.mesh.p,2);
voxels(3) = size(fem008.mesh.p,2);
voxels(4) = size(fem006.mesh.p,2);

% Analytic solution
rtarget = 0.1;
tr = (2*log(2/(1-cos(rtarget)))/(1+cos(rtarget))-1);

% Numeric solution
meanht = mean(mean_hitting,1);
errorht = 2*std(mean_hitting,1)/sqrt(N);

figure(1);
hold;
plot(voxels,tr*ones(1,nmsh),'r--');
errorbar(voxels,meanht,errorht);
xlabel('# Voxels');
ylabel('Mean hitting time (s)');
legend('Analytic solution','URDME');

figure(2);
hold;
errorbar(voxels,mean(comptime,1),2*std(comptime,1)/sqrt(N));
xlabel('# Voxels');
ylabel('Time');



