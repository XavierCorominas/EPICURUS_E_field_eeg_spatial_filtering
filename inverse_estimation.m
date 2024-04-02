clear;
clc;
close all;

% Set up the paths to needed tools:

% Simnibs for segmenting head:
addpath('/Users/tuomasmutanen/Applications/SimNIBS-4.0/simnibs_env/lib/python3.9/site-packages/simnibs/matlab_tools/')


% Add path to ISO2Mesh
addpath('/Users/tuomasmutanen/Documents/Paris_Trip/ExampleData_Finland/Analysis/ISO2MESH/')

% Add path to EEGLAB
eeglab_path = '/Users/tuomasmutanen/Documents/Paris_Trip/ExampleData_Finland/Analysis/eeglab2022.1';
addpath(eeglab_path);
eeglab;

% Add path to helsinki bem framework
addpath('/Users/tuomasmutanen/Documents/Paris_Trip/ExampleData_Finland/Analysis/hbf_lc_p-master/');
hbf_SetPaths

% Get the path of the currently running script
script_path = fileparts(mfilename('fullpath'));

% Get path to the DTI tracks
addpath('/Users/tuomasmutanen/Documents/Paris_Trip/cRETMS/DTI/matlab')

filePath = matlab.desktop.editor.getActiveFilename;
mainpath = filePath(1:(end-21));
cd(mainpath);

%%

%load data:
EEG = pop_loadset('filename','all_sp_cer_real.set','filepath','/Users/tuomasmutanen/Documents/Paris_Trip/cRETMS/DTI/');


%%

%load lead fields and cortex surf:
load("LFM_Edir_ave.mat")
LFM_Edir = LFM_Edir_ave;

load("LFM_Exyz_ave.mat")
LFM_Exyz = LFM_Exyz_ave;

load("cortex_surf.mat")

cortex.p = meshcentroid(cortex_surf.p,cortex_surf.e);



%% Load tha DTI tracks:

tracks = read_mrtrix_tracks('tck_right_cerebellum_and_left_thalamus_ends_only_nospinal.tck')
DTI_data = dlmread('tck_right_cerebellum_and_left_thalamus_ends_only_nospinal.txt');

figure;
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),'EdgeAlpha',0,'FaceAlpha',1);
view(0,75);
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
material shiny;
axis off;
hold on;
termination_coords = zeros(2*length(DTI_data),3);

for i = 1:length(DTI_data)
    termination_coords(i,:) = tracks.data{i}(1,:);
    termination_coords(length(DTI_data)+ i,:) = tracks.data{i}(end,:);
end

plot3(termination_coords(1:2:end,1),termination_coords(1:2:end,2),termination_coords(1:2:end,3),'.b')
plot3(termination_coords(2:2:end,1),termination_coords(2:2:end,2),termination_coords(2:2:end,3),'.r')

source_weighs = zeros(length(cortex.p),1);
% Let's weigh each source:
tic 
disp('Computing source weighs based on DTI')
for i = 1:length(cortex.p)
    all_distances = pdist2(cortex.p(i,:),termination_coords);
    distance_weights = 1./all_distances/sum(1./all_distances);
    source_weighs(i) = distance_weights*[DTI_data';DTI_data'];
end

source_weighs = (source_weighs - min(source_weighs));
source_weighs = source_weighs/max(source_weighs);
toc

%%
figure;
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),source_weighs,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75);
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
material shiny;
axis off;
hold on;
colorbar;
cmp = colormap('cool')
colormap(flipud(cmp))
%%

%%

ROI = source_weighs;
ROI(ROI < 0.9) = 0;

figure;
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),ROI,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75);
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
material shiny;
axis off;
hold on;
colorbar;
cmp = colormap('cool')
colormap(flipud(cmp))




%%
%close all;

figure;
plot(EEG.times,mean(EEG.data,3)','k');
xlim([-50 250])
for i = [20 40 60 80 100 200]

[~, respInd] = min(abs(EEG.times - i))
EEG_signal =mean(EEG.data(:,respInd,:),3);

LFM = LFM_Edir;%*diag(source_weighs);
LL = LFM*LFM' ;
trL = trace(LL);

lambda = 0.001;

M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

%tmp = LFM'*((LL + lambda*trace(LL)*eye(63))\EEG_signal);

tmp = LFM'*(LL_inv*EEG_signal);

MNE = abs(tmp);

figure;
subplot(1,2,1);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)
axis off;

lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
title([num2str(i), 'ms'])
subplot(1,2,2);
topoplot(EEG_signal,EEG.chanlocs)

end





%%
close all;

figure;
plot(EEG.times,mean(EEG.data,3)','k');
xlim([-50 250])
for i = [20 40 60 80 100 200]

[~, respInd] = min(abs(EEG.times - i))
EEG_signal =mean(EEG.data(:,respInd,:),3);


LL = LFM_Exyz*LFM_Exyz' ;
%trL = trace(LL);

%lambda = 1;
%tmp = LFM_Exyz'*((LL + lambda*trace(LL)*eye(63))\EEG_signal);
%MNE = tmp(1:3:end).^2 + tmp(2:3:end).^2 +tmp(3:3:end).^2;

M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

%tmp = LFM'*((LL + lambda*trace(LL)*eye(63))\EEG_signal);
%tmp = pinv(LFM_Exyz)* EEG_signal;
tmp  = LFM_Exyz'*(LL_inv*EEG_signal);
MNE = tmp(1:3:end).^2 + tmp(2:3:end).^2 +tmp(3:3:end).^2;

figure;
subplot(1,2,1);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)
axis off;

lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
title([num2str(i), 'ms'])
subplot(1,2,2);
topoplot(EEG_signal,EEG.chanlocs)

end



