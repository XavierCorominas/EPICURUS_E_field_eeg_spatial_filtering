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
mainpath = filePath(1:(end-30));
cd(mainpath);



%%

%load data:
EEG = pop_loadset('filename','all_sp_cer_real.set','filepath','/Users/tuomasmutanen/Documents/Paris_Trip/cRETMS/DTI/');


%% Load the lead fields:

%load lead fields and cortex surf:
load("LFM_Edir_ave.mat")
%LFM_Edir = LFM_Edir_ave;

load("LFM_Exyz_ave.mat")
%LFM_Exyz = LFM_Exyz_ave;

load("cortex_surf.mat")
cortex.p = meshcentroid(cortex_surf.p,cortex_surf.e);

%% Evalueate the head models: 1. Plot the overall sensitivity to all the sources. Does any of the sources stand out?

close all;

% LFM_dir

source_sentivity = sum(LFM_Edir_ave.^2,1);

figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),source_sentivity,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
axis off
title('Global LF-power of Constraint dipole model')

% LFM_Exyz

source_sentivity = sum(LFM_Exyz_ave(:,1:3:end).^2 + LFM_Exyz_ave(:,2:3:end).^2 + LFM_Exyz_ave(:,3:3:end).^2,1);

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),source_sentivity,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar
axis off
title('Global LF-power of Fre dipole model')

%% Evalueate the head models: 2. Plot the peak-localization error of point-spread-function 
% (How accurately would ideal individual sources be localized?)

close all;

LFM = LFM_Edir_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;
R_MNE = MNE_operator*LFM;

[~,inds] = max(abs(R_MNE),[],1);
localizations = cortex.p(inds,:);
PLE = sqrt(sum((localizations - cortex.p).^2,2));

figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
axis off
title('PLE of PSF of Constraint dipole model')


LFM = LFM_Exyz_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;

[~,inds] = max(abs(R_MNE),[],1);
localizations = cortex.p(inds,:);
PLE = sqrt(sum((localizations - cortex.p).^2,2));

for i = 1:length(cortex.p)
    MNE = MNE_operator*LFM_Edir_ave(:,i);
    
    [~,maxI] = max(MNE(1:3:end).^2 + MNE(2:3:end).^2 +MNE(3:3:end).^2 );
    localizations(i,:) = cortex.p(maxI,:);
end

PLE = sqrt(sum((localizations - cortex.p).^2,2));

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
%caxis([0 50])
axis off
title('PLE of PSF of Free dipole model')


%% Plot spatial deviation of MNE CTF

close all;

LFM = LFM_Edir_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;
R_MNE = MNE_operator*LFM;

S_D = zeros(length(cortex.p),1);
for i = 1:length(cortex.p)
    [~,maxInd] = max(abs(R_MNE(i,:)));
    distances = pdist2(cortex.p(maxInd,:), cortex.p)';
    S_D(i) = abs(R_MNE(i,:))*distances/sum(abs(R_MNE(i,:)));
end


figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),S_D,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
axis off
title('SD of CTF of Constraint dipole model')


LFM = LFM_Exyz_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;

S_D = zeros(length(cortex.p),1);

for i = 1:length(cortex.p)
    R_MNE_i = MNE_operator(i,:)*LFM_Edir_ave;
    [~,maxInd] = max(abs(R_MNE_i));
    distances = pdist2(cortex.p(i,:), cortex.p)';
    S_D(i) = abs(R_MNE_i)*distances/sum(abs(R_MNE_i));
end

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
%caxis([0 50])
axis off
title('SD of CTF of Free dipole model')

%% Strip brainstem:

inds = []
area1 = 2000;
area2 = 500;
[inds] = select_sources_from_surface(cortex_surf,area1,area2,inds);

%% Correct forward model:

coords_included = setdiff(1:length(cortex.p),inds);

LFM_Edir_ave = LFM_Edir_ave(:,coords_included);
cortex.p = cortex.p(coords_included,:);
cortex_surf.e = cortex_surf.e(coords_included,:);

k = 1;
for i = coords_included
    tmp(:,k:(k+2)) = LFM_Exyz_ave(:,((i - 1)*3 + 1):((i - 1)*3 + 3));
    k = k+3;
end

LFM_Exyz_ave = tmp; clear tmp;

%% Evalueate the head models: 1. Plot the overall sensitivity to all the sources. Does any of the sources stand out?

close all;

% LFM_dir

source_sentivity = sum(LFM_Edir_ave.^2,1);

figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),source_sentivity,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
axis off
title('Global LF-power of Constraint dipole model')

% LFM_Exyz

source_sentivity = sum(LFM_Exyz_ave(:,1:3:end).^2 + LFM_Exyz_ave(:,2:3:end).^2 + LFM_Exyz_ave(:,3:3:end).^2,1);

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),source_sentivity,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar
axis off
title('Global LF-power of Fre dipole model')

%% Evalueate the head models: 2. Plot the peak-localization error of point-spread-function 
% (How accurately would ideal individual sources be localized?)

close all;

LFM = LFM_Edir_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;
R_MNE = MNE_operator*LFM;

[~,inds] = max(abs(R_MNE),[],1);
localizations = cortex.p(inds,:);
PLE = sqrt(sum((localizations - cortex.p).^2,2));

figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
axis off
title('PLE of PSF of Constraint dipole model')


LFM = LFM_Exyz_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;

[~,inds] = max(abs(R_MNE),[],1);
localizations = cortex.p(inds,:);
PLE = sqrt(sum((localizations - cortex.p).^2,2));

for i = 1:length(cortex.p)
    MNE = MNE_operator*LFM_Edir_ave(:,i);
    
    [~,maxI] = max(MNE(1:3:end).^2 + MNE(2:3:end).^2 +MNE(3:3:end).^2 );
    localizations(i,:) = cortex.p(maxI,:);
end

PLE = sqrt(sum((localizations - cortex.p).^2,2));

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
%caxis([0 50])
axis off
title('PLE of PSF of Free dipole model')


%% Plot spatial deviation of MNE CTF

close all;

LFM = LFM_Edir_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;
R_MNE = MNE_operator*LFM;

S_D = zeros(length(cortex.p),1);
for i = 1:length(cortex.p)
    [~,maxInd] = max(abs(R_MNE(i,:)));
    distances = pdist2(cortex.p(maxInd,:), cortex.p)';
    S_D(i) = abs(R_MNE(i,:))*distances/sum(abs(R_MNE(i,:)));
end


figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),S_D,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
axis off
title('SD of CTF of Constraint dipole model')


LFM = LFM_Exyz_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;

S_D = zeros(length(cortex.p),1);

for i = 1:length(cortex.p)
    R_MNE_i = MNE_operator(i,:)*LFM_Edir_ave;
    [~,maxInd] = max(abs(R_MNE_i));
    distances = pdist2(cortex.p(i,:), cortex.p)';
    S_D(i) = abs(R_MNE_i)*distances/sum(abs(R_MNE_i));
end

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
%caxis([0 50])
axis off
title('SD of CTF of Free dipole model')

%% Correct individual bad sources in the cortex:

 [LFM_Edir_ave] = correct_individual_sources(cortex_surf, LFM_Edir_ave);

 %%

 [LFM_Exyz_ave] = correct_individual_sources(cortex_surf, LFM_Exyz_ave);


%% Evalueate the head models: 1. Plot the overall sensitivity to all the sources. Does any of the sources stand out?

close all;

% LFM_dir

source_sentivity = sum(LFM_Edir_ave.^2,1);

figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),source_sentivity,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
axis off
title('Global LF-power of Constraint dipole model')

% LFM_Exyz

source_sentivity = sum(LFM_Exyz_ave(:,1:3:end).^2 + LFM_Exyz_ave(:,2:3:end).^2 + LFM_Exyz_ave(:,3:3:end).^2,1);

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),source_sentivity,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar
axis off
title('Global LF-power of Fre dipole model')

%% Evalueate the head models: 2. Plot the peak-localization error of point-spread-function 
% (How accurately would ideal individual sources be localized?)

close all;

LFM = LFM_Edir_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;
R_MNE = MNE_operator*LFM;

[~,inds] = max(abs(R_MNE),[],1);
localizations = cortex.p(inds,:);
PLE = sqrt(sum((localizations - cortex.p).^2,2));

figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
axis off
title('PLE of PSF of Constraint dipole model')


LFM = LFM_Exyz_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;

[~,inds] = max(abs(R_MNE),[],1);
localizations = cortex.p(inds,:);
PLE = sqrt(sum((localizations - cortex.p).^2,2));

for i = 1:length(cortex.p)
    MNE = MNE_operator*LFM_Edir_ave(:,i);
    
    [~,maxI] = max(MNE(1:3:end).^2 + MNE(2:3:end).^2 +MNE(3:3:end).^2 );
    localizations(i,:) = cortex.p(maxI,:);
end

PLE = sqrt(sum((localizations - cortex.p).^2,2));

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
%caxis([0 50])
axis off
title('PLE of PSF of Free dipole model')


%% Plot spatial deviation of MNE CTF

close all;

LFM = LFM_Edir_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;
R_MNE = MNE_operator*LFM;

S_D = zeros(length(cortex.p),1);
for i = 1:length(cortex.p)
    [~,maxInd] = max(abs(R_MNE(i,:)));
    distances = pdist2(cortex.p(maxInd,:), cortex.p)';
    S_D(i) = abs(R_MNE(i,:))*distances/sum(abs(R_MNE(i,:)));
end


figure('Position',[100, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),S_D,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
axis off
title('SD of CTF of Constraint dipole model')


LFM = LFM_Exyz_ave;
LL = LFM*LFM';
M = 30;
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';

MNE_operator = LFM'*LL_inv;

S_D = zeros(length(cortex.p),1);

for i = 1:length(cortex.p)
    R_MNE_i = MNE_operator(i,:)*LFM_Edir_ave;
    [~,maxInd] = max(abs(R_MNE_i));
    distances = pdist2(cortex.p(i,:), cortex.p)';
    S_D(i) = abs(R_MNE_i)*distances/sum(abs(R_MNE_i));
end

figure('Position',[700, 100, 600, 400]);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),PLE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

lightangle(180,-60)

h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
colorbar;
material shiny; 
%caxis([0 50])
axis off
title('SD of CTF of Free dipole model')

