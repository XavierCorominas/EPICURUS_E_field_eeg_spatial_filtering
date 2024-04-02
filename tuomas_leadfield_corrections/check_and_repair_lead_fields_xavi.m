clear;
clc;
close all;

% Set up the paths to needed tools:

% Simnibs for segmenting head:
addpath('C:\Users\XAVIER\SimNIBS-4.0\simnibs_env\Lib\site-packages\simnibs\matlab_tools\')

% Add path to dcm2nii tools
%addpath('/Users/tpmutane/Documents/Double-coil/Dicom_2_nifti')

% Add path to ISO2Mesh
addpath('C:\Users\XAVIER\iso2mesh-1.9.6\')

% Add path to EEGLAB
eeglab_path = 'D:\SOFTWARE\eeglab2022.0';
addpath(eeglab_path);
eeglab;

% Add path to helsinki bem framework
addpath('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions\hbf_lc_p-master');
hbf_SetPaths

%
external_functions_path = ('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions');
addpath(external_functions_path);

% 
leadfield_corrections = ('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\Preprocessig_xavi_tuomas\tuomas_leadfield_corrections')
addpath(leadfield_corrections);

% Get the path of the currently running script
script_path = fileparts(mfilename('fullpath'));

% Get path to the DTI tracks
%addpath('/Users/tuomasmutanen/Documents/Paris_Trip/cRETMS/DTI/matlab')

filePath = matlab.desktop.editor.getActiveFilename;
mainpath = filePath(1:end-35);
cd(mainpath);



%% Define subejct 

%sub id
subject_id = '003';

% Add paths to leafields, passband and MRI
%addpath(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS'])
path_to_leadfields = (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS']);
addpath(path_to_leadfields);
cd(path_to_leadfields)


% Load EEG data to test the leadfields:
% BLOCK NUMBER
block = '1';   % '1' , '2'
% STIMULATION intensity
stim ='90';  %  10 , 30 , 50 , 70 , 90 , 100
% TARGET
condition = 'real'; %sham real
EEG = pop_loadset('filename',([block,'_',stim,'_M1_',condition,'_Cleaned.set']),'filepath',(['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_',subject_id,'\']));
eeglab redraw


%% Load the lead fields:

%load lead fields and cortex surf: WE HERE ONLY LOAD THE EDIR_AVE AS WE
%WORK WITH THIS NOW IN XAVIS PIPEPLINE (I.E.,WITH CORTICALLY CONSTRAINED DIPOLS)
load("LFM_Edir_ave.mat")
%LFM_Edir = LFM_Edir_ave;

%load("LFM_Exyz_ave.mat")
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


%% Strip brainstem:

inds = []
area1 = 10;
area2 = 5;
[inds] = select_sources_from_surface(cortex_surf,area1,area2,inds);

%% Correct forward model:

coords_included = setdiff(1:length(cortex.p),inds);

LFM_Edir_ave = LFM_Edir_ave(:,coords_included);
cortex.p = cortex.p(coords_included,:);
cortex_surf.e = cortex_surf.e(coords_included,:);

% k = 1;
% for i = coords_included
%     tmp(:,k:(k+2)) = LFM_Exyz_ave(:,((i - 1)*3 + 1):((i - 1)*3 + 3));
%     k = k+3;
% end
% 
% LFM_Exyz_ave = tmp; clear tmp;

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



%% Correct individual bad sources in the cortex:

 [LFM_Edir_ave] = correct_individual_sources(cortex_surf, LFM_Edir_ave);



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

%% Save new leadfield and new cortex_surf 

save_path = (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\']);


save([save_path,'LFM_Edir_ave_corrected.mat'],'LFM_Edir_ave');
save([save_path,'cortex_surf_corrected.mat'],'cortex_surf')
%save([save_path,'t.mat'],'t_corrected');
%save([save_folder,'t.mat'],'t');

close all


