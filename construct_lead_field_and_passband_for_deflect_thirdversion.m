

%************ CONSTRUCT HEAD MODEL , LEADFIELD  AND PASSBAND FOR  A INDIVIDUAL SUBJECT ***************

%% The currect script builds the LEADFIELD and the PASSBAND (passband based on the TMS-electric
% field source projection) to be used  with the Deflect framework and spatially filter a EEG dataset. 

% The script will be applyed to every subject of our cohort to obtain the
% subject-individualized leadfields and passband. Then , all TMS-EEG eeg dataset of
% the subject will be fitlered.

% XAVIER COROMINAS 01/04/2024, Paris, France.


% The  script is organized in 3 blocks:

   % - 1 construct MRI-based headmodels and compute leadfields
   % - 2 correct manually leadfields
   % - 3 compute E-field simulations and extract ROI to be further used as a passband for the spatial filters
   % - 4* apply spatial filter to a eeg dataset with the  deflect_xavi_script
 
%% **************** 1rst part of the script: build leadfield for DeFLECT ****************

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


%% Define parameters
% Add MRI data info of the subject
subject_id = '003';
path_to_mri_nifti = ['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\MRI'];

%eeg file info
% BLOCK NUMBER
block = '1';   % '1' , '2'
% STIMULATION intensity
stim ='10';  %  10 , 30 , 50 , 70 , 90 , 100
% TARGET
%target = 'M1';  





%% Create mesh model 
% Names acording to original name of mri ( they will be different if the mri is coming from other project)

subject_id_mri = '025';
sequence_mri = 'S5';

%
filename = ['v_CRETMS_03_',subject_id_mri,'_',sequence_mri,'_T1w_mp2rage_UNI_Images_MPRAGEised_biascorrected.nii'];
filename_2 = ['v_CRETMS_03_',subject_id_mri,'_S2_T1w_mp2rage_INV1.nii'];

% Change to correct path
cd(path_to_mri_nifti)
clc;
[status,cmdout] = system(['cd ',path_to_mri_nifti],'-echo');
pwd

% Run the SIMNIBS charm protocol in the terminal to segment the head tissues. Please read the SIMNIBS website for further information:
[status,cmdout] = system(['C:\Users\XAVIER\SimNIBS-4.0\bin\charm',' ',subject_id,' ',filename,' ',filename_2,  ' --forceqform --forcerun'], '-echo');


%% Crate TDCS leadfild

cd([path_to_mri_nifti,'\m2m_',subject_id,'\'])
tdcs_lf = sim_struct('TDCSLEADFIELD');
tdcs_lf.fnamehead = [subject_id,'.msh'];
tdcs_lf.pathfem = 'tdcs_leadfield\'
%tdcs_lf.field =  'J';

 run_simnibs(tdcs_lf);

%% Load mesh electrodes

%cd([path_to_mri_nifti,'\m2m_',subject_id,'\tdcs_leadfield\'])

%mesh_electrodes = mesh_load_gmsh4([subject_id,'_electrodes_EEG10-10_UI_Jurak_2007.msh']);


cd(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_00XX\MRI\m2m_XX\tdcs_leadfield\'])
mesh_electrodes = mesh_load_gmsh4(['XX_electrodes_EEG10-10_UI_Jurak_2007.msh']);

%% Load and visualize the surfaces generated with SimNIBS

% For the three-layer model we need the pial surface, skull and scalp

cd([path_to_mri_nifti,'\m2m_',subject_id])
mesh = mesh_load_gmsh4([path_to_mri_nifti ,'\m2m_',subject_id,'\',subject_id,'.msh']);


close all;
k = 1;
for i = [1,3,7,5]
surface_inds = find(mesh.triangle_regions == (1000+i));

meshes{k}.e = mesh.triangles(surface_inds,:);
meshes{k}.p = mesh.nodes;

k = k+1;
end

figure;
subplot(1,4,1);
surf_idx = 1; % pial surface
h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,2);
surf_idx = 2; % Inner Skull (CSF boundary)

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,3);

subplot(1,4,3);
surf_idx = 3; %  Outer Skull

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,4);

surf_idx = 4; % scalp

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [254, 227, 212]/254);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';


%% Preorocessing of the head surfaces:
% 
% Check and repair the surfaces with ISO2mesh and HBF routines
%
% Resample to lower resolution for computation feasibility.


for i = 1:1:4

disp(['Preprocessing and testing meshes ',num2str(i),'...'])

% disp('Smoothen the surface')
% conn = meshconn(meshes{i}.e,length(meshes{i}.p));
% meshes{i}.p = smoothsurf(meshes{i}.p,[],conn);

disp('ISO2MESH meshcheckrepair')
%[node,elem] = meshcheckrepair(meshes{i}.p,meshes{i}.e,'dup');
%[node,elem] = meshcheckrepair(meshes{i}.p,meshes{i}.e,'isolated');
%[node,elem] = meshcheckrepair(meshes{i}.p,meshes{i}.e,'deep');
[node,elem] = meshcheckrepair(meshes{i}.p,meshes{i}.e,'meshfix');


meshes{i}.e = elem; meshes{i}.p = node;

%opt.gridsize = 1;
%[node,elem] = remeshsurf(meshes{i}.p, meshes{i}.e, opt);
[node,elem] = meshresample(meshes{i}.p, meshes{i}.e, 0.08);
meshes{i}.e = elem; meshes{i}.p = node;
% 

% disp('ISO2MESH removedupelem')
% meshes{i}.e = removedupelem(meshes{i}.e);
% 
% disp('ISO2MESH surfreorient')
% 
% [node,elem] = surfreorient(meshes{i}.p,meshes{i}.e);
% meshes{i}.e = elem; meshes{i}.p = node;

disp('HBF triangle test')
[meshes{i}] = hbf_CorrectTriangleOrientation(meshes{i});

end

figure;
subplot(1,4,1);
surf_idx = 1; % pial surface
h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,2);
surf_idx = 2; % Inner Skull

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,3);

subplot(1,4,3);
surf_idx = 3; %  Outer Skull

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,4);

surf_idx = 4; % scalp

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [254, 227, 212]/254);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

%% Load the theoretical electrode locations and fiducials:

% Better, if you have the digitized locations! This is just an example for
% a head model.

% Initialize variables.
filename = 'C:\Users\XAVIER\SimNIBS-4.0\simnibs_env\Lib\site-packages\simnibs\resources\ElectrodeCaps_MNI\EEG10-10_UI_Jurak_2007.csv';
delimiter = ',';
% Format for each line of text:
formatSpec = '%s%f%f%f%s%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Create output variable
dataArray([2, 3, 4]) = cellfun(@(x) num2cell(x), dataArray([2, 3, 4]), 'UniformOutput', false);
dataArray([1, 5]) = cellfun(@(x) mat2cell(x, ones(length(x), 1)), dataArray([1, 5]), 'UniformOutput', false);
EEG1010 = [dataArray{1:end-1}];
% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

% Initialize variables.
filename = 'C:\Users\XAVIER\SimNIBS-4.0\simnibs_env\Lib\site-packages\simnibs\resources\ElectrodeCaps_MNI\Fiducials.csv';
delimiter = ',';
% Format for each line of text:
formatSpec = '%s%f%f%f%s%[^\n\r]';
% Open the text file.
fileID = fopen(filename,'r');
% Read columns of data according to the format.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string',  'ReturnOnError', false);
% Close the text file.
fclose(fileID);
% Create output variable
dataArray([2, 3, 4]) = cellfun(@(x) num2cell(x), dataArray([2, 3, 4]), 'UniformOutput', false);
dataArray([1, 5]) = cellfun(@(x) mat2cell(x, ones(length(x), 1)), dataArray([1, 5]), 'UniformOutput', false);
fiducials = [dataArray{1:end-1}];
% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;

fiducials = cell2mat(fiducials(:,2:4));

%%
mesh = mesh_load_gmsh4([path_to_mri_nifti ,'\m2m_',subject_id,'\',subject_id,'.msh']);

% load the corresponding EEG data of the subject to find the channel locations:
eeglab;
%EEG = pop_loadset('filename','example_data_S02.set','filepath','/Users/tuomasmutanen/Documents/Paris_Trip/ExampleData_Finland/Analysis/');
EEG = pop_loadset('filename',([block,'_',stim,'_M1_real_Cleaned.set']),'filepath',(['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_',subject_id,]));

%
k = 0;
for i = 1: length(EEG.chanlocs)
  k = 0;
  for j = 1:length(EEG1010)
    if  strcmpi(EEG.chanlocs(i).labels,EEG1010{j,5})
        k = k+1;
        EEG_locs(i,:) = [EEG1010{j,2},EEG1010{j,3},EEG1010{j,4}];
    end
  end
  if ~k
      EEG.chanlocs(i).labels;
  end
end

%% Add the reference

for j = 1:length(EEG1010)
    if  strcmpi('Fz',EEG1010{j,5})
        k = k+1;
        EEG_locs = [EEG_locs;[EEG1010{j,2},EEG1010{j,3},EEG1010{j,4}]];
    end
end


%%
% The channel locations are defined in MNI coordinates. Transform to
% subject_coordinates



EEG_locs_subj =  mni2subject_coords(EEG_locs,(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\MRI\m2m_',subject_id]));

fiducials =  mni2subject_coords(fiducials,(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\MRI\m2m_',subject_id]));

% Make sure the electrodes are projected to scalp:
e_loc_proj = hbf_ProjectElectrodesToScalp(EEG_locs_subj, meshes);


%% Let's check that channel locations are sensible with respect to the surfaces:

figure;
subplot(1,4,1);
surf_idx = 1; % pial surface
h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,2);
surf_idx = 2; % Skull

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

subplot(1,4,3);
surf_idx = 3; % Skull

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';


subplot(1,4,4);

surf_idx = 4; % scalp

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [254, 227, 212]/254);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';
hold on;
plot3(EEG_locs(:,1),EEG_locs(:,2),EEG_locs(:,3),'.k','MarkerSize',20);
hold on;
plot3(fiducials(:,1), fiducials(:,2),fiducials(:,3),'.r','MarkerSize',20);
hold off;


figure;

surf_idx = 1; % pial surface
h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';
hold on;

surf_idx = 2; % Skull

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',0.5, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';


surf_idx = 3; % Skull

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',0.3, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';


surf_idx = 4; % scalp

h = trisurf(meshes{surf_idx}.e, meshes{surf_idx}.p(:,1),  meshes{surf_idx}.p(:,2),...
    meshes{surf_idx}.p(:,3),'EdgeAlpha',0,'FaceAlpha',0.1, 'FaceColor', [254, 227, 212]/254);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';
hold on;
plot3(e_loc_proj.p_proj(:,1),e_loc_proj.p_proj(:,2),e_loc_proj.p_proj(:,3),'.k','MarkerSize',20);
hold on;
plot3(fiducials(:,1), fiducials(:,2),fiducials(:,3),'.r','MarkerSize',20);
hold off;



%% Compute the leadfields
cortex_surf = meshes{1};

bmeshes = meshes(1,2:4);
cortex.p = meshcentroid(cortex_surf.p,cortex_surf.e); % eliminar aixo i ficar centres a default mesh?
cortex.nn = surfacenorm(cortex_surf.p,cortex_surf.e);


contrast=100; % Omitting CSF has been taken into account. This means that we assume in reality 1/50 contrast between scalp and skull

%for contrast = 10:10:100

ci = [0.33 0.33/contrast 0.33];

co = [0.33/contrast 0.33 0];

% make BEM model
D = hbf_BEMOperatorsPhi_LC(bmeshes);


Tphi_full = hbf_TM_Phi_LC_ISA2(D,ci,co,1);
Tphi_elecs = hbf_InterpolateTfullToElectrodes(Tphi_full,bmeshes,e_loc_proj);


% make lead field matrices
LFM_Edir = hbf_LFM_Phi_LC(bmeshes,Tphi_elecs,cortex.p,cortex.nn);
LFM_Exyz = hbf_LFM_Phi_LC(bmeshes,Tphi_elecs,cortex.p);

% to original reference:

N_chan = size(LFM_Exyz,1) - 1; % Number of channels excluding the reference

LFM_Edir = LFM_Edir(1:N_chan,:) - repmat(LFM_Edir((N_chan + 1),:),[N_chan 1]);
LFM_Exyz = LFM_Exyz(1:N_chan,:) - repmat(LFM_Exyz((N_chan + 1),:),[N_chan 1]);


% To ave reference:
LFM_Edir_ave = LFM_Edir - repmat(mean(LFM_Edir ,1),[size(LFM_Edir,1) 1]);
LFM_Exyz_ave = LFM_Exyz - repmat(mean(LFM_Exyz ,1),[size(LFM_Exyz,1) 1]);



%%
close all;


random_source = randi(length(cortex.p));

EEG_signal = LFM_Edir(:,random_source);

tmp = pinv(LFM_Exyz)*EEG_signal;
MNE = tmp(1:3:end).^2 + tmp(2:3:end).^2 +tmp(3:3:end).^2;

figure;
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)

hold on;
plot3(cortex.p(random_source ,1),cortex.p(random_source ,2),cortex.p(random_source ,3),'r.','MarkerSize',50)

lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';

figure;
topoplot(EEG_signal,EEG.chanlocs)

%%

close all;

for i = [20 40 60 80 100 200 400]

EEG_signal =mean(EEG.data(:,(1000+i),:),3);


LL = LFM_Exyz_ave*LFM_Exyz_ave' ;
trL = trace(LL);



lambda = 0.1;
tmp = LFM_Exyz_ave'*((LL + lambda*trace(LL)*eye(63))\EEG_signal);
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

%%

% Save the lead-field for cortically constrained dipoles in the original
% reference:
save_folder = ['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\'];

% save_folder = [save_path,'LEADFIELDS\']
% if ~exist(save_folder)
%       mkdir([save_folder]);
% end

save([save_folder,'LFM_Edir_orig.mat'],'LFM_Edir');

% Save the lead-field for cortically free dipoles in the original
% reference:
save([save_folder,'LFM_Exyz_orig.mat'],'LFM_Exyz');

% Save the lead-field for cortically constrained dipoles in the average
% reference:
save([save_folder,'LFM_Edir_ave.mat'],'LFM_Edir_ave');

% Save the lead-field for cortically free dipoles in the average
% reference:
save([save_folder,'LFM_Exyz_ave.mat'],'LFM_Exyz_ave');

% Save the cortical surface for visualization:
save([save_folder,'cortex_surf.mat'],'cortex_surf');

% Save the boundary meshes and the channel locations for the head model
head_model.e_loc_proj = e_loc_proj;
haed_model.bmeshes = bmeshes;
save([save_folder,'head_model.mat'],'head_model');

%% ***************2nd part of the script: correct leadfields and cortex***********

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




%% ***************3nd part of the script: contruct the  passband***********

close all;
clc;
clearvars -except subject_id subject_id_mri sequence sequence_mri path_to_mri_nifti;


% Set up the paths to needed tools:

% Simnibs for segmenting head:
addpath('C:\Users\XAVIER\SimNIBS-4.0\simnibs_env\Lib\site-packages\simnibs\matlab_tools\')

% Add path to dcm2nii tools
%addpath('/Users/tpmutane/Documents/Double-coil/Dicom_2_nifti')

% Add path to ISO2Mesh
addpath('C:\Users\XAVIER\iso2mesh-1.9.6\')

% Add path to helsinki bem framework
addpath('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions\hbf_lc_p-master\');
hbf_SetPaths

% Add path to DeFleCt
addpath('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions\DeFleCT_8Nov13')
addpath('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions')

% Add eeglab path
eeglab_path = 'D:\SOFTWARE\eeglab2022.0';
addpath(eeglab_path);
eeglab;

%
  %  filename = ['v_CRETMS_03_',subject_id_mri,'_',sequence_mri,'_T1w_mp2rage_UNI_Images_MPRAGEised_biascorrected.nii'];

%
%filename = ['clean_v_CRETMS_01_',subject_id,sequence,'_T1w_mp2rage_UNI_Images.nii'];


%%
% Go to leadfields folder
cd(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\']);

%load cortex surf corrected after leadfield correction
load cortex_surf_corrected.mat;

close all;
cortex_surf.e = cortex_surf.e;
cortex_surf.p = cortex_surf.p;

figure;
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
    cortex_surf.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';




% load and correct WM from original mesh to get mesh structure
cd([path_to_mri_nifti,'\m2m_',subject_id])
mesh = mesh_load_gmsh4([path_to_mri_nifti ,'\m2m_',subject_id,'\',subject_id,'.msh']);

close all;
for i = [1] %  1 for WM
surface_inds = find(mesh.triangle_regions == 1001);
meshes.e = mesh.triangles(surface_inds,:);
meshes.p = mesh.nodes;
end

figure;
h = trisurf(meshes.e, meshes.p(:,1),  meshes.p(:,2),...
    meshes.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);

axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';

%Correct triangles
disp(['Preprocessing and testing meshes ',num2str(i),'...'])
disp('ISO2MESH meshcheckrepair')
[node,elem] = meshcheckrepair(meshes.p,meshes.e,'meshfix');
meshes.e = elem; meshes.p = node;
[node,elem] = meshresample(meshes.p, meshes.e, 0.08);
meshes.e = elem; meshes.p = node;
disp('HBF triangle test')
[meshes] = hbf_CorrectTriangleOrientation(meshes);

figure;
h = trisurf(meshes.e, meshes.p(:,1),  meshes.p(:,2),...
    meshes.p(:,3),'EdgeAlpha',0,'FaceAlpha',1, 'FaceColor', [0.9,0.9,0.9]);
axis off
view(-135,15)
%shading interp
lightangle(55,-35)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.7;
h.DiffuseStrength = 0.4;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'lit';




%%
% generate file mesh with triangles and tetrahedrom
mesh_corrected.node_data = mesh.node_data;
mesh_corrected.element_data = mesh.element_data;
mesh_corrected.nodes = cortex_surf.p;
mesh_corrected.triangles = cortex_surf.e;

triangle_regions = mesh.triangle_regions(mesh.triangle_regions==1001);
triangle_regions = triangle_regions(1:length(mesh_corrected.triangles),1);
mesh_corrected.triangle_regions = triangle_regions;


% Triangulate
tetrahedra = delaunayTriangulation(mesh_corrected.nodes); % connect nodes to report tetrahedron limits
mesh_corrected.tetrahedra = tetrahedra.ConnectivityList;
mesh_corrected.tetrahedron_regions = mesh.tetrahedron_regions(mesh.tetrahedron_regions==1001);
u = repelem(1,length(mesh_corrected.tetrahedra));
u = u';
mesh_corrected.tetrahedron_regions = u;

% 
[newelem, evol]=meshreorient(mesh_corrected.nodes,mesh_corrected.tetrahedra);
mesh_corrected.tetrahedra = newelem;
%

% Save mesh 
close all
cd ([path_to_mri_nifti ,'\m2m_',subject_id,'\']);
mesh_save_gmsh4(mesh_corrected,[subject_id,'_corrected.msh']);


%% TMS Simulation for spatial filter
% Import brainsight data and transform x axys (*be aware that X axis are inverted between
% brainsight and Simnibs space by default*)
importdata (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\NEURONAVIGATOR\Exported_Targets_s1.txt'])

TMS_target = ans;

TMS_target.dataC(:,4) = TMS_target.data(1,1:3);
TMS_target.dataC(:,1) = TMS_target.data(1,4:6);
TMS_target.dataC(1,1) = TMS_target.dataC(1,1)*-1;
TMS_target.dataC(2,1) = TMS_target.dataC(2,1)*-1;
TMS_target.dataC(3,1) = TMS_target.dataC(3,1)*-1;

TMS_target.dataC(:,2) = TMS_target.data(1,7:9);
TMS_target.dataC(:,3) = TMS_target.data(1,10:12);
TMS_target.dataC(4,1:3) = 0;
TMS_target.dataC(4,4) = 1;


%
save_path = ([path_to_mri_nifti ,'\']);

save_folder = [save_path,'simulation\']
if ~exist(save_folder)
      mkdir([save_folder]);
end

% RUN SIMULATION
%Initialize  simulation session
s = sim_struct('SESSION');
s.fnamehead = ([path_to_mri_nifti ,'\m2m_',subject_id,'\',subject_id,'_corrected.msh']);
% Output folder
s.pathfem = ([save_folder]);
s.poslist{1} = sim_struct('TMSLIST');
% Select coil
s.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');

%  4 mm distance between coil and head
s.poslist{1}.pos(1).distance = 0;
s.poslist{1}.pos.didt = 57e6;  % select desired intensity according to you dI/dt values
% Select coil centre
s.poslist{1}.pos.matsimnibs = TMS_target.dataC;

run_simnibs(s)


%% Extract passband for the deflect
mesh_load_gmsh4(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\MRI\simulation\',subject_id,'_corrected_TMS_1-0001_Magstim_70mm_Fig8_scalar.msh']);

% Create the binary mask (t) to create the spatial filter containing targeted structures
head = ans;

cortex.e =  head.triangles(head.triangle_regions == 1001,:);
cortex.p =  head.nodes;

% Make shure head.elements_data {2,1} corresponds to norme_E
E_field = head.element_data{2,1}.tridata(head.triangle_regions==1001);
percentile99_e_field = prctile(E_field,99);
delate = find(E_field>percentile99_e_field);
E_field(delate) = percentile99_e_field;

figure;
trisurf(cortex.e,cortex.p(:,1),cortex.p(:,2),cortex.p(:,3),E_field,'EdgeAlpha',0)
colorbar
%
t = E_field > 0.8*percentile99_e_field; % THIS A CRITICAL POINT, 0.1 is a arbitrary value to select the extension of the passband filter, MODIFY IT TO YOUR NEEDS
figure;
trisurf(cortex.e,cortex.p(:,1),cortex.p(:,2),cortex.p(:,3),double(t),'EdgeAlpha',0)
colorbar
%
E_field_copy = E_field;
E_field_copy(~t) = 0;
figure;
trisurf(cortex.e,cortex.p(:,1),cortex.p(:,2),cortex.p(:,3),E_field_copy,'EdgeAlpha',0)
colorbar
% Always check visually the filter and correct threshold if necessary



%%
%In case you want to correct outlayers with an apriory assumtion of the source reconstructed area (e.g., motor cortex) 
% use the following code

%select the region of interest
load(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\cortex_surf_corrected.mat'])

inds = [];
[inds] = select_sources_from_surface(cortex_surf,10,5,inds);
%inds_in_LFM = sort([3*inds - 2; 3*inds - 1; 3*inds]);
F = find(t);

% Select only overlapping elements between your region of interest and the
% E-field elements.
C = intersect(inds,F);
FG = zeros(size(t));
FG(C,1)= 1;

% 
figure;
trisurf(cortex.e,cortex.p(:,1),cortex.p(:,2),cortex.p(:,3),double(FG),'EdgeAlpha',0)
colorbar

E_field_copy_2 = E_field;
E_field_copy_2(~FG) = 0;
figure;
trisurf(cortex.e,cortex.p(:,1),cortex.p(:,2),cortex.p(:,3),E_field_copy_2,'EdgeAlpha',0)
colorbar

t_corrected = FG

%%

%% save passband
save_path = (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\']);
save_folder = [save_path,'passband\']
if ~exist(save_folder)
      mkdir([save_folder]);
end

%%
save([save_folder,'E_field.mat'],'E_field');
save([save_folder,'t_original.mat'],'t')
save([save_folder,'t.mat'],'t_corrected');
%save([save_folder,'t.mat'],'t');

close all
%% 
%% Simulation full mesh

% Import brainsight data and transform x axys (*be aware that X axis are inverted between
% brainsight and Simnibs space by default*)
importdata (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\NEURONAVIGATOR\Exported_Targets_s1.txt'])

TMS_target = ans;

TMS_target.dataC(:,4) = TMS_target.data(1,1:3);
TMS_target.dataC(:,1) = TMS_target.data(1,4:6);
TMS_target.dataC(1,1) = TMS_target.dataC(1,1)*-1;
TMS_target.dataC(2,1) = TMS_target.dataC(2,1)*-1;
TMS_target.dataC(3,1) = TMS_target.dataC(3,1)*-1;

TMS_target.dataC(:,2) = TMS_target.data(1,7:9);
TMS_target.dataC(:,3) = TMS_target.data(1,10:12);
TMS_target.dataC(4,1:3) = 0;
TMS_target.dataC(4,4) = 1;

%
% corrept depth coil:
TMS_target.dataC(3,4) = TMS_target.dataC(3,4)+6; % add 6mm distace coil-to scalp

%
save_path = ([path_to_mri_nifti ,'\']);

save_folder = [save_path,'simulationFull90\']
if ~exist(save_folder)
      mkdir([save_folder]);
end

% RUN SIMULATION
%Initialize  simulation session
s = sim_struct('SESSION');
s.fnamehead = ([path_to_mri_nifti ,'\m2m_',subject_id,'\',subject_id,'.msh']);
% Output folder
s.pathfem = ([save_folder]);

% Fields to maps
s.fields = 'eEjJ'; % Save the following results:
                   %  e: Electric field magnitude
                   %  E: Electric field vector
                   %  j: Current density magnitude
                   %  J: Current density vector
% set different output maps                   
s.map_to_surf = true;   %  Map to subject's middle gray matter surface
s.map_to_fsavg = true;  %  Map to FreeSurfer's FSAverage group template
s.map_to_vol = true;    %  Save as nifti volume
s.map_to_MNI = true;    %  Save in MNI space

% Get fields everywhere: 
s.tissues_in_niftis = 'all'
            %s.tissues_in_niftis = [1,2,3]; % Results in the niftis will be masked 
                                            % to only show WM (1), GM (2), CSF(3)
                                            % (standarE: only GM)

s.poslist{1} = sim_struct('TMSLIST');
% Select coil
s.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');

%  4 mm distance between coil and head
s.poslist{1}.pos(1).distance = 6;
s.poslist{1}.pos.didt = 102,6e6;  % select desired intensity according to you dI/dt values  11,4 = 10% / 34,2 = 30% / 57 = 50% / 79,8 = 70% / 102,6 = 90% / 114 = 100%
% Select coil centre
s.poslist{1}.pos.matsimnibs = TMS_target.dataC;

run_simnibs(s)

% % Import brainsight data and transform x axys (*be aware that X axis are inverted between
% % brainsight and Simnibs space by default*)
% importdata (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_00',subject_id,'\NEURONAVIGATOR\m1_position.txt'])
% 
% TMS_target = ans;
% 
% TMS_target.dataC(:,4) = TMS_target.data(1,1:3);
% TMS_target.dataC(:,1) = TMS_target.data(1,4:6);
% TMS_target.dataC(1,1) = TMS_target.dataC(1,1)*-1;
% TMS_target.dataC(2,1) = TMS_target.dataC(2,1)*-1;
% TMS_target.dataC(3,1) = TMS_target.dataC(3,1)*-1;
% 
% TMS_target.dataC(:,2) = TMS_target.data(1,7:9);
% TMS_target.dataC(:,3) = TMS_target.data(1,10:12);
% TMS_target.dataC(4,1:3) = 0;
% TMS_target.dataC(4,4) = 1;
% 
% TMS_target.centre = TMS_target.dataC(1:3,4)';
% TMS_target.centrey = TMS_target.dataC(1:3,2)';
% 
% % Corrept depth coil:
% % TMS_target.dataC(3,4) = TMS_target.dataC(3,4)+6; % add 6mm distace coil-to scalp
% 
% %
% save_path = ([path_to_mri_nifti ,'\']);
% 
% save_folder = [save_path,'simulationFull10b\']
% if ~exist(save_folder)
%       mkdir([save_folder]);
% end
% 
% % RUN SIMULATION
% %Initialize  simulation session
% s = sim_struct('SESSION');
% s.fnamehead = ([path_to_mri_nifti ,'\m2m_',subject_id,'\',subject_id,'.msh']);
% % Output folder
% s.pathfem = ([save_folder]);
% 
% % % Fields to maps
% s.fields = 'eEjJ'; % Save the following results:
%                    %  e: Electric field magnitude
%                    %  E: Electric field vector
%                    %  j: Current density magnitude
%                    %  J: Current density vector
% % % set different output maps                   
% s.map_to_surf = true;   %  Map to subject's middle gray matter surface
% s.map_to_fsavg = true;  %  Map to FreeSurfer's FSAverage group template
% s.map_to_vol = true;    %  Save as nifti volume
% s.map_to_MNI = true;    %  Save in MNI space
% % 
% % % Get fields everywhere: 
% s.tissues_in_niftis = 'all'
%             %s.tissues_in_niftis = [1,2,3]; % Results in the niftis will be masked 
%                                             % to only show WM (1), GM (2), CSF(3)
%                                             % (standarE: only GM)
% 
% s.poslist{1} = sim_struct('TMSLIST');
% % Select coil
% s.poslist{1}.fnamecoil = fullfile('legacy_and_other','Magstim_70mm_Fig8.ccd');
% 
% % 4 mm distance between coil and head
% s.poslist{1}.pos.didt = 11.4e6;  % select desired intensity according to you dI/dt values  11,4 = 10% / 34,2 = 30% / 57 = 50% / 79,8 = 70% / 102,6 = 90% / 114 = 100%
% % Select coil centre
% %s.poslist{1}.pos.matsimnibs = TMS_target.dataC;
% s.poslist{1}.pos.centre = TMS_target.centre;
% s.poslist{1}.pos.pos_ydir = TMS_target.centrey;
% s.poslist{1}.pos.distance = 6;
% 
% 
% run_simnibs(s)


%% Visualize results
%load headmesh corrected (with the electrode)
%head = mesh_load_gmsh4(['D:\MRI\',subject_id,'\m2m_MNI152\simulationtDCScer\MNI152_TDCS_1_scalar.msh']);

intensity = '10';
head = mesh_load_gmsh4(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\MRI\simulationFull',intensity,'\',subject_id,'_TMS_1-0001_Magstim_70mm_Fig8_scalar.msh']);


% Isolate electrode mask
electrode.e =  head.triangles(head.triangle_regions == 1002,:); %1002 for GM, 1001 for WM
electrode.p =  head.nodes;

% Make shure head.elements_data {2,1} corresponds to norme_E
E_field = head.element_data{2,1}.tridata(head.triangle_regions==1002);

% Reference values of Efield impact to the electrode.
Results.E_field_99pct = prctile(E_field,99);
Results.E_field_50pct = prctile(E_field,50);
Results.E_field_mean = mean(E_field);
Results

%visualize electrode alone with E-fields
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),E_field,'EdgeAlpha',0.5);
colorbar

% Make shure head.elements_data {2,1} corresponds to norme_E
t = E_field > 0.8*Results.E_field_99pct; % THIS A CRITICAL POINT, 0.1 is a arbitrary value to select the extension of the passband filter, MODIFY IT TO YOUR NEEDS
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),double(t),'EdgeAlpha',0.5);
colorbar

%
E_field_copy = E_field;
E_field_copy(~t) = 0;
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),E_field_copy,'EdgeAlpha',0.5)
colorbar

%% 



%% FAST E-field analyses

clear all
close all
clc, 

subject_id = '002';
intensity = '10';

head = mesh_load_gmsh4(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\MRI\simulationFull',intensity,'\',subject_id,'_TMS_1-0001_Magstim_70mm_Fig8_scalar.msh']);

% Isolate electrode mask
electrode.e =  head.triangles(head.triangle_regions == 1002,:); %1002 for GM, 1001 for WM
electrode.p =  head.nodes;

% Make shure head.elements_data {2,1} corresponds to norme_E
E_field = head.element_data{2,1}.tridata(head.triangle_regions==1002);

% Reference values of Efield impact to the electrode.
Results.E_field_max  = max(E_field);
Results.E_field_99pct = prctile(E_field,99);
Results.E_field_50pct = prctile(E_field,50);
Results.E_field_mean = mean(E_field);
Results


%visualize electrode alone with E-fields
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),E_field,'EdgeAlpha',0.5);
caxis([0, 250]);
c = colorbar;
c.Label.String = 'E-field Strength V/m';
c.FontSize =20;



% Make shure head.elements_data {2,1} corresponds to norme_E
t = E_field > 0.8*Results.E_field_99pct; % THIS A CRITICAL POINT, 0.1 is a arbitrary value to select the extension of the passband filter, MODIFY IT TO YOUR NEEDS
% figure;
 trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),double(t),'EdgeAlpha',0.5);
% c = colorbar;
% c.Label.String = 'E-field Strength V/m';
% c.FontSize =20;
title('Affected area >80%max','FontSize', 20)
%


E_field_copy = E_field;
E_field_copy(~t) = 0;
figure;
trisurf(electrode.e,electrode.p(:,1),electrode.p(:,2),electrode.p(:,3),E_field_copy,'EdgeAlpha',0.5)
caxis([0, 250]);
c = colorbar;
c.Label.String = 'E-field Strength V/m';
c.FontSize =20;
title('Affected area >80%max','FontSize', 20)


