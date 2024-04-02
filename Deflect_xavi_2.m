
%******************** Apply deflect for every dataset to filter  noise  *************


% The currect script apply the deflect filter based on a spatial filter (
% passband) of the source individual source projected E-field. Please run
% "construct_lead_field_and_passband_for_deflect.at" before this one. 

% This script needs to be applyed for every dataset in isolate.
% To run this script, modify the subject specs and the datasets specs and
% you can run all cells together.



% Xavier COROMINAS, 11/03/2023, Paris, France.


%%
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

% Add path to helsinki bem framework
addpath('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions\hbf_lc_p-master\');
hbf_SetPaths

% Add path to DeFleCt
addpath('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions\DeFleCT_8Nov13')
addpath('F:\LI_rTMS_project\eeg_analysis\Analysis\External_functions')

% Add eeglab path
eeglab_path = 'D:\SOFTWARE\eeglab2022.0';
addpath(eeglab_path);
eeglab;

% Set subject preferences

% Subject preferneces
subject_id = '002';
sequence = '_S5';

% BLOCK NUMBER
block = '1';   % '1' , '2'
% STIMULATION intensity
stim ='90';  %  10 , 30 , 50 , 70 , 90 , 100
%condition
condition = 'sham'; %sham real

        dataset =[block,'_',stim,'_M1'];


% Add paths to leafields, passband and MRI
addpath(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS'])
addpath(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\passband'])
path_to_mri_nifti = ['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\MRI\'];


%
subject_id_mri = '025';
sequence_mri = 'S5';
%
filename = ['v_CRETMS_03_',subject_id_mri,'_',sequence_mri,'_T1w_mp2rage_UNI_Images_MPRAGEised_biascorrected.nii'];
    


%% Now let's apply the deflect to YOUR dataset 

% Load the leadfield in average reference:
load(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\cortex_surf_corrected.mat'])

%Load passband
load(['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS\passband\t.mat'])

% Load the data of interest (should be in average reference):
EEG = pop_loadset('filename',([block,'_',stim,'_M1_',condition,'_Cleaned.set']),'filepath',(['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_',subject_id,'\']));
eeglab redraw


%% Find indicies of elements in mesh simulation

% localize inds of ROI (t) mask
t = t_corrected; 
[inds] = find(t);
%inds_in_LFM = sort([3*inds - 2; 3*inds - 1; 3*inds]);

%localize vertices (noces) of ROI ( t ) mask
cortex.p = meshcentroid(cortex_surf.p,cortex_surf.e);
cortex.nn = surfacenorm(cortex_surf.p,cortex_surf.e);

%first plot filter 
figure;
h = trisurf(cortex_surf.e,cortex_surf.p(:,1),cortex_surf.p(:,2),cortex_surf.p(:,3),double(t),'EdgeAlpha',0,'FaceAlpha',1)
view(-35,50)
camzoom(2.5);
axis off;
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
title(['ROI leackage'])
colorbar
%cmap=jet(20);
%cmap=cmap(11:20,:);
%colormap(cmap)
hold on 
inds_points = find(t)
points = cortex.p(inds_points,:)
plot3(points(:,1),points(:,2),points(:,3),'.k','MarkerSize',15);
view(-23,90)


%%
% Visualize crosstalk sensitivity with just only a passband
cd (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject_id,'\LEADFIELDS']);
load LFM_Edir_ave_corrected.mat
load cortex_surf_corrected.mat
load LFM_Exyz_ave.mat



% In case there are NAN rows, correct them;
A = LFM_Edir_ave; 
% Define the size of the neighborhood
neighborhood_size = 50; % Change this value to adjust the neighborhood size
% Iterate through each element
for i = 1:size(A, 1)
    for j = 1:size(A, 2)
        % Check if the element is NaN
        if isnan(A(i, j))
            % Define the bounds for surrounding elements
            row_lower = max(1, i - floor(neighborhood_size/2));
            row_upper = min(size(A, 1), i + floor(neighborhood_size/2));
            col_lower = max(1, j - floor(neighborhood_size/2));
            col_upper = min(size(A, 2), j + floor(neighborhood_size/2));
            
            % Extract surrounding elements
            surrounding_values = A(row_lower:row_upper, col_lower:col_upper);
            
            % Calculate the mean of surrounding values (excluding NaNs)
            surrounding_values = surrounding_values(~isnan(surrounding_values));
            mean_value = mean(surrounding_values);
            
            % Replace NaN with the mean of the surrounding values
            A(i, j) = mean_value;
        end
    end
end
LFM_Edir_ave = A;


%Mne estimator for cortically constraiend dipols
SNR2= 10;
L = LFM_Edir_ave;
truncval = 4;
N_chan = 63;
C = trace(L*L')*eye(N_chan);
Linv=MNestimator(L,C,SNR2,truncval);

ct_M1tip=Linv(inds,:)*L(:,:);
plotdata1=abs(sum(ct_M1tip,1));%

% Plot crosstalk leachage of the ROI.
figure;
h = trisurf(cortex_surf.e,cortex_surf.p(:,1),cortex_surf.p(:,2),cortex_surf.p(:,3),double(plotdata1),'EdgeAlpha',0,'FaceAlpha',1)
view(-35,50)
camzoom(2.5);

axis off;
cmap=jet(20);
cmap=cmap(11:20,:);

lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
%h.Marker = ".";
title(['ROI leackage'])
colorbar
%colormap(cmap)
hold on 
inds_points = find(t)
points = cortex.p(inds_points,:)
plot3(points(:,1),points(:,2),points(:,3),'.k','MarkerSize',20);
%

%% run deflect to observe change on the corstalk filter by just adding a passband 
% Makes a DeFleCT spatial filter with given passband and stopband. 
% passband:     indices to sources for which the targeted output is 1
% SVDpassband:  how many components represent the passband (optional)
% force_passband_flag:  forces the output for all passband components to 1
% stopband:             indices to sources for which the output is 0
% SVDstopband:  how many components represent the stopband (optional)
% L:  forward model
% C:    noise covariance matrix (or measurement covariance matrix)
% SNR2: assumed signal-to-noise ratio (for setting regularization)
% Csvdtrunc: number of truncated singular values for the inversion of C.
% Whitener:  the whitener that is applied to the leadfields (optional; if C
%            is given, the whitener is built from C)


% We will USE THIS ONE FOR OUR DATA. NAMEd as: w.deflect

passband = inds; 
SVDpassband = 5;
force_passband_flag = [];
stopband = [];%stop_inds_in_LFM ;% [];
SVDstopband = [];
C = trace(LFM_Edir_ave*LFM_Edir_ave')*eye(N_chan);
SNR2 = 30;  % higher SNR better final results, this is a critical value
Csvdtrunc = 4;
Whitener = [];

w_deflect=DeFleCT(passband,SVDpassband,force_passband_flag,stopband,SVDstopband,...
    LFM_Edir_ave,C,SNR2,Csvdtrunc); % if there is an error normaly is because the number of channels of the leadfield matrix do not match the covarance matrix, better redo leadfield creation to be shure that there are 63

disp('Deflect filter ready...')
ct_deflect=w_deflect*LFM_Edir_ave(:,:);
plotdata2=abs(ct_deflect);

figure;
h = trisurf(cortex_surf.e,cortex_surf.p(:,1),cortex_surf.p(:,2),cortex_surf.p(:,3),double(plotdata2),'EdgeAlpha',0,'FaceAlpha',1)
view(-35,50)
camzoom(2.5);
axis off;
cmap=jet(20);
cmap=cmap(11:20,:);
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
title(['FILTER'])
colorbar
%colormap(cmap)
hold on 
inds_points = find(t)
points = cortex.p(inds_points,:)
plot3(points(:,1),points(:,2),points(:,3),'.k','MarkerSize',20);
%TEP = mean(EEG.data,3); % 
%ROI_time_course = w_deflect*TEP; % mean data time-course along all trials of ROI



%% As a test.... Run deflect to observe change on the corstalk filter by adding a passband and the stopband
%{




%This spatial filter is going to be named as: w


% In case you want a particular ROI as a stopband use this code:
%{
      stopband_roi = [];
      [stopband_roi] = select_sources_from_surface(cortex_surf,10,5,stopband_roi);
%}

% visualize crosstalk leackage with passband and stopband: here we select
% all elements out of our ROI target as elements for the stopaband as
% example
stopband_roi = find(~t);
stop_inds_in_LFM = sort([3*stopband_roi - 2; 3*stopband_roi - 1; 3*stopband_roi]);


% Makes a DeFleCT spatial filter with given passband and stopband. 
% passband:     indices to sources for which the targeted output is 1
% SVDpassband:  how many components represent the passband (optional)
% force_passband_flag:  forces the output for all passband components to 1
% stopband:             indices to sources for which the output is 0
% SVDstopband:  how many components represent the stopband (optional)
% L:  forward model
% C:    noise covariance matrix (or measurement covariance matrix)
% SNR2: assumed signal-to-noise ratio (for setting regularization)
% Csvdtrunc: number of truncated singular values for the inversion of C.
% Whitener:  the whitener that is applied to the leadfields (optional; if C
%            is given, the whitener is built from C)



passband = inds; % !!!maybe this values should not be inds in the elafield, should only be inds normals!!!!
SVDpassband = 1;
force_passband_flag = [];
stopband = stopband_roi;%stop_inds_in_LFM ;% [];
SVDstopband = [] ;

C = trace(LFM_Edir_ave*LFM_Edir_ave')*eye(N_chan);
SNR2 = 30;  % higher SNR better final results, this is a critical value
Csvdtrunc = 4;
Whitener = [];

w=DeFleCT(passband,SVDpassband,force_passband_flag,stopband,SVDstopband,...
    LFM_Edir_ave,C,SNR2,Csvdtrunc); % if there is an error normaly is because the number of channels of the leadfield matrix do not match the covarance matrix, better redo leadfield creation to be shure that there are 63

disp('Deflect filter ready...')
ct_deflect=w*LFM_Edir_ave(:,:);
plotdata2=abs(ct_deflect);
figure;
h = trisurf(cortex_surf.e,cortex_surf.p(:,1),cortex_surf.p(:,2),cortex_surf.p(:,3),double(plotdata2),'EdgeAlpha',0,'FaceAlpha',1)
view(-35,50)
camzoom(2.5);
axis off;
cmap=jet(20);
cmap=cmap(11:20,:);
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
title(['ROI leackage'])
colorbar
colormap(cmap)
hold on 
inds_points = find(t)
points = cortex.p(inds_points,:)
plot3(points(:,1),points(:,2),points(:,3),'.k','MarkerSize',20);





%}
%% Plot source reconstructed tep of original data unisng classic MNE
TEP = mean(EEG.data,3); % 
ROI_time_course = w_deflect*TEP; % mean data time-course along all trials of ROI

% make a copy of original data
EEG_copy = EEG;
EEG_tmp = EEG;

% Plot tep of originial data
for i = [20 50 100 150 250]
EEG_signal_org =mean(EEG_tmp.data(:,(1000+i),:),3);

%set leadfields and covariance matrix
LFM = LFM_Edir_ave;
LL = LFM*LFM' ;

% Calculate MNE
M = 15; % show only 15 largest components
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';
tmp = LFM'*(LL_inv*EEG_signal_org);
MNE = tmp;
MNE = normalize(MNE, 'range', [-1.5 1.5]);

figure;
subplot(1,2,1);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)
camzoom(1.5);
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
topoplot(EEG_signal_org,EEG_tmp.chanlocs)

end



% PASS ALL DATA FOR THE  FILTER. 
% We store the data in the C3 position (column 7) ( electrode to facilitate further analyses os
% ITC and corticocortical sincronization. The data is in microvols, the
% data remaining is in miivols, hence, the row 7 contains data in
% microvolts and the rest in
% milivols.*******************************************************************


    for i = 1:size(EEG.data,3)
        EEG.data(7,:,i) = w_deflect*EEG.data(:,:,i); 
    end

% Remeber, data is a diferent scale, hence apply this step here.
% In case you want to transform the data at that point to vols jus have to
% multiple the data by 0.001. after the deflect application so: 
% 
  %  EEG.data(7,:,:) = EEG.data(7,:,:)*0.001;

eeglab redraw

% Plot the data again
for i = [20 50 100 150 250]

EEG_signal =mean(EEG.data(:,(1000+i),:),3);

% Calculate MNE
%set leadfields and covariance matrix
LFM = LFM_Edir_ave;
LL = LFM*LFM' ;

% Calculate MNE
M = 15; % show only 15 largest components
[U,S,V] = svd(LL,'econ');
S_diag = diag(S);
S_inv = [1./S_diag(1:M);zeros(size(S,1)-M,1)];
S_inv = diag(S_inv);
LL_inv = V*S_inv*U';
tmp = LFM'*(LL_inv*EEG_signal);
MNE = tmp;

MNE = normalize(MNE, 'range', [-1.5 1.5]);

figure;
subplot(1,2,1);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)
axis off;
%hold on
%plot3(points(:,1),points(:,2),points(:,3),'.k','MarkerSize',20);
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


%% Plot TF and ITC of C3 = ROI passband.
close all
figure; pop_newtimef( EEG, 1, 7, [-2000   1999], [3         0.8] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C3', 'baseline',[-500 -10], 'erspmax', [[-8 8]], 'plotphase', 'off', 'padratio', 1, 'winsize', 1000);
figure; pop_newtimef( EEG_tmp, 1, 7, [-2000   1999], [3         0.8] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C3', 'baseline',[-500 -10], 'erspmax', [[-8 8]], 'plotphase', 'off', 'padratio', 1, 'winsize', 1000);


%% save eeg dataset
%EEG = EEG_tmp;
save_path =['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_',subject_id];

save_folder = [save_path,'\deflect\'];

if ~exist(save_folder)
      mkdir([save_folder]);
end 

EEG = pop_saveset( EEG, 'filename',[dataset,'_',condition,'_Cleaned_deflect.set'],'filepath',[save_folder]);

