
%******************** Apply spatial filter for every individual dataset  *************


% The currect script apply the a MNE-like filter based of the  individual source projected E-field. 
% Please preprocess your TMS-EEG data and run "Leadfield_mesh_passband.mat" before this one. 

% This script needs to be applyed for every dataset.

% To run this script, modify the subject specs and the datasets specs.

% Xavier COROMINAS, Paris, France.
% Contributors: Xavier COROMINAS, Martina BRACCO, Tuomas P MUTANEN.

% Modify paths <> along the script.

%% Add paths
clear;
clc;
close all;

% Set up the paths to needed tools:
% Simnibs for segmenting head:
addpath('<>\SimNIBS-4.0\simnibs_env\Lib\site-packages\simnibs\matlab_tools\')

% Add path to ISO2Mesh
addpath('<>\iso2mesh-1.9.6\')

% Add path to helsinki bem framework
addpath('<>\External_functions\hbf_lc_p-master\');
hbf_SetPaths

% Add path to DeFleCt
addpath('<>\External_functions\DeFleCT_8Nov13')
addpath('<>\External_functions')

% Add eeglab path
eeglab_path = '<>\eeglab2022.0';
addpath(eeglab_path);
eeglab;

%% Set subject preferences

% Subject preferneces
subject_id = 'XX';

%eeg dataset
dataset =['<your_tms_eeg_dataset>'];


% Add paths to subject leafields, passband and MRI
addpath(['<>\LEADFIELDS'])
addpath(['<>\LEADFIELDS\passband'])
path_to_mri_nifti = ['<>\MRI\'];


%% Now let's apply the filter to YOUR dataset 

% Load the leadfield in average reference:
load(['<>\LEADFIELDS\cortex_surf_corrected.mat'])

%Load passband
load(['<>\LEADFIELDS\passband\t.mat'])

% Load the data of interest (should be in average reference):
EEG = pop_loadset('filename',(['<your_tms_eeg_dataset>.set']),'filepath',(['<path_to_your_tms_eeg_dataset>\']));
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


%% Visualize crosstalk sensitivity with just only a passband

cd (['<>\LEADFIELDS']);
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


%MNE estimator for cortically constraiend dipols
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

%% Run deflect framework to observe changes on the corstalk filter by just adding a passband 

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


% We will USE THIS ONE FOR OUR DATA. NAMED as: w.deflect

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



%% As a test.... you can try to tun deflect to observe change on the corstalk filter by adding a stopband
%{


%This spatial filter is going to be named as: w

% In case you want a particular cortical regions as a stopband use this code:
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



%% PASS ALL DATA FOR THE  FILTER. 
% IMPORTANT*** We store the data in the C3 position (column 7) ( electrode to facilitate further analyses os
% ITC and corticocortical sincronization. The data is in microvols, the
% data remaining is in milivols, hence, the row 7 contains data in
% microvolts and the rest in
% milivols.*********


    for i = 1:size(EEG.data,3)
        EEG.data(7,:,i) = w_deflect*EEG.data(:,:,i); 
    end

% Remeber,filtered data is a diferent scale, hence apply this step here.
% In case you want to transform the data at that point to vols. Just have to
% multiply the data by 0.001. after the filter, so: 

  %  EEG.data(7,:,:) = EEG.data(7,:,:)*0.001;


eeglab redraw

% Plot the data again
for i = [20 50 100 150 250]

EEG_signal =mean(EEG.data(:,(1000+i),:),3);

% Calculate MNE
%set leadfields and covariance matrix
LFM = LFM_Edir_ave;
LL = LFM*LFM' ;

% Calculate MNE. Use a MNE with a  value decomposition keeping the largest
% components

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


%% Plot ERP, ERSP and ITC of C3 = ROI passband.

close all

% Analyses of fitlered data
figure; pop_newtimef( EEG, 1, 7, [-2000   1999], [3         0.8] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C3', 'baseline',[-500 -10], 'erspmax', [[-8 8]], 'plotphase', 'off', 'padratio', 1, 'winsize', 1000);

% Analyses of original data
figure; pop_newtimef( EEG_tmp, 1, 7, [-2000   1999], [3         0.8] , 'topovec', 1, 'elocs', EEG.chanlocs, 'chaninfo', EEG.chaninfo, 'caption', 'C3', 'baseline',[-500 -10], 'erspmax', [[-8 8]], 'plotphase', 'off', 'padratio', 1, 'winsize', 1000);


%% save eeg dataset
%EEG = EEG_tmp;
save_path =['<>'];

save_folder = [save_path,'\deflect\'];

if ~exist(save_folder)
      mkdir([save_folder]); 
end 

EEG = pop_saveset( EEG, 'filename',[dataset,'_filtered.set'],'filepath',[save_folder]);

%% END

