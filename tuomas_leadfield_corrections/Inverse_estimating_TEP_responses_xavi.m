
close all; clear;

eeglab_path = '/Users/tuomasmutanen/Documents/Paris_Trip/ExampleData_Finland/Analysis/eeglab2022.1/';
addpath(eeglab_path);
eeglab;

% Load the leadfield in average reference:

load('/Users/tuomasmutanen/Documents/Paris_Trip/cRETMS/MRI/LFM_Exyz_ave.mat')
load('/Users/tuomasmutanen/Documents/Paris_Trip/cRETMS/MRI/cortex_surf.mat')

% Load the data of interest (should be in average reference):

EEG = pop_loadset('filename','all_sp_cer_real.set','filepath','/Users/tuomasmutanen/Documents/Paris_Trip/cRETMS/Processed/cRETMS_01_004/');

figure; pop_timtopo(EEG, [-20  300], [20   60   80  100 120 160], 'ERP data and scalp maps of Merged datasets');


%% Let's first test very simple minimum-norm estimate for estimating the current distributions
%
% The method:
%
% Hämäläinen, Matti S., and Risto J. Ilmoniemi. "Interpreting magnetic 
% fields of the brain: minimum norm estimates." Medical & biological 
% engineering & computing 32 (1994): 35-42.
%
% Has been used for instance in this classic Massimini paper:
%
% Massimini, Marcello, et al. "Breakdown of cortical effective connectivity 
% during sleep." Science 309.5744 (2005): 2228-2232.



for i = [20 60 80 100 120 160]

EEG_signal =mean(EEG.data(:,(1000+i),:),3);
N_chan = size(EEG.data,1);

L = LFM_Exyz_ave;
LL = LFM_Exyz_ave*LFM_Exyz_ave' ;
trL = trace(LL);

lambda = 0.01; % Regularization set to assume SNR=100, lambda = 1/SNR

G_MNE = L'/(LL + lambda*trL*eye(N_chan));

%R_MNE = G_MNE*L;

%mne_comp = LFM_Exyz_ave'*((LL + lambda*trL*eye(63))\EEG_signal);
mne_comp = G_MNE*EEG_signal;
MNE_power = mne_comp(1:3:end).^2 + mne_comp(2:3:end).^2 +mne_comp(3:3:end).^2;

figure;
subplot(1,2,1);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE_power,'EdgeAlpha',0,'FaceAlpha',1);
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

%% You can also estimate time courses with MNE

% First select your ROI
inds = [];
[inds] = select_sources_from_surface(cortex_surf,10,5,inds);

inds_in_LFM = sort([3*inds - 2; 3*inds - 1; 3*inds]);

%%
TEP = mean(EEG.data,3);
mne_time_comp = G_MNE*TEP;
MNE_time_power = mne_time_comp(inds_in_LFM(1:3:end),:).^2 + mne_time_comp(inds_in_LFM(2:3:end),:).^2 +mne_time_comp(inds_in_LFM(3:3:end),:).^2;

ROI_time_course = sum(MNE_time_power,1);
figure;
plot(EEG.times,ROI_time_course);

%% Next we test SLORETA, which is a MNE-like source estimator...

% Let's run the LORETA correction to MNE
N_source = size(G_MNE,1);
G_SLOR = G_MNE;
for i = 1:N_source
    W_slor_i = G_MNE(i,:)*L(:,i);
    G_SLOR(i,:) = 1/sqrt(W_slor_i)*G_MNE(i,:);
end


% Then we inverse estimate the source distributions:

for i = [20 60 80 100 120 160]

EEG_signal =mean(EEG.data(:,(1000+i),:),3);


SLORETA_comp = G_SLOR*EEG_signal;
SLORETA_power = SLORETA_comp(1:3:end).^2 + SLORETA_comp(2:3:end).^2 +SLORETA_comp(3:3:end).^2;

figure;
subplot(1,2,1);
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),SLORETA_power,'EdgeAlpha',0,'FaceAlpha',1);
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

%% Let's visualize the time course of SLORETA


TEP = mean(EEG.data,3);
SLOR_time_comp = G_SLOR*TEP;
SLOR_time_power = SLOR_time_comp(inds_in_LFM(1:3:end),:).^2 +...
    SLOR_time_comp(inds_in_LFM(2:3:end),:).^2 + SLOR_time_comp(inds_in_LFM(3:3:end),:).^2;

ROI_time_course = sum(SLOR_time_power,1);
figure;
plot(EEG.times,ROI_time_course);



