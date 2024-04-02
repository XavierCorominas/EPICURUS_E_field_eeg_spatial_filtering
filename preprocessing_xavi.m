

%% Start assigning paths
clear all
close all
clc

%
FieldTrip_path = ['D:\SOFTWARE\fieldtrip-20220625'];
addpath(FieldTrip_path);
ft_defaults;

%
eeglab_path = ['D:\SOFTWARE\eeglab2022.0'];
addpath(eeglab_path);
eeglab;

%
Tesa_path = ['D:\SOFTWARE\eeglab2022.0\plugins\TESA1.1.1'];
addpath(Tesa_path);

%
External_functions = ['F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions'];
addpath(External_functions);

save_path =['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\'];


%% SET CONDITION

% SUBJECT ID
subject = '002';   


% BLOCK NUMBER
block = '1';   % '1' , '2'
% STIMULATION intensity
stim ='50';  %  10 , 30 , 30 , 70 , 90 , 100
% TARGET
%target = 'M1';  


% TRIGGERS -S1 real , S2 sham
% triggercode = ['S  1'];    %    'S  1'   ,   'S  2'
% condition = 'real'; %sham

% 
triggercode = ['S  1'];    %    'S  1'   ,   'S  2'
condition = 'real'; %sham


%%

mainpath = (['F:\BRAINMAG_Li-rTMS\DATA\BrainNMAG_05_',subject,'\EEG\']);
%mainpath = (['F:\LI_rTMS_project\eeg_analysis\DATA\',subject,'\EEG\']);

addpath(mainpath);
cd(mainpath);

% % File path:
raw_data_path = [mainpath,'1session_threshold_estimation\'];

% CONDITION
if triggercode == 'S  1'
    cond = 'real';
elseif triggercode == 'S  2'
    cond = 'sham';
end 

% DATASET
dataset =[block,'_',stim,'_','M1'];

%SAVE FOLDER
save_folder = [save_path,'\BraiNMAG_05_',subject,'\'];
if ~exist(save_folder)
      mkdir([save_folder]);
end


%% OPEN FILE

% This reads the raw data to EEGLAB format
EEG = pop_loadbv([raw_data_path], [dataset,'.vhdr'], [], []); 


% If the file is not found localize the header file (.vhdr of the dataset
% on your data dolfer) and open it in notebook. Then check that the
% markerfile and datafile names are correct.


EEG= pop_chanedit(EEG, 'lookup',[eeglab_path,'\plugins\dipfit\standard_BESA\standard-10-5-cap385.elp'],'changefield',{1 'type' 'EEG'},'changefield',{2 'type' 'EEG'},'changefield',{3 'type' 'EEG'},'changefield',{4 'type' 'EEG'},'settype',{'1:63' 'EEG'});

 eeglab redraw 
 
%% This section checks whether the dataset has triggers and corrects if necessary 

% Check manually the EEG.event to be shure that there is nothing strange

% when Tesa Findpulse is activated: set the rate of charge of TMS artifact
% to 20e4 (or 10e4), no need to modify any other parameter. You can modify
% the electrode reference if you need it, here we selected C3 as the coil
% has been placed above it and the artefacts are going to be stronger.

   if length(EEG.event) < 260
    disp('The data lacks triggers! Checking automatically the onsets of the trials.')
    EEG_tmp = pop_tesa_findpulse( EEG, 'CP3', 'refract', 3, 'rate',20e4); % if 240 pulses are not found, decrease the rate

    load('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\Preprocessig_xavi_tuomas\Old\rhythmic_markers.mat')
    rhytmic_markers = load('F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\Preprocessig_xavi_tuomas\Old\rhythmic_markers.mat')
   end

EEG_org = EEG;

%Correct for non stimulus codes
copy_array = EEG.event;
% Find indices where code differs from 'Stimulus'
indices_to_remove = find(~strcmp({copy_array.code}, 'Stimulus'));
% Remove rows from the matrix
copy_array(indices_to_remove) = [];
EEG.event = copy_array;


% Find If there si any strange trigger
s4Indices = find(strcmp({EEG.event.type}, triggercode));
EEG_m = EEG.event(s4Indices);
if strcmp(triggercode, 'S  1') && length(EEG_m) == 160 && all(strcmp({EEG_m.code}, 'Stimulus'))
    disp('No artefacts missing for S1');
elseif strcmp(triggercode, 'S  2') && length(EEG_m) == 80 && all(strcmp({EEG_m.code}, 'Stimulus'))
    disp('No artefacts missing for S2');
else
    disp('Artefacts missing');
end
%


%% If artefacts missing, Check manually EEG.event and delate rows of non real events. If after delating rows there are less than 140 events in the 
% rhythmic/arrhythmic or less 60 in the SP, then localize the missing
% events and copy paste the correspongind rows from the EEG_tmp.event to
% EEG.event



%% Change the name of the markers for repetitive stimulation

if  ~contains(stim,'sp');

x = length({EEG.event.latency});

for n = 2:4:x
EEG.event(1,n).type = 'TMS2';
end
for n = 3:4:x
EEG.event(1,n).type = 'TMS3';
end
for n = 4:4:x
EEG.event(1,n).type = 'TMS4';
end


end
%% 

% Do first save (only for real datasets to not duplicate unnecessarly
% saving space)
% if contains(cond,'real')
% EEG = pop_saveset( EEG, 'filename',[ dataset,'_' cond, '_events.set'],'filepath',[save_path,'LITMS_01_',subject]);
% 
% EEG_1= EEG;
% end
% EEG = EEG_1;
% EEG = pop_loadset('filename',[ dataset,'_real' , '_events.set'],'filepath',[save_path,'LITMS_01_',subject]);
% 

%%

% EPOCH DATA (-1.2 ms to 1.2 ms)
% and select the trials of interest:

EEG = pop_epoch( EEG, {  triggercode  }, [-1.2         1.2], 'newname', [dataset, ' epochs'], 'epochinfo', 'yes');

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,30000:30400)')
%EEG = pop_rmbase( EEG, [-1200 -100] ,[]);
EEG_plot = pop_reref( EEG, []);
figure; pop_timtopo(EEG_plot, [-130  230],  [], 'ERP data and scalp maps of  resampled');
figure; pop_timtopo(EEG_plot, [-5  10],  [], 'ERP data and scalp maps of  resampled');

figure;
plot(EEG.times,mean(EEG.data,3)','k'); xlim([-5 10]); ylim([-100 100])


% HERE IN CASE THAT THERE IS A ELECTRODE VERY BAD WITHA  VERY BIG DECAY WE
% CAN INTERPOLATE IT WITH EEGLAB. iN THE gui GO TO TOOLS-->INTERPOLATE
% ELECTRODE AND SELECT MANUALLY YHE ELECTRODE

%% check TMS latency if repetitive stimulation
if  ~contains(stim,'sp');

    % check TMS latency if repetitive stimulation
for x = 1:length(EEG.epoch);
TMS_lat(x,:) = cell2mat([EEG.epoch(x).eventlatency])
end 


TMS2 = mean(TMS_lat(:,2))
TMS2_sd = std(TMS_lat(:,2))
TMS3 = mean(TMS_lat(:,3) - TMS_lat(:,2))
TMS3_sd = std(TMS_lat(:,3) - TMS_lat(:,2))
TMS4 = mean(TMS_lat(:,4) - TMS_lat(:,3))
TMS4_sd = std(TMS_lat(:,4) - TMS_lat(:,3))

end



%% Demean data (-1200 ms to 1200 ms)
EEG = pop_rmbase( EEG, [-500  -10]);

% Visualize data
EEGLAB_plot_EEG_wrapper(EEG, [-1200 1200], 20, 1:63)

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,29000:40400)')
figure; pop_timtopo(EEG, [-100  500],  [], 'ERP data and scalp maps of  resampled');
figure; pop_timtopo(EEG, [0 30],  [], 'ERP data and scalp maps of  resampled');


%% Remove TMS pulse artifact and peaks of TMS-evoked muscle activity (-2 to 10 ms)... 
% and interpolate missing data around TMS pulse

 %EEG = pop_tesa_removedata( EEG, [-2 4]); Check for every subject and
 %modify if necessary for everyone. 
 % 
% For 10% [-2 4]
%  30% [-2 6]
%  30% [-2 8]
%  70% [-2 10]
%  90% [-2 12]
%  100% [-2 14]

 EEG = pop_tesa_removedata( EEG, [-2 14], [], {triggercode, 'TMS2', 'TMS3', 'TMS4'});


%% Visualize data
EEGLAB_plot_EEG_wrapper(EEG, [-1200 1200], 20, 1:63, 20)

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,30000:30400)')

plot(data(:,29000:40430)')
figure; pop_timtopo(EEG, [-230  300],  [], 'ERP data and scalp maps of  resampled');


%% Remove the recharging artifact
% First zoom to better visialize the artifact
% Press ENTER when ready
% Select the correct time window
% Press ENTER when done



EEG = EEGLAB_remove_and_interpolate_recharging_artifact(EEG, 2); % 1 + single pulse, 2 = repetitive pulse


%% Visualize data
EEGLAB_plot_EEG_wrapper(EEG, [-1200 1200], 20, 1:63)

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,30000:30400)')




%% Downsample data (23000 Hz to 1000 Hz), Strange downsample error sometime ocurre, first delate signal path and readd it 
rmpath('C:\Program Files\MATLAB\R2022a\toolbox\signal\signal');
addpath('C:\Program Files\MATLAB\R2022a\toolbox\signal\signal');
%% 
EEG = pop_resample( EEG, 1000);

%% Prepare the data for ICA: First remove for a while the very worst channels!

EEG_evo = mean(EEG.data,3);
[~, sigmas] = DDWiener(EEG_evo);
figure;
plot(sigmas,'*')


%% labeling the very worst channels to not affect the ICA run
badC = find(sigmas > (median(sigmas) + 5*std(sigmas)));
goodC = setdiff(1:length(sigmas),badC);

EEG2ICA = pop_select( EEG, 'nochannel', badC);


%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition


tmpdata=reshape(EEG2ICA.data,[size(EEG2ICA.data,1) size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]); 
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG2ICA.data,2)*size(EEG2ICA.data,3)]);
tempCpca=svd(tmpdata);
th=0.01;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th));
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca


%% 3. ICA - First round

Cpca = 60;

EEG.ICA_dims1 = Cpca;

EEG2ICA=ICA_analysis(EEG2ICA,EEG2ICA.data,Cpca);% run ICA using the runica method from EEGlab

EEGpart.ICA=EEG2ICA;
EEG2ICA.comp2remove=[];
'ICA has been computed'



%% select and remove artefact-contaminated independent components (ICs)
% this function opens a GUI that allows you to inspect single independent
% components and to decide whether to remove or to keep them, looking at
% the topography of ICs weights (bottom right), at the time course of the
% ICs (bottom left), at the power spectrum of the ICs (center right), at
% the single-trial time course of the ICs (center left), at the original
% average TMS-evoked potentials (top left), at the average TMS-evoked
% potentials after subtraction of the ICs that represent artefacts


ICA_selectandremove(EEG2ICA)
% click on button LOAD


%% Update the data after choosing the bad ICs

EEG.ICA2rem1 = length(data_GUI.comp2remove);
EEG2ICA.data = data_GUI.compproj;

EEG.data(goodC,:,:) = EEG2ICA.data;
EEG.goodC = goodC;
EEG.badC = badC;


%% Visualize data
EEGLAB_plot_EEG_wrapper(EEG, [-1200  1200], 20, 1:size(EEG.data,1))

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,1151:1701)')
figure; pop_timtopo(EEG, [-30  300],  [], 'ERP data and scalp maps of  resampled');

%% Perform the baseline correction again after ICA:
EEG = pop_rmbase( EEG, [-500 -10] ,[]);

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,1151:1701)')
figure; pop_timtopo(EEG, [-30  300],  [], 'ERP data and scalp maps of  resampled');


%% Basic high-pass filtering from 2 Hz

EEG = tesa_filtbutter( EEG, 2, [], 4, 'highpass' );

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,1151:1701)')

%% Save the current form of dataset

% EEG = pop_saveset( EEG, 'filename',[ dataset,'_' cond, '_ICA.set'],'filepath',[save_path,'LITMS_01_',subject]);
% 
% close all

%%
%%%
%%%
% dataset
% cond
% EEG = pop_loadset('filename',[dataset,'_' cond,'_ICA.set'],'filepath',['F:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);


%%
% Sound step

EEG = pop_tesa_sound(EEG, 'lambdaValue', 0.1, 'iter', 10 );



%% Visualize data
EEGLAB_plot_EEG_wrapper(EEG, [-1200  1200], 10, 1:size(EEG.data,1))

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,1151:1701)')


%% SSPSIR step

[EEG] = pop_tesa_sspsir(EEG, 'artScale', 'manual','timeRange', [0,210]);

EEGLAB_plot_EEG_wrapper(EEG, [-1300  1300], 10, 1:size(EEG.data,1));

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,1151:1701)')




%% Extend data removal to 12 ms (-2 to 12 ms)

% 
% if contains(target, 'li' );
%     if contains(stim, 'rith') || contains(stim,'arith');
%     EEG = pop_tesa_removedata( EEG, [-2 6], [], {triggercode, 'TMS2', 'TMS3', 'TMS4'});
%     elseif contains(stim, 'sp');
%     EEG = pop_tesa_removedata( EEG, [-2 6], [], {triggercode});
%     end 
% end 
% 
% 
% 
% if contains(target, 'hi' );
%     if contains(stim, 'rith') || contains(stim,'arith');
 %    EEG = pop_tesa_removedata( EEG, [-2 14], [], {triggercode, 'TMS2', 'TMS3', 'TMS4'});
%     elseif contains(stim, 'sp');
%     EEG = pop_tesa_removedata( EEG, [-2 14], [], {triggercode});
%     end 
% end 
% 


 EEG = pop_tesa_removedata( EEG, [-2 14], [], {triggercode, 'TMS2', 'TMS3', 'TMS4'});

%% 
clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,1151:1701)')


%% Interpolate missing data around TMS pulse
EEG = pop_tesa_interpdata( EEG, 'cubic', [1,1] );

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,930:1901)')
%}
%% Doing the final band-pass filtering

EEG = pop_tesa_filtbutter( EEG, 2, 100, 4, 'bandpass' );

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,1151:1701)')

%% Trials rejection - 1
close all
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed for first round of trials rejection'

 
%% Trials rejection - 2

if ~isempty(TMPREJ)
    [trialrej elecrej]=eegplot2trial(TMPREJ,size(EEG.data,2),size(EEG.data,3));
else
    trialrej=[];
end

EEG.BadTr =find(trialrej==1);
EEG = pop_rejepoch( EEG, EEG.BadTr ,0); 



%% Before applying ICA to reduce artifacts, check which is the maximum number of independent components via singular value decomposition

tmpdata=reshape(EEG.data,[size(EEG.data,1) size(EEG.data,2)*size(EEG.data,3)]); 
tmpdata=tmpdata-repmat(mean(tmpdata,2),[1,size(EEG.data,2)*size(EEG.data,3)]);
tempCpca=svd(tmpdata);
th=0.01;% this threshold is arbitrary but it works in most cases
figure,semilogy(tempCpca,'.-')% plot of the eigenvalues of the data matrix
a=ylim;
hold on
Cpca=max(find(tempCpca>th))
line([Cpca Cpca],[a(1) a(2)],'Color','k','LineStyle','--');
hold off
'check whether the vertical dashed line correctly separates the eigenvalues close to zero from the ones much higher than zero, otherwise change the threshold th'

% if the vertical dotted line correctly separates the eigenvalues close to
% zero from the ones much higher than zero, Cpca correctly corresponds to
% the total number of independent components that can be estimated from the
% data matrix. if this is not the case, you have to change the threshold above or to manually correct Cpca

%% Final ICA for residual artifacts

Cpca = 50;

EEG.ICA_dims1 = Cpca;

EEG2ICA=ICA_analysis(EEG,EEG.data,Cpca);% run ICA using the runica method from EEGlab

EEGpart.ICA=EEG2ICA;
EEG2ICA.comp2remove=[];
'ICA has been computed'



%% select and remove artefact-contaminated independent components (ICs)
% this function opens a GUI that allows you to inspect single independent
% components and to decide whether to remove or to keep them, looking at
% the topography of ICs weights (bottom right), at the time course of the
% ICs (bottom left), at the power spectrum of the ICs (center right), at
% the single-trial time course of the ICs (center left), at the original
% average TMS-evoked potentials (top left), at the average TMS-evoked
% potentials after subtraction of the ICs that represent artefacts
ICA_selectandremove(EEG2ICA)
% click on button LOAD


%% Update the data after choosing the bad ICs

EEG.ICA2rem1 = length(data_GUI.comp2remove);

EEG2ICA.data = data_GUI.compproj;


%% Plot eeg layout
EEGLAB_plot_EEG_wrapper(EEG2ICA, [-1200  1200], 20, 1:size(EEG2ICA.data,1))

clear data

data = squeeze(mean(EEG2ICA.data,3));
figure()
plot(data(:,930:1701)')

%Actualize eeg data
if ~isempty(data_GUI.comp2remove)
    EEG.data = EEG2ICA.data(:, 1:2400, :);
end


%% Doing the final band-pass filtering

EEG = pop_tesa_filtbutter( EEG, 2, 50, 4, 'bandpass' );

clear data

data = squeeze(mean(EEG.data,3));
figure()
plot(data(:,930:1701)')

%% Trial rejection 2
TMPREJ=[];
eegplot(EEG.data,'winlength',5,'command','pippo','srate',EEG.srate,'limits',[EEG.times(1) EEG.times(end)]);
'data have been displayed'


%% Select a shorter epoch to remove possible artifacts
 
 EEG = pop_resample( EEG, 1000);

 EEG = pop_select(EEG, 'time', [-1 1]);


 % Visualize the data:
figure; pop_timtopo(EEG, [-600  600],  [], 'ERP data and scalp maps of  resampled');



%% If additionaly finally you want to eliminate and interpolate a specific electrode you can run:
%(Take into account that during ICA bad channel have been eliminated)

%{
1 option: manual rejection. To'_' add the label of the cannel (ex.'FP9')youwant to delate.
EEG.allchan = EEG.chanlocs;
EEG = pop_select( EEG,'nochannel', {'C3' 'C5'});
%}

%{
1.2 option: automatic rejection
EEG.allchan = EEG.chanlocs;
EEG = pop_rejchan(EEG, 'elec',[1:size(EEG.data,1)] ,'threshold',5,'norm','on','measure','kurt');
%}

%{
2. interpolate msising channels: 
EEG = pop_interp(EEG, EEG.allchan, 'spherical');
%}

%{
3. Rereference to the average: 
EEG = pop_reref( EEG, []);
%}

%4.Visualize the data
%EEGLAB_plot_EEG_wrapper(EEG, [-1000  1000], 20, 1:size(EEG.data))

%% Save the current form of dataset

EEG = pop_saveset( EEG, 'filename',[ dataset,'_',condition,'_Cleaned.set'],'filepath',[save_path,'\BraiNMAG_05_',subject,'\']);

close all

%% FIN