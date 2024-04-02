
%%
clear all
close all
clc

%mainpath = 'E:\LI_rTMS_project\eeg_analysis\DATA\';
%cd(mainpath);

eeglab_path = ['D:\SOFTWARE\eeglab2022.0'];
addpath(eeglab_path);
eeglab;

%
External_functions = ['F:\BRAINMAG_Li-rTMS\ANALYSES\Analysis\External_functions'];
addpath(External_functions);

data_path = ['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\'];

save_path = ['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\Merged_all\'];
addpath(save_path)

%%  Merge Blocks SP,real,hi


% Subject preferneces
subject_id = 'XX';
% BLOCK NUMBER
block = '1';   % '1' , '2'
% STIMULATION intensity
stim ='30';  %  10 , 30 , 50 , 70 , 90 , 100
% TARGET
%target = 'M1';  
condition = 'real'; %sham real



dataset =[block,'_',stim,'_M1_',condition]; % modify this and delate block

% add number of subjects
for nsub = [1 2 3 4 5 6 7 8 10]; % 

    close all
      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned_deflect.set'],'filepath',['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_00',subject_id,'\deflect']);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned_deflect.set'],'filepath',['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_00',subject_id,'\deflect']);
EEG1 = pop_select(EEG1, 'time', [-0.99 0.99]);
EEG2 = pop_select(EEG2, 'time', [-0.99 0.99]);
EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[subject,'_',dataset '.set'], 'filepath',[save_path]);
pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%%



% Subject preferneces
subject_id = 'XX';
% BLOCK NUMBER
block = '1';   % '1' , '2'
% STIMULATION intensity
stim ='50';  %  10 , 30 , 50 , 70 , 90 , 100
% TARGET
%target = 'M1';  
condition = 'real'; %sham real


dataset =[block,'_',stim,'_M1_',condition];
% % add number of subjects
% for nsub = [1]; % 
% 
%     close all
%       if nsub <= 100
 %subject = ['00', num2str(nsub)]
%       end 
% 
%         for i_site = 1%:size(site)
% 
%             for i_stim = 1%:size(stim)
% 
%                 i_cond = 1
% 
%                 dataset= [block{i_site},'_',stim{i_stim},'_M1_',condition{i_cond}]

     
EEG1 = pop_loadset('filename',[dataset,'_Cleaned_deflect.set'],'filepath',['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_00',subject_id,'\deflect']);
%EEG2 = pop_loadset('filename',[dataset,'_Cleaned_deflect.set'],'filepath',['F:\BRAINMAG_Li-rTMS\ANALYSES\Analyzed_data\EEG\Preprocessed\BraiNMAG_05_00',subject_id,'\deflect']);
EEG1 = pop_select(EEG1, 'time', [-0.99 0.99]);
%EEG2 = pop_select(EEG2, 'time', [-0.99 0.99]);
%EEG = pop_mergeset(EEG1);
EEG = EEG1;
subject = '001'

EEG = pop_saveset( EEG, 'filename',[subject,'_',dataset '.set'], 'filepath',[save_path]);
pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');

%




%% Merge blocks Rith,real,hi


site = {'hi'};
stim = {'rith'};
cond = {'real'};

target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 
    close all
      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG1 = pop_select(EEG1, 'time', [-0.99 0.99]);
EEG2 = pop_select(EEG2, 'time', [-0.99 0.99]);
EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ subject,'_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%% Merge blocks Arrith,real,hi 

cond = {'real'};
site = {'hi'};
stim = {'arith'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%% %%  Merge Blocks SP,real,li

cond = {'real'};
site = {'li'};
stim = {'sp'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%% Merge blocks Rith,real,li

cond = {'real'};
site = {'li'};
stim = {'rith'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%% Merge blocks Arrith,real,li

cond = {'real'};
site = {'li'};
stim = {'arith'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 



%%  Merge Blocks SP,sham,hi

cond = {'sham'};
site = {'hi'};
stim = {'sp'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%% Merge blocks Rith,sham,hi

cond = {'sham'};
site = {'hi'};
stim = {'rith'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 
%% Merge blocks Arrith,sham,hi 

cond = {'sham'};
site = {'hi'};
stim = {'arith'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 
%% %%  Merge Blocks SP,sham,li

cond = {'sham'};
site = {'li'};
stim = {'sp'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 
%% Merge blocks Rith,sham,li

cond = {'sham'};
site = {'li'};
stim = {['rith']};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%% Merge blocks Arrith,sham,li


cond = {'sham'};
site = {'li'};
stim = {'arith'};
target = 'M1';

% add number of subjects
for nsub = [1 4 5 6 7 8 10 11 13 14 15 16 17 18 21 22 23]; % 

      if nsub <= 100
subject = ['00', num2str(nsub)]
      end 

        for i_site = 1%:size(site)

            for i_stim = 1%:size(stim)

                i_cond = 1

                dataset= [site{i_site},'_' stim{i_stim},'_' target,'_' cond{i_cond}]

     
EEG1 = pop_loadset('filename',['1_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);
EEG2 = pop_loadset('filename',['2_' dataset,'_Cleaned.set'],'filepath',['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\LITMS_01_',subject]);

EEG = pop_mergeset(EEG1,EEG2);



EEG = pop_saveset( EEG, 'filename',[ 'all_',dataset '.set'], 'filepath',[save_path]);

pop_tesa_plot(EEG, 'xlim', [-200,600], 'ylim', [-25,25])
figure; pop_timtopo(EEG, [-100  600], [20 50 70 100 160 220 290 400], 'ERP data and scalp maps of  resampled');


clear EEG1
clear EEG2
clear EEG

    end
    end 
end 

%%




