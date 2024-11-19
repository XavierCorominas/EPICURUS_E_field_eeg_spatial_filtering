%% Analyses and resutls of spatially vs non-spatially fitlered data.

% Xavier COROMINAS, 2024.

%% Start assigning paths

clear all
close all
clc

%
FieldTrip_path = ['D:\SOFTWARE\fieldtrip-20230118'];
addpath(FieldTrip_path);
ft_defaults;
%
External_functions = ['E:\LI_rTMS_project\eeg_analysis\Analysis\External_functions'];
addpath(External_functions);

%% LOAD UNFILTERED DATA

data_path = ['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\Merged_all_org\'];
save_path =['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\Analysis\Data_org\TEP'];


% Load sensor  DATA and filter it from 2 to 45Hz
for nsub = [1, 4, 5, 6, 7, 8, 10, 11, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24];

    subject = ['00', num2str(nsub)]

site = 'hi';
stim = 'sp';
cond = 'real';
target = 'M1';
     cfg = []; 
     cfg.lpfilter        = 'yes';
     cfg.lpfreq          = 45;
     cfg.dataset = [data_path, subject,'_',site,'_', stim,'_','M1_',cond,'.set'];
     hi_sp_real{nsub} = ft_preprocessing(cfg);

end 



%%  Load filtered data

data_path = ['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\Merged_all_deflect\'];
save_path =['E:\LI_rTMS_project\eeg_analysis\Analyzed_data\Analysis\Data_deflect\TEP'];


for nsub = [1, 4, 5, 6, 7, 8, 10, 11, 13, 14, 15, 16, 17, 18, 21, 22, 23, 24];

    subject = ['00', num2str(nsub)]
    site = 'hi';
    stim = 'sp';
    cond = 'real';
    target = 'M1';
     cfg = []; 
     cfg.lpfilter        = 'yes';
     cfg.lpfreq          = 45;
     cfg.dataset = [data_path, subject,'_',site,'_', stim,'_','M1_',cond,'.set'];
     hi_sp_real_filtered{nsub} = ft_preprocessing(cfg);

end 


% Loop through the outer cell array
for i = 1:length(hi_sp_real_filtered)

    % Check if the current cell is non-empty and contains the 'trial' field
    if ~isempty(hi_sp_real_filtered{1,i})

        % Loop through the inner cell array
        for j = 1:numel(hi_sp_real_filtered{1, i}.trial)

            % Check if the current inner cell is not empty
            if ~isempty(hi_sp_real_filtered{1, i}.trial{j})

                % Multiply the value inside the inner cell by 10
                hi_sp_real_filtered{1, i}.trial{j}(7,:) = hi_sp_real_filtered{1, i}.trial{j}(7,:) * 0.001;

            end
        end
    end
end


%% TEPs unfiltered data

% cluster of electrodes around stimualted left M1 region. 
row_unique =[7 ; 40;  39 ; 42 ; 38 ]; % correspond tot he electrodes under the stimualted region in left motor cortex C1,C3,C5,CP3,FC3

%
for nsub = [001, 004, 005, 006, 007, 008, 0010, 0011, 0013, 0014, 0015, 0016, 0017, 0018, 0021, 0022,0023, 0024];
cfg = [];
cfg.channel = {'C1','C3','CP1','CP3'} ; %C1,C3,C5,CP3,FC3
cfg.latency  = [-0.5 0.6];
cfg.keeptrials = 'no'; 
tl_hi_sp_real{nsub} = ft_timelockanalysis(cfg,hi_sp_real{nsub});
end 

% baseline correction
for nsub = [001, 004, 005, 006, 007, 008, 0010, 0011, 0013, 0014, 0015, 0016, 0017, 0018, 0021, 0022, 0023, 0024];

%cfg = [];
%cfg.baseline = [-0.5 0.6];
%cfg.baselinetype = 'zscore'; % 'absolute'
%tl_hi_sp_real{1,nsub} = ft_timelockbaseline(cfg, tl_hi_sp_real{1,nsub});

cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'absolute'; % 'absolute'
bc2_tl_hi_sp_real{1,nsub} = ft_timelockbaseline(cfg, tl_hi_sp_real{1,nsub});
end 



% GAV
nsub = [001, 005, 006, 007, 008, 0010, 0011, 0013, 0014, 0015, 0016, 0017, 0018, 0021, 0022,0023, 0024];
cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
gav_tl_hi_sp_real       = ft_timelockgrandaverage(cfg, bc2_tl_hi_sp_real{nsub});  

%  mean of all individuals with rows of intereset according to channels 
copy2(:,1,:) = mean(gav_tl_hi_sp_real.individual,2);
gav_tl_hi_sp_real.individual = copy2;



%% TEP filtered

% timelock of original data
for nsub = [001, 004, 005, 006, 007, 008, 0010, 0011, 0013, 0014, 0015, 0016, 0017, 0018, 0021, 0022,0023, 0024];

cfg = [];
cfg.channel = 'C3';
cfg.latency  = [-0.5 0.6]
cfg.keeptrials = 'no';  
tl_hi_sp_real_filtered{nsub} = ft_timelockanalysis(cfg,hi_sp_real_filtered{nsub});

end 

%baseline correction
% zscore
for nsub = [001, 004, 005, 006, 007, 008, 0010, 0011, 0013, 0014, 0015, 0016, 0017, 0018, 0021, 0022, 0023, 0024];
%cfg = [];
%cfg.baseline = [-0.5 0.6];
%cfg.baselinetype = 'zscore'; % 'absolute'
%tl_hi_sp_real_filtered{1,nsub} = ft_timelockbaseline(cfg, tl_hi_sp_real_filtered{1,nsub});

cfg = [];
cfg.baseline = [-0.5 -0.2];
cfg.baselinetype = 'absolute'; % 'absolute'
bc2_tl_hi_sp_real_filtered{1,nsub} = ft_timelockbaseline(cfg, tl_hi_sp_real_filtered{1,nsub});
end 

% GAV
nsub = [001, 005, 006, 007, 008, 0010, 0011, 0013, 0014, 0015, 0016, 0017, 0018, 0021, 0022,0023, 0024];
cfg = [];
cfg.channel   = 'C3';
cfg.latency   = 'all';
cfg.parameter = 'avg';
cfg.keepindividual = 'yes';
gav_tl_hi_sp_real_filtered       = ft_timelockgrandaverage(cfg, bc2_tl_hi_sp_real_filtered{nsub});  

%
gav_tl_hi_sp_real.label =gav_tl_hi_sp_real_filtered.label
gav_tl_hi_sp_real.elec =  gav_tl_hi_sp_real_filtered.elec


%% STAT hi sp real vs filtered 

cfg = [];
cfg.latency          = [0 0.6];
cfg.parameter        = 'individual';
cfg.channel          = {'C3'};
%cfg.latency          = [0 0.3];
cfg.method           = 'montecarlo';
cfg.statistic        = 'depsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.025;
cfg.clusterstatistic = 'maxsum';
%cfg.minnbchan        = 1;
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.correcttail         = 'alpha';
cfg.alpha               = 0.05;
cfg.numrandomization =  5000 

design = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17; 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 ]
cfg.design = design;
cfg.uvar   = 1;
cfg.ivar   = 2;
% Real vs Sham for rith TMS
stat_hi_sp_real_sham = ft_timelockstatistics(cfg, gav_tl_hi_sp_real_filtered, gav_tl_hi_sp_real);

% correct pvalues of positive clusters
test = struct2cell(stat_hi_sp_real_sham.posclusters);
test = squeeze(test(1,1,:));
test = cell2mat(test);

for i = 1:length(test)
    if test(i,1) < 0.5
      test(i,1) = test(i,1)*2; 
    elseif test(i,1) > 0.5 
      test(i,1) = 2*(1-test(i,1));
    end 
end 
for i= 1:length(test)
    test(i,2) = i;
end 
stat_hi_sp_real_sham.pvaluesposclusterscorrected = test;


% correct pvalues of negative clusters
test_neg = struct2cell(stat_hi_sp_real_sham.negclusters);
test_neg = squeeze(test_neg(1,1,:));
test_neg = cell2mat(test_neg);

for i = 1:length(test_neg)
    if test_neg(i,1) < 0.5
      test_neg(i,1) = test_neg(i,1)*2; 
    elseif test_neg(i,1) > 0.5 
      test_neg(i,1) = 2*(1-test_neg(i,1));
    end 
end 
for i= 1:length(test_neg)
    test_neg(i,2) = i;
end 
stat_hi_sp_real_sham.pvaluesnegclusterscorrected = test_neg;
%

%% Time series plot filtred vs unfiltered

mean_hi_sp_real = squeeze(nanmean(gav_tl_hi_sp_real.individual(:,1,:),1))';
std_hi_sp_real = squeeze(std(gav_tl_hi_sp_real.individual(:,1,:),0,1))';
error_hi_sp_real =  (std_hi_sp_real./sqrt(size(gav_tl_hi_sp_real.individual(:,1,:),1))).*1.96;

mean_hi_sp_sham = squeeze(nanmean(gav_tl_hi_sp_real_filtered.individual(:,1,:),1))';
std_hi_sp_sham= squeeze(std(gav_tl_hi_sp_real_filtered.individual(:,1,:),0,1))';
error_hi_sp_sham =  (std_hi_sp_sham./sqrt(size(gav_tl_hi_sp_real_filtered.individual(:,1,:),1))).*1.96;


% PLOT sp conditions
x = (1:1101);
figure() 
% Set the font for the entire figure
set(gca, 'FontName', 'Times New Roman');  % Set axes text to Times New Roman
set(gcf, 'DefaultTextFontName', 'Times New Roman');  % Set general text to Times New Roman

h3 = shadedErrorBar(x, mean_hi_sp_sham, error_hi_sp_sham);
h3.mainLine.LineWidth = 2;
h3.mainLine.Color = [1 0 0]; % red
h3.patch.FaceColor = h3.mainLine.Color;

h4 = shadedErrorBar(x, mean_hi_sp_real, error_hi_sp_real);
h4.mainLine.LineWidth = 2;
h4.mainLine.Color = [0 0 1]; % blue
h4.patch.FaceColor = h4.mainLine.Color;

% Assign x label, y label, and other plot elements
ylim([-3 3]);
xlim([400 1101]);
xticks([400 500 600 700 800 900 1000 1100]);
xticklabels({'-100','0','100','200','300','400','500','600'});
yticks([-3 -2 -1 0 1 2 3]);
yticklabels({'-3','-2','-1','0','1','2','3'});
ylabel('Amplitude (µV)', 'FontSize', 20,'FontName', 'Times New Roman');
xlabel('Time (ms)', 'FontSize', 20,'FontName', 'Times New Roman');
title([]);

% Customize legend
legend({'Filtered', 'Unfiltered'}, 'EdgeColor', 'none', 'FontSize', 25,'FontName', 'Times New Roman');

% Increase tick label font size for both axes
set(gca, 'FontSize', 20);

% Plot additional elements and annotations
line([500 500], get(gca, 'ylim'), 'Linestyle', '-', 'Color', 'k', 'DisplayName', 'TMSp1', 'linewidth', 1, 'HandleVisibility', 'off');
yline(0, '-.', 'LineWidth', 1.5, 'HandleVisibility', 'off');

% Plot rectangles for clusters
pos1 = stat_hi_sp_real_sham.posclusterslabelmat == 1;
cluster_pos = find(pos1);
min_cluster_pos = min(cluster_pos);
max_cluster_pos = max(cluster_pos);
distance_pos = max_cluster_pos - min_cluster_pos;
r = rectangle('Position', [min_cluster_pos + 500, -3, distance_pos, 6], 'FaceColor', [0.5 0.5 0.5 0.4], 'EdgeColor', 'none');

neg1 = stat_hi_sp_real_sham.negclusterslabelmat == 1;
cluster_neg = find(neg1);
min_cluster_neg = min(cluster_neg);
max_cluster_neg = max(cluster_neg);
distance_neg = max_cluster_neg - min_cluster_neg;
r = rectangle('Position', [min_cluster_neg + 500, -3, distance_neg, 6], 'FaceColor', [0.5 0.5 0.5 0.4], 'EdgeColor', 'none');




%% ************ FIGURES

%% PLOTS FOR UNFILTERED DATA 

close all 
clear all
clc;

% Data
cd('E:\LI_rTMS_project\eeg_analysis\Analyzed_data\Analysis\Data_org\TEP');
load("gav_teps_plots_supl_figure_hi_sp_real.mat");

% MNI152 mesh and leadfields
cd ('D:\MRI\MNI152\leadfields')
load('cortex_surf.mat');
load('LFM_Edir_ave.mat');
load('LFM_Exyz_ave.mat');
cd ('D:\MRI\MNI152\leadfields')


%% Visualize
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.showcomment   =  'no';
cfg.ylim = [-3 3];
cfg.xlim = [-0.05 0.4]
figure; ft_multiplotER(cfg, gav_tl_hi_sp_real)

%% Generate sensor level topographies

close all
for i=[0.015, 0.040, 0.1, 0.175, 0.230, 0.300]
cfg = [];

cfg.xlim = [i i];
cfg.channel = {'all'};
cfg.zlim = [-1.5 1.5];
cfg.layout = 'biosemi64.lay';
cfg.comment   =  'no';
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    cmap = (flipud(brewermap(64,'RdBu')))
    cfg.colormap = cmap;
    cfg.colormap            = cmap;

ft_topoplotER(cfg,gav_tl_hi_sp_real); 
%c = colorbar
%c.FontSize = 15
%c.Label.String ='Amplitude (uV)','FontSize',20;
set(gca, 'FontName', 'Times New Roman');  % Set axes text to Times New Roman
set(gcf, 'DefaultTextFontName', 'Times New Roman');  % Set general text to Times New Roman
%title([num2str(i*1000), 'ms'],'FontSize',20)
end



%% Souce based topographies. Classic MNE with cortically constrained dipoles ( SVD 15 components)

close all
for i = [15,40,100,175,230,300]
    i = i;
    EEG_signal_org = gav_tl_hi_sp_real.avg(:,(500+i)); % slect dataset of interest
    EEG_signal_org = mean(EEG_signal_org,2);

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

% Normalize 
MNE = normalize(MNE, 'range', [-1.5 1.5]);

figure;
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    cmap = (flipud(brewermap(64,'RdBu')))
    cfg.colormap = cmap;
    cfg.colormap            = cmap

h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'FaceAlpha',1,'EdgeAlpha',0);
view(0,90)
axis off;
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
%title([num2str(i), 'ms'])
%c =colorbar
caxis([-1.5 1.5])
colormap(cmap)

end 


%%  PLOTS FOR FILTERED DATA

cd('E:\LI_rTMS_project\eeg_analysis\Analyzed_data\Analysis\Data_org\TEP');
load("gav_tep_ohb_stats_filtered.mat");

EEG_filtered = gav_tl_hi_sp_real;
EEG_filtered.individual = squeeze(mean(EEG_filtered.individual));
EEG_filtered.individual = EEG_filtered.individual';
%

%
cd('E:\LI_rTMS_project\eeg_analysis\Analyzed_data\Analysis\Data_org\TEP');
load("gav_teps_plots_supl_figure_hi_sp_real.mat");
gav_tl_hi_sp_real.avg(7,:) = EEG_filtered.individual;


% Visualize

cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.showcomment   =  'no';
cfg.ylim = [-3 3];
cfg.xlim = [-0.05 0.4]
figure; ft_multiplotER(cfg, gav_tl_hi_sp_real)

%% Normalize data filtered

gav_tl_hi_sp_real.avg(1:6,:) = 0;
gav_tl_hi_sp_real.avg(8:63,:) = 0;

MNE = gav_tl_hi_sp_real.avg(7,:);
MNE = normalize(MNE, 'range', [-1.5 1.5]);
gav_tl_hi_sp_real.avg(7,:) = MNE;


% Visualize
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.showcomment   =  'no';
cfg.ylim = [-3 3];
cfg.xlim = [-0.05 0.4]
figure; ft_multiplotER(cfg, gav_tl_hi_sp_real)

%% Simulate chanel level topographies

close all
for i=[0.020, 0.040, 0.1, 0.175, 0.230, 0.300]
cfg = [];
cfg.xlim = [i i];
cfg.channel = {'all'};
cfg.zlim = [-1.5 1.5];
cfg.layout = 'biosemi64.lay';
cfg.comment   =  'no';
    ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
    cmap = (flipud(brewermap(64,'RdBu')))
    cfg.colormap = cmap;
    cfg.colormap            = cmap;

ft_topoplotER(cfg,gav_tl_hi_sp_real); 
set(gca, 'FontName', 'Times New Roman');  % Set axes text to Times New Roman
set(gcf, 'DefaultTextFontName', 'Times New Roman');  % Set general text to Times New Roman
%title([num2str(i*1000), 'ms'],'FontSize',20)
end


%% Simulate source level topographies

close all

for i = [20,40,100,175,230,300]

    EEG_signal_org = gav_tl_hi_sp_real.avg(:,(500+i)); % slect dataset of interest
    
% Calculate MNE
LL = LFM_Edir_ave*LFM_Edir_ave' ;
trL = trace(LL);
lambda = 0.1;
MNE = LFM_Edir_ave'*((LL + lambda*trace(LL)*eye(63))\EEG_signal_org);


figure;
set(gca, 'FontName', 'Times New Roman');  % Set axes text to Times New Roman
set(gcf, 'DefaultTextFontName', 'Times New Roman');  % Set general text to Times New Roman

h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'EdgeAlpha',0,'FaceAlpha',1);
    cmap = (flipud(brewermap(64,'RdBu')))
    cfg.colormap = cmap;
    cfg.colormap            = cmap
view(0,75)
axis off;
%colorbar
view(0,90)
axis off;
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
caxis([-1.5 1.5])
colormap(cmap)
end 

%%


%%