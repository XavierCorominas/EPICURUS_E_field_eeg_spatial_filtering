
% Plots TEP suplementary figure paper. suplementari figure 2

%%
close all 
clear all
clc;


cd('F:\LI_rTMS_project\eeg_analysis\Analyzed_data\Analysis\Data_org\TEP');
load("gav_teps_plots_supl_figure_hi_sp_real.mat");
cd ('D:\MRI\MNI152\leadfields')
load('cortex_surf.mat');
load('LFM_Edir_ave.mat');
load('LFM_Exyz_ave.mat');
cd ('D:\MRI\MNI152\leadfields')
load('CustomColormap.mat');
load('CustomColormap1.mat');




%%
cfg = [];
cfg.layout = 'biosemi64.lay';
cfg.showcomment   =  'no';
cfg.ylim = [-3 3];
cfg.xlim = [-0.05 0.4]
figure; ft_multiplotER(cfg, gav_tl_hi_sp_real)

%%

close all
for i=[0.015, 0.045, 0.1, 0.175, 0.230, 0.300]
cfg = [];
cfg.xlim = [i i];
cfg.channel = {'all'};
cfg.zlim = [-1.5 1.5];
cfg.layout = 'biosemi64.lay';
cfg.comment   =  'no';
ft_topoplotER(cfg,gav_tl_hi_sp_real); 
c = colorbar
c.FontSize = 15
c.Label.String ='Amplitude (uV)','FontSize',20;
title([num2str(i*1000), 'ms'],'FontSize',20)
end



%%
close all

 for i = [20,40,100,175,230,300]
%for i = [100]
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

%plot

% % extract only tails of the data to be plotted
% pct10 = prctile(MNE,10);
% pct90 = prctile(MNE,90);
% MNE_C = MNE;
% idx = find(MNE > pct10 & MNE < pct90);
% MNE(idx) = 0;
% normalize 
MNE = normalize(MNE, 'range', [-1.5 1.5]);



figure;
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
c =colorbar
caxis([-1.5 1.5])
colormap parula
%colormap(CustomColormap1)
end 


%% 
%Do the same for the spatially filtered data


close all

for i = [20,40,100,175,230,300]
    i = i-5;
    EEG_signal_org = gav_tl_hi_sp_real.avg(:,(500+i):(500+i+10)); % slect dataset of interest

    EEG_signal_org = mean(EEG_signal_org,2);
    EEG_signal_org(2:63,1)=0;
    EEG_signal_org(7,1)= EEG_signal_org(1,1);
   EEG_signal_org(1,1) = 0;

   % Calculate MNE
LL = LFM_Edir_ave*LFM_Edir_ave' ;
trL = trace(LL);
lambda = 0.1;
MNE = LFM_Edir_ave'*((LL + lambda*trace(LL)*eye(63))\EEG_signal_org);

%{
% Another way of computing the classic MNE without assuming a lambda value
% but assuming SNR

SNR = 100;
lambda = trace(LFM_Edir_ave*LFM_Edir_ave')/SNR;
G = LFM_Edir_ave'/(LFM_Edir_ave*LFM_Edir_ave' + lambda*eye(63));
MNE =G*EEG_signal_org;
%}

% extract only tails of the data to be plotted
%  pct10 = prctile(MNE,10);
%  pct90 = prctile(MNE,90);
%  idx = find(MNE > pct10 & MNE < pct90);
%  MNE(idx) = 0;


figure;
h = trisurf(cortex_surf.e, cortex_surf.p(:,1),  cortex_surf.p(:,2),...
   cortex_surf.p(:,3),MNE,'EdgeAlpha',0,'FaceAlpha',1);
view(0,75)
axis off;
colorbar
view(0,90)
axis off;
lightangle(180,-60)
h.FaceLighting = 'gouraud';
h.AmbientStrength = 0.8;
h.DiffuseStrength = 0.8;
h.SpecularStrength = 0.9;
h.SpecularExponent = 25;
h.BackFaceLighting = 'unlit';
title([num2str(i+5), 'ms'])
%colorbar
caxis([-1.5 1.5])
%colormap(CustomColormap)
end 



%% 
%Tranform data to vols in case of the spatially fitlered data;

gav_tl_hi_sp_real.individual = gav_tl_hi_sp_real.individual*0.001;
gav_tl_hi_sp_sham.individual = gav_tl_hi_sp_sham.individual*0.001;


mean_hi_sp_real = squeeze(nanmean(gav_tl_hi_sp_real.individual(:,1,:),1))';
std_hi_sp_real = squeeze(std(gav_tl_hi_sp_real.individual(:,1,:),0,1))';
error_hi_sp_real =  (std_hi_sp_real./sqrt(size(gav_tl_hi_sp_real.individual(:,1,:),1))).*1.96;

mean_hi_sp_sham = squeeze(nanmean(gav_tl_hi_sp_sham.individual(:,1,:),1))';
std_hi_sp_sham= squeeze(std(gav_tl_hi_sp_sham.individual(:,1,:),0,1))';
error_hi_sp_sham =  (std_hi_sp_sham./sqrt(size(gav_tl_hi_sp_sham.individual(:,1,:),1))).*1.96;


% PLOT sp conditions
x = (1:1101);
figure() 
h3 = shadedErrorBar (x, mean_hi_sp_real, error_hi_sp_real)
h3.mainLine.LineWidth = 2;
h3.mainLine.Color = [1.0000    0.8398         0]; % gold
h3.patch.FaceColor = h3.mainLine.Color;

h4 = shadedErrorBar (x, mean_hi_sp_sham, error_hi_sp_sham)
h4.mainLine.LineWidth = 2;
h4.mainLine.Color = [0.7383    0.7148    0.4180]; % kaki
h4.patch.FaceColor = h4.mainLine.Color;

% asing x label
ylim([-3 3])
xlim([400 1101])
xticks([400 450 500 550 600 650 700 750 800 850 900 950 1000 1050 1100])
xticklabels({'-100','-50','0','50','100','150','200','250','300','350','400','450','500','550','600'})
yticks([-3 -2.5 -2 -1.5 -1 -0.5 0 0.5 1 1.5 2 2.5 3])
yticklabels({'-3','-2.5','-2','-1.5','-1','-0.5','0','0.5','1','1.5','2','2.5','3'})
%asign title and make tms pulse lines
%ylabel('Amplitude (uV)','FontSize',14);
%xlabel('Time (ms)','FontSize',14);
title('Single Pulse','FontSize',20);
legend({'Hi-Sp' 'Sham Hi-Sp'}, 'EdgeColor','none');
fontsize(legend,14,'points')
line([500 500],get(gca,'ylim'),'Linestyle','-','Color','k','DisplayName','TMSp1','linewidth',1,'HandleVisibility','off')
yline(0,'-.','LineWidth',1.5,'HandleVisibility','off')

%%

test_data = ans.gav_tl_hi_sp_real;
test_data.avg = test_data.avg*0;
test_data.var = test_data.var*0;

data1 = normalize(gav_tl_hi_sp_real.avg, 'range', [-1.5 1.5]);
test_data.avg(7,:)= data1;

data2 = normalize(gav_tl_hi_sp_real.var, 'range', [-1.5 1.5]);
test_data.var(7,:)= data2;


close all
for i=[0.02, 0.04, 0.1, 0.175, 0.230, 0.300]
cfg = [];
cfg.xlim = [i i];
cfg.channel = {'all'};
cfg.zlim = [-1.5 1.5];
cfg.layout = 'biosemi64.lay';
cfg.comment   =  'no';
ft_topoplotER(cfg,test_data); 
c = colorbar
c.FontSize = 15
c.Label.String ='Amplitude (uV)','FontSize',20;
title([num2str(i*1000), 'ms'],'FontSize',20)
end


%%