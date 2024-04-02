function [] = EEGLAB_plot_EEG_wrapper(EEG, timeLim, amp_limit,  Chans, trial)

[~, Tind1] = min(abs(EEG.times - timeLim(1)))
[~, Tind2] = min(abs(EEG.times - timeLim(2)))

if nargin <5
    data_plot = squeeze(mean(EEG.data(Chans,Tind1:Tind2,:),3));
else
    data_plot = squeeze(EEG.data(Chans,Tind1:Tind2,trial));
end
    
figure,
plottopo(data_plot,EEG.chanlocs(Chans),0,[timeLim(1) timeLim(2) -amp_limit amp_limit],'',0,[.07 .07],{},1,[0])


end

