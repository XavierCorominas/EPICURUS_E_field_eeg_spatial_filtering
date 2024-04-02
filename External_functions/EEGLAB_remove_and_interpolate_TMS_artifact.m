function [EEG_out] = EEGLAB_remove_and_interpolate_recharging_artifact(EEG_in, x1, x2, stim)

% This function allows the user visually inspect the data and point the
% onset of the recharging artifact

[~, onset] = min(abs(EEG_in.times + 20));
[~, end_of_interval] = min(abs(EEG_in.times - 250));

evoked_response = mean(EEG_in.data,3);

x1 = find(EEG_in.times == x1);
x2 = find(EEG_in.times == x2);

y1 = (setdiff(1:size(EEG_in.data,2), x1:x2 ));

if stim == 1    
    for i = 1:size(EEG_in.data,1)
    disp(['Interpolating the TMS artifact... Channel ',num2str(i)])
   for k = 1:size(EEG_in.data,3)
       
        EEG_in.data(i,x1:x2,k) = spline(y1, EEG_in.data(i,y1,k), x1:x2);
        
         
   end
    end

else 
    
for x = 1:length(EEG_in.epoch);
TMS_lat(x,:) = cell2mat([EEG_in.epoch(x).eventlatency]);
end 

for i = 1:size(EEG_in.data,1)
    disp(['Interpolating the TMS artifact... Channel ',num2str(i)])
   for k = 1:size(EEG_in.data,3)
        
        x3 = find(EEG_in.times == round(EEG_in.times(x1) + (TMS_lat(k,2))));
        x4 = find(EEG_in.times == round(EEG_in.times(x2) + (TMS_lat(k,2))));

        x5 = find(EEG_in.times == round(EEG_in.times(x1) + (TMS_lat(k,3))));
        x6 = find(EEG_in.times == round(EEG_in.times(x2) + (TMS_lat(k,3))));
        
        x7 = find(EEG_in.times == round(EEG_in.times(x1) + (TMS_lat(k,4))));
        x8 = find(EEG_in.times == round(EEG_in.times(x2) + (TMS_lat(k,4))));
         
        y1 = (setdiff(1:size(EEG_in.data,2), x1:x2 ));
        y2 = (setdiff(1:size(EEG_in.data,2), x3:x4 ));
        y3 = (setdiff(1:size(EEG_in.data,2), x5:x6 ));
        y4 = (setdiff(1:size(EEG_in.data,2), x7:x8 ));
                
        EEG_in.data(i,x1:x2,k) = spline(y1, EEG_in.data(i,y1,k), x1:x2);
        EEG_in.data(i,x3:x4,k) = spline(y2, EEG_in.data(i,y2,k), x3:x4);
        EEG_in.data(i,x5:x6,k) = spline(y3, EEG_in.data(i,y3,k), x5:x6);
        EEG_in.data(i,x7:x8,k) = spline(y4, EEG_in.data(i,y4,k), x7:x8);
         
   end
end
end 

tmpfig = figure;

subplot(2,1,1)
evoked_response_clean = mean(EEG_in.data,3);
plot(evoked_response', 'r');
hold on;
plot(evoked_response_clean', 'k');

xlim([onset, end_of_interval]);

title('Date after interpolation. Press any key to continue.')
subplot(2,1,2)
plot(evoked_response'-evoked_response_clean', 'k');
xlim([onset, end_of_interval]);
title('The rejected recharging artifact.')
pause();

close(tmpfig);

EEG_out = EEG_in;

end
