function [LFM_correct] = correct_individual_sources_automatic(mesh, LFM, percentile_threshold)




% Generate new submatrix
source_sentivity = sum(LFM.^2,1);

% Find the 90th percentile value
threshold_value = prctile(source_sentivity(:), percentile_threshold);
outliers_index = source_sentivity > threshold_value;
inds = find(outliers_index);

%interpolate
edgenb=edgeneighbors(mesh.e);
LFM_correct = LFM;

for i = inds
        
        % Interpolate each bad source topography based on good neighbours
        good_neighbours = setdiff(edgenb(i,:),i);
        LFM_correct(:,i) = nanmedian(LFM(:,good_neighbours),2);

%            LFM_correct(:,i) = nanmedian(prctile(LFM,percentile_interpolation));
    
end

end
