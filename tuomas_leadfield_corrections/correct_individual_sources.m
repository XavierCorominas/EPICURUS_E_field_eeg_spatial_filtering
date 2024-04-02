function [LFM_correct] = correct_individual_sources(mesh, LFM)

% This functions allows to click and select sources from the cortical mesh
% surface

% Tuomas Mutanen, 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input:
%
% mesh = the cortical mesh
% LFM = the lead field matrix
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N_sources = size(mesh.e, 1);
LFM_correct = LFM;
inds = [];

if size(LFM,2) == N_sources
    source_sentivity = sum(LFM.^2,1);
else
    source_sentivity = sum(LFM(:,1:3:end).^2 + LFM(:,2:3:end).^2 + LFM(:,3:3:end).^2,1);

end
    
plot_vec = source_sentivity;

tri_areas = zeros(size(mesh.e, 1), 1);
%[conn,~,~]=meshconn(mesh.e,length(mesh.p));
edgenb=edgeneighbors(mesh.e);

for i = 1:length(mesh.e)
        A = mesh.p(mesh.e(i,1), :);
        B = mesh.p(mesh.e(i,2), :);
        C = mesh.p(mesh.e(i,3), :);
        
        ABC = [A; B; C]';
              
        tri_areas(i) = 0.5*sqrt( det( [ABC(1,:); ABC(2,:); 1 1 1] )^2 ...
        + det( [ABC(2,:); ABC(3,:); 1 1 1] )^2 ...
        + det( [ABC(3,:); ABC(1,:); 1 1 1] )^2 );

end

% Reading the cortical surface from the average head

cortex.vc = mesh.p;
cortex.tri = mesh.e;

tri_midpoint = Tri_Midpoints(cortex.vc, cortex.tri);

N_cortex_points = length(cortex.tri);

% Implementing the visualization and the GUI

    fig = figure(999);
    
    % Get the screen size
    screenSize = get(0, 'ScreenSize');
    
    % Define the position and size for the figure
    figurePosition = [screenSize(1), screenSize(4) * 0.25, screenSize(3) * 0.75, screenSize(4) * 0.75];
    
    % Set the figure's position and size
    set(gcf, 'Position', figurePosition);

    plot_vec(inds) = 1;
    
    h = trisurf(cortex.tri, cortex.vc(:,1), cortex.vc(:,2), cortex.vc(:,3),plot_vec,'EdgeAlpha',0);
    grid off;
    axis off;
   % caxis([-1 1])
    view([0, 90]);
    colormap('jet')
    colorbar;   
    lightangle(90,-90)
     h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.5;
    h.DiffuseStrength = 0.9;
    h.SpecularStrength = 0.2;
    h.SpecularExponent = 50;
    h.BackFaceLighting = 'unlit';
    [az,el] = view;
    
   hold on; 
    plot3(tri_midpoint(:,1),tri_midpoint(:,2),tri_midpoint(:,3),'k.','MarkerSize',20)
    hold off;

continue_choosing = 1;

while continue_choosing

    figure(999);
    plot_vec(inds) = 0;
    h = trisurf(cortex.tri, cortex.vc(:,1), cortex.vc(:,2), cortex.vc(:,3),plot_vec,'EdgeAlpha',0);
    grid off;
    axis off;
%    caxis([-1 1])
    view([0, 90]);
    colormap('jet')
    colorbar;   
    lightangle(90,-90)
     h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.5;
    h.DiffuseStrength = 0.9;
    h.SpecularStrength = 0.2;
    h.SpecularExponent = 50;
    h.BackFaceLighting = 'unlit';
    view(az,el)
    
    hold on; 
    plot3(tri_midpoint(:,1),tri_midpoint(:,2),tri_midpoint(:,3),'k.','MarkerSize',20)
    hold off;
    
    datacursormode on;

    title('Choose a point reflecting ROI, press a KEY to SAVE the point')
    pause();
    dcm_obj = datacursormode(fig);
    c_info = getCursorInfo(dcm_obj);
    %try
        %current_indx = find(sqrt(sum((tri_midpoint - repmat(c_info.Position, [N_cortex_points, 1])).^2,2)) < radius1);
        
        [current_indx] = return_triangle(mesh,c_info.Position);
        inds = [inds, current_indx];


    %catch
     %   break;
    %end

    title('Press any KEY to CONTINUE, CLICK to END selecting ROI points')
    continue_choosing = waitforbuttonpress;
    [az,el] = view;

end

continue_choosing = 1;

while continue_choosing

    figure(999);
    plot_vec = source_sentivity;
    plot_vec(inds) = 0;
    h = trisurf(cortex.tri, cortex.vc(:,1), cortex.vc(:,2), cortex.vc(:,3),plot_vec,'EdgeAlpha',0);
    grid off;
    axis off;
  %  caxis([-1 1])
    view([0, 90]);
    colormap('jet')
    colorbar;   
    lightangle(90,-90)
     h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.5;
    h.DiffuseStrength = 0.9;
    h.SpecularStrength = 0.2;
    h.SpecularExponent = 50;
    h.BackFaceLighting = 'unlit';
    view(az,el);
    
    %if ~isempty(reference_points)
    hold on; 
    plot3(tri_midpoint(:,1),tri_midpoint(:,2),tri_midpoint(:,3),'k.','MarkerSize',20)
    hold off;
    %end 
    datacursormode on;

    title('Trim accidental extra points, press a KEY to DELETE the point')
    pause();
    dcm_obj = datacursormode(fig);
    c_info = getCursorInfo(dcm_obj);
  %  try
   %     current_indx = find(sqrt(sum((tri_midpoint - repmat(c_info.Position, [N_cortex_points, 1])).^2,2)) < radius2);
        [current_indx] = return_triangle(mesh,c_info.Position);
        inds = setdiff(inds, current_indx);
   % catch
    %    break;
    %end

    title('Press any KEY to CONTINUE, CLICK to END deleting ROI points')
    continue_choosing = waitforbuttonpress;
    [az,el] = view; 
    
end

if size(LFM,2) == N_sources
    % The cortically constrained case:

    for i = inds
        
        % Interpolate each bad source topography based on good neighbours
        good_neighbours = setdiff(edgenb(i,:),inds);
        LFM_correct(:,i) = mean(LFM(:,good_neighbours),2);
    
    end


else
    % The free dipole case:

    for i = inds
        
        % Interpolate each bad source topography based on good neighbours
        good_neighbours = setdiff(edgenb(i,:),inds);
        LFM_correct(:,(i-1)*3+1) = mean(LFM(:,(good_neighbours-1)*3+1),2);
        LFM_correct(:,(i-1)*3+2) = mean(LFM(:,(good_neighbours-1)*3+2),2);
        LFM_correct(:,(i-1)*3+3) = mean(LFM(:,(good_neighbours-1)*3+3),2);
    end
end

    close(999)

end

function midpoints=Tri_Midpoints(nodes,elements)
% function midpoints=TriangleMidpoints(nodes,elements)
% Calculates midpoints of the mesh triangles.

p1=nodes(elements(:,1),:);
p2=nodes(elements(:,2),:);
p3=nodes(elements(:,3),:);
midpoints=(p1+p2+p3)/3;
    
end

function [patch_indices] = return_triangle(mesh,patch_midpoint)

% This function returns a cortical surface patch centered at
% patch_midpoint


% Reading the cortical surface from the average head

cortex.vc = mesh.p;
cortex.tri = mesh.e;

tri_midpoint =  Tri_Midpoints(cortex.vc, cortex.tri);

% Look for the stargint triangle:

[~, start_tri] = min( sum((tri_midpoint - patch_midpoint).^2,2 ) );

patch_indices = start_tri;


end



