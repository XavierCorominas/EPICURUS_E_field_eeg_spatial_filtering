function [inds] = select_sources_from_surface(mesh,area1,area2,inds)

% This functions allows to click and select sources from the cortical mesh
% surface

% Tuomas Mutanen, 2019

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% input:
%
% sa = the head model
% radius1 = the radius used to picking up the sources of interest
% radius2 = optional, the radius used to deleting the extra sources of interest
% ind = optional, already previously defined indices to enable to continue
% from the previously saved point
% reference_points = Anatomical landmarks in MNI coordinates to help to identify ROIs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 3
    area2 = area1;
end

if nargin < 4
    inds = [];
end

if nargin < 5
    reference_points = [];
end

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

    fig = figure(999)
    set(fig,'Position',[0 0 1600 800]);

    
    plot_vec = zeros(1,N_cortex_points);
    plot_vec(inds) = 1;
    
    h = trisurf(cortex.tri, cortex.vc(:,1), cortex.vc(:,2), cortex.vc(:,3),plot_vec,'EdgeAlpha',0);
    grid off;
    axis off;
    caxis([-1 1])
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
    plot3(tri_midpoint(:,1),tri_midpoint(:,2),tri_midpoint(:,3),'k.','MarkerSize',1)
    hold off;

continue_choosing = 1;

while continue_choosing

    figure(999);
    plot_vec(inds) = 1;
    h = trisurf(cortex.tri, cortex.vc(:,1), cortex.vc(:,2), cortex.vc(:,3),plot_vec,'EdgeAlpha',0);
    grid off;
    axis off;
    caxis([-1 1])
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
    plot3(tri_midpoint(:,1),tri_midpoint(:,2),tri_midpoint(:,3),'k.','MarkerSize',1)
    hold off;
    
    datacursormode on;

    title('Choose a point reflecting ROI, press a KEY to SAVE the point')
    pause();
    dcm_obj = datacursormode(fig);
    c_info = getCursorInfo(dcm_obj);
    %try
        %current_indx = find(sqrt(sum((tri_midpoint - repmat(c_info.Position, [N_cortex_points, 1])).^2,2)) < radius1);
        
        [current_indx] = return_cortical_patch(mesh, c_info.Position, area1, edgenb, tri_areas);
        inds = [inds; current_indx];


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
    plot_vec = zeros(1,N_cortex_points);
    plot_vec(inds) = 1;
    h = trisurf(cortex.tri, cortex.vc(:,1), cortex.vc(:,2), cortex.vc(:,3),plot_vec,'EdgeAlpha',0);
    grid off;
    axis off;
    caxis([-1 1])
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
    plot3(tri_midpoint(:,1),tri_midpoint(:,2),tri_midpoint(:,3),'k.','MarkerSize',1)
    hold off;
    %end 
    datacursormode on;

    title('Trim accidental extra points, press a KEY to DELETE the point')
    pause();
    dcm_obj = datacursormode(fig);
    c_info = getCursorInfo(dcm_obj);
  %  try
   %     current_indx = find(sqrt(sum((tri_midpoint - repmat(c_info.Position, [N_cortex_points, 1])).^2,2)) < radius2);
        [current_indx] = return_cortical_patch(mesh, c_info.Position, area2, edgenb, tri_areas);
        inds = setdiff(inds, current_indx);
   % catch
    %    break;
    %end

    title('Press any KEY to CONTINUE, CLICK to END deleting ROI points')
    continue_choosing = waitforbuttonpress;
    [az,el] = view; 
    
end

end

function midpoints=Tri_Midpoints(nodes,elements)
% function midpoints=TriangleMidpoints(nodes,elements)
% Calculates midpoints of the mesh triangles.

p1=nodes(elements(:,1),:);
p2=nodes(elements(:,2),:);
p3=nodes(elements(:,3),:);
midpoints=(p1+p2+p3)/3;
    
end

function [patch_indices] = return_cortical_patch(mesh, patch_midpoint, patch_surface_area, edgenb, tri_areas)

% This function returns a cortical surface patch centered at
% patch_midpoint


% Reading the cortical surface from the average head

cortex.vc = mesh.p;
cortex.tri = mesh.e;

tri_midpoint =  Tri_Midpoints(cortex.vc, cortex.tri);

% Look for the stargint triangle:

[~, start_tri] = min( sum((tri_midpoint - patch_midpoint).^2,2 ) );

patch_indices = start_tri;
current_tri = start_tri;
tot_A = tri_areas(start_tri);
k = 1;
while tot_A < patch_surface_area
    for i = 1:length(edgenb(current_tri,:))
        if find(patch_indices == edgenb(current_tri,i))
        else
            patch_indices = [patch_indices ; edgenb(current_tri,i)];
            
            tot_A = tot_A + tri_areas(edgenb(current_tri,i));
        end
    end
    current_tri = patch_indices(k+1);
    k = k + 1;
end

end



