function hp=PlotMeshData(mesh,data,varargin)
% function hp=PlotDataOnMesh(mesh,data,varargin)
% mesh: a struct with fields "p" for vertices and "e" for faces
% data: scalar data to be plotted on the mesh, [number-of-vertices x 1]
% varargin = options for and settings for the plot
%   - a struct with optional fields:
%       - 'colormap', <colormap array> or <number of colors for jet colormap>
%       - 'caxis', <colormap axis, either [min, max] or [absmax]>
%       - 'inflated', <use inflated mesh? 1 or 0 (needs mesh.pinf)>
%       - 'colorbar', <plot colorbar? 1 or 0>
%       - 'view',   <view angles in a format taken by view()-function>
%    - Or, if no struct given, parameter-value pairs as described above.
% If any struct is given as input, no other arguments are parsed.
%
% Copyright (c) 2011--2013 Matti Stenroos (matti.stenroos@aalto.fi)
%  -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
%                !! There is no warranty of any kind !!
% -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
% 17 Oct 2013

%if arguments are passed as varargin from another function, they need to be
%converted:
if iscell(varargin) && length(varargin)==1 && iscell(varargin{1})
    varargin=varargin{1};
end
plotoptions=ParseOptions(varargin);

%make sure data is given the correct way
if size(data,1)==1;
    data=data';
end
%check that the data is of correct size...
if size(data,1)~=size(mesh.p,1), % default . p
    error('PlotDataOnMesh: The sizes of mesh and datavector do not match.');
end
%OK, let's start...

%if plotoptions.inflated,p=mesh.pinf;else
    p=mesh.p;
%end

%this alloww good-quality isocolor surface visualization:
set(gcf,'renderer','zbuffer');if ~ishold,cla;end
%...make the plot
hp = patch('faces',mesh.e,'vertices',p,'facevertexcdata',data,'facecolor','interp','edgecolor','none');

%set colormap
if size(plotoptions.colormap,2)==3
    colormap(plotoptions.colormap)
else
    colormap(jet(plotoptions.colormap));
end
caxis(SetPlotScale(data,plotoptions.caxis));%set scale

if plotoptions.colorbar%add colorbar
    colorbar;
end
view(plotoptions.view);%set view

%and make everything look nice...
axis tight equal off;material dull;lighting gouraud;
if isempty(findobj(gca,'Type','light'));
    camlight
end

function plotoptions=ParseOptions(varargin)
if iscell(varargin) && length(varargin)==1 && iscell(varargin{1})
    varargin=varargin{1};
end
%defaults
plotdef.caxis=[];
plotdef.colorbar=1;
cmap=jet(20);
cmap=cmap(11:20,:);
plotdef.colormap=cmap;
plotdef.inflated=1;
plotdef.view=[-90 0];

%parse options
plotoptions=GetStruct(varargin);%try, if there is any struct...?
if ~isstruct(plotoptions)%if not, check the parameters...
    plotoptions=struct('kind','plotoptions');
    if IsParameter(varargin,'caxis'),plotoptions.caxis=GetValue(varargin,'caxis');end
    if IsParameter(varargin,'colorbar'),plotoptions.colorbar=GetValue(varargin,'colorbar');end
    if IsParameter(varargin,'colormap'),plotoptions.colormap=GetValue(varargin,'colormap');end
    if IsParameter(varargin,'inflated'),plotoptions.inflated=GetValue(varargin,'inflated');end
    if IsParameter(varargin,'view'),plotoptions.view=GetValue(varargin,'view');end
end
%fill missing option values with defaults.
if ~isfield(plotoptions,'caxis'),plotoptions.caxis=plotdef.caxis;end
if ~isfield(plotoptions,'colorbar'),plotoptions.colorbar=plotdef.colorbar;end
if ~isfield(plotoptions,'colormap'),plotoptions.colormap=plotdef.colormap;end
if ~isfield(plotoptions,'inflated'),plotoptions.inflated=plotdef.inflated;end
if ~isfield(plotoptions,'view'),plotoptions.view=plotdef.view;end

function cscale=SetPlotScale(data,cscale)
% function cscale=SetPlotScale(data,cscale)
% Compute scale for data to be plotted.
% cscale (optional): if empty or omitted, autoscale to the abs(max(plotdata))

if ~isempty(cscale) && length(cscale)==1,
    if all(data>=0)
        cscale=[0 1]*abs(cscale);
    elseif all(data<=0)
        cscale=[-1 0]*abs(cscale);
    else
        cscale=[-1 1]*abs(cscale);
    end
elseif all(data>=0)
    cscale=[0 1]*max(data);
elseif all(data<=0)
    cscale=[1 0]*min(data);
else
    cscale=[-1 1]*max(abs(data));
end

function [res,index]=GetStruct(paralist)
% function [res,index]=GetStruct(paralist)
% Checks, whether there is a struct in the input list and
% returns it and the index to the struct (or 0 and []). If there are many
% structs, this returns only the first one.
N=length(paralist);
res=0;index=[];
for I=1:N,
    if isstruct(paralist{I}),
        res=paralist{I};
        index=I;
        break
    end
end

function res=IsParameter(input,parameter)
% function res=IsParameter(input,parameter)
% Checks, whether a value exists in a parameter/value list
% input: parameter/value list; X1,Y1,X2,Y2,...
% parameter: the parameter to find
% res: TRUE or FALSE

% paralist=input(1:2:end);
paralist=input;
test=strcmpi(paralist,parameter);
res=any(test);
