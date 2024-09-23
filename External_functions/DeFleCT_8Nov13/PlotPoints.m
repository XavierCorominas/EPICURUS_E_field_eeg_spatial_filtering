function h=PlotPoints(points,style,pointsize)
%function PlotPoints(points,style,pointsize)
%Plots a set of 3D points;this is just a wrapper for the Matlab plot3
%function 
%
%points: N x 3-vector
%pointsize and style are optional; see the help of plot3 function
%h: plot handle
%
%Matti Stenroos (matti.stenroos@aalto.fi)
%version 130610
if nargin==1 || isempty(style),style='k.';end
if nargin<3 || isempty(pointsize),pointsize=20;end
h=plot3(points(:,1),points(:,2),points(:,3),style);
set(h,'MarkerSize',pointsize);
