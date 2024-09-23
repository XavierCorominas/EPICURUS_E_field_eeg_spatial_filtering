% DeFleCT example 2: Targeting a ROI and blocking some crosstalk
%
% This script shows how to focus a spatial filter to a ROI and block
% crosstalk from some source regions using the DeFleCT method --- and what
% happens, if we try to block "too much". The script implements the
% examples presented in Figures 2 and 3 of 
% Hauk and Stenroos: A Framework for the Design of Flexible Cross-Talk
% Functions for Spatial Filtering of EEG/MEG Data: DeFleCT. HBM 2013.
%
% Matti Stenroos (matti.stenroos@aalto.fi) 
% MRC CBU, Cambridge, UK & Aalto University, Espoo, Finland 
% version 130610

clear
%load model
%set here your path to the modeldata:
datapath='/proj/bioem/exampledata/deflectmodel';
load(datapath);
% EMEG data from MNE sample dataset, imported from MNE software:
%   C: noise covariance matrix for 364 "good" sensors
%   L: lead field matrix / forward model for 364 sensors and 8192 sources
%   sources: source meshes for which the L was computed.
% ROIs used in this example: ROI_V1tip, ROI_V1tipstop, ROI_LOstop, 
%   ROI_V1tipstop_addition 
% Parameters for regularization: truncval, SNR2
%%
%choose the left hemisphere mesh
plotmesh=sources.smeshes{1};
plotind=1:plotmesh.nop;

%make the MN estimator for all sources in the model
SNR2=10;
Linv=MNestimator(L,C,SNR2,truncval);
%let's look at the left hemisphere here only
Linv=Linv(plotind,:);
% We could make a resolution kernel for all left hemisphere sources...
% R=Linv*L(:,plotind);
% but let's focus directly on our target...

%compute MN crosstalk for all 10 source points in V1tip
ct_V1tip=Linv(ROI_V1tip,:)*L(:,plotind);
ct_V1tip=Linv(ROI_V1tip,:)*L(:,:);

%compute crosstalk for the whole ROI by summing; the absolute value is just a
%convention.
plotdata1=abs(sum(ct_V1tip,1));%

%make a plot...
plotname1=sprintf('Sum of all CTFs for V1tip');
cmap=jet(20);
cmap=cmap(11:20,:);
scale=.16;
viewangle=[10 -10];
set(figure(1),'position',[0 100 600 350],'toolbar','none','menubar','none','name',plotname1);clf;hold on
%PlotMeshData(plotmesh,plotdata1,'inflated',1,'caxis',scale,'colorbar',1,'colormap',cmap,'view',viewangle);
PlotMeshData(plotmesh,plotdata1);
camzoom(1.45);
PlotPoints(plotmesh.pinf(ROI_V1tip,:),'b.',10);%target
PlotPoints(plotmesh.pinf(ROI_LOstop,:),'r.',10);%sources at LO area that we want to block
PlotPoints(plotmesh.pinf(ROI_V1tipstop,:),'k.',10);%other sources to be blocked
%%

%now filter away the crosstalk 
ROI_stop=union(ROI_V1tipstop,ROI_LOstop);%make stopband
%make deflect filter with
%passband: first component of ROI_V1tip, force response to 1
%stopband: ROI_stop, force all to zero

w=DeFleCT(ROI_V1tip,1,1,ROI_stop,[],L,C,SNR2,truncval);
ct_deflect=w*L(:,plotind);
plotdata2=abs(ct_deflect);
scale=.5;
plotname2=sprintf('Optimized crosstalk for V1tip');
set(figure(2),'position',[600 100 600 350],'toolbar','none','menubar','none','name',plotname2);clf;hold on
%PlotMeshData(plotmesh,plotdata2,'inflated',1,'caxis',scale,'colorbar',1,'colormap',cmap,'view',viewangle);
PlotMeshData(plotmesh,plotdata2)
camzoom(1.3);
PlotPoints(plotmesh.pinf(ROI_V1tip,:),'b.',10);
PlotPoints(plotmesh.pinf(ROI_LOstop,:),'r.',10);
PlotPoints(plotmesh.pinf(ROI_V1tipstop,:),'k.',10);
%%
%Now try too much...
%add two sources to stopband
ROI_stop2=union(ROI_V1tipstop,ROI_V1tipstop_addition);
%make deflect filter with
%passband: first component of ROI_V1tip, force response to 1
%stopband: ROI_stop2, force all to zero
w=DeFleCT(ROI_V1tip,1,1,ROI_stop2,[],L,C,SNR2,truncval);
ct_deflect=w*L(:,plotind);
plotdata3=abs(ct_deflect);
scale=.5;
plotname3=sprintf('Crosstalk optimization gone wrong');
set(figure(3),'position',[600 450 600 350],'toolbar','none','menubar','none','name',plotname3);clf;hold on
%PlotMeshData(plotmesh,plotdata3,'inflated',1,'caxis',scale,'colorbar',1,'colormap',cmap,'view',viewangle);
PlotMeshData(plotmesh,plotdata3);
camzoom(1.3);
PlotPoints(plotmesh.pinf(ROI_V1tip,:),'b.',10);
PlotPoints(plotmesh.pinf(ROI_LOstop,:),'r.',10);
PlotPoints(plotmesh.pinf(ROI_V1tipstop,:),'k.',10);
PlotPoints(plotmesh.pinf(ROI_V1tipstop_addition,:),'y.',10);

