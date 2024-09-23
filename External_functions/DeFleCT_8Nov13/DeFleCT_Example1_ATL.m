% DeFleCT example 1: Combining Discrete and Distributed Source Constraints
%
% This script shows, how to to block a ROI from a minimum-norm-type spatial
% filter using the DeFleCT method. The script implements the example
% presented in Figure 1 of 
% Hauk and Stenroos: A Framework for the Design of Flexible Cross-Talk
% Functions for Spatial Filtering of EEG/MEG Data: DeFleCT. HBM 2013.
%
% Matti Stenroos (matti.stenroos@aalto.fi) 
% MRC CBU, Cambridge, UK & Aalto University, Espoo, Finland 
% 17 Oct 2013

clear
%load model:
%set here your path to the modeldata:
% datapath='/imaging/ms08/exampledata/deflectmodel';
datapath='/proj/bioem/exampledata/deflectmodel';
load(datapath);
% EMEG data from MNE sample dataset, imported from MNE software:
%   C: noise covariance matrix for 364 "good" sensors
%   L: lead field matrix / forward model for 364 sensors and 8192 sources
%   sources: source meshes for which the L was computed.
% ROIs used in this example: ROI_ATL, ROI_IFG, 
% Parameters for regularization: truncval, SNR2

%choose the left hemisphere mesh
plotmesh=sources.smeshes{1};
plotind=1:plotmesh.nop;

%make the MN estimator for all sources in the model
Linv=MNestimator(L,C,SNR2,truncval);
%let's look at the left hemisphere here only
Linv=Linv(plotind,:);
% We could make a resolution kernel for all left hemisphere sources...
% R=Linv*L(:,plotind);

%...but let's go straight to the crosstalk of an example source (at ATL)
anaind=34;
ct_MNE=Linv(ROI_ATL(anaind),:)*L(:,plotind);%crosstalk for the example source position
plotdata1=abs(ct_MNE);%following the common convention, omit the sign...

%make a plot...
plotname1=sprintf('Sample CTF (ATL)');
cmap=jet(20);
cmap=cmap(11:20,:);
scale=10e-3;
set(figure(1),'position',[0 100 600 350],'toolbar','none','menubar','none','name',plotname1);clf;hold on
PlotMeshData(plotmesh,plotdata1,'inflated',1,'caxis',scale,'colorbar',1,'colormap',cmap);
camzoom(1.3);
PlotPoints(plotmesh.pinf(ROI_ATL,:),'k.',5);%add the ATL, IFG, and target
PlotPoints(plotmesh.pinf(ROI_IFG,:),'r.',5);
PlotPoints(plotmesh.pinf(ROI_ATL(anaind),:),'b.',20);

%Now filter away the crosstalk from the IFG region:
%Make deflect filter with
%passband: ROI_ATL(anaind), do not force response to 1.
%stopband: ROI_IFG, force 6 first components to 0
w=DeFleCT(ROI_ATL(anaind),[],0,ROI_IFG,6,L,C,SNR2,truncval);
ct_deflect=w*L(:,plotind);%compute crosstalk
plotdata2=abs(ct_deflect);

plotname2=sprintf('Sample CTF (optimized ATL)');
set(figure(2),'position',[600 100 600 350],'toolbar','none','menubar','none','name',plotname2);clf;hold on
PlotMeshData(plotmesh,plotdata2,'inflated',1,'caxis',scale,'colorbar',1,'colormap',cmap);
camzoom(1.3);
PlotPoints(plotmesh.pinf(ROI_ATL,:),'k.',5);
PlotPoints(plotmesh.pinf(ROI_IFG,:),'r.',5);
PlotPoints(plotmesh.pinf(ROI_ATL(anaind),:),'b.',20);

