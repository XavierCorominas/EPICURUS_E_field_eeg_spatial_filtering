function varargout = ICA_selectandremove(varargin)
% ICA_SELECTANDREMOVE MATLAB code for ICA_selectandremove.fig
%      ICA_SELECTANDREMOVE, by itself, creates a new ICA_SELECTANDREMOVE or raises the existing
%      singleton*.
%
%      H = ICA_SELECTANDREMOVE returns the handle to a new ICA_SELECTANDREMOVE or the handle to
%      the existing singleton*.
%
%      ICA_SELECTANDREMOVE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ICA_SELECTANDREMOVE.M with the given input arguments.
%
%      ICA_SELECTANDREMOVE('Property','Value',...) creates a new ICA_SELECTANDREMOVE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ICA_selectandremove_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ICA_selectandremove_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ICA_selectandremove

% Last Modified by GUIDE v2.5 01-Mar-2017 12:15:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ICA_selectandremove_OpeningFcn, ...
                   'gui_OutputFcn',  @ICA_selectandremove_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ICA_selectandremove is made visible.
function ICA_selectandremove_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ICA_selectandremove (see VARARGIN)

% Choose default command line output for ICA_selectandremove
handles.output = hObject;

handles.EEG= varargin{1};
handles.index_ch2remove=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ICA_selectandremove wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ICA_selectandremove_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_remICs.
function listbox_remICs_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_remICs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_remICs contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_remICs


% --- Executes during object creation, after setting all properties.
function listbox_remICs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_remICs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_allIC.
function listbox_allIC_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_allIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_allIC contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_allIC


% --- Executes during object creation, after setting all properties.
function listbox_allIC_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_allIC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in rem.
function rem_Callback(hObject, eventdata, handles)
% hObject    handle to rem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% rej_butt_Callback(hObject, eventdata, handles)
tmp=get(handles.listbox_allIC,'Value')';
handles.index_ch2remove=str2num(get(handles.listbox_remICs,'String'));
handles.index_ch2remove=unique([handles.index_ch2remove;tmp]);%Indici dei canali bad,ordinati e senza ripetizioni
set(handles.listbox_remICs,'String',handles.index_ch2remove); %Aggiorno la lista dei canali che ho rimosso
[compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff(1:size(handles.EEG.icaweights,1), handles.index_ch2remove));
compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);             
axes(handles.current_bfly);
cla(handles.current_bfly);
plot(handles.EEG.times,mean(compproj,3)','b'),xlim([handles.EEG.times(1) handles.EEG.times(end)])
%variance
% tmp=str2num(get(handles.listbox_remICs,'String'));
set(handles.current_var,'string',strcat('variance:',num2str(varegg*100)));
ICs_list_Callback(hObject, eventdata, handles)
% tmpcurrent=get(handles.ICs_list,'Value');
% [compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff(1:size(handles.EEG.icaweights,1), unique(union(handles.index_ch2remove,tmpcurrent))));
% compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);             
% axes(handles.proposed_bfly);
% cla(handles.proposed_bfly);
% plot(handles.EEG.times,mean(compproj,3)','b'),xlim([-600 600])
% %variance
% % tmp=str2num(get(handles.listbox_remICs,'String'));
% set(handles.residual_var,'string',strcat('variance:',num2str(100-round(sum(handles.EEG.compvars(unique(union(handles.index_ch2remove,tmpcurrent))))*10)/10)));
% 
% % handles.goodchannels=setdiff(1:63,handles.index_ch2remove);
% % setappdata(handles.rej_butt,'index_ch2remove',handles.index_ch2remove);
% % setappdata(handles.rej_butt,'goodchannels',handles.goodchannels);
% % setappdata(handles.rej_butt,'flag_remove',0);
% set(handles.listbox_allIC,'Value',1);
guidata(hObject, handles);


% --- Executes on button press in keep.
function keep_Callback(hObject, eventdata, handles)
% hObject    handle to keep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%add_butt_Callback(hObject, eventdata, handles)
handles.index_ch2remove=str2num(get(handles.listbox_remICs,'String'));
if ~isempty(get(handles.listbox_remICs,'String'))
    handles.index_ch2insert=get(handles.listbox_remICs,'Value');%Indice dei canale da reinserire
    handles.index_ch2remove(handles.index_ch2insert)=[];
    if ~isempty(handles.index_ch2remove)
        set(handles.listbox_remICs,'Value',1);
        set(handles.listbox_remICs,'String',num2str(handles.index_ch2remove));
    else set(handles.listbox_remICs,'Value',[]);
        set(handles.listbox_remICs,'String','');
    end
    
    [compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff(1:size(handles.EEG.icaweights,1), handles.index_ch2remove));
    compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);
    axes(handles.current_bfly);
    cla(handles.current_bfly);
    plot(handles.EEG.times,mean(compproj,3)','b'),xlim([handles.EEG.times(1) handles.EEG.times(end)])
    %variance
%     tmp=str2num(get(handles.listbox_remICs,'String'));
    set(handles.current_var,'String',strcat('variance:',num2str(varegg*100)));
    
%     tmpcurrent=get(handles.ICs_list,'Value');
%     [compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff(1:size(handles.EEG.icaweights,1), unique(union(handles.index_ch2remove,tmpcurrent))));
%     compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);
%     axes(handles.proposed_bfly);
%     cla(handles.proposed_bfly);
%     plot(handles.EEG.times,mean(compproj,3)','b'),xlim([-600 600])
%     %variance
%     % tmp=str2num(get(handles.listbox_remICs,'String'));
%     set(handles.residual_var,'string',strcat('variance:',num2str(100-round(sum(handles.EEG.compvars(unique(union(handles.index_ch2remove,tmpcurrent))))*10)/10)));
%     
%     
%     %     handles.goodchannels=setdiff(1:63,handles.index_ch2remove);
% %     setappdata(handles.rej_butt,'goodchannels',handles.goodchannels);
% %     setappdata(handles.rej_butt,'index_ch2remove',handles.index_ch2remove);
%     set(handles.listbox_remICs,'Value',1);
ICs_list_Callback(hObject, eventdata, handles)
end
guidata(hObject, handles);


% --- Executes on selection change in ICs_list.
function ICs_list_Callback(hObject, eventdata, handles)
% hObject    handle to ICs_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.ICs_list,'String',[1:size(handles.EEG.icaact,1)]');

selcomp=get(handles.ICs_list,'Value');
tmp=get(handles.listbox_remICs,'String');
alreadyrem=str2num(tmp);

    [compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff(1:size(handles.EEG.icaweights,1), union(alreadyrem,selcomp)));
    compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);
    axes(handles.proposed_bfly);
    cla(handles.proposed_bfly);
    tmpy=get(handles.current_bfly,'YLim');
    plot(handles.EEG.times,mean(compproj,3)','b'),ylim(tmpy),xlim([handles.EEG.times(1) handles.EEG.times(end)]);
    %variance
    % tmp=str2num(get(handles.listbox_remICs,'String'));
    set(handles.residual_var,'string',strcat('variance:',num2str(varegg*100)));

    
    
%     camp_butterfly=get(handles.ICs_list,'Value');
data=handles.EEG.icaact(selcomp,:,:);
[spectra,freqs,speccomp,contrib,specstd] = spectopo( data, handles.EEG.pnts, handles.EEG.srate, 'mapnorm', handles.EEG.icawinv(:,selcomp),'plot','off');
freqslim = 50;
a = max(find(freqs <= freqslim));
%
%                             axes(handles.Display.CurAnalysis.AxSpec)
axes(handles.ICspectrum)
plot(freqs(1:a),spectra(1:a),'-r');
% set(he,'LineWidth',2)
% set(gca, 'xlim', [0 min(freqslim, EEG.srate/2)]);
% set(gca, 'ylim', [min(spectra(1:a)) max(spectra(1:a))]);
% set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)');
% set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)');
% title('Power Spectrum and Scalp Map Power')

%topo
axes(handles.ICtopo)
cla(handles.ICtopo)
topoplot((handles.EEG.icawinv(:,selcomp)),handles.EEG.chanlocs,'electrodes','on','maplimits','maxmin'),caxis([-max(abs(handles.EEG.icawinv(:,selcomp))) max(abs(handles.EEG.icawinv(:,selcomp)))]) ;



%imagesc
axes(handles.component_imagesc);
tmptoplot=permute(handles.EEG.icaact(selcomp,:,:),[3,2,1]);
% imagesc(handles.EEG.times,[1:size(handles.EEG.icaact,3)],tmptoplot),colormap('jet'),colorbar,caxis([-max(max(abs(tmptoplot)))*1/3 max(max(abs(tmptoplot)))*1/3])
% [outdata,outvar,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx,erpsig]=ssp_erpimage( tmptoplot, ones(1,size(handles.EEG.icaact,3)), handles.EEG.times , '', 3, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','noshow','on');
cla(handles.component_imagesc);
imagesc(tmptoplot),colormap('jet'), caxis([-max(max(abs(tmptoplot)))*1/3 max(max(abs(tmptoplot)))*1/3])%,colorbar 
set(handles.component_imagesc,'XTickLabel',[]);
%  axis([handles.EEG.times(1) handles.EEG.times(end) 1 handles.EEG.trials]);
% erpimage(permute(handles.EEG.icaact(tmp,:,:),[3,2,1]),[1:size(handles.EEG.icaact,3)],handles.EEG.times,'',3,1,'caxis',2/3,'cbar','erp', 'yerplabel','noshow','on');

% [outdata,outvar,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx,erpsig] =ssp_erpimage( data-offset, ones(1,EEG.trials)*10000, EEG.times , '', 3, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','noshow','on');
% [spectra,freqs,speccomp,contrib,specstd] = spectopo( data, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,v),'plot','off');
% freqslim = 50;
% a = max(find(freqs <= freqslim));
% pause(1)
% 
% 
%averageIC
axes(handles.averageIC);
tmptoplot=permute(handles.EEG.icaact(selcomp,:,:),[3,2,1]);
plot(handles.EEG.times,mean(tmptoplot,1),'b'),xlim([handles.EEG.times(1) handles.EEG.times(end)]);



guidata(hObject, handles);

% Hints: contents = cellstr(get(hObject,'String')) returns ICs_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ICs_list


% --- Executes during object creation, after setting all properties.
function ICs_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ICs_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in up.
function up_Callback(hObject, eventdata, handles)
% hObject    handle to up (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.ICs_list,'String',[1:size(handles.EEG.icaact,1)]');

guidata(hObject, handles);

% --- Executes on button press in down.
function down_Callback(hObject, eventdata, handles)
% hObject    handle to down (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function current_var_Callback(hObject, eventdata, handles)
% hObject    handle to current_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


guidata(hObject, handles);

% Hints: get(hObject,'String') returns contents of current_var as text
%        str2double(get(hObject,'String')) returns contents of current_var as a double


% --- Executes during object creation, after setting all properties.
function current_var_CreateFcn(hObject, eventdata, handles)
% hObject    handle to current_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function residual_var_Callback(hObject, eventdata, handles)
% hObject    handle to residual_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of residual_var as text
%        str2double(get(hObject,'String')) returns contents of residual_var as a double


% --- Executes during object creation, after setting all properties.
function residual_var_CreateFcn(hObject, eventdata, handles)
% hObject    handle to residual_var (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.listbox_allIC,'String',[1:size(handles.EEG.icaact,1)]');
set(handles.ICs_list,'String',[1:size(handles.EEG.icaact,1)]');
handles.index_ch2remove=[];
%butterfly
for i=1:size(handles.EEG.data,1)
    l=plot(handles.current_bfly,handles.EEG.times,zeros(1,length(handles.EEG.times)),'DisplayName','ch','LineWidth',1,'Color','b');
    set(handles.current_bfly,'XLim',[handles.EEG.times(1) handles.EEG.times(end)])
    handles.butterfly_current(i)=l; %Puntatori alle righe
    hold(handles.current_bfly,'on');
end
if ~isempty(handles.EEG.comp2remove)
    set(handles.listbox_remICs,'String',handles.EEG.comp2remove);
[compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff(1:size(handles.EEG.icaweights,1), handles.EEG.comp2remove));
compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);             
   
    tmp_butterfly=mat2cell(mean(compproj,3),ones(size(handles.EEG.data,1),1),length(handles.EEG.times));
else
    tmp_butterfly=mat2cell(mean(handles.EEG.data,3),ones(size(handles.EEG.data,1),1),length(handles.EEG.times));
       
end
       set(handles.butterfly_current,{'YData'},tmp_butterfly,'LineWidth',1,'Color','b');



for i=1:size(handles.EEG.data,1)
    l=plot(handles.proposed_bfly,handles.EEG.times,zeros(1,length(handles.EEG.times)),'DisplayName','ch','LineWidth',1,'Color','b');
        set(handles.proposed_bfly,'XLim',[handles.EEG.times(1) handles.EEG.times(end)])
    handles.butterfly_proposed(i)=l; %Puntatori alle righe
    hold(handles.proposed_bfly,'on');
end
tmp1=str2num(get(handles.listbox_remICs,'String'));
tmp2=get(handles.ICs_list,'Value');
tmp3=union(tmp1,tmp2);
[compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff([1:size(handles.EEG.icaact,1)],tmp3));
compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);             
tmp_butterfly=mat2cell(mean(compproj,3),ones(size(handles.EEG.data,1),1),length(handles.EEG.times));
set(handles.butterfly_proposed,{'YData'},tmp_butterfly,'LineWidth',1,'Color','b');

%variance
set(handles.residual_var,'string',strcat('variance:',num2str(varegg*100)));

%IC spectrum
camp_butterfly=get(handles.ICs_list,'Value');
data=handles.EEG.icaact(camp_butterfly,:,:);
[spectra,freqs,speccomp,contrib,specstd] = spectopo( data, handles.EEG.pnts, handles.EEG.srate, 'mapnorm', handles.EEG.icawinv(:,camp_butterfly),'plot','off');
freqslim = 50;
a = max(find(freqs <= freqslim));
%
%                             axes(handles.Display.CurAnalysis.AxSpec)
axes(handles.ICspectrum)
plot(freqs(1:a),spectra(1:a),'-r');
% set(he,'LineWidth',2)
% set(gca, 'xlim', [0 min(freqslim, EEG.srate/2)]);
% set(gca, 'ylim', [min(spectra(1:a)) max(spectra(1:a))]);
% set( get(gca, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)');
% set( get(gca, 'xlabel'), 'string', 'Frequency (Hz)');
% title('Power Spectrum and Scalp Map Power')

%topo
camp_butterfly=get(handles.ICs_list,'Value');
axes(handles.ICtopo)
topoplot((handles.EEG.icawinv(:,camp_butterfly)),handles.EEG.chanlocs,'electrodes','on','maplimits','maxmin'),caxis([-max(abs(handles.EEG.icawinv(:,camp_butterfly))) max(abs(handles.EEG.icawinv(:,camp_butterfly)))]);%,colorbar ;


%variance
tmp=str2num(get(handles.listbox_remICs,'String'));
[compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff([1:size(handles.EEG.icaact,1)],tmp'));
set(handles.current_var,'string',strcat('variance:',num2str(varegg*100)));


%imagesc
axes(handles.component_imagesc)
tmp=get(handles.ICs_list,'Value');
tmptoplot=permute(handles.EEG.icaact(tmp,:,:),[3,2,1]);
% imagesc(handles.EEG.times,[1:size(handles.EEG.icaact,3)],tmptoplot),colormap('jet'),colorbar,caxis([-max(max(abs(tmptoplot)))*1/3 max(max(abs(tmptoplot)))*1/3])
% [outdata,outvar,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx,erpsig]=ssp_erpimage( tmptoplot, ones(1,size(handles.EEG.icaact,3)), handles.EEG.times , '', 3, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','noshow','on');
imagesc(tmptoplot),colormap('jet'), caxis([-max(max(abs(tmptoplot)))*1/3 max(max(abs(tmptoplot)))*1/3])%,colorbar 
set(handles.component_imagesc,'XTickLabel',[]);
%  axis([handles.EEG.times(1) handles.EEG.times(end) 1 handles.EEG.trials]);
% erpimage(permute(handles.EEG.icaact(tmp,:,:),[3,2,1]),[1:size(handles.EEG.icaact,3)],handles.EEG.times,'',3,1,'caxis',2/3,'cbar','erp', 'yerplabel','noshow','on');

% [outdata,outvar,outtrials,limits,axhndls,erp,amps,cohers,cohsig,ampsig,outamps,phsangls,phsamp,sortidx,erpsig] =ssp_erpimage( data-offset, ones(1,EEG.trials)*10000, EEG.times , '', 3, 1, 'caxis', 2/3, 'cbar','erp', 'yerplabel', '','noshow','on');
% [spectra,freqs,speccomp,contrib,specstd] = spectopo( data, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,v),'plot','off');
% freqslim = 50;
% a = max(find(freqs <= freqslim));
% pause(1)
% 
% 
%averageIC
axes(handles.averageIC);
tmp=get(handles.ICs_list,'Value');
tmptoplot=permute(handles.EEG.icaact(tmp,:,:),[3,2,1]);
plot(handles.EEG.times,mean(tmptoplot,1),'b'),xlim([handles.EEG.times(1) handles.EEG.times(end)]);



guidata(hObject, handles);


% --- Executes on button press in ICscroll.
function ICscroll_Callback(hObject, eventdata, handles)
% hObject    handle to ICscroll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
eegplot( handles.EEG.icaact, 'srate', handles.EEG.srate, 'title', 'Scroll component activities -- eegplot()', ...
    'limits', [handles.EEG.xmin handles.EEG.xmax]*1000 , 'command', []);



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in pruned_ICA.
function pruned_ICA_Callback(hObject, eventdata, handles)
% hObject    handle to pruned_ICA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.index_comp2remove=str2num(get(handles.listbox_remICs,'String'));
[compproj, varegg ] = SC_eeglab_compvar( handles.EEG.data, { handles.EEG.icasphere handles.EEG.icaweights }, handles.EEG.icawinv, setdiff(1:size(handles.EEG.icaweights,1), handles.index_comp2remove));
compproj = reshape(compproj, handles.EEG.nbchan, handles.EEG.pnts, handles.EEG.trials);             
data_GUI=struct('compproj',compproj,'comp2remove',handles.index_comp2remove,'var_compproj',varegg);
assignin('base','data_GUI',data_GUI)


% --- Executes on button press in pushbutton8_back.
function pushbutton8_back_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8_back (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selcomp=get(handles.ICs_list,'Value');
if selcomp>1
    set(handles.ICs_list,'Value',selcomp-1);
    ICs_list_Callback(hObject, eventdata, handles);
end

guidata(hObject, handles);


% --- Executes on button press in pushbutton9_forward.
function pushbutton9_forward_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9_forward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

selcomp=get(handles.ICs_list,'Value');
if selcomp<size(handles.EEG.icaact,1);
    set(handles.ICs_list,'Value',selcomp+1);
    ICs_list_Callback(hObject, eventdata, handles);
end

guidata(hObject, handles);
