function varargout = AMGLab_GUI(varargin)
amg_globals;
% AMGLAB_GUI M-file for AMGLab_GUI.fig
%      AMGLAB_GUI, by itself, creates a new AMGLAB_GUI or raises the existing
%      singleton*.
%
%      H = AMGLAB_GUI returns the handle to a new AMGLAB_GUI or the handle to
%      the existing singleton*.
%
%      AMGLAB_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AMGLAB_GUI.M with the given input arguments.
%
%      AMGLAB_GUI('Property','Value',...) creates a new AMGLAB_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AMGLab_GUI_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AMGLab_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Copyright 2002-2003 The MathWorks, Inc.

% Edit the above text to modify the response to help AMGLab_GUI

% Last Modified by GUIDE v2.5 12-Jul-2006 16:51:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @AMGLab_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @AMGLab_GUI_OutputFcn, ...
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


% --- Executes just before AMGLab_GUI is made visible.
function AMGLab_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
amg_globals;
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AMGLab_GUI (see VARARGIN)

% Choose default command line output for AMGLab_GUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AMGLab_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);
viewSettings_Callback(hObject, eventdata, handles);

% --- Outputs from this function are returned to the command line.
function varargout = AMGLab_GUI_OutputFcn(hObject, eventdata, handles) 
amg_globals;
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in viewSettings.
function viewSettings_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to viewSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panelText, 'String', '');
set(handles.panelGraph, 'Visible', 'Off');
displaytext=get_settings();
set(handles.panelText, 'String', displaytext);
set(handles.panelText, 'Visible', 'On');
drawnow;

% --- Executes on button press in changeSettings.
function changeSettings_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to changeSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.panelGraph, 'Visible', 'Off');
set(handles.panelText, 'String', 'Changing AMGLab Settings...');
set(handles.panelText, 'Visible', 'On');
settingsFig = settingsGUI();
uiwait( settingsFig );

set(handles.panelText, 'String', '');
displaytext=get_settings();
set(handles.panelText, 'String', displaytext);
set(handles.panelText, 'Visible', 'On');
drawnow;

% --- Executes on button press in runButton.
function runButton_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to runButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
str = 'Executing Specified AMG Solve...                          ';
set(handles.panelGraph, 'Visible', 'Off');
set(handles.panelText, 'String', str);
set(handles.panelText, 'Visible', 'On');
drawnow;

% Derrick Cerwinsky
IH = zeros(COARSEST,1); % Added for smoother itteration count.

SOLN = amg_solve(handles);


% --- Executes on button press in defaultSettings.
function defaultSettings_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to defaultSettings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panelText, 'String', '');
amg_defaults;
viewSettings_Callback(hObject, eventdata, handles);

% --- Executes on button press in viewSolution.
function viewSolution_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to viewSolution (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panelText, 'String', '');
set(handles.panelText, 'Visible', 'Off');

% here's how to plot a 2d result
% xvals = -5:0.5:5;
% t=size(xvals);
% for i=1:t(2),
%    for j=1:t(2),
%        ZMAT(i,j) = sin(xvals(i)) + cos(xvals(j));
%    end
% end
% axes(handles.panelGraph); surf(xvals, xvals, ZMAT);
temp = size(SOLN);
xvals = zeros( temp(1), 1 );
for i=1:temp(1)
    xvals( i , 1 ) = i;
end
axes(handles.panelGraph);

if PLOT_ACTUAL == YES
    check = A(1).matrix \ RHS;
    plot(xvals, SOLN,'o-' ,xvals, check , xvals, abs(check-SOLN));
    legend('Computed Solution', 'Actual Solution', 'Difference');
else
    plot(xvals, SOLN);
    legend('Computed Solution');
end

set(handles.panelGraph, 'Visible', 'On');
zoom on;
drawnow;

% --- Executes on button press in viewResid.
function viewResid_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to viewResid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panelText, 'String', '');
set(handles.panelText, 'Visible', 'Off');

temp = size(RH1);
Lstrings = '';
colorvals = ['g';'b';'r';'c';'m';'k'];
colorvals2 = ['g-.';'b-.';'r-.';'c-.';'m-.';'k-.'];

axes(handles.panelGraph);
    
for i=1:temp(2) % for each level
    xvals(i) = i; % get the x value for the graph
end
    
for i=1:temp(1) % for each cycle
    % plot data for this cycle's initial residuals using alternating colors
    legendstring = sprintf('Cycle %d Initial         ', i);
    Lstrings=[Lstrings;legendstring];
    legendstring = sprintf('Cycle %d Corrected       ', i);
    Lstrings=[Lstrings;legendstring];
    plot( xvals, RH1(i,:) , colorvals( (mod(i,6)+1),: ) ),hold on;
    plot( xvals, RH2(i,:) , colorvals2( (mod(i,6)+1),: ) ),hold on;
    legend(Lstrings), hold on;
end
title('Residual Behavior'),hold on;
xlabel('Coursening Level'),hold on
ylabel('Residual 2-Norm'),hold on;
zoom on;

hold off;

set(handles.panelGraph, 'Visible', 'On');
drawnow;

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.panelText, 'String', '');
set(handles.panelText, 'Visible', 'Off');

temp = size(IH);
Lstrings = '';
colorvals = ['g';'b';'r';'c';'m';'k'];
colorvals2 = ['g-.';'b-.';'r-.';'c-.';'m-.';'k-.'];

axes(handles.panelGraph);
    
for i=1:temp(2) % for each level
    xvals(i) = i; % get the x value for the graph
end
    
for i=1:temp(1) % for each cycle
    % plot data for this cycle's initial residuals using alternating colors
    legendstring = sprintf('Cycle %d Pre-Smooth ', i);
    Lstrings=[Lstrings;legendstring];
    plot( xvals, IH(:) , colorvals( (mod(i,6)+1),: ) ),hold on;
    %if POST_SMOOTH % if we are post-smoothing at each level
    %    plot( xvals, IH2(i,:) , colorvals2( (mod(i,6)+1),: ) ),hold on;
    %    legendstring = sprintf('Cycle %d Post-Smooth', i);
    %    Lstrings=[Lstrings;legendstring];
    %end
    legend(Lstrings), hold on;
end
title('Iteration Behavior'),hold on;
xlabel('Coursening Level'),hold on
ylabel('Iteration Count'),hold on;
zoom on;

hold off;

set(handles.panelGraph, 'Visible', 'On');
drawnow;

% --- Executes on button press in exitButton.
function exitButton_Callback(hObject, eventdata, handles)
amg_globals;
% hObject    handle to exitButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clear all;
close('all');

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
amg_globals;
if get(hObject, 'Value') == 1
    PLOT_ACTUAL = YES;
else
    PLOT_ACTUAL = NO;
end


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
amg_globals;
if get(hObject, 'Value') == 1
    SHOW_PROFILE = YES;
else
    SHOW_PROFILE = NO;
end


