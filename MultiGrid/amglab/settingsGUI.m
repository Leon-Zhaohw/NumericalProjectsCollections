function varargout = settingsgui(varargin)
amg_globals;

% SETTINGSGUI Application M-file for settingsgui.fig
%    FIG = SETTINGSGUI launch settingsgui GUI.
%    SETTINGSGUI('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 25-Feb-2010 10:41:08
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
    end
    
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end


% | ABOUT CALLBACKS:
% | GUIDE automatically appends subfunction prototypes to this file, and 
% | sets objects' callback properties to call them through the FEVAL 
% | switchyard above. This comment describes that mechanism.
% |
% | Each callback subfunction declaration has the following form:
% | <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
% |
% | The subfunction name is composed using the object's Tag and the 
% | callback type separated by '_', e.g. 'slider2_Callback',
% | 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
% |
% | H is the callback object's handle (obtained using GCBO).
% |
% | EVENTDATA is empty, but reserved for future use.
% |
% | HANDLES is a structure containing handles of components in GUI using
% | tags as fieldnames, e.g. handles.figure1, handles.slider2. This
% | structure is created at GUI startup using GUIHANDLES and stored in
% | the figure's application data using GUIDATA. A copy of the structure
% | is passed to each callback.  You can store additional information in
% | this structure at GUI startup, and you can change the structure
% | during callbacks.  Call guidata(h, handles) after changing your
% | copy to replace the stored original so that subsequent callbacks see
% | the updates. Type "help guihandles" and "help guidata" for more
% | information.
% |
% | VARARGIN contains any extra arguments you have passed to the
% | callback. Specify the extra arguments by editing the callback
% | property in the inspector. By default, GUIDE sets the property to:
% | <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
% | Add any extra arguments after the last argument, before the final
% | closing parenthesis.
% --------------------------------------------------------------------

% --------------------------------------------------------------------
function src_user_Callback(hObject, eventdata, handles)
% hObject    handle to src_user (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=USER_SPEC_FD;

% --------------------------------------------------------------------
function src_usr2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5src_usr2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=USER_SPEC_FEM;

% --------------------------------------------------------------------
function src_usr3_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_5src_usr3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=USER_SPEC_UF;

% --------------------------------------------------------------------
function setup_rugestueben_Callback(hObject, eventdata, handles)
% hObject    handle to setup_rugestueben (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SETUP_ALG=RUGE_STUEBEN;

% --------------------------------------------------------------------
function smooth_jac_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_jac (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SMOOTHER=JACOBI;

% --------------------------------------------------------------------
function smooth_SOR_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_SOR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SMOOTHER=SOR;

% --------------------------------------------------------------------
function smooth_rich_Callback(hObject, eventdata, handles)
% hObject    handle to smooth_rich (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SMOOTHER=RICHARDSON;

% --------------------------------------------------------------------
function slove_smooth_Callback(hObject, eventdata, handles)
% hObject    handle to slove_smooth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
C_SOLVER=THE_SMOOTHER;

% --------------------------------------------------------------------
function solve_elim_Callback(hObject, eventdata, handles)
% hObject    handle to solve_elim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
C_SOLVER=DIRECT_ELIM;

% --------------------------------------------------------------------
function ps_on_Callback(hObject, eventdata, handles)
% hObject    handle to ps_on (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
POST_SMOOTH=YES;

% --------------------------------------------------------------------
function ps_off_Callback(hObject, eventdata, handles)
% hObject    handle to ps_off (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
POST_SMOOTH=NO;

% --------------------------------------------------------------------
function TwoLevels_Callback(hObject, eventdata, handles)
% hObject    handle to TwoLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COARSEST;
COARSEST=2;

% --------------------------------------------------------------------
function ThreeLevels_Callback(hObject, eventdata, handles)
% hObject    handle to ThreeLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COARSEST;
COARSEST=3;

% --------------------------------------------------------------------
function FourLevels_Callback(hObject, eventdata, handles)
% hObject    handle to FourLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COARSEST;
COARSEST=4;

% --------------------------------------------------------------------
function FiveLevels_Callback(hObject, eventdata, handles)
% hObject    handle to FiveLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COARSEST;
COARSEST=5;

% --------------------------------------------------------------------
function TenLevels_Callback(hObject, eventdata, handles)
% hObject    handle to TenLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COARSEST;
COARSEST=10;

% --------------------------------------------------------------------
function FifteenLevels_Callback(hObject, eventdata, handles)
% hObject    handle to FifteenLevels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global COARSEST;
COARSEST=15;

% --------------------------------------------------------------------
function OneCycle_Callback(hObject, eventdata, handles)
% hObject    handle to OneCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CYCLES;
CYCLES=1;

% --------------------------------------------------------------------
function TwoCycle_Callback(hObject, eventdata, handles)
% hObject    handle to TwoCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CYCLES;
CYCLES=2;

% --------------------------------------------------------------------
function ThreeCycle_Callback(hObject, eventdata, handles)
% hObject    handle to ThreeCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CYCLES;
CYCLES=3;

% --------------------------------------------------------------------
function FourCycle_Callback(hObject, eventdata, handles)
% hObject    handle to FourCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CYCLES;
CYCLES=4;

% --------------------------------------------------------------------
function FiveCycle_Callback(hObject, eventdata, handles)
% hObject    handle to FiveCycle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global CYCLES;
CYCLES=5;

% --------------------------------------------------------------------
function threshold3_Callback(hObject, eventdata, handles)
% hObject    handle to threshold5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_THRESHOLD;
STOP_VALUE=0.001;

% --------------------------------------------------------------------
function threshold5_Callback(hObject, eventdata, handles)
% hObject    handle to threshold5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_THRESHOLD;
STOP_VALUE=0.00001;

% --------------------------------------------------------------------
function threshold7_Callback(hObject, eventdata, handles)
% hObject    handle to threshold7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_THRESHOLD;
STOP_VALUE=0.0000001;

% --------------------------------------------------------------------
function threshold9_Callback(hObject, eventdata, handles)
% hObject    handle to threshold9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_THRESHOLD;
STOP_VALUE=0.000000001;

% --------------------------------------------------------------------
function reduction10_Callback(hObject, eventdata, handles)
% hObject    handle to reduction10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_REDUCE;
STOP_VALUE=0.1;

% --------------------------------------------------------------------
function reduction1_Callback(hObject, eventdata, handles)
% hObject    handle to reduction1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_REDUCE;
STOP_VALUE=0.01;

% --------------------------------------------------------------------
function reduction01_Callback(hObject, eventdata, handles)
% hObject    handle to reduction01 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_REDUCE;
STOP_VALUE=0.001;

% --------------------------------------------------------------------
function reduction001_Callback(hObject, eventdata, handles)
% hObject    handle to reduction001 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE=RESID_REDUCE;
STOP_VALUE=0.0001;


% --- Executes on button press in Finished_Settings.
function Finished_Settings_Callback(hObject, eventdata, handles)
% hObject    handle to Finished_Settings (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pause off;
close;


% --------------------------------------------------------------------
function setup_atonce_Callback(hObject, eventdata, handles)
% hObject    handle to setup_atonce (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SETUP_OPT = AT_ONCE;

% --------------------------------------------------------------------
function setup_ateach_Callback(hObject, eventdata, handles)
% hObject    handle to setup_ateach (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SETUP_OPT = AT_EACH;

% --------------------------------------------------------------------
function setup_amgm_Callback(hObject, eventdata, handles)
% hObject    handle to setup_amgm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SETUP_ALG = AMGm;

% --------------------------------------------------------------------
function setup_smoothed_Callback(hObject, eventdata, handles)
% hObject    handle to setup_smoothed (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SETUP_ALG = SMOOTH_AGGREGATE;

% --------------------------------------------------------------------
function setup_amge_Callback(hObject, eventdata, handles)
% hObject    handle to setup_amge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SETUP_ALG = AMGe;


% --------------------------------------------------------------------
function setup_identity_Callback(hObject, eventdata, handles)
% hObject    handle to setup_identity (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
SETUP_ALG = IDENTITY;


% --------------------------------------------------------------------
function cycletype_v_Callback(hObject, eventdata, handles)
% hObject    handle to cycletype_v (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cycletype_w_Callback(hObject, eventdata, handles)
% hObject    handle to cycletype_w (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function cycletype_f_Callback(hObject, eventdata, handles)
% hObject    handle to cycletype_f (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function NinePoints_Callback(hObject, eventdata, handles)
% hObject    handle to NinePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = STIFFNESS;
FINEPOINTS = 9;

% --------------------------------------------------------------------
function SixteenPoints_Callback(hObject, eventdata, handles)
% hObject    handle to SixteenPoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = STIFFNESS;
FINEPOINTS = 16;

% --------------------------------------------------------------------
function TwentyfivePoints_Callback(hObject, eventdata, handles)
% hObject    handle to TwentyfivePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = STIFFNESS;
FINEPOINTS = 25;

% --------------------------------------------------------------------
function FortyninePoints_Callback(hObject, eventdata, handles)
% hObject    handle to FortyninePoints (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = STIFFNESS;
FINEPOINTS = 49;


% --------------------------------------------------------------------
function EighteenElement_Callback(hObject, eventdata, handles)
% hObject    handle to EighteenElement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = ELEMENTS;
FINEPOINTS=18;

% --------------------------------------------------------------------
function ThirtyTwoElement_Callback(hObject, eventdata, handles)
% hObject    handle to ThirtyTwoElement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = ELEMENTS;
FINEPOINTS=32;

% --------------------------------------------------------------------
function FiftyElement_Callback(hObject, eventdata, handles)
% hObject    handle to FiftyElement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = ELEMENTS;
FINEPOINTS=50;

% --------------------------------------------------------------------
function NinetyEightElement_Callback(hObject, eventdata, handles)
% hObject    handle to NinetyEightElement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
PROB_SRC=EXAMPLE1;
PROB_TYPE = ELEMENTS;
FINEPOINTS=98;


% --------------------------------------------------------------------
function maxitrs_1_Callback(hObject, eventdata, handles)
% hObject    handle to maxitrs_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE = MAX_ITRS;
STOP_VALUE = 1;

% --------------------------------------------------------------------
function maxitrs_10_Callback(hObject, eventdata, handles)
% hObject    handle to maxitrs_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE = MAX_ITRS;
STOP_VALUE = 10;

% --------------------------------------------------------------------
function maxitr_20_Callback(hObject, eventdata, handles)
% hObject    handle to maxitr_20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE = MAX_ITRS;
STOP_VALUE = 20;

% --------------------------------------------------------------------
function maxitr_50_Callback(hObject, eventdata, handles)
% hObject    handle to maxitr_50 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE = MAX_ITRS;
STOP_VALUE = 50;

% --------------------------------------------------------------------
function maxitr_100_Callback(hObject, eventdata, handles)
% hObject    handle to maxitr_100 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;
STOP_TYPE = MAX_ITRS;
STOP_VALUE = 100;

% --------------------------------------------------------------------
function setAbsMaxItrs_Callback(hObject, eventdata, handles)
% hObject    handle to setAbsMaxItrs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

amg_globals;

ABSOLUTE_MAX_ITRS = input('\nPlease enter the desired iteration limit.\n(Press Return for the Default 100):');
if isempty(ABSOLUTE_MAX_ITRS)
    ABSOLUTE_MAX_ITRS = 100;
end


% --------------------------------------------------------------------
function setSOROmega_Callback(hObject, eventdata, handles)
% hObject    handle to setSOROmega (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

amg_globals;

OMEGA = input('\nSOR omega.\n(Default 1.0):');
if isempty(OMEGA)
    OMEGA = 1.0;
end


% --------------------------------------------------------------------
function OneLevel_Callback(hObject, eventdata, handles)
% hObject    handle to OneLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% This was added to have just one level of multigrid.
global COARSEST;
COARSEST=1;


% --------------------------------------------------------------------
function iter_yes_Callback(hObject, eventdata, handles)
% hObject    handle to iter_yes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;

SHOW_IH = YES;

% --------------------------------------------------------------------
function iter_no_Callback(hObject, eventdata, handles)
% hObject    handle to iter_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;

SHOW_IH = NO;


% --------------------------------------------------------------------
function tg_beck_Callback(hObject, eventdata, handles)
% hObject    handle to tg_beck (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;

SETUP_ALG = BECK;


% --------------------------------------------------------------------
function Set_theta_Callback(hObject, eventdata, handles)
% hObject    handle to Set_theta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;

THETA = input('\nRuge-Stueben theta for strong connections.\n(Default 0.25):');
if isempty(THETA)
    THETA = 0.25;
end


% --------------------------------------------------------------------
function Untitled_9_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_10_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;

DEBUG = 1;


% --------------------------------------------------------------------
function Untitled_11_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
amg_globals;

DEBUG = 0;
