% AMG_SOLVE
% 
%        x = AMG_SOLVE(GUI_HANDLES)
%            GUI_HANDLES allows this function to access the disply in the
%            AMGLab GUI

%  Ryan McKenzie
%  Department of Computational Sciences
%  University of Kentucky

function x = amg_solve(GUI_HANDLES);

amg_globals;
str=get(GUI_HANDLES.panelText,'String');

% Check the source of the problem specification
if PROB_SRC==USER_SPEC_FD % if the user finite difference matrix is to be used
    str=strcat(str,sprintf('\n    Initializing Problem to User Specifications              '));
    str=strcat(str,sprintf('\n    USER: Please enter your input file name at the console...'));
    set(GUI_HANDLES.panelText,'String',str);
    drawnow;
    amg_usersetFD; % set up globals as the user has specified them
elseif PROB_SRC==USER_SPEC_FEM % if the user finite element matrix is to be used
    str=strcat(str,sprintf('\n    Initializing Problem to User Specifications              '));
    str=strcat(str,sprintf('\n    USER: Please enter your input file name at the console...'));
    set(GUI_HANDLES.panelText,'String',str);
    drawnow;
    amg_usersetFEM; % set up globals as the user has specified them
elseif PROB_SRC==USER_SPEC_UF % if the user finite difference matrix is to be used
    str=strcat(str,sprintf('\n    Initializing Problem to User Specifications              '));
    str=strcat(str,sprintf('\n    USER: Please enter your input file name at the console...'));
    set(GUI_HANDLES.panelText,'String',str);
    drawnow;
    amg_usersetUF; % set up globals as the user has specified them
elseif PROB_SRC==EXAMPLE1 % if example 1 is called for
    str=strcat(str,sprintf('\n    Initializing Problem to Poisson Example             '));
    set(GUI_HANDLES.panelText,'String',str);
    drawnow;
    amg_example1; % set up globals for example 1
end

profile clear;

if SETUP_OPT==AT_ONCE % if user wishes to find all coarse grids up front
    %  see "restrict.m" for what happens otherwise
    str=strcat(str,sprintf('\n    Calling Setup Algorithm to Determine Coarse Grids   '));
    set(GUI_HANDLES.panelText,'String',str);
    drawnow;
    if SHOW_PROFILE==YES
        profile on;
    end
    amg_setup(1); % call setup algorithm
    if SHOW_PROFILE==YES
        profile off;
    end
end

INIT_RESID = norm(RHS - (A(1).matrix * X_Guess)); % get initial residual
str=strcat(str,sprintf('\n    2-Norm of Initial Residual is %d                     ', INIT_RESID));
set(GUI_HANDLES.panelText,'String',str);
drawnow;

x = X_Guess;
cycleCount = 1; % count the number of AMG cycles executed
while (cycleCount<=CYCLES)
    str=strcat(str,sprintf('\n    Starting AMG Cycle Number %i...                     ', cycleCount));
    set(GUI_HANDLES.panelText,'String',str);
    drawnow;
    if SHOW_PROFILE==YES
        profile resume;
    end
    x = amg_cycle(cycleCount,1,RHS,x);  %  execute as many AMG cycles if necessary
    if SHOW_PROFILE==YES
        profile off;
    end
    str=strcat(str,sprintf('\n    2-Norm of Residual after AMG Cycle Number %i is %d  ', cycleCount, norm(RH2(cycleCount,1))));
    set(GUI_HANDLES.panelText,'String',str);
    drawnow;
    cycleCount=cycleCount+1; % incriment cycle count
end

%  Divya
%  Display the result on screen for each cycle----

x
RH1
RH2
IH
%IH2

% ------------------------------------------------

FINAL_RESID = norm(RHS - (A(1).matrix * x)); % get final residual
str=strcat(str,sprintf('\n    Total Smoothing Iterations: %d                       ', sum(IH) ));

str=strcat(str,sprintf('\nFinished...                                              '));
set(GUI_HANDLES.panelText,'String',str);
drawnow;

if SHOW_PROFILE==YES
    profile viewer;
end
