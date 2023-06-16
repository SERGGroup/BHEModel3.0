function varargout = PreosPropsMenu(varargin)
% PREOSPROPSMENU M-file for PreosPropsMenu.fig ver 1,0
% see end of file for summary of major changes.
% This file should be used with ProeosPropsMenu.fig to call and display results
% for PreosPropsMenu.m

% Most 'items' of the menu have two functions below. One of the functions is to 
% create the 'item' and the other 'callback' is executes when the item is changed by typing or clicking.  

% The 'handles' variables are held in the Menu figure workspace.
% Variables in the Menu workspace are initialized using the 'set' function.
% Variables are passed from the 'figure' workspace to the 'base' (default)  workspace using
% the 'assignin' function or the 'evalin' function. 
% Array elements cannot be assigned in the 'base' workspace using the 'assignin' function. Rather
% they must be assigned using an explicit 'evalin' call. Also, when the assigned value is a variable, 
% a temporary variable is created in the 'base' workspace, assigned to the array element,  and then destroyed. 
% See the 'Value' function to see how the value 'match(3)' is passed to the 'base' workspace.
% To view the items used in the Menu, open the PreosPropsMenu.fig in the GUIDE editor, then right-click on a 
% menu item and select the 'properties inspector'. Look for the 'TAG' or 'String'. Most other properties were 
% usually left as defaults.

% The 'handles' structure used in the functions below holds properties for all the 'items' on the Menu. 
% The 'handles' structure holds a child variable for each 'item' TAG on the menu.
% Each of the 'items' has properties.  For example, the entry box for the reference temperature has 
% a 'Tag' name 'Tref'. Futher, 'handles.Tref' has many properties including a 'Value' and a 'String'.
% Most variables are displayed/entered as 'Strings' and are converted to 'Values' or numbers when passed to
% the 'base' workspace.

% The functions below were created by using the GUIDE GUI development feature in Matlab. 
% The functions were given TAG names in the GUIDE environment that were meaningful to their purpose.
% Then this file was edited to process the variables from the menu and communicate with the 'base' workspace.
% For example, the 'Calculate' button was edited to run ProeosProps.
% Comments are included below.

% Many of the comments below were genereated automatically by the GUIDE during menu creation.

%      PREOSPROPSMENU, by itself, creates a new PREOSPROPSMENU or raises the existing
%      singleton*.
%
%      H = PREOSPROPSMENU returns the handle to a new PREOSPROPSMENU or the handle to
%      the existing singleton*.
%
%      PREOSPROPSMENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PREOSPROPSMENU.M with the given input arguments.
%
%      PREOSPROPSMENU('Property','Value',...) creates a new PREOSPROPSMENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PreosPropsMenu_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PreosPropsMenu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Last Modified by GUIDE v2.5 15-Mar-2012 00:04:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PreosPropsMenu_OpeningFcn, ...
                   'gui_OutputFcn',  @PreosPropsMenu_OutputFcn, ...
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


% --- Executes just before PreosPropsMenu is made visible.
function PreosPropsMenu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PreosPropsMenu (see VARARGIN)

% Choose default command line output for PreosPropsMenu
handles.output = hObject;

%Provide default values in both the GUI and base workspace
set(handles.Tref,'String',298.15);
assignin('base','Tref',298.15);
set(handles.Pref,'String',0.1);
assignin('base','Pref',0.1);
set(handles.T,'String',298.15);
assignin('base','T',298.15);
set(handles.P,'String',0.1);
assignin('base','P',0.1);

%controls that are 'selectors' must be initialized.
set(handles.RunType,'SelectedObject',handles.noMatch);
set(handles.rootToUse,'SelectedObject',handles.largeZ);
% set a default value for the 'match' array
assignin('base','match',[0 1 0 1]);

% Update handles structure
guidata(hObject, handles);

% run the program once to initialize
calculate_Callback(hObject, eventdata, handles);

% UIWAIT makes PreosPropsMenu wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PreosPropsMenu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% currently not changed from default. CTL

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in calculate.
function calculate_Callback(hObject, eventdata, handles)
% hObject    handle to calculate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%upon click of 'Calculate' button
R = 8.3143; %MPa.cm^3/mol.K

% get current state of byAdjusting
switch get(handles.byAdjusting,'Value')   % Get selected object index
    case 1
		evalin('base','match(2)=1;');
    case 2
		evalin('base','match(2)=2;');
    otherwise
       sprintf('Error - Unable to interpret byAdjusting input from dropdown box')
end
% the 'evalin' function runs using input form the 'base' workspace, but returns values to the Menu workspace
[Z, H, S, U, phi, info] = evalin('base','PreosProps(Tref, Pref, T, P, match)');

% Transfer the results to the 'base' workspace.
assignin('base','Z',Z);
assignin('base','H',H);
assignin('base','U',U);
assignin('base','phi',phi);
assignin('base','info',info);

%update Menu strings to display results on the menu
set(handles.resultZ,'String',sprintf('%g\t  %g',Z));
% T and P resulting from the calculation are stored in info cell array. The following statement uses Z*R*T/P to calculate volume.
set(handles.resultV,'String',sprintf('%g\t  %g',Z*R*info{1}(1)/info{2}(1)));
set(handles.resultH,'String',sprintf('%g\t  %g', H));
set(handles.resultS,'String',sprintf('%g\t  %g', S));
set(handles.resultU,'String',sprintf('%g\t  %g', U));
% the fugacity is calculated by phi*P where P is extracted from the 'info' cell array.
set(handles.resultFug,'String',sprintf('%g\t  %g', phi*info{2}(1)))
set(handles.resultUserObj,'String', sprintf('%g\t %g',info{5}))
set(handles.resultT,'String',info{1}(1));
set(handles.resultP,'String',info{2}(1));
set(handles.name,'String',info(4));
set(handles.name2,'String',info(4));
exitflag = info{3}(1);

% I'm not sure this is working correctly. to display warning if the solver fails.
if exitflag < 0
   set(handles.textWarn, 'String','Abnormal termination. See command window for error.');
   else
   set(handles.textWarn, 'String','');
end
  

% --- Executes on selection change in byAdjusting.
function byAdjusting_Callback(hObject, eventdata, handles)
% hObject    handle to byAdjusting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns byAdjusting contents as cell array
%        contents{get(hObject,'Value')} returns selected item from byAdjusting
	
% The first  option is to adjust T, second is to adjust P. Set 'match' accordingly.
switch get(hObject,'Value')   % Get selected object index
    case 1
		evalin('base','match(2)=1;');
    case 2
		evalin('base','match(2)=2;');
    otherwise
       sprintf('Error - Unable to interpret byAdjusting input from dropdown box')
end


% --- Executes during object creation, after setting all properties.
function byAdjusting_CreateFcn(hObject, eventdata, handles)
% hObject    handle to byAdjusting (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'Value',1);
% Save the new  value in gui for later processing when
guidata(hObject,handles)



function value_Callback(hObject, eventdata, handles)
% hObject    handle to value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of value as text
%        str2double(get(hObject,'String')) returns contents of value as a double
value = str2double(get(hObject, 'String'));
if isnan(value)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

assignin('base','temp',value);
evalin('base','match(3) = temp;');
%evalin('base','clear temp');
% Save the new  value in gui for later processing when
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function value_CreateFcn(hObject, eventdata, handles)
% this function creates the 'value' textbox
% hObject    handle to value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function RunType_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to RunType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in the RunType radio button panel.
function RunType_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in RunType 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% For a set of radio buttons the selected button is returned in the tag of the button  panel.
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'noMatch'
        evalin('base', 'match(1)=0;');
    case 'matchU'
        evalin('base', 'match(1)=1;');
    case 'matchH'
        evalin('base', 'match(1)=2;');
    case 'matchS'
        evalin('base', 'match(1)=3;');
    case 'matchV'
        evalin('base', 'match(1)=4;');
    case 'findSat'
        evalin('base', 'match(1)=5;');
    case 'userObj'
        evalin('base', 'match(1)=6;');
        evalin('base', 'match(2)=3;');
    otherwise
       sprintf('Error - Unable to interpret RunType input from button box')
end
%store results for later processing
%guidata(hObject,handles)



% --------------------------------------------------------------------
function rootToUse_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to rootToUse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected object is changed in button panel rootToUse.
function rootToUse_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in rootToUse 
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% For radio buttons, the selected buttons is returned as the tag of the button panel.
switch get(hObject,'Tag')   % Get Tag of selected object
    case 'largeZ'
		evalin('base','match(4)=1;');
    case 'smallZ'
		evalin('base','match(4)=2;');
    otherwise
       sprintf('Error - Unable to interpret rootToUse input from button box')
end
%store results for later processing


% ----Executes when the text is changed in textbox T
function T_Callback(hObject, eventdata, handles)
% hObject    handle to T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of T as text
%        str2double(get(hObject,'String')) returns contents of T as a double
T = str2double(get(hObject, 'String'));
if isnan(T)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new T value
assignin('base','T',T);


% --- Executes during object creation, after setting all properties for textbos T.
function T_CreateFcn(hObject, eventdata, handles)
% hObject    handle to T (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ----Executes when the text is changed in textbox P
function P_Callback(hObject, eventdata, handles)
% hObject    handle to P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P as text
%        str2double(get(hObject,'String')) returns contents of P as a double
P = str2double(get(hObject, 'String'));
if isnan(P)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end

% Save the new P value
assignin('base','P',P);

% --- Executes during object creation, after setting all properties for textbos P.
function P_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ----Executes when the text is changed in textbox Tref
function Tref_Callback(hObject, eventdata, handles)
% hObject    handle to Tref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Tref as text
%        str2double(get(hObject,'String')) returns contents of Tref as a double
Tref = str2double(get(hObject, 'String'));
if isnan(Tref)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
assignin('base','Tref',Tref);



% --- Executes during object creation, after setting all properties for textbos Tref.
function Tref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Tref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ----Executes when the text is changed in textbox Pref
function Pref_Callback(hObject, eventdata, handles)
% hObject    handle to Pref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Pref as text
%        str2double(get(hObject,'String')) returns contents of Pref as a double
Pref = str2double(get(hObject, 'String'));
if isnan(Pref)
    set(hObject, 'String', 0);
    errordlg('Input must be a number','Error');
end
assignin('base','Pref',Pref);


% --- Executes during object creation, after setting all properties for textbox Pref.
function Pref_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Pref (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ver 1.02 - 3/7/14 fix bug where 'byAdjusting' did not register when restarting.
% ver 1.01 - 12/12/11 fix typo in comments. Added text above root box. 
% ver 1.0 - 2/22/09 initial release
