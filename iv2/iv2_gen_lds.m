function    iv2_gen_lds(fun,i2); 
% To prepare local adaptoer files (dxetc4iv2_gen_lds.m) for your site:
%       
%       usage:      iv2_gen_lds('set','dxetc4xxx')
%       
% 
% (cL)2022    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
if exist(which(['local_',lower(fun)]));
                                feval(['local_',lower(fun)], i2);
else;                           disp(['.unknown local function: ',fun]);                            end;
return;
%%

function                        local_set(i2);
%%
if ~isempty(findobj(groot,'Name','IDAE local adapter Generator'));               	
  	figure(findobj(groot,'Name','IDAE local adapter Generator'));                   return;         end;
%
%
[qdx, qnm, qex]                 = fileparts(i2);
ok                              = 1;
if ~isempty(which(qnm));        disp(['.problem! function: ',qnm,'.m already present.']);
                                ok                          = 0;                                    end;
if ~isempty(qex);               disp(['.problem! no extention (',qex,') is allowed']);
                                ok                          = 0;                                    end;
if ~isempty(qdx);               disp(['.problem! no directory (',qdx,filesep,') is allowed']);
                                ok                          = 0;                                    end;
if ok<1;                                                                            return;         end;
%
%
ud.lds                          = deblank(i2);
ud.idae.code_path               = fileparts(which(mfilename));
cd(ud.idae.code_path);
resume                          = 0;
if exist(fullfile(fileparts(which(mfilename)), [ud.lds,'.mat']),'file');
                                load(fullfile(fileparts(which(mfilename)), [ud.lds,'.mat']));
                                resume                      = 1;                                    end;
%
ns                              = 22;
tstr                            = 'iProject generator';
[fH, bpos]                      = bWindow([], ...
                                'nrs',                      ns+2,               ...
                                'bwd',                      600,                ...
                                'ttl',                      'IDAE local adapter Generator');
set(fH,                         'Toolbar',                  'none',             ...
                                'Menubar',                  'none',             ...
                                'Tag',                      'iv2_gen_lds');
%
bHs                             = postJBs(fH,               'B',bpos(1,:),[8,2;1,1]);
c1bHs(1)                     	= bHs(1);
set(c1bHs(1),   'Tag','iv2_gen_lds_1_1',    'BackgroundColor',iv2_bgcs(12), 'FontWeight','bold');
c2bHs(1)                        = bHs(2);
set(c2bHs(1),   'Tag','iv2_gen_lds_1_2');
%
%
for i=2:1:size(bpos,1)-5;
    bHs                        	= postJBs(fH,               'B',bpos(i,:),[8,2;1,1]);
    c1bHs(i)                   	= bHs(1);
    c2bHs(i)                  	= bHs(2);
    set(bHs(1), 'String',' ',   'Tag',['iv2_gen_lds_',int2str(i),'_1']);
    set(bHs(2), 'String',' ',   'Tag',['iv2_gen_lds_',int2str(i),'_2']);                            end;
%
% information board:
bHs                             = postJBs(fH, 'B',[bpos(end-1, 1:3)+[3,0,-6], bpos(1,4).*4+2],[1;1]);
set(bHs(1), 'String','info',    'Style','text', 'Tag','iv2_gen_lds_infoB',  ...
                                'FontName','Consolas',  'FontSize',10,  'HorizontalAlignment','left');
% utility GUIs @bottom row:
bHs                             = postJBs(fH,               'B',bpos(end,:), [1,1,1,1;1,1,1,1]);
set(bHs(1), 'String','Help',    'Tag','iv2_gen_lds_end_1',	'CallBack','iv2_gen_lds(''help'',[]);');
set(bHs(2), 'String','Display', 'Tag','iv2_gen_lds_end_2',  'CallBack','iv2_gen_lds(''disp'',[]);');
set(bHs(3), 'String','Save',    'Tag','iv2_gen_lds_end_3',  'CallBack','iv2_gen_lds(''save'',[]);');
set(bHs(4), 'String','Quit',    'Tag','iv2_gen_lds_end_4',  'CallBack','delete(gcf);');
%
ud.c1bHs                        = c1bHs;
ud.c2bHs                        = c2bHs;
% 
% when resuming previously saved variables:
if resume>0;                    set(gcf,    'UserData',ud);
                                local_resume(ud);                                   return;         end;
                                
ud.stage                        = 'idae_0';
set(gcf,    'UserData',ud);
%
set(c1bHs(1),   'String','Definitions of computer species for this session');
set(c2bHs(1),   'String','Move on', 'CallBack','iv2_gen_lds(''done'',[]);');
set(c1bHs(2),   'String','Workstations to perform data-analysis with IDAE');
set(c1bHs(3),   'String','Data Server (n=1) to store analysis inputs/outputs');
set(c1bHs(4),   'String','Image Server (n=1) to provide reconstructed PET & MRI files'); 
set(c1bHs(5),   'String','(Data and Image servers could be physically identical)');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Entertained combinations of Workstations/Data Server/Image Server:',     ...
    '1. Windows/Windows/Windows',   ...
    '2. Linux/Linux/Linx, and', ...
    '3. Windows using Linux terminal emulator/Linux/Linux', ...
    '* Not ready for other combinations for now (Quit)'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%

function                        local_pet_is(to_resume);
%%
ud                              = get(gcf,  'UserData');
set(gco,    'Enable','off');
%
if to_resume>0;                 path_x                      = ud.pet.is.real;
else;
    [fname, path_x]          	= uigetfile(fullfile(pwd,'*'));
    if ~ischar(fname);         	set(gco,    'Enable','on');                         return;         end;
    %
    ud.pet.is.real            	= path_x;                                                           end;
%
if path_x(1)==filesep;          s0                          = filesep;
else;                           s0                          = '';                                   end;
%
path_x(path_x==' ')             = '*';
path_x(path_x==filesep)         = ' ';
sc                              = getLseg(path_x,   [0,2]);
%
ud.stage                        = 'pet_is_1';
%
% construction of source PET path in a cell array:
for i=1:1:numel(sc);            sc{i}(sc{i}=='*')           = ' ';                                  end;
sc{1}                           = [s0,sc{1}];
%
ud.pet_is.real_c                = sc;
%
set(ud.c1bHs(1), 'String','PET Path in Image server. Select one each from right column. Then, hit >');
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
set(ud.c1bHs(numel(sc)+2),  'String','Want to re-select path? If so, restart >');
set(ud.c2bHs(numel(sc)+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''idae_done''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);']);
%
s2x                             = {'Select','fixed (= left column)','variable','Subject_ID','PET_ID'};
ud.sc_pet_is                  	= s2x;
s2e                             = {' ','Common to all subjects/scans',     ...
                                    'To vary across subjects/scans',                    ...
                                    'Subject folder, if any','Segment to pin down PET scans'};
s1{1}                           = '* Selections for path segments:';
s1c1                            = char(s2x(2:end));
s1c2                            = char(s2e(2:end));
for i=1:1:size(s1c1,1);         s1{i+1}                     = [' ',s1c1(i, :),'  : ',s1c2(i, :)];  	end;
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',s1,    'FontName','Consolas',  ...
                                'FontSize',10,  'HorizontalAlignment','left')
%
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'),    'UserData',s1);
% dispCharArrays(1,char(s2x(2:end)),3,char(s2e(2:end)));
% disp('< end list');
for i=1:1:numel(sc);
    set(ud.c1bHs(i+1),  'String',sc{i});
    set(ud.c2bHs(i+1),  'Value',1,  'Style','popupmenu',    'String',s2x);                          end;
%
set(gcf,    'UserData',ud);
return;
%%

function                        local_done(i2);
%%
ud                              = get(gcf,  'UserData');
feval(['local_done_',ud.stage], ud);
return;
%%

function                      	local_done_pet_is_1(ud);
%%
v                               = zeros(numel(ud.pet_is.real_c),    1);
v(:)                            = local_check_res(v,    ud.c2bHs);
if any(v<2);                                                                        return;         end;
%
clear ss;
for i=1:1:size(v, 1);           ss{i}                       = ' ';                                  end;
im1                             = umo_cstrs(char(ud.sc_pet_is), ['fix';'var';'Sub';'PET'],  'im1');
% converting segment selections into symbolic segments:
for i=find(v'==im1(1));         ss{i}                       = ud.pet_is.real_c{i};                  end;
for i=find(v'==im1(2));         ss{i}                       = '*';                                  end;
for i=find(v'==im1(3));         ss{i}                       = '$Subject_ID$';                       end;
ic                              = 0;
for i=find(v'==im1(4));         ic                          = ic + 1;
                                ss{i}                       = ['$PET_ID_',int2str(ic),'$'];     	end;
%
set(ud.c1bHs(1), 'String','Done!. See infoBord. Hit > to move on');
set(ud.c2bHs(1), 'String','Move on')
set(ud.c2bHs(2:end),    'Enable','off');
out                             = ss{1};
for i=2:1:numel(ss);            out                         = fullfile(out, ss{i});                 end;
ud.pet.is.symbolic              = out;
ud.pet_is.symbolic_c            = ss;
ud.stage                        = 'pet_is_2';
set(gcf,    'UserData',ud);                                                     
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
	{'* Generalized PET path for image server: ',   ...
    [' ',ud.pet.is.symbolic],'> hit ''Display'' to view more'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
return;
%%

function    out                 = local_check_res(v,c2bHs);
%%
out                             = zeros(size(v));
for i=1:1:size(v,1);
    out(i, :)                   = get(c2bHs(i+1),  'Value');
    if out(i)<2;                set(c2bHs(i+1),  'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(c2bHs(i+1),  'BackgroundColor',iv2_bgcs(0));            end;    end;
%    
if any(out<2);                  
    set(findobj(gcf, 'Tag','iv2_gen_lds_1_1'), 'BackgroundColor',iv2_bgcs(11));
 	pause(0.5);
    set(findobj(gcf, 'Tag','iv2_gen_lds_1_1'), 'BackgroundColor',iv2_bgcs(12));                     end;
return;
%%

function                     	local_done_pet_is_2(ud);
%%
for i=umo_cstrs(char(ud.pet_is.symbolic_c), '$PET', 'im1');
 	set(ud.c1bHs(i+1),  'Style','edit',  'BackgroundColor',iv2_bgcs(6));
    set(ud.c2bHs(i+1),  'Value',1,  'Style','pushbutton',   'String',ud.pet_is.real_c{i});          end;
%    
set(ud.c1bHs(1), 'String','Work on date/time elements as instructed below, if any. Hit > when done');
set(ud.c2bHs(1), 'String','Done',   'Enable','on');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
	{'* Identify PET-ID segments with date/time elements (in integers). ',   	...
    '> If identified, hit the GUI and replace integers with  yyyy or yy (year), ', 	...
    '    mm (month), dd (day), HH (hour), MM (minute), or SS (second).',            ...
    '  Leave non-date/time elements unchanged. (e.g., pet-yyyymmdd)',               ...
   	'> Leave PET_ID segments without date/time elements unchanged.',                ...
    '> Original inputs in right column for your convenience'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
ud.stage                        = 'pet_is_3';
set(gcf,    'UserData',ud);                                                     
return;
%%

function                     	local_done_pet_is_3(ud);
%%
set(findobj(gcf, 'Style','edit'),   'Style','pushbutton');
%
clear ss;
for i=1:1:numel(ud.pet_is.real_c);   
                                ss{i}                       = '*';                                  end;
%
for i=umo_cstrs(char(ud.pet_is.symbolic_c), '$PET', 'im1');
    if ~strcmpi(ud.pet_is.real_c{i}, get(ud.c1bHs(i+1), 'String'));
       	ss{i}                  	= deblank(get(ud.c1bHs(i+1), 'String'));                            end;
    set(ud.c1bHs(i+1),  'CallBack','iv2_gen_lds(''toggle'',''mri'');');
    set(ud.c2bHs(i+1),  'String',ud.pet_is.symbolic_c{i}(ud.pet_is.symbolic_c{i}~='$'));            end;
%    
ud.pet_is.format_c              = ss;
%
out                             = ss{1};
for i=2:1:numel(ss);            out                         = fullfile(out, ss{i});                 end;
%
ud.pet.is.format                = out;
ud.stage                        = 'pet_is_4';
set(ud.c1bHs(1), 'String','Hit/highlight path segment(s) as instructed below. Hit > when done');
set(ud.c2bHs(1), 'String','Approve',    'CallBack','iv2_gen_lds(''done'',[]);');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',                           ...
	{'* IDAE would like to know if it is possible to narrow down the search for',   ...
    '  PET source files with one or a few PET_ID segments (light green)',   	...
    '> Hit a left column GUI of PET-related (mostlikely the one @PET_ID_1)'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
% marking / setting segments with $PET_ID_i$: 
set(ud.c1bHs(umo_cstrs(char(ud.pet_is.symbolic_c),'$PET', 'im1')+1),    ...
    'BackgroundColor',iv2_bgcs(6),  'CallBack','iv2_gen_lds(''toggle'',''pet'');');
set(gcf,    'UserData',ud);
return;
%%

function                          	local_done_pet_is_4(ud);
%%
% just checking if 'search' paths are saved:
if ~isfield(ud.pet_is,'search_c') || ~isfield(ud.pet.is,'search');                 	return;         end;
%
set(ud.c1bHs(1), 'String','Done for PET image server! Starting PET data server');  
set(ud.c2bHs(1), 'String','Start',  'CallBack','iv2_gen_lds(''pet_ds'',0)');
% reset left & right columns:
local_reset_guis(ud);
%
ud.stage                        = 'pet_end';
set(gcf,    'UserData',ud);
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',                           ...
    {'* Instructions: ',    ...
    '> Using file navigator that pops-up, select a PET file in Data server',        ...
    '  to import a representative PET path in Data server',                         ...
    '- Easier if the same scan (chosen for image server) is selected. ',            ...
    ['  (= ',ud.pet.is.real,')']});

return;
%%

function                        local_reset_guis(ud);
%%
set(ud.c1bHs(2:end),    'String',' ',   'Value',1,  'Style','pushbutton',   ...
                                'BackgroundColor',iv2_bgcs(0),  'CallBack',' ', 'Enable','on');
set(ud.c2bHs(2:end),    'String',' ',   'Value',1,  'Style','pushbutton',   ...
                                'BackgroundColor',iv2_bgcs(0),  'CallBack',' ', 'Enable','on');
return;
%%

function                        local_pet_ds(to_resume);
%%
ud                              = get(gcf,  'UserData');
set(gco,    'Enable','off');
%
if to_resume>0;                 path_x                      = ud.pet.ds.real;
else;                           [fname, path_x]           	= uigetfile(fullfile(pwd,'*'));
    if ~ischar(fname);       	set(gco,    'Enable','on');                         return;         end;
                                ud.pet.ds.real            	= path_x;                               end;
%
if path_x(1)==filesep;          s0                          = filesep;
else;                           s0                          = '';                                   end;
%                           
path_x(path_x==' ')             = '*';
path_x(path_x==filesep)         = ' ';
sc                              = getLseg(path_x,   [0,2]);
%
% construction of source PET path in a cell array:
for i=1:1:numel(sc);            sc{i}(sc{i}=='*')           = ' ';                                  end;
sc{1}                           = [s0,sc{1}];
%
ud.pet_ds.real_c                = sc;
%
ud.stage                        = 'pet_ds_1';
%
set(ud.c1bHs(1), 'String','Assign roles / rules to individual PET path segments from left column');
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
set(ud.c1bHs(numel(sc)+2),  'String','Want to re-select path? If so, restart >');
set(ud.c2bHs(numel(sc)+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''pet_is_4''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'], 'Enable','on');
%
s2x                             = {'Select','fixed (=left column)'};
cc                              = char(ud.pet_is.symbolic_c);
q                               = ones(1, numel(ud.pet_is.symbolic_c));
q(:, cc(:,1)=='$')              = 0;
ic                              = numel(s2x);
for i=find(q>0);                  
    ic                          = ic + 1;
  	s2x{ic}                     = ['= segment_',int2str(i),' (',ud.pet_is.real_c{i},')'];         	end;
%
pet                             = find(cc(:,1)=='$' & cc(:,2)=='P');
for i=find(cc(:,1)'=='$');   
    ic                          = ic + 1;
    s2x{ic}                     = ['= ',ud.pet_is.symbolic_c{i}(ud.pet_is.symbolic_c{i}~='$'),   ...
                                                            ' (',ud.pet_is.real_c{i},')'];
    ic                          = ic + 1;
    s2x{ic}                     = ['= ',ud.pet_is.symbolic_c{i}(ud.pet_is.symbolic_c{i}~='$'),   ...
                                                            ' - modify'];                           end;
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Selections for path segments (image server examples in ()):',   ...
    ' fixed (=left column)  :  Common to all subjects/scans',           ...
    ' = segment_i           :  Copy segment_i of image server path ',           	... 
    ' = segment_ID          :  Copy Subject_ID / PET_ID_i of image server path ',  	...
    ' = segment_ID - modify :  Copy Segment_ID with user-defined modifications'},   ...
                                'FontName','Consolas',  'FontSize',10,  'HorizontalAlignment','left');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
ud.sc_pet_ds                 	= s2x;
for i=1:1:numel(sc);
    set(ud.c1bHs(i+1),  'String',sc{i});
    set(ud.c2bHs(i+1),  'Value',1,  'Style','popupmenu',    'String',s2x,   'Enable','on');         end;
%
set(gcf,    'UserData',ud);
return;
%%

function                     	local_done_pet_ds_1(ud);
%%
% unmarking problematic segments, if any:
set(ud.c2bHs,   'BackgroundColor',iv2_bgcs(0));     
%
v                               = zeros(numel(ud.pet_ds.real_c), 1);
v(:)                            = local_check_res(v, ud.c2bHs);
if any(v<2);                                                                        return;         end;
% expected PET_ID_i from the image server:
ic                              = 0;
for i=umo_cstrs(char(ud.sc_pet_ds), '= PET', 'im1');
    ic                          = ic + 1;
    pet_seg{ic}                 = ['= ',getLseg(ud.sc_pet_ds{i}, 2)];                            	end;
% PET_ID_i may not be applied to more than one PET_ID segment
cm1                             = umo_cstrs(char(pet_seg), [], 'cm1');
im1                             = umo_cstrs(char(ud.sc_pet_ds(v)), char(pet_seg(cm1(:,2)>0)), 'im1');
if size(im1,2)>1;   
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
        {'* Problem! PET_ID_i may not be applied to more than one PET_ID segment (pink)',  	...
        '> Correct them and hit ''Done'' one more time'});                          
    for i=find(im1(:,2)>0)';
        set(ud.c2bHs(im1(i, im1(i,:)>0)+1),  'BackgroundColor',iv2_bgcs(11));                       end;
                                                                                    return;         end;
% no PET_ID_i was selected:
if sum(im1>0)<1;
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
        {'* Problem! Select at least one PET_ID_i for PET_ID segments (currently none)',  	...
        '> Correct them and hit ''Done'' one more time'});                          return;         end;
%

% the last string for '= PET' is '= PET_ID_last - fit': 
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
    {'* IDAE requests to make PET paths unique to IDAE',           	...
    '  to avoid over-writing of files by IDAE by chance ',        	...
    '> Assign ''- modify'' to at least one PET_ID segment',         ...
    ['  Preferably, to the last PET_ID segment (i.e., ''',getLseg(ud.sc_pet_ds{end},2),''')']});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
set(ud.c1bHs(1), 'String','Review / modify selections per instructions below as needed. Then hit > ');
set(ud.c2bHs(1), 'String','Done');
ud.stage                        = 'pet_ds_2';
set(gcf,    'UserData',ud);
return;
%%

function                     	local_done_pet_ds_2(ud);
%%
v                               = zeros(numel(ud.pet_ds.real_c), 1);
v(:)                            = local_check_res(v, ud.c2bHs);
if any(v<2);                                                                        return;         end;
% % expected PET_ID_i from the image server:
% ic                              = 0;
% for i=umo_cstrs(char(ud.sc_pet_ds), '= PET', 'im1');
%     ic                          = ic + 1;
%     pet_seg{ic}                 = ['= ',getLseg(ud.sc_pet_ds{i}, 2)];                            	end;
% % 
% cm1                             = umo_cstrs(char(pet_seg), [], 'cm1');
set(ud.c1bHs([1:1:size(v,1)]+1),    'Style','pushbutton',   'BackgroundColor',iv2_bgcs(0));

mmm                             = zeros(size(v));
for i=umo_cstrs(char(ud.sc_pet_ds(v)), '= PET', 'im1');
    if size(getLseg(ud.sc_pet_ds{v(i)},3),2)==1;
        mmm(i, :)             	= getLseg(ud.sc_pet_ds{v(i)},3)=='-';                       end;    end; 
%
if sum(mmm)<1;
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
        {'* Not right! ',  '> Assign ''- modify'' to at least one PET_ID segment',         ...
        ['  Preferably, to the last PET_ID segment (i.e., ''',getLseg(ud.sc_pet_ds{end},2),''')']});
                                                                                    return;         end;
%
set(ud.c2bHs([1:1:size(v,1)]+1),    'Enable','off');
drawnow;
%
s1                              = {' ','Fix it to a new name across subjects / scans', 	...
                                    ' ','More complex - define rules later',' '};
                                    
for i=find(mmm'>0);
    s1{1}                       = ['Current: ',ud.pet_ds.real_c{i}];
    s1{3}                       = ['Append a string to ',getLseg(ud.sc_pet_ds{v(i)},2), ...
                                                            ' from image server'];
    s1{5}                       = ['Copy ',getLseg(ud.sc_pet_ds{v(i)},2),' (loophole)'];
    set(ud.c1bHs(i+1),  'Value',1,  'Style','popupmenu',    'String',s1,    ...
                                'CallBack','iv2_gen_lds(''repl_seg'',[]);');                        end;
%
ud.stage                        = 'pet_ds_3';
set(gcf,    'UserData',ud);
%
set(ud.c1bHs(1), 'String','Work on popup-menu GUIs on left column as instructed below.');
set(ud.c2bHs(1), 'String','Done');
%
set(ud.c1bHs(size(v,1)+3), 'String','Want to re-select segments to modify? if yes, hit >');
set(ud.c2bHs(size(v,1)+3), 'String','Re-select',   'CallBack',['ud=get(gcf,''UserData'');', ...
    'ud.stage=''pet_ds_2''; set(gcf,''UserData'',ud); set(ud.c2bHs,''Enable'',''on'');']);
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Define modification rules', '> Select a desired rule from pulldown menu',   ...
    '> Then, enter a string of your chice (no spaces)', '  When all are done hit ''Done'' @top right'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%

function                     	local_done_pet_ds_3(ud);
%%
%
set(findobj(gcf, 'Style','edit'),   'Style','pushbutton');
im1                             = umo_cstrs(char(get(ud.c1bHs, 'Style')), 'popupmenu', 'im1');
if im1(1)>0;                    
    set(ud.c1bHs(im1),  'BackgroundColor',iv2_bgcs(11));
   	pause(0.5);
 	set(ud.c1bHs(im1),  'BackgroundColor',iv2_bgcs(0));                             return;         end;
%
v                               = zeros(numel(ud.pet_ds.real_c), 1);
v(:)                            = local_check_res(v, ud.c2bHs);
if any(v<2);                                                                        return;         end;
set(ud.c2bHs([1:size(v,1)+2]+1),  	'enable','off');

clear ss;
for i=1:1:size(v,1);            ss{i}                       = ' ';                                  end;
%
imx                             = umo_cstrs(char(ud.sc_pet_ds(v)), ...
                                  	['fixed'; '= seg'; '= Sub'; '= PET'], 'im1')
% coping with the cases where ( is still left on the GUI:
ok                              = 1;
for i=imx(4, imx(4, :)>0); 
    sc1                         = getLseg(get(ud.c1bHs(i+1), 'String'), 1);
    if sc1(1)=='(';             set(ud.c1bHs(i+1),  'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(ud.c1bHs(i+1),  'BackgroundColor',iv2_bgcs(6));
                                ok                          = 0;                            end;    end;
if ok<1;                                                                            return;         end;                                
% fixed segments:
for i=imx(1, imx(1, :)>0);      ss{i}                       = ud.pet_ds.real_c{i};                 	end;
% segment_i:
for i=imx(2, imx(2, :)>0);      
    ss{i}                       = ['$',getLseg(ud.sc_pet_ds{v(i)},2),'$'];                          end;
% subject_ID segment:
for i=imx(3, imx(3, :)>0); 
    sc3                         = getLseg(ud.sc_pet_ds{v(i)}, 3);
    if sc3(1)=='-';
        ss{i}                   = ['$',getLseg(ud.sc_pet_ds{v(i)}, 2),'-modify$'];
    else;
        ss{i}                   = ['$',getLseg(ud.sc_pet_ds{v(i)}, 2),'$'];                 end;    end;
% PET_ID_i segment:
for i=imx(4, imx(4, :)>0); 
    sc3                         = getLseg(ud.sc_pet_ds{v(i)}, 3);
    if sc3(1)=='-';
        if get(ud.c1bHs(i+1), 'UserData')==2;
            sc1                 = get(ud.c1bHs(i+1), 'String');
            ss{i}               = sc1(sc1~=' ');
        elseif get(ud.c1bHs(i+1), 'UserData')==3;
            sc1                 = get(ud.c1bHs(i+1), 'String');
            ss{i}               = ['$',getLseg(ud.sc_pet_ds{v(i)}, 2),'-modify$',sc1(sc1~=' ')];
        elseif get(ud.c1bHs(i+1), 'UserData')==4;
            ss{i}               = ['$',getLseg(ud.sc_pet_ds{v(i)}, 2),'-modify$']; 
        else;
            ss{i}               = ['$',getLseg(ud.sc_pet_ds{v(i)}, 2),'$'];                         end;
    else;
         ss{i}                  = ['$',getLseg(ud.sc_pet_ds{v(i)}, 2),'$'];                 end;    end;
%
ud.pet_ds.symbolic_c            = ss;
out                             = ss{1};
for i=2:1:numel(ss);            out                         = fullfile(out, ss{i});                 end;
ud.pet.ds.symbolic              = out;
ud.stage                        = 'pet_done';
set(gcf,    'UserData',ud);
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Generalized PET path in data server: ',    ...
    ['   ',ud.pet.ds.symbolic],  ...
    ' Segments without $   :   fixed to them across subjects / scans',     ...
    ' $whatever_ID$        :   to copy $whatever_ID$ from image server',   ...
    ' $whatever_ID-modify$ :   to copy $whatever_ID$ after modifications', ...
    ' $whatever_ID$xxx     :   to copy $whatever_ID$ with xxx appended'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
set(ud.c1bHs(1), 'String','Done for PET paths. Observe generalized path in Data server below');
set(ud.c2bHs(1), 'String','Start MRI');
return;
%%

function                     	local_done_pet_done(ud);
%%
% reset left & right columns:
local_reset_guis(ud);
%
set(ud.c1bHs(1), 'String','Answer below question and Hit > to select a MRI souce file');
set(ud.c2bHs(1), 'String','Start',  'CallBack','iv2_gen_lds(''mri_is'',0);');
%
set(ud.c1bHs(2), 'String','Confirm if MRI source files in DICOM or Phillips PET/REC format');
set(ud.c2bHs(2),   'Value',1,  'Style','popupmenu',    'Enable','on',  ...
                                'String',{'Select','Yes','No, not always'});
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Reassurance:',  ' Image server: to provide reconstructed MRI files',         ...
    ' Data server: to hold input / output files for analyses'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%%

function                            local_disp(i2);
%%
ud                              = get(gcf,  'UserData');
disp(['> output dxetc4xxx: ',ud.lds]);
if isfield(ud,'pet');   
    if isfield(ud.pet,'is');    disp('> PET Path variables - Image server:');
                                disp(ud.pet.is);                                                  	end;
    if isfield(ud.pet,'ds');    disp('> PET Path variables - Data server:');
                                disp(ud.pet.ds);                                          	end;    end;
%
if isfield(ud,'mri');   
    if isfield(ud.mri,'is');    disp('> MRI Path variables - Image server:');
                                disp(ud.mri.is);                                                  	end;
    if isfield(ud.mri,'ds');    disp('> MRI Path variables - Data server:');
                                disp(ud.mri.ds);                                          	end;    end;
%
if isfield(ud,'idae');          disp('> IDAE-related Paths');
                                disp(ud.idae);                                                      end;
if isfield(ud,'etc');          	disp('> Other critical settings');
                                disp(ud.etc);                                                       end;
return;
%%

function                     	local_toggle(i2);
%%
set(gco,    'Enable','off');
drawnow;
if sum((get(gco, 'BackgroundColor') - iv2_bgcs(6)).^2)<10^-6;
                            	set(gco, 'BackgroundColor',iv2_bgcs(18));
else;                        	set(gco, 'BackgroundColor',iv2_bgcs(6), 'Enable','on'); 
                                                                                    return;         end;
%
ud                              = get(gcf, 'UserData');
bgcs                            = cell2mat(get(ud.c1bHs, 'BackgroundColor'));
bg18                            = iv2_bgcs(18);
ii                              = find(sum((bgcs - bg18(ones(size(bgcs,1),1), :)).^2,2)<10^-6)-1;
n                               = length(ii);
if n<1;                                                                             return;         end;
if n>2;                         set(gco, 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(gco, 'BackgroundColor',iv2_bgcs(6));            return;         end;
%                            
eval(['symbolic_c               = ud.',i2,'_is.symbolic_c;']);
eval(['real_c                   = ud.',i2,'_is.real_c;']);
%
jj                              = zeros(2,  max(ii));
jj(1,   ii)                     = 1;
jj(2,   umo_cstrs(char(symbolic_c(1:max(ii))), '$', 'im1'))  	= 1;
ss                            	= symbolic_c(1:max(ii));
s2                              = ss;
for i=find(jj(1,:)<1 & jj(2,:)>0);
                                s2{i}                       = '*';
                                ss{i}                       = '*';                                  end;
for i=find(jj(1,:)>0 & jj(2,:)>0);
                                s2{i}                       = symbolic_c{i};
                                ss{i}                       = real_c{i};                            end;
%
eval(['ud.',i2,'_is.search_c    = s2;']);
% ud.pet_is.search_c            	= s2;
ssc                             = ss{1};
s2c                             = s2{1};
for i=2:1:numel(s2);            ssc                         = fullfile(ssc, ss{i});
                                s2c                         = fullfile(s2c, s2{i});                 end;
% ud.pet.is.search                = s2c;
eval(['ud.',i2,'.is.search      = s2c;']);
%
d                               = dir(fullfile(ssc, '*'));
cm1                             = umo_cstrs(char(d.folder), [], 'cm1');
%
s1                              = {['* current search criteria: ',ssc],     ...    
                                    ['> ',int2str(sum(cm1(:,2)>0)),' directories to narrow down ',  ...
                                    'in search of ',upper(i2),' source files']};
set(gco,    'Enable','on');
if sum(cm1(:,2)>0)<10;          s1{end+1}                   = '< Good enough! Ok to Approve';
else;  
    if n<2;                     s1{end+1}                   = '< not good. Add a segment';
    else;                       
        s1{end+1}            	= '< not good. Erase one and add another, if any';
        s1{end+1}               = '  Or accept current criteria as is';                     end;    end;
%
set(gcf,    'UserData',ud);
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',s1); 
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%

function                            local_help(i2);
%%
coud                                = get(gco,  'Userdata');
if isempty(get(gco, 'Userdata'));                                                 	return;         end;
disp(char(get(gco, 'Userdata')));
return;
%%


function                        local_save(i2);
%%
ud                              = get(gcf,  'UserData');
save(fullfile(fileparts(which(mfilename)), [ud.lds,'.mat']),    'ud');
disp(['> gatherd information saved .. ',10, ...
                                ' file : ',fullfile(fileparts(which(mfilename)), [ud.lds,'.mat'])]);
return;
%%

function                        local_resume(ud);
%%

s1{1}                           = 'Resume from (select one):';
cb{1}                           = ' ';
%
im1                             = umo_cstrs(char(fieldnames(ud)), ...
                                    char('idae','pet_is','pet_ds','mri_is','mri_ds','etc'), 'im1');
% idae has to be present - otherwise start from scratch:
if im1(1)<1;                    local_done_idae_0(ud);                              return;         end;
%
s1{end+1}                       = '- Selection of your home directory in Data Server';
cb{end+1}                       = ['ud.stage=''idae_1''; ud.reuse.idae=0;',     ...
                                    'set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'];
%
s1{end+1}                       = '- Recoded your home directory in Data Server';
cb{end+1}                       = ['ud.stage=''idae_1''; ud.reuse.idae=1;',     ...
                                    'set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'];
% pet_is is present:
if im1(2)>0;
    s1{end+1}                 	= '- Selection of a representative PET path in Image Server';
    cb{end+1}                	= 'iv2_gen_lds(''pet_is'',0);';                                     
    %
  	s1{end+1}                 	= '- Recorded PET path in Image Server';
    cb{end+1}                	= 'iv2_gen_lds(''pet_is'',1);';                                     end;
% pet_ds is present:
if im1(3)>0;
    s1{end+1}                 	= '- Selection of a representative PET path in Data Server';
    cb{end+1}                	= 'iv2_gen_lds(''pet_ds'',0);';                                     
    %
  	s1{end+1}                 	= '- Recorded PET path in Data Server';
    cb{end+1}                	= 'iv2_gen_lds(''pet_ds'',1);';                                     end;
% mri_is is present:
if im1(4)>0;
    s1{end+1}                 	= '- Selection of a representative MRI path in Image Server';
    cb{end+1}                	= ['ud.stage=''pet_done''; ',    ...
                                    'set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'];
    %
  	s1{end+1}                 	= '- Recorded MRI path in Image Server';
    cb{end+1}                	= 'iv2_gen_lds(''mri_is'',1);';                                     end;
% mri_ds is present:
if im1(5)>0;
    s1{end+1}                 	= '- Selection of a representative MRI path in Data Server';
    cb{end+1}                	= 'iv2_gen_lds(''mri_ds'',0);';                                     
    %
  	s1{end+1}                 	= '- Recorded MRI path in Data Server';
    cb{end+1}                	= 'iv2_gen_lds(''mri_ds'',1);';                                     end;
% etc is recorded:
if im1(6)>0;
    s1{end+1}                 	= '- The last section (SPM12 etc)';
    cb{end+1}                	= ['ud.stage=''mri_ds_done''; ',    ...
                                    'set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'];         end;
%
set(ud.c1bHs(end),  'Value',1,  'Style','popupmenu',    'String',s1,    'UserData',cb,  ...
                                'CallBack','iv2_gen_lds(''resume_doit'',[]);');
return;
%%

function                        local_resume_doit(i2);
%%
v                               = get(gco,  'Value');
if v<2;                                                                             return;         end;
%
ud                              = get(gcf,  'UserData');
cb                              = get(gco,  'UserData');
eval(cb{v});
return;
%%

%% MRI_is

function                        local_mri_is(to_resume);
%%
ud                              = get(gcf,  'UserData');
set(gco,    'Enable','off');
%
if to_resume>0;                 path_x                      = ud.mri.is.real;
else;
    % making sure that MRI source files are in DICOM or PAR/REC format:   
    v                         	= 0;
    v(:)                        = local_check_res(v, ud.c2bHs);
    if any(v<2);                                                                    return;         end;
    if v(1)>2;
        set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
            {'* Problem! IDAE cannot hundle other formats than DICOM or PAR/REC format', 	...
            '> Save current variables (Hit ''Save'') and exit',     ...
            '> Report the problem to IDAE team (See the IDAE site in GitHub)'});  return;         end;
    %
    [fname, path_x]           	= uigetfile(fullfile(ud.pet.is.real,'*'));
    if ~ischar(fname);       	set(gco,    'Enable','on');                         return;         end;
                                ud.mri.is.real            	= path_x;                               end;
%
if path_x(1)==filesep;          s0                          = filesep;
else;                           s0                          = '';                                   end;
%                           
path_x(path_x==' ')             = '*';
path_x(path_x==filesep)         = ' ';
sc                              = getLseg(path_x,   [0,2]);
%
% construction of source PET path in a cell array:
for i=1:1:numel(sc);            sc{i}(sc{i}=='*')           = ' ';                                  end;
sc{1}                           = [s0,sc{1}];
%
ud.mri_is.real_c                = sc;
%
ud.stage                        = 'mri_is_1';
%
set(ud.c1bHs(1), 'String','Assign roles / rules to individual MRI path setments from right menu');
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
set(ud.c1bHs(numel(sc)+2),  'String','Want to re-select example MRI path? If yes, hit >');
set(ud.c2bHs(numel(sc)+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''pet_done''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);']);
%
s2x                             = {'Select','fixed (=left column)','variable',  ...
                                                            'Subject_ID','MRI_ID','Series_ID'};
ud.sc_mri_is                  	= s2x;
s2e                             = {' ','Common to all subjects/scans',     ...
                                    'To vary across subjects/scans',                    ...
                                    'Subject folder, if any','Segments specific to MRI scans',  ...
                                    'MRI series folder (special MRI folder), if any'};
s1{1}                           = '* Selections for path segments:';
s1c1                            = char(s2x(2:end));
s1c2                            = char(s2e(2:end));
for i=1:1:size(s1c1,1);         s1{i+1}                     = [' ',s1c1(i, :),'  : ',s1c2(i, :)];  	end;
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',s1,    'FontName','Consolas',  ...
                                'FontSize',10,  'HorizontalAlignment','left')
%
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'),    'UserData',s1);
% dispCharArrays(1,char(s2x(2:end)),3,char(s2e(2:end)));
% disp('< end list');
for i=1:1:numel(sc);
    set(ud.c1bHs(i+1),  'String',sc{i});
    set(ud.c2bHs(i+1),  'Value',1,  'Style','popupmenu',    'String',s2x,   'Enable','on');      	end;
%
set(gcf,    'UserData',ud);
return;
%%

function                        local_done_mri_is_1(ud);
%%
v                               = zeros(numel(ud.mri_is.real_c),    1);
v(:)                            = local_check_res(v, ud.c2bHs);
if any(v<2);                                                                        return;         end;
%
im1                             = umo_cstrs(char(ud.sc_mri_is), ['fix';'var';'Sub';'MRI';'Ser'], 'im1');
%
if ~any(v'==im1(4));
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
        {'* Problem! IDAE cannot hundle cases without MRI_ID segments',         ...
        '  Review segments (or folders listed in left column) carefully',       ...
        '  See the manual for further information'});                            	return;         end;
%
clear ss;
for i=1:1:size(v,1);            ss{i}                       = ' ';                                  end;
%
for i=find(v'==im1(1));         ss{i}                       = ud.mri_is.real_c{i};                  end;
for i=find(v'==im1(2));         ss{i}                       = '*';                                  end;
for i=find(v'==im1(3));         ss{i}                       = '$Subject_ID$';                     	end;
ic                              = 0;
for i=find(v'==im1(4));         ic                          = ic + 1;
                                ss{i}                       = ['$MRI_ID_',int2str(ic),'$'];         end;
for i=find(v'==im1(5));         ss{i}                       = '$Series_ID$';                       	end;                   
%
ud.mri_is.symbolic_c            = ss;
%
ddd                             = ss{1};
for i=2:1:numel(ss);            ddd                         = fullfile(ddd, ss{i});                 end;
ud.mri.is.symbolic              = ddd;
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
	{'* Generalized MRI path for image server:',[' ',ddd],'> hit ''Display'' to view more'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
ud.stage                        = 'mri_is_2';
set(gcf,    'UserData',ud);
set(ud.c1bHs(1), 'String','Done! Review info board below. Hit > to move on');
set(ud.c2bHs(1), 'String','Move on')
set(ud.c2bHs(2:end),    'Enable','off');
return;
%%

function                        local_done_mri_is_2(ud);
%
im1                             = umo_cstrs(char(ud.mri_is.symbolic_c), ['$MRI';'$Ser'], 'im1');
%
for i=sort(im1(im1>0))';
 	set(ud.c1bHs(i+1),  'Style','edit',  'BackgroundColor',iv2_bgcs(6));
    set(ud.c2bHs(i+1),  'Value',1,  'Style','pushbutton',   'String',ud.mri_is.real_c{i});          end;
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
	{'* Identify MRI-related segments with date/time elements (in integers). ',   	...
    '> If identified, hit the GUI and replace integers with  yyyy or yy (year), ', 	...
    '    mm (month), dd (day), HH (hour), MM (minute), or SS (second).',            ...
    '    Leave non-date/time elements unchanged. (e.g., yyyy-mm-dd_MRI)',           ...
   	'  Leave PET_ID segments without date/time elements unchanged.'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
% 
set(ud.c1bHs(1), 'String','Work on date/time elements as instructed below, if any. Then, hit >');
set(ud.c2bHs(1), 'String','Done');
ud.stage                        = 'mri_is_3';
set(gcf,    'UserData',ud);
return;
%%

function                        local_done_mri_is_3(ud);
%%
set(findobj(gcf, 'Style','edit'),   'Style','pushbutton');
%
clear ss;
for i=1:1:numel(ud.mri_is.real_c);   
                                ss{i}                       = '*';                                  end;
%
im1                             = umo_cstrs(char(ud.mri_is.symbolic_c), ['$MRI';'$Ser'], 'im1');
for i=sort(im1(im1>0))';
    if ~strcmpi(ud.mri_is.real_c{i}, get(ud.c1bHs(i+1), 'String'));
       	ss{i}                  	= deblank(get(ud.c1bHs(i+1), 'String'));                            end;
    set(ud.c1bHs(i+1),  'CallBack','iv2_gen_lds(''toggle'',''mri'');');
    set(ud.c2bHs(i+1),  'String',ud.mri_is.symbolic_c{i}(ud.mri_is.symbolic_c{i}~='$'));            end;
%    
ud.mri_is.format_c              = ss;
%
out                             = ss{1};
for i=2:1:numel(ss);            out                         = fullfile(out, ss{i});                 end;
%
ud.mri.is.format                = out;
ud.stage                        = 'mri_is_done';
set(ud.c1bHs(1), 'String','Hit/highlight path segment(s) as instructed below. Hit > when done');
set(ud.c2bHs(1), 'String','Approve',    'CallBack','iv2_gen_lds(''done'',[]);');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',                           ...
	{'* IDAE would like to know if it is possible to narrow down the search for',   ...
    '  MRI source files with one or a few MRI-related segments (light green)',   	...
    '> Hit a left column GUI of MRI-related (start from @MRI_ID_1)'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
% 
set(gcf,    'UserData',ud);
return;
%%

function                        local_done_mri_is_done(ud);
%%
% just checking if 'search' paths are saved:
if ~isfield(ud.mri_is,'search_c') || ~isfield(ud.mri.is,'search');                 	return;         end;
%
set(ud.c1bHs(1), 'String','Done for MRI image server! Let''s work on MRI data server');  
set(ud.c2bHs(1), 'String','Start',  'CallBack','iv2_gen_lds(''mri_ds'',0)');
% reset left & right columns:
local_reset_guis(ud);
%
ud.stage                        = 'mri_is_end';
set(gcf,    'UserData',ud);
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',                                   ...
    {'* Instructions: ',    ...
    '> Using file navigator that pops-up, select a MRI file in data server',    ...
    '  to get a representative MRI path to work with', 	...
  	'- Easier if the same scan (chosen for image server) is selected. ',               ...
  	['  (',ud.mri.is.real,')']});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%


function                        local_mri_ds(to_resume);
%%
ud                              = get(gcf,  'UserData');
set(gco,    'Enable','off');
%
if to_resume<1;                 [fname, ud.mri.ds.real]   	= uigetfile(fullfile(ud.pet.ds.real,'*'));
    if ~ischar(fname);       	set(gco,    'Enable','on');                         return;         end;
                                                                                                    end;
%
sc                              = local_path2sc(ud.mri.ds.real);
ud.mri_ds.real_c                = sc;
%
set(ud.c1bHs(1), 'String','Assign definition rules to individual MRI Path segments from right column') 
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on');
%
set(ud.c1bHs(numel(sc)+2),  'String','Want to reselect example path? If yes, restart >');
set(ud.c2bHs(numel(sc)+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''mri_is_done''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'], 'Enable','on');
%
s2x                             = {'Select','fixed (=left column)'};
cc                              = char(ud.mri_is.symbolic_c);
q                               = ones(1, numel(ud.mri_is.symbolic_c));
q(:, cc(:,1)=='$')              = 0;
ic                              = numel(s2x);
for i=find(q>0);                  
    ic                          = ic + 1;
  	s2x{ic}                     = ['= segment_',int2str(i),' (',ud.mri_is.real_c{i},')'];         	end;
%
for i=find(cc(:,1)'=='$');   
    ic                          = ic + 1;
    s2x{ic}                     = ['= ',ud.mri_is.symbolic_c{i}(ud.mri_is.symbolic_c{i}~='$'),   ...
                                                            ' (',ud.mri_is.real_c{i},')'];
    ic                          = ic + 1;
    s2x{ic}                     = ['= ',ud.mri_is.symbolic_c{i}(ud.mri_is.symbolic_c{i}~='$'),   ...
                                                            ' - modify'];                           end;
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Definitions of the selections (image server examples in ()):',   ...
    ' fixed (=left column)  :  Common to all subjects / scans',        	...
    ' = segment_i           :  Copy segment_i of image server path ',               ... 
    ' = segment_ID          :  Copy Subject/MRI/Series_ID of image server path ',  	...
    ' = segment_ID - modify :  Copy Segment_ID after user-defined modifications'},  ...
                                'FontName','Consolas',  'FontSize',10,  'HorizontalAlignment','left');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
ud.sc_mri_ds                 	= s2x;
for i=1:1:numel(sc);
    set(ud.c1bHs(i+1),  'String',sc{i});
    set(ud.c2bHs(i+1),  'Value',1,  'Style','popupmenu',    'String',s2x,   'Enable','on');         end;
%
ud.stage                        = 'mri_ds_1';
set(gcf,    'UserData',ud);
return;
%%

function                     	local_done_mri_ds_1(ud);
%%
% unmarking problematic segments, if any:
set(ud.c2bHs,   'BackgroundColor',iv2_bgcs(0));     
%
v                               = zeros(numel(ud.mri_ds.real_c), 1);
v(:)                            = local_check_res(v, ud.c2bHs);
if any(v<2);                                                                        return;         end;
% checking if MRI_ID_i / Series_ID are properly selected:
im1                             = local_check_mri_ids(ud,v);
if isempty(im1);                                                                    return;         end;
% 
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
    {'* IDAE requests to make MRI paths unique to IDAE',           	...
    '  to avoid over-writing of files by IDAE by chance ',        	...
    '> Assign ''- modify'' to at least one of MRI_ID_i or Series_ID', ...
    '  to go with user-defined rules (to specify later)'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
set(ud.c1bHs(1), 'String','Review / modify selections per instructions below as needed. Then, hit > ');
set(ud.c2bHs(1), 'String','Done');
ud.stage                        = 'mri_ds_2';
set(gcf,    'UserData',ud);
return;
%%

function    im1                 = local_check_mri_ids(ud,v);
%%
% MRI_ID_i may not be applied to more than one MRI_ID segment
im1                             = umo_cstrs(char(ud.sc_mri_ds(v)), ['= MRI';'= Ser'], 'im1');
for i=1:1:size(v,1);            mri_seg{i}                  = int2str(i);                           end;
for i=im1(im1>0)';              mri_seg{i}                  = getLseg(ud.sc_mri_ds{v(i)}, 2);    	end;
cm1                             = umo_cstrs(char(mri_seg), [], 'cm1');
if any(cm1(:,2)<1);   
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
        {'* MRI_ID_i / Series_ID may be assigned to one and only one MRI_ID segment',   ...
        '> Correct them and hit ''Done'' one more time'});                          
    for i=find(cm1(:,2)>1)';
        for j=find(cm1(:,1)==cm1(i,1))';
            set(ud.c2bHs(j+1),  'BackgroundColor',iv2_bgcs(11));
            pause(0.5);
            set(ud.c2bHs(j+1),  'BackgroundColor',iv2_bgcs(0));                             end;    end;
    im1                         = [];                                               return;         end;
%    
if im1(1,1)<1;
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
        {'* Select at least one MRI_ID_i for MRI_ID segments (currently none)',         ...
        '> Correctselections and hit ''Done'' (top right) one more time'});                          
    im1                         = [];                                               return;         end;
return;
%%

function                     	local_done_mri_ds_2(ud);
%%
v                               = zeros(numel(ud.mri_ds.real_c), 1);
v(:)                            = local_check_res(v, ud.c2bHs);
if any(v<2);                                                                        return;         end;
% 
% checking if MRI_ID_i / Series_ID are properly selected:
im1                             = umo_cstrs(char(ud.sc_mri_ds(v)), ['= MRI';'= Ser'], 'im1');
v2                              = zeros(size(v));
v2(im1(im1>0))                  = 1;
for i=im1(im1>0)';              s3                          = getLseg(ud.sc_mri_ds{v(i)},3);
                                v2(i,   :)                  = s3(1)=='-';                           end;
% - modify not found @MRI_ID_i / Series_ID:
if sum(v2)<1;
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',               ...
        {'* Assign ''- modify'' to at least one MRI_ID_i or Series_ID',    	...
        '  to make PET paths in Data server unique to IDAE',                ...
        '> Correct selections and hit ''Done'' (top right) one more time'});
    % blinking:
    for i=im1(im1>0)';          set(ud.c2bHs(i+1),  'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(ud.c2bHs(i+1),  'BackgroundColor',iv2_bgcs(0));                 end;
                                                                                    return;         end;
%
for i=1:1:size(v,1);
    set(ud.c1bHs(i+1),  'Value',1,  'Style','pushbutton',   'String',ud.mri_ds.real_c{i},   ...
                                                            'BackgroundColor',iv2_bgcs(0));         end;
set(ud.c2bHs([1:1:size(v,1)]+1),    'Enable','off');
drawnow;
%
s1                              = {' ','Fix it to a new name across subjects / scans', 	...
                                    ' ','More complex - define rules later',' '};
for i=find(v2'>0);
    s1{1}                       = ['Current: ',ud.mri_ds.real_c{i}];
    s1{3}                       = ['Append a string to ',getLseg(ud.sc_mri_ds{v(i)},2), ...
                                                            ' from image server'];
    s1{5}                       = ['Copy ',getLseg(ud.sc_mri_ds{v(i)},2),               ...
                                                            ' from image server (loophole)'];
    set(ud.c1bHs(i+1),  'Value',1,  'Style','popupmenu',    'String',s1,    ...
                                'CallBack','iv2_gen_lds(''repl_seg'',[]);');                        end;
%
set(ud.c1bHs(1), 'String','Work on MRI segments with pop-up menu (left column)');
set(ud.c2bHs(1), 'String','Done',   'CallBack','iv2_gen_lds(''done'',[]);');
%
set(ud.c1bHs(size(v,1)+3), 'String','Want to re-select segments to modify? if yes, hit >');
set(ud.c2bHs(size(v,1)+3), 'String','Re-select',   'CallBack',['ud=get(gcf,''UserData'');', ...
    'ud.stage=''mri_ds_1''; set(gcf,''UserData'',ud); set(ud.c2bHs,''Enable'',''on'');']);
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Define modification rules', '> Select a desired rule from pulldown menu',   ...
    '> Then, enter the input', '  When all are done hit ''Done'' @top right'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
ud.stage                        = 'mri_ds_3';
set(gcf,    'UserData',ud);
return;
%%

function                        local_repl_seg(i2);
%%
if get(gco, 'Value')<2;                                                             return;         end;
set(gco,    'UserData', get(gco, 'Value'));
if get(gco, 'Value')==2;
    set(gco,    'Value',1,  'Style','edit',     'String','(enter new folder name)',     ...
                                'CallBack',' ',     'BackgroundColor',iv2_bgcs(6));
elseif get(gco, 'Value')==3;
    set(gco,    'Value',1,  'Style','edit',     'String','(enter string to append)',    ...
                                'CallBack',' ',     'BackgroundColor',iv2_bgcs(6));
elseif get(gco, 'Value')==4;
    set(gco,    'Value',1,  'Style','pushbutton',   'String','[set rules later]',       ...
                                'CallBack',' ',     'BackgroundColor',iv2_bgcs(6));  
else;
    set(gco,    'Value',1,  'Style','pushbutton',   'String','[= right column; no -modify]',    ...
                                'CallBack',' ',     'BackgroundColor',iv2_bgcs(6));                 end;
return;
%%

function                     	local_done_mri_ds_3(ud);
%%
v                               = zeros(numel(ud.mri_ds.real_c), 1);
v(:)                            = local_check_res(v, ud.c2bHs);
if any(v<2);                                                                        return;         end;
%
% disp(char(ud.sc_mri_ds(v(im1(im1>0)))))
for i=1:1:size(v,1);            ss{i}                       = ' ';                                  end;
imx                             = umo_cstrs(char(ud.sc_mri_ds(v)),    ...
                                  	['fixed';'= seg';'= Sub';'= MRI';'= Ser'],  'im1');
%
for i=1:1:size(v,1);            ss{i}                       = ' ';                                  end;
% fixed:
for i=imx(1, imx(1, :)>0);      ss{i}                       = ud.mri_ds.real_c{i};                  end;
% segments from image server:
for i=imx(2, imx(2, :)>0);      
    ss{i}                       = ['$',getLseg(ud.sc_mri_ds{v(i)},2),'$'];                          end;
% Subject_ID field:
for i=imx(3, imx(3, :)>0);
    s3                          = getLseg(ud.sc_mri_ds{v(i)}, 3);
    if s3(1)=='-';              
        ss{i}                 	= ['$',getLseg(ud.sc_mri_ds{v(i)},2),'-modify$'];   
    else;
         ss{i}                 	= ['$',getLseg(ud.sc_mri_ds{v(i)},2),'$'];                  end;    end;
%
imy                             = imx(4:end,    :);
for i=imy(imy(:)>0)';
    s3                          = getLseg(ud.sc_mri_ds{v(i)}, 3);
    if s3(1)=='-';
        if get(ud.c1bHs(i+1), 'UserData')==2;
            sc1                 = get(ud.c1bHs(i+1), 'String');
            ss{i}               = sc1(sc1~=' ');
        elseif get(ud.c1bHs(i+1), 'UserData')==3;
            sc1                 = get(ud.c1bHs(i+1), 'String');
            ss{i}               = ['$',getLseg(ud.sc_mri_ds{v(i)},2),'$',sc1(sc1~=' ')];
        else;       
            ss{i}             	= ['$',getLseg(ud.sc_mri_ds{v(i)},2),'-modify$'];                   end;     
    else;                   
         ss{i}                 	= ['$',getLseg(ud.sc_mri_ds{v(i)},2),'$'];                  end;    end;
%
ud.mri_ds.symbolic_c            = ss;
out                             = ss{1};
for i=2:1:size(v,1);            out                         = fullfile(out, ss{i});                 end;
%
ud.mri.ds.symbolic              = out;
ud.stage                        = 'mri_ds_done';
set(gcf,    'UserData',ud);
%
set(ud.c1bHs(1), 'String','Review generalized MRI path in Data server @info-board');
set(ud.c2bHs(1), 'String','Move on');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Generalized MRI path in Data server:',  ['   ',ud.mri.ds.symbolic],  ...     
    ' Segments without $   :   fixed to them across subjects / scans',     ...
    ' $whatever_ID$        :   to copy $whatever_ID$ from image server',   ...
    ' $whatever_ID-modify$ :   to copy $whatever_ID$ after modifications', ...
    ' $whatever_ID$xxx     :   to copy $whatever_ID$ with xxx appended'});
%
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%

function                     	local_done_idae_0(ud);
%%
%
% reset left & right columns:
local_reset_guis(ud);
%
set(ud.c1bHs(1), 'String','Follow instructions below.  Then, hit > ');
set(ud.c2bHs(1), 'String','Start',   'CallBack','iv2_gen_lds(''done'',[]);');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',       ...
    {'* Let''s create a folder to place files to run IDAE (analysis packages, etc.)', ...
    '> Select a file in your home directory in Data Server',    ...
    '  (your home directory usually ends with your username)'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
ud.stage                        = 'idae_1';
set(gcf,    'UserData',ud);
return
%%

function    sc                  = local_path2sc(path_x);
%%
%                           
if path_x(1)==filesep;          s0                          = filesep;
else;                           s0                          = '';                                   end;
% just in case speces are allowed in the system:
path_x(path_x==' ')             = '*';
path_x(path_x==filesep)         = ' ';
sc                              = getLseg(path_x,   [0,2]);
%
% construction of source PET path in a cell array:
for i=1:1:numel(sc);            sc{i}(sc{i}=='*')           = ' ';                                  end;
sc{1}                           = [s0,sc{1}];
return;
%%

function                        local_idae_toggle(i2);
%%
ud                              = get(gcf,  'UserData');
if sum((get(gco, 'BackgroundColor') - iv2_bgcs(0)).^2)>10^-6;                       return;         end;
%
set(ud.c1bHs(2:numel(ud.idae.real_c)),  'BackgroundColor',iv2_bgcs(0));
set(ud.c2bHs(2:numel(ud.idae.real_c)),  'String',' ');
set(gco, 'BackgroundColor',iv2_bgcs(18));
c_tag                           = get(gco,  'Tag');
set(findobj(gcf, 'Tag',[c_tag(1, 1:end-1),'2']),   'String','User_ID');
return;
%%

function                     	local_done_idae_1(ud);
%%
set(gco,    'Enable','off');
if isfield(ud,'reuse') && isfield(ud.reuse,'idae') && ud.reuse.idae>0;
                                path_x                      = fileparts(ud.idae.real);
else;                           [fname, path_x]             = uigetfile(fullfile(pwd,'*'));        
    if ~ischar(fname);         	set(gco,    'Enable','on');                         return;         end;
                                                                                                    end;
%
ud.idae.real                    = fullfile(path_x, '(enter a folder name without spaces)');
ud.idae.real_c                  = local_path2sc(ud.idae.real);
%                           
%
for i=1:1:numel(ud.idae.real_c)-1;
   	set(ud.c1bHs(i+1), 'String',ud.idae.real_c{i}, 'CallBack','iv2_gen_lds(''idae_toggle'',[]);');	end;
%
set(ud.c1bHs(i+2), 'Style','edit',  'String',ud.idae.real_c{end}, 'BackgroundColor',iv2_bgcs(6));
set(ud.c2bHs(i+2), 'String','IDAE-ID'); 
%
set(ud.c1bHs(i+3), 'String','Do not see your user name on left column GUIs? If No, hit >');
set(ud.c2bHs(i+3), 'String','Restart',  'CallBack',['ud=get(gcf,''UserData''); ',  ...
                    'ud.stage=''idae_0''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);']); 
%
set(ud.c1bHs(1), 'String','Your home path segments. Follow below instructions. Then, hit >');
set(ud.c2bHs(1), 'String','Done',   'CallBack','iv2_gen_lds(''done'',[]);');
%
ud.stage                        = 'idae_2';
set(gcf,    'UserData',ud);
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Complete below tasks: ','1. Hit left column GUI showing your user name',  ...
    '2. Replace light green GUI with a non-existing folder name (no spaces please)',    ...
    '   (See MATLAB command window for existing folders)'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
%
d                               = dir(fullfile(fileparts(ud.idae.real),'*'));
d_name                          = char(d.name);
disp(['.directories in: ',fileparts(ud.idae.real)]);
dispCharArrays(1, char(d(d_name(:,1)~='.' & [d.isdir]'>0).name));
disp('> do not enter above existing folders in light green GUI.');
set(gco,    'Enable','on');
return;
%%

function                     	local_done_idae_2(ud);
%%
im1                             = umo_cstrs(char(get(ud.c2bHs, 'String')), ['User';'IDAE'], 'im1');
h0                              = findobj(gcf, 'Style','edit');
set(h0,     'Style','pushbutton');
% checking if User_ID is present:
if im1(1)<1;    
    for i=2:1:im1(2,1)-1;       set(ud.c1bHs(i),  'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(ud.c1bHs(i),  'BackgroundColor',iv2_bgcs(0));                   end;
    set(h0, 'Style','edit');                                                    	return;         end;
%
% constructing ud.idae.symbolic_c:
ss_end                          = get(ud.c1bHs(im1(2,1)),    'String');
ud.idae.real_c{end}             = ss_end(ss_end~=' ');
ud.idae.real                  	= fullfile(fileparts(ud.idae.real), ud.idae.real_c{end});
if exist(ud.idae.real,'dir');
    set(ud.c1bHs(im1(2,1)),     'String',[ud.idae.real_c{end},' - exists. replace it'],     ...
                                'Style','edit', 'BackgroundColor',iv2_bgcs(11));
 	pause(0.5);
    set(ud.c1bHs(im1(2,1)),     'BackgroundColor',iv2_bgcs(6));                     return;         end;
%
s                               = mkdir(ud.idae.real);
if s<1;
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
        {'* Unable to create your IDAE home directory ',    ...
        '  (No full-access to your home directory in Data Server?)', 	...
        '> Gathered information saved. ''Quit'' and Consult your local manager'});   
    local_save([]);                                                                 return;         end;
%
ss                              = ud.idae.real_c;
ss{im1(1,1)-1}                  = '$User_ID$';
ss{end}                         = ss_end(ss_end~=' ');
ud.idae.symbolic_c              = ss;
ddd                             = ss{1};
for i=2:1:numel(ss);            ddd                         = fullfile(ddd, ss{i});                 end;
ud.idae.symbolic                = ddd;
%
ud.stage                        = 'idae_done';
set(gcf,    'UserData',ud);
set(ud.c1bHs(1), 'String','Done for IDAE home directory! Observe info-Board. Then, hit >')
set(ud.c2bHs(1), 'String','Move on',   'CallBack','iv2_gen_lds(''done'',[]);');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
    {'* Your IDAE home directory: ',[' ',ud.idae.real],             ...
    '* Generalized IDAE home directory, shared by users with the same conventions:',	...
    [' ',ud.idae.symbolic],    '> Hit ''Display'' for more info'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%

function                     	local_done_idae_done(ud);
%%
% reset left & right columns:
local_reset_guis(ud);

set(ud.c1bHs(1), 'String','Select a PET source file in Image Server');
set(ud.c2bHs(1), 'String','Start',   'CallBack','iv2_gen_lds(''pet_is'',0);');
% 
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
    {'* Select a PET source file in Image Server to learn the directory conventions'});
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%

function                        local_done_mri_ds_done(ud);
% reset left & right columns:
local_reset_guis(ud);
%
s_count                         = 0;
%% SPM12
spm_code                        = which('spm_vol.m');
ic                              = 2;
disp('.sorting out SPM12'' home directory:');
if isempty(spm_code);       
    set(ud.c1bHs(ic),   'String','Find SPM12''s home directory (select spm_vol.m)');
    set(ud.c2bHs(ic),   'String','Start',   'CallBack','iv2_gen_lds(''search_spm'',[]);');
    s_count                     = s_count + 1;
else;                           
    ud.etc.spm_home             = fileparts(spm_code);
  	disp(['* SPM12''s home directory: ',ud.etc.spm_home]);
  	set(ud.c1bHs(ic),   'String','SPM12'' home directory');
 	set(ud.c2bHs(ic),   'String','Done');                                                           end;
%% Freesurfer
disp('.sorting out directory to place scripts for FreeSurfer:');
s                               = 1;
if ~exist(fullfile(ud.idae.real,'fs_scripts'),'dir');
    s                          	= mkdir(fullfile(ud.idae.real,'fs_scripts'));                       end;
ic                              = ic + 1;
set(ud.c1bHs(ic),   'String','Folder for Freesuefer scripts');
if s>0;                         
    disp('> created/present');
 	disp([' ',fullfile(ud.idae.real,'fs_scripts')]);
   	set(ud.c2bHs(ic),   'String','Done');
    ud.etc.fs_script_ws_real    = fullfile(ud.idae.real,'fs_scripts');
    ud.etc.fs_script_ws         = fullfile(ud.idae.symbolic,'fs_scripts');
    ud.etc.fs_script_ws_c       = ud.idae.symbolic_c;
    ud.etc.fs_script_ws_c{end+1}                            = 'fs_scripts';     
    qqq                         = {'FS_done.txt', 'FS_submitted.txt'};
    for i=1:1:numel(qqq);
        write2ptf(fullfile(ud.idae.real,'fs_scripts',qqq{i}), 'symbolic');                      	end;
else;                           disp('> unable to create folder for Freesuefer scripts');
                                disp([' you have no full-control of: ',ud.idae.real]);
                                set(ud.c2bHs(ic),   'String','Failed');                             end;
% window:
if any(ud.pet.ds.real=='\') && s>0;
    ic                          = ic + 1;
    set(ud.c1bHs(ic),   'String','Enter Freesurfer''s script folder seen from Linux side');
    set(ud.c2bHs(ic),   'String','Start',   'CallBack','iv2_gen_lds(''fs_script'',''w'');');
    s_count                     = s_count + 1;                                                      end;
if ~any(ud.pet.ds.real=='\') && s>0;
    ud.etc.fs_script_ds         = ud.etc.fs_script_ws;
    ud.etc.fs_script_ds_c       = ud.etc.fs_script_ws_c;                                            end;
% 
ic                              = ic + 1;
disp('.sorting out Freesurfer''s ''subject'' directory:');
if any(ud.pet.ds.real=='\');
    set(ud.c1bHs(ic),   'String','Enter Freesurfer''s subject folder as instructed');
    set(ud.c2bHs(ic),   'String','Start',   'CallBack','iv2_gen_lds(''fs_subject'',''w'');');
  	s_count                     = s_count + 1;
else;
    set(ud.c1bHs(ic),   'String','Select sample-001.mgz in Freesurfer''s subject folder');
    set(ud.c2bHs(ic),   'String','Start',   'CallBack','iv2_gen_lds(''fs_subject'',''l'');');
    s_count                     = s_count + 1;                                                      end;
%    
% ''dump'' folder
disp('.sorting out ''dump'' folder to place temporary files:');
s2                              = 1;
if ~exist(fullfile(ud.idae.real,'tmp'),'dir');  
    s2                          = mkdir(fullfile(ud.idae.real,'tmp'));                              end;
ic                              = ic + 1;
set(ud.c1bHs(ic),   'String','Dump folder to place temporary files');
if s2>0;                      	disp('> created/present');
                                disp([' ',fullfile(ud.idae.real,'tmp')]);
                                set(ud.c2bHs(ic),   'String','Done');
    ud.etc.dump                 = fullfile(ud.idae.symbolic,'tmp');
    write2ptf(fullfile(ud.idae.real,'tmp','scratch.m'), 'presence only');
else;                           disp('- unable to create ''dump'' folder');
                                disp([' you have no full-control of: ',ud.idae.real]);
                                set(ud.c2bHs(ic),   'String','Failed');                             end;
%% data unit:
disp('.sorting out radioactivity units:');
ic                              = ic + 1;
set(ud.c1bHs(ic),   'String','Select the radioativity unit to use');
set(ud.c2bHs(ic),   'Value',1,  'Style','popupmenu',    'String',{'Select','nCi/mL','Bq/mL'},   ...
                                'CallBack','iv2_gen_lds(''unit_done'',[]);');

ud.stage                        = 'all_done';

ud.matlab.startup               = which('startup.m');
set(gcf,    'UserData',ud);

set(ud.c1bHs(1), 'String','A few more to go. Follow the instructions shown in info-board below')
set(ud.c2bHs(1), 'String','Done',   'CallBack','iv2_gen_lds(''done'',[]);');
%
if s_count>0;
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
        {['* ',int2str(s_count),' ''Start'' GUIs (right column) & radioactivity unit to go!'],   ...
        '> Visit each ''Start'' and complete the task',   '> Slect the radioactivity unit'});
else;
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
        {'* One more to go', '> Slect the radioactivity unit'});
    local_save([]);                                                                                 end;
set(findobj(gcf, 'Tag','iv2_gen_lds_end_1'), ...
                                'Userdata',get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String'));
return;
%%

function                        local_get_spm_path(i2);
%%
set(gco,    'Enable','off');
ud                              = get(gcf, 'UserData');
[fname, ud.idae.code_path]      = uigetfile(fullfile(pwd,'*'));
if ~ischar(fname);                                                                  return;         end;
if ~exist(fullfile(ud.idae.code_path, 'spm_vol.m'),'file');
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
        {'* Wrong folder for SPM12? Unable to find spm_vol.m there', '> Try one more time'});
    set(gco,    'Enable','on');                                                   	return;         end;
%
set(gcf,    'UserData',ud);
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',{'* SPM12 folder: ',[' ',ud.idae.code_path]});
return;
%%

function                        local_fs_script(i2);
%%
ud                              = get(gcf,      'UserData');
set(gco,    'Enable','off');
drawnow;
if i2=='w';
    disp(['> Freesurfer script folder (Windowds): ',ud.etc.fs_script_ws_real]);
    disp('  Find it in you the Linux system to run Freesurfer (typically starts with /mnt)');
    disp('  There should be FS_done.txt  FS_submitted.txt in the flder by now');
    disp('> Copy & past it below');
    s                           = input(' Enter the folder name seen from Linux (Ret=cancel): ','s');
    if isempty(s);                                                                  return;         end;
    im1                         = umo_cstrs(char(ud.idae.symbolic_c),'$User_ID', 'im1');
    if im1(1)<1 || size(im1,2)>1;                                                   return;         end;
    %
    s(s=='/')                   = ' ';
    sc                          = getLseg(s,    [0,2]);
    sc{1}                       = ['/',sc{1}];
    im2                         = umo_cstrs(char(sc),  ud.idae.real_c{im1}, 'im1');
    if im2(1)<1 || size(im2,2)>1;                                                   return;         end;
    %
    sc{im2(1)}                  = '$User_ID$';
    ud.etc.fs_script_ds_c       = sc;
    out                         = sc{1};
    for i=2:1:numel(sc);        out                         = [out,'/',sc{i}];                  	end;
    ud.etc.fs_script_ds         = out;
    
    set(gcf,    'UserData',ud);
    set(gco,    'String','Done',    'CallBack',' ', 'Enable','on');                 return;         end;
return;
%%

function                        local_fs_subject(i2);
%%
ud                              = get(gcf,      'UserData');
set(gco,    'Enable','off');
drawnow;
if i2=='w';
    disp('> Inquiring Freesurfer''s subject folder');
    disp('  Go to Freesurfer''s subject folder in your Linux system');
    disp('  (Tyoically try [your linux]: cd freesurfer/subject)');
    disp('  Then, type [your linux]: pwd')
    disp('  Copy & past the output below');
  	s                           = input(' Freesurfer''s subject folder in Linux) (Ret=cancel): ','s');
    if isempty(s);                                                                  return;         end;
    im1                         = umo_cstrs(char(ud.idae.symbolic_c),'$User_ID', 'im1');
    if im1(1)<1 || size(im1,2)>1;                                                   return;         end;
    %
    s(s=='/')                   = ' ';
    sc                          = getLseg(s,    [0,2]);
    sc{1}                       = ['/',sc{1}];
    im2                         = umo_cstrs(char(sc),  ud.idae.real_c{im1}, 'im1');
    if im2(1)<1 || size(im2,2)>1;                                                   return;         end;
    %
    sc{im2(1)}                  = '$User_ID$';
    ud.etc.fs_subject_ds_c      = sc;
    %
    out                         = sc{1};
    for i=2:1:numel(sc);        out                         = [out,'/',sc{i}];                  	end;
    ud.etc.fs_subject_ds        = out;
    set(gcf,    'UserData',ud);
    set(gco,    'String','Done',    'CallBack',' ', 'Enable','on');                 return;         end;
%

  
return;
%%

function                        local_unit_done(i2);
%%
if get(gco, 'Value')<2;                                                             return;         end;
ud                              = get(gcf,      'UserData');
s                               = get(gco,  'String');
ud.etc.unit                     = s{get(gco, 'Value')};
set(gcf,    'UserData',ud);
return;
%%

function                        local_done_all_done(ud);
%%
ud                              = get(gcf,      'UserData');
%
ok                              = 1;
%
fck                             = {'idae','pet_is','pet_ds','mri_is','mri_ds','etc'};
fn2ck{1}                        = {'code_path','symbolic_c','symbolic'};
fn2ck{2}                        = {'symbolic_c','format_c','search_c'};
fn2ck{3}                        = {'symbolic_c'};
fn2ck{4}                        = {'symbolic_c','format_c','search_c'};
fn2ck{5}                        = {'symbolic_c'};
fn2ck{6}                        = {'spm_home','fs_script_ws','fs_subject_ds','fs_script_ds',    ...
                                                        	'dump','unit'};
%
for i=1:1:numel(fck);
    if isfield(ud,fck{i});
        clear fnms;
        eval(['fnms            	= fieldnames(ud.',fck{i},');']);
        im1                   	= umo_cstrs(char(fnms), char(fn2ck{i}), 'im1');
        if any(im1<1);        	disp(['.missing critical fields @ud.',fck{i}]);
                                dispCharArrays(1, char(fnms(im1<1)));
                                ok                          = 0;                                    end;
    else;                     	disp(['.problem! ''',fck{i},''' field missing']);
                                ok                          = 0;                        	end;    end;
%
if ok<1;                        set(ud.c1bHs(1), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(ud.c1bHs(1), 'BackgroundColor',iv2_bgcs(12));   return;         end;
%
iv2_gen_lds_2(ud);
% 
% save(fullfile(lds.iv2_home,['aid_',ud.lds,'.mat']),     'lds');
% disp('.done! (information for your local adapeter saved)');
% disp([' output: ',fullfile(lds.iv2_home,['aid_',ud.lds,'.mat'])]);    
return;
%%


