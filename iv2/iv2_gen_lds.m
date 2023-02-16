function    iv2_gen_lds(fun,i2); 
% To prepare local adaptoer files (dxetc4iv2_gen_lds.m) for your site:
%       
%       usage:      iv2_gen_lds('set','dxetc4xxx')
%       
% 
% (cL)2022    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;
if exist(which(['local_',lower(fun)]))
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
set(c1bHs(1),   'Tag','iv2_gen_lds_1_1',    'BackgroundColor',iv2_bgcs(6), 'FontWeight','bold');
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
local_sort_idae(ud);
return;
%%

function                        local_sort_idae(ud);
%%
%
set(ud.c1bHs(1),   'String','Functional definitions of computer species for this session');
set(ud.c2bHs(1),   'String','Move on', 'CallBack','iv2_gen_lds(''done'',[]);');
set(ud.c1bHs(2),   'String','Workstations (WS) perform data-analysis with IDAE');
set(ud.c1bHs(3),   'String','Image Server (IS; n=1) provides reconstructed PET & MRI files'); 
set(ud.c1bHs(4),   'String','Data Server (DS; n=1) stores analysis inputs/outputs (n=1), and');
set(ud.c1bHs(5),   'String','has a folder housing user folders for all users, including you');
%
set(ud.c1bHs(6),   'String','Questionaries about your home directory in Data server', ...
                                                            'BackgroundColor',iv2_bgcs(16));
%
ic                              = 7;
set(ud.c1bHs(ic), 'String','Enter your username that was set in Data server');
if isfield(ud.idae,'user_name');
    set(ud.c2bHs(ic), 'String','Done', 'CallBack','iv2_gen_lds(''copy_paste'',[]);');
    set(ud.c1bHs(ic+1), 'Style','text', 'UserData','user_name',     ...
                            'String',ud.idae.user_name,  'HorizontalAlignment','center');
else
    set(ud.c2bHs(ic), 'String','Start', 'CallBack','iv2_gen_lds(''copy_paste'',[]);');
    set(ud.c1bHs(ic+1), 'Style','text', 'UserData','user_name',     ...
                            'String',' ',  'HorizontalAlignment','center');                         end
%
% user's home directory
ic                              = 9;
set(ud.c1bHs(ic), 'String','Select a file from your home directory in Data server');
if isfield(ud.idae,'user_home')
    set(ud.c2bHs(ic), 'String','Done', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');
    set(ud.c1bHs(ic+1), 'String',ud.idae.user_home, 'Style','text', ...
                                'UserData','user_home', 'HorizontalAlignment','center');
else;
    set(ud.c2bHs(ic), 'String','Start', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');    
    set(ud.c1bHs(ic+1), 'String',' ', 'Style','text', ...
                                'UserData','user_home', 'HorizontalAlignment','center');            end
%
% user's IDAE home directory
ic                              = 11;
set(ud.c1bHs(ic), 'String','Generate your IDAE home directory in Data server');
if isfield(ud.idae,'user_home')
    set(ud.c2bHs(ic), 'String','Done', 'CallBack','iv2_gen_lds(''idae_home'',[]);');
    set(ud.c1bHs(ic+1), 'String',fullfile(ud.idae.user_home,'iv2'), ...
                                'Style','text', 'HorizontalAlignment','center');
else;
    set(ud.c2bHs(ic), 'String','Do it', 'CallBack','iv2_gen_lds(''idae_home'',[]);');  
    set(ud.c1bHs(ic+1), 'String',' ', 'Style','text', 'HorizontalAlignment','center');              end

%
if isfield(ud.idae,'new_user')
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
        {'* Setting IDAE for a new user'
        '> Enter your username (Done > type your username > Enter)'
        '  (skip Questionaries 2 & 3 if no change to the conventions)'
        '> Hit ''Move on'' @top-right'});
else
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
        {'* Observe functional definitions of WS/IS/DS for this session'
        ' - WS/IS/DS may be phisically identical or separated in any ways'
        ' - IS is usullay managed by the PET/MRI Centers'
        '> Complete questionaries about your home directory in Data server '
        ' - The last segment of your home directory must agree with your username'
        ' - This module creates a folder to place files/folders for performance of IDAE'});         end
%
ud.stage                        = 'set_idae';
set(gcf, 'UserData',ud);
return;
%%

function                        local_pet_is(ud);
%%
%
if isempty(ud);                 ud                          = get(gcf, 'UserData');                 end;
local_reset_guis(ud);
%
ud.pet_is.real_c                = local_str2cell(ud.pet.is.real);
ns                              = numel(ud.pet_is.real_c);
%
set(ud.c1bHs(1), 'String','Starting the section of PET path in Image server');
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
set(ud.c1bHs(ns+2),  'String','Want to re-select path? If so, restart >');
set(ud.c2bHs(ns+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''idae_done''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);']);
%
for i=1:1:ns;
    set(ud.c1bHs(i+1),  'String',ud.pet_is.real_c{i});
    set(ud.c2bHs(i+1),  'String',['seg_',int2str(i)], 'CallBack','iv2_gen_lds(''toggle_pet'',0);'); end
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'FontName','Consolas',  'FontSize',10, ...
    'HorizontalAlignment','left', 'String',     ...
    {'* PET path segments (left) and segment #s (right) are shown'
    '> Identify at least one each of (using right GUIs which toggles):'
    ' 1. Fixed segments (i.e., consistent across all subjects/scans), and'
    ' 2. the first PET-proper segment'
    '> Hit ''Done'' when all are done'});
%
ud.stage                        = 'pet_is_1';
set(gcf,    'UserData',ud);
return;
%%

function                        local_toggle_pet(i2);
%%
set(findobj(gcf, 'Tag','iv2_gen_lds_1_2'), 'Enable','on')
if strncmpi(get(gco,'String'),'seg',3);
    set(gco, 'String','fixed', 'BackgroundColor',iv2_bgcs(18));
elseif strncmpi(get(gco,'String'),'fix',3);
    set(gco, 'String','PET folder', 'BackgroundColor',iv2_bgcs(16)); 
else;
    s0                          = get(gco, 'Tag');
    s1                          = find(s0=='_',2, 'last');
    set(gco, 'String',['seg_',int2str(str2num(s0(s1(1)+1:s1(2)-1))-1)], ...
                                                            'BackgroundColor',iv2_bgcs(0));         end
return;
%%

function                        local_done(i2);
%%
ud                              = get(gcf,  'UserData');
% ud.stage
feval(['local_done_',ud.stage], ud);
return;
%%

function                      	local_done_pet_is_1(ud);
%%
[symbolic_c, search_c]          = local_done_xxx_is_1(ud.pet_is.real_c,ud.c1bHs,ud.c2bHs);
%
if isempty(symbolic_c);
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String',  ...
        {'* Incomplete responses'
        '> Identify at least one each of (using right GUIs which toggles):'
        ' 1. Fixed segments (i.e., consistent across all subjects/scans), and'
        ' 2. the first PET-proper segment'});                                       return;         end
%
% removing callback from right column GUIs:
set(ud.c2bHs(2:numel(symbolic_c)+1), 'CallBack',' ');
%
ud.pet_is.symbolic_c            = symbolic_c;
ud.pet_is.search_c              = search_c;
ud.pet.is.symbolic              = local_str2cell(symbolic_c);

local_done_pet_is_2(ud);
%
return
%% 

function    [o1, o2]            = local_done_xxx_is_1(real_c,c1bHs,c2bHs);
%%
fixSer                          = umo_cstrs(['fix';'Ser';'PET'],  ...
                                    char(get(c2bHs(2:numel(real_c)+1), 'String')), 'im1');
%
o1                              = [];       % = symbolic_c
o2                              = [];       % = search_c
ok                              = 1;
if ~any(fixSer==1);             ok                          = 0;                                    end
if ~any(fixSer==2) && ~any(fixSer==3)
                                ok                          = 0;                                    end
if ok<1;                                                                            return;         end
%
o1                              = real_c;
for i=1:1:size(fixSer, 1)
    if fixSer(i)==2;            o1{i}                       = '$Series_ID$';
    elseif fixSer(i)==3;        o1{i}                       = '$PET_ID$';
    elseif fixSer(i)<1;         o1{i}                       = '#';                          end;    end
%
o2i                             = o1;
for i=find(fixSer'>1);          o2i{i}                      = '#';                                  end
%
% defining segments for source file search:
fixSer(find(fixSer>1,1):end)    = max(fixSer);
% 
% to include 'pet folder' for PET (not including 'Series folder' for MRI:
if max(fixSer)==3;              fixSer(find(fixSer==3,1))   = 0;                                    end;
%
set(c1bHs(find(fixSer<1)+1),    'BackgroundColor',iv2_bgcs(18));
%
o2                              = o2i(1:find(fixSer>1,1)-1);
%     fixSer(find(fixSer==2,1):end)                           = 2;
%     set(c1bHs())
% 
% % 
% % setting callback for non-fixed left-column GUIs:
% set(ud.c1bHs(find(fixSer~=1)+1), 'CallBack',['iv2_gen_lds(''toggle'',''',pom,''');']);
% 
% % 
% set(ud.c1bHs(1), 'String','Follow the instructions below. Hit > when done');
% set(ud.c2bHs(1), 'String','Accept',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','off');
% 
% %
% set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',                           ...
% 	{'* Try if one segment can narrow down the search for '
%     ['  ',upper(pom),' source files sufficiently (to use when IDAE is up)']
%     '> Hit a path segment (left) you think the best (i.e., the least error-prone)'
%     '  in your system and observed/respond the results shown here'});
return;
%%

function                        local_done_pet_is_2(ud);
%%
local_done_xxx_is_2(ud.c1bHs,'pet')
ud.stage                        = 'pet_is_3';
%
set(ud.c2bHs(1), 'String','Done');
%
set(gcf,    'UserData',ud);
return;
%%

function                        local_done_xxx_is_2(c1bHs,pom);
%%
% extracting marked segments (=the search criteria):
c18                             = iv2_bgcs(18);
mkd                             = find( sum( abs(cell2mat(get(c1bHs, 'BackgroundColor')) - ...
                                    c18(ones(numel(c1bHs),1), :)), 2) < 10.^-6) -1;
%
if isempty(mkd);                set(c1bHs(1), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5);
                                set(c1bHs(1), 'BackgroundColor',iv2_bgcs(6));   return;         end
% 
% making path segments of the search criteria editable:
set(c1bHs(mkd+1), 'Style','edit', 'Callback',' ');
%
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
	{['* Any orange segment(s) (left; editable) include date of ',upper(pom),'?']
    '> If yes, replace date integers with yyyy or yy (year), mm (month), and dd (day)'
    ' - Leave fixed elements unchanged (e.g., ''abc-'' of abc-yyyymmdd), and'
    ' - Replace inconsistent elements with * in the segment(s)'
    '> Hit ''Done'' @top-right when done'})
% 
set(c1bHs(1), 'String',['Work on ',upper(pom),' date elements as instructed below, if any.']);
return;
%%

function                     	local_done_pet_is_3(ud);
%%
%
ud.pet_is.search_c              = local_done_xxx_is_3(ud.pet_is.real_c,ud.pet_is.search_c,ud.c1bHs,'pet');
%
set(ud.c1bHs(1), 'String','For information alone. Review the info board below');
set(ud.c2bHs(1), 'String','Next',    'CallBack','iv2_gen_lds(''pet_ds'',[]);');
%
ud.stage                        = 'end_pet_is';
set(gcf, 'UserData',ud)
%
return
%%

function    search_c            = local_done_xxx_is_3(real_c,search_c,c1bHs,pom);
%%
%
% identifying modified segments (for date of service):
im1                             = umo_cstrs(char(get(c1bHs(2:numel(real_c)+1), 'String')), ...
                                                            char(real_c), 'im1');
%
im1(numel(search_c)+1:end, :)   = 5;
for i=find(im1'<1);             s1                          = get(c1bHs(i+1), 'String');
                                s1                          = s1(s1~=' ');
    if any(s1=='*');            s1c                         = getLseg(replace(s1,'*',' '), [0,2]);
                                s1x                         = '';
        for j=1:1:numel(s1c);   s1x                         = [s1x,s1c{j},'*'];                     end
        if s1(1)=='*';          s1x                         = ['*',s1x];                            end
        if s1(end)~='*';        s1x                         = s1x(1, 1:end-1);                      end
                                s1                          = s1x;                                  end
    search_c{i}                 = ['$',s1,'$'];                                                     end
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',       ...
	{['* Starting search path for ',upper(pom),' source files (to use when IDAE is up):']
    ['   ',local_str2cell(search_c)]
    '  where segments with # and $ will be supplied from the study log file'
    ['< This is the end of the section of ',upper(pom),' paths in Image server']
    '* Hit ''Next'' @top-right to start the Data server section'});


function                        local_reset_guis(ud);
%%
set(ud.c1bHs(2:end),    'String',' ',   'Value',1,  'Style','pushbutton',  'UserData',[],    ...
                                'BackgroundColor',iv2_bgcs(0),  'CallBack',' ', 'Enable','on');
set(ud.c2bHs(2:end),    'String',' ',   'Value',1,  'Style','pushbutton',  'UserData',[],    ...
                                'BackgroundColor',iv2_bgcs(0),  'CallBack',' ', 'Enable','on');
return;
%%

function                        local_pet_ds(ud);
%%
if isempty(ud);                 ud                          = get(gcf, 'UserData');                 end;
local_reset_guis(ud)
%
ud.pet_ds.real_c                = local_str2cell(ud.pet.ds.real);
ns                              = numel(ud.pet_ds.real_c);
% %
ud.stage                        = 'pet_ds_1';
%
% constructing selections
ss                              = {'fixed', 'Subject ID','Study ID','copy from study Log file'};
ic                              = numel(ss);
for i=1:1:numel(ud.pet_is.symbolic_c)
    ic                          = ic + 1;
    ss{ic}                      = ['is_seg_',int2str(i),' (',ud.pet_is.real_c{i},')'];              end
% 
for i=1:1:ns;
    set(ud.c1bHs(i+1), 'String',ud.pet_ds.real_c{i});
    set(ud.c2bHs(i+1), 'Value',1, 'Style','popupmenu', 'String',ss);                                end
%
set(ud.c1bHs(1), 'String','Starting the section of PET paths in Data server')
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
set(ud.c1bHs(ns+2),  'String','Want to re-select path? If so, restart >');
set(ud.c2bHs(ns+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''pet_is_4''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'], 'Enable','on');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Segments of user-selected PET-path are shown on left GUIs'
    '> Select one each from right pulldown menu according to the local path conventions'
    '  - ''fixed'' to apply the string of left column to all scans' 
    '  - Subject ID, Study ID, and other items to copy from the study log file'
    '  - ''is_seg_i'' to make use of Segment i of the Image server path (in parentheses)'
    '> Review the selections. Hit ''Done'' if all look good'});
set(gcf,    'UserData',ud);
return;
%%

function                     	local_done_pet_ds_1(ud);
%%
% unmarking problematic segments, if any:
% set(ud.c2bHs,   'BackgroundColor',iv2_bgcs(0));     
%
for i=1:1:numel(ud.pet_ds.real_c);
    v                           = get(ud.c2bHs(i+1), 'Value');
    if i==1;                    s0                          = get(ud.c2bHs(i+1), 'String');         end
    s2{i}                       = s0{v};
    ss{i}                       = getLseg(s0{v}, 1); 
    if strncmpi(ss{i},'is_seg_',7);
        set(ud.c2bHs(i+1), 'Value',1, 'Style','edit', 'String',s0{v}); 
    else;
        set(ud.c2bHs(i+1), 'Value',1, 'Style','pushbutton', 'String',s0{v});                end;    end 
%
% disabling 'restart' GUIs:
set(ud.c1bHs(numel(ud.pet_ds.real_c)+2), 'Enable','off');
set(ud.c2bHs(numel(ud.pet_ds.real_c)+2), 'Enable','off');
%
%
set(ud.c1bHs(1), 'String','Work on right GUIs as instructed below')
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',       ...
    {'* Modify eligible right GUIs (= editable) per local path conventions'
    '  or to make resulting paths specific to IDAE, as needed/desired'
    '  - Prepend/append a string (or both) to is_seg_i'
    '  - Add * anywhere to specify the modifications later'
    '  - Replace is_seg_i with a string to use as a fixed segment'
    '* Do not erase any (i.e., add only) except replacing is_seg_i'});
%
ud.pet_ds.is_segs_c             = ss;
ud.pet_ds.c2_str_c              = s2;
ud.stage                        = 'pet_ds_2';
set(gcf,    'UserData',ud);
return;
%%

function                     	local_done_pet_ds_2(ud);
%%
%
[ud.pet_ds.symbolic_c, ud.pet_ds.fixSer]                    = local_done_xxx_ds_2(ud,'pet',ud.pet_ds);
%
ud.stage                        = 'pet_ds_3';
set(gcf,    'UserData',ud);
%
return;
%%

function   [symbolic_c, fixSer] = local_done_xxx_ds_2(ud,pom,xxx_ds);
%%
%
symbolic_c                      = xxx_ds.is_segs_c;
fixSer                          = umo_cstrs(['fixed ';'is_seg';'Series'],char(xxx_ds.is_segs_c), 'im1');
% working on 'fixed' segments:
for i=find(fixSer'==1);
    set(ud.c1bHs(i+1), 'String',[xxx_ds.real_c{i},' (fixed)']);
    symbolic_c{i}               = xxx_ds.real_c{i};                                                 end
%
% working on study ID, Study subject ID, and others from study log file:
for i=find(fixSer'<1);
    sR                          = get(ud.c2bHs(i+1), 'String');
    if strcmpi(getLseg(sR,1),'copy');
        set(ud.c1bHs(i+1), 'String',['(to ',sR,')']);
        symbolic_c{i}           = '#study_log_file#';
    else;
        sR(sR==' ')             = '_';
        set(ud.c1bHs(i+1), 'String', [sR, ' (from study log file)']);
        symbolic_c{i}           = ['#',sR,'#'];                                             end;    end
%
% 
for i=find(fixSer'==2);
    sR                          = get(ud.c2bHs(i+1), 'String');
    xxx_is_str                  = xxx_ds.c2_str_c{i}(1, find(xxx_ds.c2_str_c{i}=='(',1)+1: ...
                                    find(xxx_ds.c2_str_c{i}==')',1,'last')-1);
    % to specify modification strings later (more complex)
    if any(sR=='*');
        set(ud.c1bHs(i+1), 'String',[xxx_is_str,' (to modify)']);
        symbolic_c{i}           = ['$',xxx_ds.is_segs_c{i},'$*'];
    % the eligible string (is_seg_i or Series ID) is replaced by a string: 
    elseif ~contains(sR,getLseg(xxx_ds.c2_str_c{i},1))
        set(ud.c1bHs(i+1), 'String',[getLseg(sR,1),' (to fix)']);
        symbolic_c{i}                             = getLseg(sR,1);
    % no changes at all:
    elseif strcmpi(sR(sR~=' '),xxx_ds.c2_str_c{i}(xxx_ds.c2_str_c{i}~=' '));
        % set(ud.c1bHs(i+1), 'String',[ud.pet_ds.is_segs_c{i},' (to copy)']);
        set(ud.c1bHs(i+1), 'String',[xxx_is_str,' (as is)']);
        symbolic_c{i}                           = ['$',xxx_ds.is_segs_c{i},'$'];
    % prepend/append option:
    else;
        if any(sR=='(');        sR(1, find(sR=='(',1):end)  = ' ';                                  end;
        sR2                     = replace(sR,xxx_ds.is_segs_c{i},' * ');
        sR2_c                   = getLseg(sR2, [0,2]);
        if sR2_c{1}(1)=='*';
            % str                 = [ud.pet_ds.is_segs_c{i},sR2_c{2},' (to modify)'];
            str                 = [xxx_is_str,sR2_c{2},' (to modify)'];
            str_2               = ['$',xxx_ds.is_segs_c{i},'$',sR2_c{2}];
        else;
            % str                 = [sR2_c{1},xxx_ds.is_segs_c{i}];
            str                 = [sR2_c{1},xxx_is_str];
            str_2               = [sR2_c{1},'$',xxx_ds.is_segs_c{i},'$'];
            if numel(sR2_c)>2;  str                         = [str, sR2_c{3}];             
                                str_2                       = [str_2, sR2_c{3}];                    end
            str                 = [str,' (modified)']; 
            set(ud.c1bHs(i+1), 'String',str)
            symbolic_c{i}       = str_2;                                            end;    end;    end
% 
for i=find(fixSer'==3);
    sR                          = get(ud.c2bHs(i+1), 'String');
    sR(sR=='*')                 = ' ';
    sR2                         = replace(sR,'Series_ID', ' * ');
    cL                          = sR2(1, 1:find(sR2=='*',1)-1);
    cR                          = sR2(1, find(sR2=='*',1)+1:end);
    set(ud.c1bHs(i+1), 'String',[cL(cL~=' '), 'Series_ID', cR(cR~=' '),' (from MRI)']);
    symbolic_c{i}               = [cL(cL~=' '), '&Series_ID&', cR(cR~=' ')];                        end

% to make revisions feasible
set(ud.c1bHs(size(fixSer,1)+2), 'String','Not quite right? Revise selections and hit >', 'Enable','on');
set(ud.c2bHs(size(fixSer,1)+2), 'String','Update', ...
    'CallBack', ['ud = get(gcf, ''UserData''); iv2_gen_lds(''done_',pom,'_ds_2'',ud);'], 'Enable','on');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Data server path (left GUIs) that was generated per local path conventions'
    '  for user-enterd Image server path. Rule smmaries (in parentheses) are:'
    '  - ''fixed''/''to fix'' to use shown strings for all scans'
    '  - ''copy'' to obtain the strings from the study log file'
    '  - ''to modify''/''modified'' to modify ''is_seg_i'' later or as shown here'
    '> Review them carefully. Hit ''Take'' if all look fine'})
%
set(ud.c1bHs(1), 'String',['Generated ',upper(pom),' path segments in Data server (rules).']);
set(ud.c2bHs(1), 'String','Take', 'Callback','iv2_gen_lds(''done'',[]);');
%
%
return;
%%

function                     	local_done_pet_ds_3(ud);
%%
%
% resetting left & right column GUIs:
% local_reset_guis(ud)
local_mri_is(ud)
%
return;
%%

function                            local_disp(i2);
%%
ud                              = get(gcf,  'UserData');
disp(['> output dxetc4xxx: ',ud.lds]);
%
if isfield(ud,'idae');          disp('> IDAE-related Paths');
                                disp(ud.idae);                                                      end
if isfield(ud,'pet_is');        disp('> PET Path variables - Image server:');
                                disp(ud.pet_is);                                                  	end
if isfield(ud,'pet_ds');        disp('> PET Path variables - Data server:');
                                disp(ud.pet_ds);                                          	        end
%
if isfield(ud,'mri_is');        disp('> MRI Path variables - Image server:');
                                disp(ud.mri_is);                                                  	end
if isfield(ud,'mri_ds');        disp('> MRI Path variables - Data server:');
                                disp(ud.mri_ds);                                          	        end
if isfield(ud,'etc');          	disp('> Other critical settings');
                                disp(ud.etc);                                                       end
return;
%%

% function                     	local_toggle(i2);
% %%
% % deselecting gco, if it is selected already:
% if sum(abs(get(gco, 'BackgroundColor') - iv2_bgcs(18)))<10.^-6;
%     set(gco, 'BackgroundColor',iv2_bgcs(0));                                        return;         end;
% %
% set(gco, 'BackgroundColor',iv2_bgcs(18));
% %
% ud                              = get(gcf, 'UserData');
% % 
% eval(['nr                       = numel(ud.',i2,'_is.real_c);']);
% % i2 is eiter pet or mri:
% i2(:)                           = lower(i2);
% % disabling C1 GUIs:
% 
% set(ud.c1bHs(2:nr+1), 'Enable','off');
% drawnow;
% %
% 
% %
% bgcs                            = cell2mat(get(ud.c1bHs, 'BackgroundColor'));
% bg18                            = iv2_bgcs(18);
% ii                              = find(sum((bgcs - bg18(ones(size(bgcs,1),1), :)).^2,2)<10^-6)-1;
% %
% im1                             = umo_cstrs('fix',char(get(ud.c2bHs(2:max(ii)+1), 'String')), 'im1');
% im1(ii, :)                      = 1;
% %
% sstr                            = '';
% for i=1:1:max(ii);
%     if im1(i)>0;                eval(['ss{i}                = ud.',i2,'_is.real_c{i};']);
%     else;                       ss{i}                       = '*';                                  end
%     sstr                        = fullfile(sstr,ss{i});                                             end
% %
% if eval(['isfield(ud.',i2,'_is,''search_res'')']);
%     eval(['sss                  = ud.',i2,'_is.search_res;']);
%     % previously tried:
%     im1                         = umo_cstrs(char(sss.str),sstr,'im1');
%     if im1>0;
%         set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',       ...
%             {'* Previousely tried search criteria:'
%                 [' - String          : ',sstr]
%                 [' - Elapsed time    : ',num2str(sss(im1).etime)]
%                 [' - # of paths found: ',int2str(sss(im1).n)]
%                 '> Accept it if # < 20'
%                 '> Add more segments (max = 2), or try another segment (deselect = hit again)'});
%         set(ud.c1bHs(2:nr+1), 'Enable','on');                                       return;         end
%     %
%     q                           = numel(sss);
% % when not tried yet:
% else;                           q                           = 0;                                    end
% %   
% %
% tic;
% dxs                             = dir(sstr);
% t1                              = toc;
% %
% sss(q+1).str                    = sstr;
% sss(q+1).etime                  = t1;
% cm1                             = umo_cstrs(char(dxs.folder),[],'cm1');
% sss(q+1).n                      = sum(cm1(:,2)>0);
% %
% set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
%     {'* Current search criteria:'
%     [' - String          : ',sstr]
%     [' - Elapsed time    : ',num2str(sss(q+1).etime,3),' (sec)']
%     [' - # of paths found: ',int2str(sss(q+1).n)]
%     '> Accept it if # < 20'
%     '> Add more segments (max = 2), or try another segment (deselect = hit again)'});
% %
% eval(['ud.',i2,'_is.search_res  = sss;']);
% %
% set(gcf, 'UserData',ud);
% %
% set(ud.c1bHs(2:nr+1), 'Enable','on')
% set(ud.c2bHs(1), 'Enable','on')
% %
% return;
% %%

function                            local_help(i2);
%%
disp(['** Step of: ',get(findobj(gcf, 'Tag','iv2_gen_lds_1_1'), 'String')])
disp(char(get(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String')))
disp('> (the end)')
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

s1{1}                           = 'Resume from?';
cb{1}                           = ' ';
%
im1                             = umo_cstrs(char(fieldnames(ud)), ...
                                    char('idae','pet_is','pet_ds','mri_is','mri_ds','etc'), 'im1');
% idae has to be present - otherwise start from scratch:
if im1(1)<1;                    local_done_set_idae(ud);                                 return;         end
%
s1{end+1}                       = '- Re-enter your home directory in Data Server';
cb{end+1}                       = 'local_sort_idae(ud);';
%
%
s1{end+1}                       = '- Set IDAE for a new user';
cb{end+1}                       = ['ud.idae.new_user = ''yes''; ',  ...
                                    'set(findobj(gcf, ''String'',''Save''),''Enable'',''off''); ',  ...
                                        'local_sort_idae(ud);'];
%
%
s1{end+1}                       = '- Reselect PET/MRI paths in Image/Data servers';
cb{end+1}                       = 'local_get_is_ds(ud)';
%    
% pet_is is present:
if im1(2)>0 
    %
  	s1{end+1}                 	= '- Resume from the section of PET paths in Image Server';
    cb{end+1}                	= 'local_pet_is(ud);';                                              end
% pet_ds is present:
if im1(3)>0;
    %
  	s1{end+1}                 	= '- Resume from the section of PET paths in Data Server';
    cb{end+1}                	= 'iv2_gen_lds(''pet_ds'',);';                                      end
% mri_is is present:
if im1(4)>0;
    %
  	s1{end+1}                 	= '- Resume from the section of MRI paths in Image Server';
    cb{end+1}                	= 'local_mri_is(ud);';                                              end
% mri_ds is present:
if im1(5)>0;
    %
  	s1{end+1}                 	= '- Resume from the section of MRI paths in Data Server';
    cb{end+1}                	= 'local_mri_ds(ud);';                                              end
% etc is recorded:
if im1(6)>0;
    s1{end+1}                 	= '- Resume from the section of SPM12 & Freesurfer';
    cb{end+1}                	= 'local_sort_etc(ud)';                                             end
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

function                        local_mri_is(ud);
%%
%
if isempty(ud);                 ud                          = get(gcf, 'UserData');                 end;
local_reset_guis(ud);
% 
ud.mri_is.real_c                = local_str2cell(ud.mri.is.real);
ns                              = numel(ud.mri_is.real_c);
%
% construction of source PET path in a cell array:
for i=1:1:ns            
    set(ud.c1bHs(i+1), 'String',ud.mri_is.real_c{i});
    set(ud.c2bHs(i+1), 'String',['seg_',int2str(i)], 'CallBack','iv2_gen_lds(''toggle_mri'',[]);'); end

%
set(ud.c1bHs(1), 'String','Starting the section of MRI paths in Image server');
set(ud.c2bHs(1), 'String','Done',  'CallBack','iv2_gen_lds(''done'',[]);')  % 'Enable','off')
%
set(ud.c1bHs(ns+2),  'String','Want to re-select example MRI path? If yes, hit >');
set(ud.c2bHs(ns+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''pet_done''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);']);
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'FontName','Consolas', 'FontSize',10,  ... 
    'HorizontalAlignment','left', 'String',     ...
    {'* MRI path segments (left) and segment #s (right) are shown.'
    '> Identify the following segments (hit right-column GUIs; toggles): '
    '  1. fixed segments that apply for all subjects/scans, and'
    '  2. MRI series folder (the one that is most relevant to MRI series)'
    '> Hit ''Done'' when all look OK'})
ud.stage                        = 'mri_is_1';
set(gcf,    'UserData',ud);
return;
%%

function                        local_done_mri_is_1(ud);
%%
[symbolic_c, search_c]          = local_done_xxx_is_1(ud.mri_is.real_c,ud.c1bHs,ud.c2bHs);
%
if isempty(symbolic_c);
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'), 'String',  ...
        {'* Incomplete responses'
    '> Identify the following segments (hit right-column GUIs; toggles): '
    '  1. fixed segments that apply for all subjects/scans, and'
    '  2. MRI series folder (the one that is most relevant to MRI series)'});       return;         end
%
% removing callback from right column GUIs:
set(ud.c2bHs(2:numel(symbolic_c)+1), 'CallBack',' ');
%
ud.mri_is.symbolic_c            = symbolic_c;
ud.mri_is.search_c              = search_c;
ud.mri.is.symbolic              = local_str2cell(symbolic_c);

local_done_mri_is_2(ud);
%
return
%%

function                        local_toggle_mri(i2);
%%
set(findobj(gcf, 'Tag','iv2_gen_lds_1_2'), 'Enable','on')
if strncmpi(get(gco,'String'),'seg',3);
    set(gco, 'String','fixed', 'BackgroundColor',iv2_bgcs(18));                     return;   
elseif strncmpi(get(gco,'String'),'fix',3);
    set(gco, 'String','Series folder', 'BackgroundColor',iv2_bgcs(16));             return;         end;
%
s0                              = get(gco, 'Tag');
s1                              = find(s0=='_',2, 'last');
set(gco, 'String',['seg_',int2str(str2num(s0(s1(1)+1:s1(2)-1))-1)], 'BackgroundColor',iv2_bgcs(0));
return;
%%

function                        local_done_mri_is_2(ud);
%% 
local_done_xxx_is_2(ud.c1bHs,'mri')
%
set(ud.c2bHs(1), 'String','Done');
%
ud.stage                        = 'mri_is_3';
set(gcf,    'UserData',ud);
return;
%%

function                        local_done_mri_is_3(ud);
%%
%
ud.mri_is.search_c              = local_done_xxx_is_3(ud.mri_is.real_c,ud.mri_is.search_c,ud.c1bHs,'mri');
%
set(ud.c1bHs(1), 'String','For information alone. Review the info board below');
set(ud.c2bHs(1), 'String','Next',    'CallBack','iv2_gen_lds(''mri_ds'',[]);');
%
ud.stage                        = 'end_mri_is';
set(gcf, 'UserData',ud)
%
return;
%%

function                        local_mri_ds(ud);
%%
if isempty(ud);                 ud                          = get(gcf, 'UserData');                 end
% set(gco,    'Enable','off');
%
local_reset_guis(ud)

ud.mri_ds.real_c                = local_path2sc(ud.mri.ds.real);
ns                              = numel(ud.mri_ds.real_c);

%
% constructing selections
ss                              = {'fixed', 'Subject ID','Study ID','copy from study Log file'};

ic                              = numel(ss);
for i=1:1:ns
    ic                          = ic + 1;
    ss{ic}                      = ['is_seg_',int2str(i),' (',ud.mri_is.real_c{i},')'];              end
% 
for i=1:1:ns-1;
    set(ud.c1bHs(i+1), 'String',ud.mri_ds.real_c{i}, ...
        'Callback','set(findobj(gcf, ''Tag'',''iv2_gen_lds_1_2''),''Enable'',''on'');');
    set(ud.c2bHs(i+1), 'Value',1, 'Style','popupmenu', 'String',ss);                                end
%
set(ud.c1bHs(ns+1), 'String',ud.mri_ds.real_c{ns});
set(ud.c2bHs(ns+1), 'String','Series_ID', 'BackgroundColor',iv2_bgcs(16));   
%
%
set(ud.c1bHs(1), 'String','Starting the section of MRI paths in Data server')
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
set(ud.c1bHs(ns+2),  'String','Want to re-select path? If so, restart >');
set(ud.c2bHs(ns+2),  'String','Restart', 'CallBack',['ud=get(gcf,''UserData''); ',   ...
    'ud.stage=''mri_is_4''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);'], 'Enable','on');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Segments of user-selected MRI path in Data server are shown on left'
    '> Select one each from right pulldown menu per your path conventions '
    '  - ''fixed'' to apply the left-column string for all scans' 
    '  - Subject ID, Study ID, and other items to copy from the study log file'
    '  - ''is_seg_i'' to make use of Segment i of the Image server path (in parentheses)'
    '* The last segment is MRI Series ID in IDAE (prepend/append later if desired)'
    '> Review the selections. Hit ''Done'' if all look good'},  ...
                                'FontName','Consolas',  'FontSize',10,  'HorizontalAlignment','left');
%
%
ud.stage                        = 'mri_ds_1';
set(gcf,    'UserData',ud);
return;
%%

function                     	local_done_mri_ds_1(ud);
%%
%
for i=1:1:numel(ud.mri_ds.real_c)-1;
    v                           = get(ud.c2bHs(i+1), 'Value');
    if i==1;                    s0                          = get(ud.c2bHs(i+1), 'String');         end
    s2{i}                       = s0{v};
    ss{i}                       = getLseg(s0{v}, 1); 
    if strncmpi(ss{i},'is_seg_',7);
        set(ud.c2bHs(i+1), 'Value',1, 'Style','edit', 'String',s0{v}); 
    else;
        set(ud.c2bHs(i+1), 'Value',1, 'Style','pushbutton', 'String',s0{v});                end;    end 
%
% series ID GUI
i                               = i + 1;
set(ud.c2bHs(i+1), 'Value',1, 'Style','edit'); 
s2{i}                           = 'Series_ID';
ss{i}                           = 'Series_ID';
%
% disabling 'restart' GUIs:
set(ud.c1bHs(numel(ud.mri_ds.real_c)+2), 'Enable','off');
set(ud.c2bHs(numel(ud.mri_ds.real_c)+2), 'Enable','off');
%
%
set(ud.c1bHs(1), 'String','Selections of right-column GUIs are confired')
set(ud.c2bHs(1), 'String','Done',    'CallBack','iv2_gen_lds(''done'',[]);', 'Enable','on')
%
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',       ...
    {'* Modify eligible right GUIs (= editable) per local path conventions '
    '  or to make resulting paths specifit to IDAE, as needed/desired'
    '- Preprend/append a string (or both) to is_seg_i or Series ID'
    '- Add * anywhere to specify modification rules later'
    '- Replace is_seg_i with a string to apply for all scans '
    '* Do not erase any (i.e., add only) except replacing is_seg_i'});
%
ud.mri_ds.is_segs_c             = ss;
ud.mri_ds.c2_str_c              = s2;
ud.stage                        = 'mri_ds_2';
set(gcf,    'UserData',ud);
% 
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
%
[ud.mri_ds.symbolic_c, ud.mri_ds.fixSer]                    = local_done_xxx_ds_2(ud,'mri',ud.mri_ds);
%
ud.stage                        = 'mri_ds_3';
set(gcf,    'UserData',ud);
return;
%%

function                     	local_done_mri_ds_3(ud);
%%
%
local_sort_etc(ud);
%
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


function                     	local_done_set_idae(ud);
%%
%
sss                             = {'user_name','user_home','real'};
bNos                            = 8:2:12;
im1                             = umo_cstrs(char(fieldnames(ud.idae)), char(sss), 'im1');
if any(im1<1);
    set(ud.c1bHs(bNos(find(im1'<1,1))), 'BackgroundColor',iv2_bgcs(11));
    pause(0.5)
    set(ud.c1bHs(bNos(find(im1'<1,1))), 'BackgroundColor',iv2_bgcs(0));
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
        {'* Incomplete entry found (blinked in pink).'
        ' - Hit associated right GUI and complete the left GUI'});                  return;         end
%
if ~isfield(ud.idae,'new_user');
    [idx, unm2]                 = fileparts(ud.idae.user_home);
    if ~strcmpi(ud.idae.user_name, unm2)
        set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
            {'* Critical Problem! ' 
            ' - The last segment of your home directory not identical to your username ' 
            '* IDAE cannot not function with this problem' 
            '> Create a folder in Data server to hold all user''s user folders for IDAE'})
                                                                                    return;         end
                                                                                                    end
%
% sorting out variables of ud.idae:
if isfield(ud.idae,'new_user');
    ud.idae.user_home           = fullfile(fileparts(ud.idae.user_home),ud.idae.user_name);
    ud.idae.real                = fullfile(ud.idae.user_home, 'iv2');                               end
%
ud.idae.real_c                  = local_str2cell(ud.idae.real);
ud.idae.symbolic                = fullfile(fileparts(ud.idae.user_home),'$User_ID$','iv2');
ud.idae.symbolic_c              = local_str2cell(ud.idae.symbolic);
%
% creating idae's home directory and subdirectories and
%   creating the files needed for performance of IDAE
%
if ~exist(ud.idae.real,'dir');  mkdir(ud.idae.real);                                                end;
%
d2mk                            = {'fs_scripts','tmp'};
for i=1:1:numel(d2mk);
    if ~exist(fullfile(ud.idae.real,d2mk{i}),'dir');
                                mkdir(fullfile(ud.idae.real,d2mk{i}));                      end;    end
%
f2mk                            = {'FS_done.txt','FS_submitted.txt'};
for i=1:1:numel(f2mk);
    if ~exist(fullfile(ud.idae.real,d2mk{1},f2mk{i}),'file');
        write2ptf(fullfile(ud.idae.real,d2mk{1},f2mk{i}), 'symbolic');                      end;    end
%
if ~exist(fullfile(ud.idae.real,d2mk{2},'scratch.m'),'file');
    write2ptf(fullfile(ud.idae.real,d2mk{2},'scratch.m'), 'symbolic');                              end
%
% start_IDAE.m to make IDAE codes available to the user
% 
if ~exist(fullfile(ud.idae.user_home,'start_IDAE.m'),'file');
    write2ptf(fullfile(ud.idae.user_home,'start_IDAE.m'),     ...
        ['% this file adds paths of IDAE codes',10, ...
        'addpath ',ud.idae.code_path,' -begin',10, ...
        'addpath ',fullfile(ud.idae.real,d2mk{2}),' -begin']);                                      end
%
set(gcf, 'UserData',ud);
%
if isfield(ud.idae,'new_user'); local_sort_etc(ud);                                 return;         end

local_get_is_ds(ud);
return;
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

% function                     	local_done_idae_1(ud);
% %%
% set(gco,    'Enable','off');
% if isfield(ud,'reuse') && isfield(ud.reuse,'idae') && ud.reuse.idae>0;
%                                 path_x                      = fileparts(ud.idae.real);
% else;                           [fname, path_x]             = uigetfile(fullfile(pwd,'*'));        
%     if ~ischar(fname);         	set(gco,    'Enable','on');                         return;         end;
%                                                                                                     end;
% %
% ud.idae.real                    = fullfile(path_x, '(enter a folder name without spaces)');
% ud.idae.real_c                  = local_path2sc(ud.idae.real);
% %                           
% %
% for i=1:1:numel(ud.idae.real_c)-1;
%    	set(ud.c1bHs(i+1), 'String',ud.idae.real_c{i}, 'CallBack','iv2_gen_lds(''idae_toggle'',[]);');	end;
% %
% set(ud.c1bHs(i+2), 'Style','edit',  'String',ud.idae.real_c{end}, 'BackgroundColor',iv2_bgcs(6));
% set(ud.c2bHs(i+2), 'String','IDAE-ID'); 
% %
% set(ud.c1bHs(i+3), 'String','Do not see your user name on left column GUIs? If No, hit >');
% set(ud.c2bHs(i+3), 'String','Restart',  'CallBack',['ud=get(gcf,''UserData''); ',  ...
%                     'ud.stage=''idae_0''; set(gcf,''UserData'',ud); iv2_gen_lds(''done'',[]);']); 
% %
% set(ud.c1bHs(1), 'String','Your home path segments. Follow below instructions. Then, hit >');
% set(ud.c2bHs(1), 'String','Done',   'CallBack','iv2_gen_lds(''done'',[]);');
% %
% ud.stage                        = 'idae_2';
% set(gcf,    'UserData',ud);
% %
% set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
%     {'* Complete below tasks: ','1. Hit left column GUI showing your user name (turns orange when done)',  ...
%     '2. Replace light green GUI with a non-existing folder name (no spaces please)',    ...
%     '   (See MATLAB command window for existing folders)'});
% %
% d                               = dir(fullfile(fileparts(ud.idae.real),'*'));
% d_name                          = char(d.name);
% disp(['.directories in: ',fileparts(ud.idae.real)]);
% dispCharArrays(1, char(d(d_name(:,1)~='.' & [d.isdir]'>0).name));
% disp('> do not enter above existing folders in light green GUI.');
% set(gco,    'Enable','on');
% return;
% %%

% function                     	local_done_idae_2(ud);
% %%
% im1                             = umo_cstrs(char(get(ud.c2bHs, 'String')), ['User';'IDAE'], 'im1');
% h0                              = findobj(gcf, 'Style','edit');
% set(h0,     'Style','pushbutton');
% % checking if User_ID is present:
% if im1(1)<1;    
%     for i=2:1:im1(2,1)-1;       set(ud.c1bHs(i),  'BackgroundColor',iv2_bgcs(11));
%                                 pause(0.5);
%                                 set(ud.c1bHs(i),  'BackgroundColor',iv2_bgcs(0));                   end
%     set(h0, 'Style','edit');                                                    	return;         end
% %
% % when spces are presetn in the IDAE segment;
% if im1(2)>0 && any(get(ud.c1bHs(im1(2)),'String')==' ');
%     set(ud.c1bHs(im1(2)), 'BackgroundColor',iv2_bgcs(11));
%     pause(0.5);
%     set(ud.c1bHs(im1(2)), 'BackgroundColor',iv2_bgcs(6));
%     set(h0, 'Style','edit');                                                    	return;         end
% 
% % constructing ud.idae.symbolic_c:
% ss_end                          = get(ud.c1bHs(im1(2,1)),    'String');
% ud.idae.real_c{end}             = ss_end(ss_end~=' ');
% ud.idae.real                  	= fullfile(fileparts(ud.idae.real), ud.idae.real_c{end});
% if exist(ud.idae.real,'dir');
%     set(ud.c1bHs(im1(2,1)),     'String',[ud.idae.real_c{end},' - exists. replace it'],     ...
%                                 'Style','edit', 'BackgroundColor',iv2_bgcs(11));
%  	pause(0.5);
%     set(ud.c1bHs(im1(2,1)),     'BackgroundColor',iv2_bgcs(6));                     return;         end;
% %
% s                               = mkdir(ud.idae.real);
% if s<1;
%     set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
%         {'* Unable to create your IDAE home directory ',    ...
%         '  (No full-access to your home directory in Data Server?)', 	...
%         '> Gathered information saved. ''Quit'' and Consult your local manager'});   
%     local_save([]);                                                                 return;         end;
% %
% ss                              = ud.idae.real_c;
% ss{im1(1,1)-1}                  = '$User_ID$';
% ss{end}                         = ss_end(ss_end~=' ');
% ud.idae.symbolic_c              = ss;
% ddd                             = ss{1};
% for i=2:1:numel(ss);            ddd                         = fullfile(ddd, ss{i});                 end;
% ud.idae.symbolic                = ddd;
% %
% ud.stage                        = 'idae_done';
% set(gcf,    'UserData',ud);
% set(ud.c1bHs(1), 'String','Done for IDAE home directory! Observe info-Board. Then, hit >')
% set(ud.c2bHs(1), 'String','Move on',   'CallBack','iv2_gen_lds(''done'',[]);');
% %
% set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
%     {'* Your IDAE home directory: ',[' ',ud.idae.real],             ...
%     '* Generalized IDAE home directory, shared by users with the same conventions:',	...
%     [' ',ud.idae.symbolic],    '> Hit ''Display'' for more info'});
% return;
% %%
% 
% function                     	local_done_idae_done(ud);
% %%
% % reset left & right columns:
% local_reset_guis(ud);
% 
% set(ud.c1bHs(1), 'String','Select a PET source file in Image Server');
% set(ud.c2bHs(1), 'String','Start',   'CallBack','iv2_gen_lds(''pet_is'',0);');
% % 
% set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',           ...
%     {'* Select a PET source file in Image Server to register it''s directory conventions',  ...
%     '  - any file is OK, so long as it is in a regular PET souce file folder'});
% return;
% %%

function                        local_search_files(i2);
%%
coUD                            = get(gco, 'UserData');
ud                              = get(gcf, 'UserData');
%
[fname, path_x]           	    = uigetfile(fullfile(pwd,coUD.fln),coUD.title);

if ~ischar(fname);         	    set(gco,    'Enable','on');                         return;         end
path_x                          = fileparts(fullfile(path_x,fname));
%
eval([coUD.var,'                = path_x;'])
set(coUD.bH, 'String','Done');
set(ud.c1bHs(coUD.bNo+1), 'String',path_x);
%
return;
%%

function                        local_sort_etc(ud);
% reset left & right columns:
local_reset_guis(ud);
%
%% SPM12
spm_home                        = which('spm_vol.m');
ic                              = 2;
if isempty(spm_home);       
    set(ud.c1bHs(ic),   'String','Find SPM12''s home directory using the file navigator');
    set(ud.c2bHs(ic),   'String','Start',   'CallBack','iv2_gen_lds(''search_files'',[]);', ...
        'UserData',struct('fln','spm_vol.m', 'var','ud.etc.spm_home', 'bH',ud.c2bHs(ic), 'bNo',ic,  ...
        'title','Identify spm_vol.m from SPM12''s main folder'));
    set(ud.c1bHs(ic+1), 'Style','text', 'String',' ', 'HorizontalAlignment','center');
else;                           
    ud.etc.spm_home             = fileparts(spm_home);
  	set(ud.c1bHs(ic),   'String','SPM12''s home directory');
 	set(ud.c2bHs(ic),   'String','Done');
    set(ud.c1bHs(ic+1), 'Style','text', 'String',ud.etc.spm_home, 'HorizontalAlignment','center');  end
%% Freesurfer
s                               = 1;
if ~exist(fullfile(ud.idae.real,'fs_scripts'),'dir');
    s                          	= mkdir(fullfile(ud.idae.real,'fs_scripts'));                       end
%
ic                              = ic + 2;
set(ud.c1bHs(ic),   'String','Folder for Freesuefer scripts');
%
if s>0;                         
 	% disp([' ',fullfile(ud.idae.real,'fs_scripts')]);
    ud.etc.fs_script_ws_real    = fullfile(ud.idae.real,'fs_scripts');
    ud.etc.fs_script_ws_symbolic= replace(ud.etc.fs_script_ws_real,ud.idae.user_name,'$User_ID$');
    % ud.etc.fs_script_ws_c       = ud.idae.symbolic_c;
    % ud.etc.fs_script_ws_c{end+1}                            = 'fs_scripts';     
    qqq                         = {'FS_done.txt', 'FS_submitted.txt'};
    for i=1:1:numel(qqq);
        write2ptf(fullfile(ud.etc.fs_script_ws_real,qqq{i}), 'symbolic');                      	    end
    set(ud.c2bHs(ic),   'String','Done');
else;                           
    % disp('> unable to create folder for Freesuefer scripts');
    % disp([' you have no full-access tp: ',ud.idae.real]);
    ud.etc.fs_script_ws_real    = ['Faile: No full-access to: ',ud.idae.real];
    set(ud.c2bHs(ic),   'String','Failed');                                                         end
%
set(ud.c1bHs(ic+1), 'Style','text', 'String',ud.etc.fs_script_ws_real, 'HorizontalAlignment','center');

% windows 10:
if any(ud.pet.ds.real=='\') && s>0;
    ic                          = ic + 2;
    set(ud.c1bHs(ic),   'String','Enter Freesurfer''s script folder seen from the Linux side');
    set(ud.c2bHs(ic),   'String','Start', 'CallBack','iv2_gen_lds(''copy_paste'',[]);');
    % 
    if isfield(ud.etc,'fs_script_ds_real');
        if isfield(ud.idae,'new_user');
            ud.etc.fs_script_ds_real    ...
                                = replace(ud.etc.fs_script_ds_symbolic,'$User_ID$', ...
                                                            ud.idae.user_name);                     end
        set(ud.c1bHs(ic+1), 'String',ud.etc.fs_script_ds_real,    ...
                                'Style','text',  'UserData','fs_script_ds');
        set(ud.c2bHs(ic),   'String','Done');
    else
        set(ud.c1bHs(ic+1), 'String','Copy & paste from the Linux machine', ...
                                'Style','text',  'UserData','fs_script_ds');                        end
else
    ud.etc.fs_script_ds_real    = ud.etc.fs_script_ws_real;
    ud.etc.fs_script_ds_symbolic= ud.etc.fs_script_ws_symbolic;                                     end

%% data unit:
% disp('.sorting out radioactivity units:');
ic                              = ic + 2;
set(ud.c1bHs(ic),   'String','Select the radioativity unit to use');
set(ud.c2bHs(ic),   'Value',1,  'Style','popupmenu',    'String',{'Select','nCi/mL','Bq/mL'},   ...
                                'CallBack','iv2_gen_lds(''unit_done'',[]);');
%
if isfield(ud.etc,'unit')
    set(ud.c2bHs(ic), 'Value',umo_cstrs(char({'Select','nCi/mL','Bq/mL'}),ud.etc.unit, 'im1'));     end
%
set(ud.c1bHs(1), 'String','A few items to go. Follow the instructions below')
set(ud.c2bHs(1), 'String','Done',   'CallBack','iv2_gen_lds(''done'',[]);');
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* Complete the items marked as ''Start'' (right) as instructed, and'
    '  select the radioactivity unit'
    '> Hit ''Done'' @top-right when all are done'});
%
ud.etc.variables                = {'spm_home','fs_script_ws_real','fs_script_ws_symbolic',  ...  
                                    'fs_script_ds_real','fs_script_ds_symbolic','unit'};
ud.stage                        = 'check_etc';
set(gcf, 'UserData',ud);
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

% function                        local_fs_script(i2);
% %%
% ud                              = get(gcf,      'UserData');
% set(gco,    'Enable','off');
% drawnow;
% if i2=='w';
%     disp(['> Freesurfer script folder (Windowds): ',ud.etc.fs_script_ws_real]);
%     disp('  Find it in you the Linux system to run Freesurfer (typically starts with /mnt)');
%     disp('  There should be FS_done.txt  FS_submitted.txt in the flder by now');
%     disp('> Copy & past it below');
%     s                           = input(' Enter the folder name seen from Linux (Ret=cancel): ','s');
%     if isempty(s);                                                                  return;         end;
%     im1                         = umo_cstrs(char(ud.idae.symbolic_c),'$User_ID', 'im1');
%     if im1(1)<1 || size(im1,2)>1;                                                   return;         end;
%     %
%     s(s=='/')                   = ' ';
%     sc                          = getLseg(s,    [0,2]);
%     sc{1}                       = ['/',sc{1}];
%     im2                         = umo_cstrs(char(sc),  ud.idae.real_c{im1}, 'im1');
%     if im2(1)<1 || size(im2,2)>1;                                                   return;         end;
%     %
%     sc{im2(1)}                  = '$User_ID$';
%     ud.etc.fs_script_ds_c       = sc;
%     out                         = sc{1};
%     for i=2:1:numel(sc);        out                         = [out,'/',sc{i}];                  	end;
%     ud.etc.fs_script_ds         = out;
%     
%     set(gcf,    'UserData',ud);
%     set(gco,    'String','Done',    'CallBack',' ', 'Enable','on');                 return;         end;
% return;
% %%
% 
% function                        local_fs_subject(i2);
% %%
% ud                              = get(gcf,      'UserData');
% set(gco,    'Enable','off');
% drawnow;
% if i2=='w';
%     disp('> Inquiring Freesurfer''s subject folder');
%     disp('  Go to Freesurfer''s subject folder in your Linux system');
%     disp('  (Tyoically try [your linux]: cd freesurfer/subject)');
%     disp('  Then, type [your linux]: pwd')
%     disp('  Copy & past the output below');
%   	s                           = input(' Freesurfer''s subject folder in Linux) (Ret=cancel): ','s');
%     if isempty(s);                                                                  return;         end;
%     im1                         = umo_cstrs(char(ud.idae.symbolic_c),'$User_ID', 'im1');
%     if im1(1)<1 || size(im1,2)>1;                                                   return;         end;
%     %
%     s(s=='/')                   = ' ';
%     sc                          = getLseg(s,    [0,2]);
%     sc{1}                       = ['/',sc{1}];
%     im2                         = umo_cstrs(char(sc),  ud.idae.real_c{im1}, 'im1');
%     if im2(1)<1 || size(im2,2)>1;                                                   return;         end;
%     %
%     sc{im2(1)}                  = '$User_ID$';
%     ud.etc.fs_subject_ds_c      = sc;
%     %
%     out                         = sc{1};
%     for i=2:1:numel(sc);        out                         = [out,'/',sc{i}];                  	end;
%     ud.etc.fs_subject_ds        = out;
%     set(gcf,    'UserData',ud);
%     set(gco,    'String','Done',    'CallBack',' ', 'Enable','on');                 return;         end;
% % 
% % Linux version
% disp(i2)
%     s                           = getenv('SUBJECTS_DIR');
%     ud.etc.fs_subject_ds        = s;
% 
%     s(s=='/')                   = ' ';
%     sc                          = getLseg(s,    [0,2]);
%     sc{1}                       = ['/',sc{1}];
%     ud.etc.fs_subject_ds_c      = sc;
%     set(gcf,    'UserData',ud);
%     set(gco,    'String','Done',    'CallBack',' ', 'Enable','on');              
% 
%   
% return;
% %%

function                        local_copy_paste(i2)
%%
ud                              = get(gcf,  'UserData');
rR                              = find(ud.c2bHs==gco);
% add to sss when a new string is added as a possibility: 
sss                             = ['fs_scri';'user_na'];
if ~isempty(rR)
    set(gco, 'String','In progress')
    set(ud.c1bHs(rR+1), 'Style','edit', 'CallBack','iv2_gen_lds(''copy_paste'',[]);');
    im1                         = umo_cstrs(sss, get(ud.c1bHs(rR+1),'UserData'),'im1');
    if im1(1)<0;                set(gco, 'String','???');                           return;         end
    if im1(1)==1
        set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
            {['* ',get(ud.c1bHs(rR), 'String')] 
            '> In the Linux machine with Freesurfer, find the directory that correspond to: '
            ['   ',get(ud.c1bHs(rR-1), 'String')]
            '  (you should be able to see FS_done.txt in the folder by now)'
            '> Copy the full directory to the designated GUI'});                                    end
                                                                                    return;         end
% done on a left-column GUI:
rL                              = find(ud.c1bHs==gco);
if isempty(rL);                                                                     return;         end
%
im1                             = umo_cstrs(sss, get(ud.c1bHs(rL),'UserData'),'im1');
if im1(1)<0;                    set(ud.c2bHs(rL-1), 'String','???');                return;         end

%
s1                              = get(gco, 'String');
% for fs_script_ds:
if im1(1)==1;
    ud.etc.fs_script_ds_real    = s1(s1~=' ');
    %
    im1                         = umo_cstrs(char(local_str2cell(ud.etc.fs_script_ws_real)), ...
                                    char(local_str2cell(ud.etc.fs_script_ds_real)), 'im1');
    % lower segments must all agree
    if ~any(im1>0) || any(im1(find(im1>0,1):end)<1) || ...
        ~contains(ud.etc.fs_script_ds_real,ud.idae.user_name)
        set(gco, 'BackgroundColor',iv2_bgcs(11));
        pause(0.5)
        set(gco, 'BackgroundColor',iv2_bgcs(0));                                    return;         end

    % replacing the username with '$User_ID$':
    ud.etc.fs_script_ds_symbolic= replace(ud.etc.fs_script_ds_real,ud.idae.user_name,'$User_ID$');
%
% for username:
elseif im1(1)==2;               ud.idae.user_name           = s1(s1~=' ');                          end;
%
set(ud.c2bHs(rL-1), 'String','Done')
set(gcf, 'UserData',ud);
return;
%%

function    out                 = local_str2cell(i2);
%%
% inserted to: iv2_setAproj.m/
%
if iscell(i2);
    out                         = i2{1};
    for i=2:1:numel(i2);        out                         = fullfile(out,i2{i});                  end
    if i2{1}(1)=='/';           out(out==filesep)           = '/';                                  end
                                                                                    return;         end
%
i2x                             = i2;
i2x(i2=='/' | i2=='\')          = ' ';
if i2(1)=='/' || i2(1)=='\';    i2x(1)                      = i2(1);                                end
out                             = getLseg(i2x, [0,2]);
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
fn2ck{6}                        = {'spm_home','fs_script_ws','fs_subject_ds','dump','unit'};
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
                                set(ud.c1bHs(1), 'BackgroundColor',iv2_bgcs(6));   return;         end;
%
iv2_gen_lds_2(ud);
% 
% save(fullfile(lds.iv2_home,['aid_',ud.lds,'.mat']),     'lds');
% disp('.done! (information for your local adapeter saved)');
% disp([' output: ',fullfile(lds.iv2_home,['aid_',ud.lds,'.mat'])]);    
return;
%%

function                        local_get_is_ds(ud)
%%
%
% ud                              = get(gcf, 'UserData');
local_reset_guis(ud);

set(ud.c1bHs(1), 'String','Complete below inquiries, after reviwing the info-board ')
set(ud.c2bHs(1), 'String','Move on',   'CallBack','iv2_gen_lds(''done'',[]);');

% PET source file in Image server:
ic                              = 2;
set(ud.c1bHs(ic), 'String','Select a representative PET source file in Image server');
if isfield(ud,'pet') && isfield(ud.pet,'is') && isfield(ud.pet.is,'real');
    set(ud.c2bHs(ic), 'String','Done', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');
    set(ud.c1bHs(ic+1), 'String',ud.pet.is.real, 'Style','text', ...
                                'UserData','ud.pet.is.real', 'HorizontalAlignment','center');
else;
    set(ud.c2bHs(ic), 'String','Start', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');
    set(ud.c1bHs(ic+1), 'String',' ', 'Style','text', ...
                                'UserData','ud.pet.is.real', 'HorizontalAlignment','center');       end
%
% PET analysis output file in Data server:
ic                              = 4;
set(ud.c1bHs(ic), 'String','Select a representative PET analysis output in Data server');
if isfield(ud,'pet') && isfield(ud.pet,'ds') && isfield(ud.pet.ds,'real');
    set(ud.c2bHs(ic), 'String','Done', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');
    set(ud.c1bHs(ic+1), 'String',ud.pet.ds.real, 'Style','text', ...
                                'UserData','ud.pet.ds.real', 'HorizontalAlignment','center');
else;
    set(ud.c2bHs(ic), 'String','Start', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');    
    set(ud.c1bHs(ic+1), 'String',' ', 'Style','text', ...
                                'UserData','ud.pet.ds.real', 'HorizontalAlignment','center');       end
%
% A MRI source file in Image server
ic                              = 6;
set(ud.c1bHs(ic), 'String','Select a representative MRI source file in Image server');
if isfield(ud,'mri') && isfield(ud.mri,'is') && isfield(ud.mri.is,'real');
    set(ud.c2bHs(ic), 'String','Done', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');
    set(ud.c1bHs(ic+1), 'String',ud.mri.is.real, 'Style','text', ...
                                'UserData','ud.mri.is.real', 'HorizontalAlignment','center');
else;
    set(ud.c2bHs(ic), 'String','Start', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');
    set(ud.c1bHs(ic+1), 'String',' ', 'Style','text', ...
                                'UserData','ud.mri.is.real', 'HorizontalAlignment','center');       end
%
% A MRI analysis output file in Data server:
ic                              = 8;
set(ud.c1bHs(ic), 'String','Select a representative MRI analysis output in Data server');
if isfield(ud,'mri') && isfield(ud.mri,'ds') && isfield(ud.mri.ds,'real');
    set(ud.c2bHs(ic), 'String','Done', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');
    set(ud.c1bHs(ic+1), 'String',ud.mri.ds.real, 'Style','text', ...
                                'UserData','ud.mri.ds.real', 'HorizontalAlignment','center');
else;
    set(ud.c2bHs(ic), 'String','Start', 'CallBack','iv2_gen_lds(''uigetfile'',[]);');    
    set(ud.c1bHs(ic+1), 'String',' ', 'Style','text', ...
                                'UserData','ud.mri.ds.real', 'HorizontalAlignment','center');       end
%
%
set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ...
    {'* The responses (= example folders) will be used to make it possible:'
    ' - To identify source files as seamlessly as possible using the study log file'
    '   which lists all scans of a study (manually maintained by the study team), and'
    ' - To follow the existing conventions when to set folders for outputs from IDAE'
    '* Acceptable formats for PET & MRI source files are:'
    '   DICOM or ECAT 7 (*.v) for PET & DICOM or Phillips PET/REC for MRI'})
%
ud.stage                        = 'get_is_ds';
set(gcf, 'UserData',ud);
return;
%%

function                        local_done_get_is_ds(ud);
%%
bNos                            = 2:2:8;
im1                             = umo_cstrs('Done', char(get(ud.c2bHs(bNos), 'String')), 'im1'); 
if any(im1<1);
    for i=find(im1'<1);         set(ud.c1bHs(bNos(i)),  'BackggroundColor',iv2_bgcs(11));
                                pause(0.5)
                                set(ud.c1bHs(bNos(i)),  'BackggroundColor',iv2_bgcs(0));            end
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ... 
        {'* Insufficient entries found (blinked in pink)'
        '> Complete them using the right-column GUIs'});                            return;         end
% 
ud.stage                        = 'pet_is_1'

local_pet_is(ud);

return;
%%

function                        local_uigetfile(i2);
%%
%
ud                              = get(gcf, 'UserData');

set(gco, 'String','In progress')
drawnow

rR                              = find(ud.c2bHs==gco,1);

[fname, path_x]          	    = uigetfile('*', get(ud.c1bHs(rR), 'String'));
if ~ischar(fname);                                                                  return;         end
%
%
path_s                          = fileparts(fullfile(path_x, fname));
set(ud.c1bHs(rR+1), 'String',path_s);
% 
sss                             = ['ud.pet.is';'ud.pet.ds';'ud.mri.is';'ud.mri.ds';'user_home'];
im1                             = umo_cstrs(sss, get(ud.c1bHs(rR+1), 'UserData'), 'im1');
if im1(1)<1;                    set(ud.c2bHs(rR), 'String','???');                  return;         end;
%
if im1(1)==1;                   ud.pet.is.real              = path_s;
elseif im1(1)==2;               ud.pet.ds.real              = path_s;
elseif im1(1)==3;               ud.mri.is.real              = path_s;
elseif im1(1)==4;               ud.mri.ds.real              = path_s; 
elseif im1(1)==5;               ud.idae.user_home           = path_s;                               
                                ud.idae.real                = fullfile(path_s,'iv2');
                                set(ud.c2bHs(rR+2), 'String','Done');
                                set(ud.c1bHs(rR+3), 'String',ud.idae.real);                         end;
%
set(ud.c2bHs(rR), 'String','Done');
set(gcf, 'UserData',ud);
return
%%

function                        local_done_check_etc(ud)
%%
im1                             = umo_cstrs(char(get(ud.c2bHs,'style')), 'popupmenu', 'im1');

ok                              = 1;
if ~isfield(ud.etc,'unit');     set(ud.c2bHs(im1), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5)
                                set(ud.c2bHs(im1), 'BackgroundColor',iv2_bgcs(0));
                                ok                          = 0;                                    end
%
im2                             = umo_cstrs('Done', char(get(ud.c2bHs(2:2:im1(1)-1),'String')), 'im1');
if any(im2<1);
    for i=find(im2'<1).*2;      set(ud.c2bHs(i), 'BackgroundColor',iv2_bgcs(11));
                                pause(0.5)
                                set(ud.c2bHs(i), 'BackgroundColor',iv2_bgcs(0));                    end
                                ok                          = 0;                                    end
%
if ok<1; 
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ... 
        {'* Insufficient entries found (blinked in pink)'
        '> Complete them using the right-column GUIs'});                            return;         end
%
%
if isfield(ud.idae,'new_user')
    set(findobj(gcf, 'Tag','iv2_gen_lds_infoB'),    'String',   ... 
        {['* IDAE is ready for user: ',ud.idae.user_name]
        [' - IDAE home directory: ',ud.idae.real]
        ' - Type as follows in Matlab''s command window when to use IDAE'
        ['  >> run ',fullfile(ud.idae.user_home,'start_IDAE.m')]});                 return;         end
    

return;
%%