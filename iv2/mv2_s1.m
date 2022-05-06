function    mv2_s1(i1,i2); 

% To set up prepMP for a project (iProj). Need to have scanDB.m and iProj.iv2. 
%       
%       usage 1:    mv2_s1('full/path/iproj.iv2','usr')
%       
%   to set-up initial GUI windows.
%
%       usage 2:    mv2_s1('fun',[]);
%       
%   called from GUIs of windows that are created by this code
%   
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if size(i1,2)==3;               feval(['local_',i1], 	double(gcf),double(gco));   return;         end;
%
[idx, inm, iex]                 = fileparts(i1);
if strcmpi(iex,'.iv2');         local_c01(i1,   i2);
else;
end;

return;
%%

function                        local_c01(i1,   i2);
%%
h                               = findobj('Tag',            'prepMP Selector');
if ~isempty(h);                 figure(h(1));                                       return;         end;
[lds, mri, spec, usr, s1] = gei(i1,	'lds4dxetc','haveMRI',  'subjSpec','IDAEuser','condNames');
% special if it is a dosimatry study:
if size(s1,1)==1 && strcmpi(s1,'dosimetry');
                                [idx, inm]                  = fileparts(i1);
                                dm2_genDMDASet('set', fullfile(idx,inm));           return;         end;
                                
if (ischar(mri) && mri~='1') || (isnumeric(mri) && mri~=1);                    
                                disp('.not applicable .. No MRI');                  return;         end;
if spec~='h';                   disp('.not applicable .. a Non-human study');       return;         end;
if isempty(lds);                disp(['Not ready for iv2 ... ',i1]);                return;         end;
idx                             = feval(lds,'idx',          usr);

[idx, ipj]                      = fileparts(i1);

% reading coms file, if present
coms                            = fullfile(idx,             [ipj,'.ev2']);
if exist(coms,'file');          
    [c1, c2]                    = umo_getptf(coms,0,        1);
    c3                          = char(zeros(size(c1,1),1)  + 32);
    c3x                         = find(c1(:,1)=='#');
    for i=1:1:length(c3x);      c3(c3x(i):end,  :)          = char(abs('A') - 1 + i);               end;
    im1                         = umo_cstrs(c1,'prepMP',    'im1');
    for i=1:1:length(im1);
        if any(c1(im1(i),:)=='_');
            c31{i}              = [deblank(c1(im1(i),:)),')'];
            c31{i}(c31{i}=='_') = '(';
        else;
            c31{i}              = [deblank(c1(im1(i),:)),' (',c3(im1(i)),')'];                      end;
        c32{i}                  = deblank(c2(im1(i),    :));                                        end;
    gini                        = char(sum(c1(:,1)=='#')+'A');
else;                           
    c1                          = [];
    c31                         = [];
    gini                        = 'A';                                                              end;
%
str1                            = {['Adding prepMP to Group ',gini],    'Quit'};  
%    
odx                             = fullfile(idx,ipj,         'iv2');
if ~exist(odx,'dir');           mkdir(odx);                                                         end;

fff                             = iv2_prepMP([],[],         'set','on');
%
nrs                             = 1 + 1 + numel(fff.set) + ~isempty(c1) + numel(c31) + 1;
[fNo, bpos]                     = bWindow([], ...
                                'nrs',                      nrs,                ...
                                'bwd',                      600,                ...
                                'ttl',                      'prepMP Selector');
set(fNo,                        'CloseRequestFcn',          ' ',                ...
                                'tag',                      'prepMP Selector',  ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');
% First row GUIs:
ic                              = 1;
jHs                             = postJBs(fNo,              'B',bpos(ic,:),[4,1;1,1]);
for i=1:1:2;
    set(jHs(i),                 'String',                   str1{i},            ...
                                'BackgroundColor',          iv2_bgcs(1));                           end;
set(jHs(2),                     'CallBack',                 'delete(gcf);'); 
% Second row GUIs:
ic                              = ic + 1;
jHs                             = postJBs(fNo,              'B',bpos(ic,:),[1;1]);
set(jHs(1),                     'BackgroundColor',          iv2_bgcs(2),        ...
                                'String',   'Available Stage-1 packages to choose from');
% displaying selections:
for i=1:1:numel(fff.set);
    ic                          = ic + 1;
    jHs                         = postJBs(fNo,              'B',bpos(ic,:), [1,4;1,1]);
    set(jHs(1),                 'String',                   fff.set{i},         ...
                                'BackgroundColor',          iv2_bgcs(2),        ...
                                'Tag',                      'mv2_s1_prepMP',    ...
                                'Style',                    'radiobutton',      ...
                                'CallBack',                 'mv2_s1(''c02'',[]);');
    set(jHs(2),                 'Tag',                      [fff.set{i},'_d'],  ...
                                'String',                   fff.descrip{i});                        end;
if ~isempty(c1);
    ic                          = ic + 1;
    jHs                         = postJBs(fNo,              'B',bpos(ic,:),[1;1]);
    set(jHs(1),                 'BackgroundColor',          iv2_bgcs(6),        ...
                                'String',                   'Existing prepMPs (group)');            end;
% copying from coms file:
for i=1:1:numel(c31);           
    ic                          = ic + 1;
    jHs                         = postJBs(fNo,              'B',bpos(ic,:), [1,4;1,1]);
    set(jHs(1),                 'String',                   c31{i},             ...
                                'BackgroundColor',          iv2_bgcs(6));
    set(jHs(2),                 'String',                   c32{i});                                end;
% the last row GUIs:
ic                              = ic + 1;
jHs                             = postJBs(fNo,              'B',bpos(ic,:), [4,1;1,1]);
strn                            = 'Hit a Column 1 GUI to select. If OK,  hit ''Done'' >';
set(jHs(1),                     'String',                   strn,           ...
                                'BackgroundColor',          iv2_bgcs(12));
set(jHs(2),                     'String',                   'Done',         ...
                                'BackgroundColor',          iv2_bgcs(12),   ...
                                'CallBack',                 'mv2_s1(''c03'',[]);');
%
pflg                            = {'nocrop',                'crop'};
cwUD                            = struct('gini',gini,       'iv2',i1,           ...
                                'odx',odx,                  'ev2',coms,         'lds',lds);
%                                 'pflg',                     pflg{hrrt(1)+1});
% cwUD.c31                        = c31;
% cwUD.c32                        = c32;
% cwUD.pio                        = feval(lds,'pio',          []);
%                                
% cwUD
set(fNo,                        'userData',                 cwUD);
return;
%%

function                        local_c02(fNo,  oNo);
%%
tstr                            = get(oNo,                  'Tag');
if get(oNo,'Value')==1;
    hs                          = findobj(fNo,  'Tag',      tstr);
    set(hs(hs~=oNo),            'Value',                    0);                                     end;

return;
%%

function                        local_c03(fNo,  oNo);
%% prepMPxxx selected > making the iPack:
h                               = findobj(gcf,              'Tag','mv2_s1_prepMP',  'Value',1);
if numel(h)~=1;                                                                   	return;         end;
cwUD                            = get(gcf,                  'userData');

ipk                             = get(h,                    'String');
%
cwUD.iPack                      = ipk;
cwUD.descrip                    = get(findobj(gcf,  'Tag',[ipk,'_d']),              'String');
if isempty(cwUD.descrip);                                                         	return;         end;
% retreiving lines for the iPac (=x)
x                               = iv2_prepMP(ipk,           []);
% retreiving selections for ???
y                               = iv2_prepMP([],[],         'sst','on', 'lds',cwUD.lds);
% ignored: pio is given in dxetc4xxx (=lds)
fnms                            = fieldnames(y);
f2u                             = ones(numel(fnms),     1);
[rs, cs]                        = find(char(fnms)=='_');
f2u(rs, :)                      = 0;
for i=find(f2u>0)';             f2u(i, :)                   = numel(eval(['y.',fnms{i}]));          end;
% 
delete(gcf);
nrs                             = 1 + sum(f2u>0) + sum(f2u) + 1;
[fNo, bpos]                     = bWindow([],               ...
                                'nrs',                      nrs,                ...
                                'bwd',                      480,                ...
                                'ttl',                      'prepMP - Final Step');
set(fNo,                        'CloseRequestFcn',          ' ',                ...
                                'tag',                      'prepMP - Final',   ...
                                'Toolbar',                  'none',             ...
                                'Menubar',                  'none');
% First row GUIs:
ic                              = 1;
str1                            = {'Select one @ 1st column for each variable',  'Quit'};
jHs                             = postJBs(fNo,              'B',bpos(ic,:),[5,1;1,1]);
for i=1:1:2;
    set(jHs(i),                 'String',                   str1{i},            ...
                                'BackgroundColor',          iv2_bgcs(1));                           end;
set(jHs(1),                     'Tag',                      'prepMP_infoB');
set(jHs(2),                     'CallBack',                 'delete(gcf);'); 
% 
for i=find(f2u>0)';
    ic                          = ic + 1;
    jHs                         = postJBs(fNo,              'B',bpos(ic,:),[1;1]);
    set(jHs(1),                 'String',   eval(['y.',fnms{i},'_desc']),               ...
                                'FontWeight',               'bold',                     ...
                                'BackgroundColor',          iv2_bgcs(6)); 
    for j=1:1:f2u(i);
        ic                      = ic + 1;
        jHs                     = postJBs(fNo,              'B',bpos(ic,:),[1,5;1,1]);
        set(jHs(1),             'String',eval(['y.',fnms{i},'{j}']),                    ...
                                'Style',                    'radiobutton',              ...
                                'Tag',                      ['prepMP_',fnms{i}],        ...
                                'BackgroundColor',          iv2_bgcs(2),                ... 
                                'CallBack',                 'mv2_s1(''c04'',[]);');
        set(jHs(2),             'String',eval(['y.',fnms{i},'_tips{j}']));                  end;    end;
%
% the last row GUIs:
ic                              = ic + 1;
jHs                             = postJBs(fNo,              'B',bpos(ic,:), [5,1;1,1]);
set(jHs(1),                     'String','Check selections. If OK  hit ''Done'' >',  	...
                                'FontWeight',               'bold',                     ...
                                'BackgroundColor',          iv2_bgcs(12));
set(jHs(2),                     'String',                   'Done',                     ...
                                'BackgroundColor',          iv2_bgcs(12),               ...
                                'CallBack',                 'mv2_s1(''c05'',[]);');
%
cwUD.iv2_lines                  = x;
cwUD.f2c                        = fnms(f2u>0);
set(gcf,    'UserData',cwUD);
return;
%%

function                        local_c04(fNo,  oNo);
%%
set(findobj(gcf,  'Tag',get(gco,'Tag')),      'Value',0);
set(gco,    'Value',1);
return;
%%

function                        local_c05(fNo,  oNo);
%% Selections on pio/avr etc are done:
if ~strcmpi(get(oNo,'String'),'Done');                                              return;         end;
set(gco,    'Enable','off');
%
cwUD                            = get(gcf,                  'userData');
h                               = findobj(gcf,  'Style','radiobutton',  'Value',1);
if numel(h)<numel(cwUD.f2c);
    hx                          = findobj(gcf,  'Tag','prepMP_infoB');
    set(hx(1),      'BackgroundColor',iv2_bgcs(11));
    pause(1);
    set(hx(1),      'BackgroundColor',iv2_bgcs(1));
    set(gco,        'Enable','on');                                                 return;         end;
% 
cropORnot                       = get(findobj(gcf, 'Tag','prepMP_crop', 'Value',1), 'String');
w2u                             = zeros(2,  2);
ic                              = 0;
for i=numel(cwUD.iv2_lines)-[1:-1:0];
    ic                          = ic + 1;
    w2u(ic, :)                  = [i, double(contains(cwUD.iv2_lines{i}{1},'No-cropping'))];        end;
%
iv2_lines                       = cwUD.iv2_lines(1:end-2);
if strcmpi(cropORnot,'crop');   k                           = w2u(w2u(:,2)<1, 1);
else;                           k                           = w2u(w2u(:,2)>0, 1);                   end;
for i=1:1:numel(cwUD.iv2_lines{k});
                                iv2_lines{end+1}            = cwUD.iv2_lines{k}{i};                 end;
for i=1:1:numel(iv2_lines);
    if contains(iv2_lines{i},'???');
        c1                      = getLseg(iv2_lines{i}, 1);
        iv2_lines{i}            = [iv2_lines{i}(1,  1:find(iv2_lines{i}=='?',1)-1),     ...
                                    get(findobj(gcf, 'Tag',['prepMP_',c1],  'Value',1), 'String')]; 
                                                                                            end;    end;
%
delete(gcf);
%
fH                              = fopen(cwUD.ev2,           'a');
fwrite(fH,                      ['#',cwUD.gini,' ',cwUD.descrip, 10],               'char');
fwrite(fH,                      [ cwUD.iPack,' ',cwUD.descrip, 10],                 'char');
fclose(fH);
%
if ~exist(cwUD.odx,'dir');      mkdir(cwUD.odx);                                                    end;
fH                              = fopen(fullfile(cwUD.odx, [cwUD.iPack,'.m']),  'w');
for i=1:1:numel(iv2_lines);     fwrite(fH,  [iv2_lines{i},10],  'char');                            end;
fclose(fH);
disp('.new analysis package is ready');
mv2_startIDAE(fullfile(cwUD.odx, [cwUD.iPack,'.m']),      	cwUD.iv2);
return;
%%

