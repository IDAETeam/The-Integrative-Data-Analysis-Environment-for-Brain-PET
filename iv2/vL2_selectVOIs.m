function    vL2_selectVOIs(i1,i2, varargin);  

% vL2_selectVOIs:  To set-up a GUI window listing VOIs (for external tasks)
%       
%       usage:      vL2_selectVOIs(ListOfVOIIDNos,info)
%            
%   VOIIDNos  	a list (n by 1) of VOIIDNos to display 
%   info       	[figure# etc to use when called from outside]
%
% Generated window:
%   top row         strval,L,R,W    handels are stored as cwUD.bHs    
%   main rows       anatomical label,L,R,W          (cwUD.gHs)
%   bottom rows     nbr (see below) rows of GUIs    (cwUD.qHs)
%
% Options:
%   ttl,val    	Figure Name and Tag     (default: 'vL2 VOI selector')
%   tjb,val     what to do when figures of the same tag exist
%            	(default: 'figure(h(1)); return;')
%   str,val    	the string to display on the 1st top-row GUI
%   nbr,val    	# of bottom rows to add (default: 1)
%   brc,val    	how to prepare bottom row GUIs 
%             	(default: [1;3] for 3 equal columns)
%   bwd,val   	to specify the width of the window (default: 360)
%   mnc,val     to set minimal number of columns (dfault: depends)
%              	'mnc',2 will make # of columns 2 (inclusive) or more
%   snm,'on'    to add short names (default: long labels alone)
%
% Notes:
%   vL2_selectVOIs('fun',[]) is also valid.
%               see inside for available functions or use get_localFun.m
% (cL)2014    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if isempty(i2);                 feval(['local_',lower(i1)],double(gcf),double(gco));return;         end;

ttlval                          = 'vL2 VOI selector';
tjbval                          = 'figure(h(1)); return;';
strval                          = 'Existing VOIs';
nbrval                          = 1;
brcval                          = [1;3];
bwdval                          = 240+20.*3;
mncval                          = 1;
snmval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

h                               = findobj('Tag',            ttlval);
if ~isempty(h);                 eval(tjbval);                                    	return;         end;

v0                              = i1(:);
v2                              = consolidVOINos(v0, []);
v3                              = zeros(size(v2,1),     4);
for i=1:1:3;                    v3(:,   [1,i+1])            = consolidVOINos(v0,v2+(i-1).*100);     end;
v3(:, 1)                        = v2;
vv2                             = VOIdef(v2);
if strcmpi(snmval,'on');
    for i=1:1:size(vv2.anm,1);  
        vv2c{i}               	= [deblank(vv2.anm(i, :)),' (',deblank(vv2.snm(i, :)),')'];         end;
    vv2.anm                     = char(vv2c);                                                       end;
%
vv2.anm(:, 1)                   = upper(vv2.anm(:,1));
[vv3, is]                       = sortrows(vv2.anm);

v2(:)                           = v2(is);
cwUD                            = struct('gHs',zeros(size(v2,1), 3), 'vnos',[v2+100, v2+200, v2],   ...
                                                            'vss',double(v3(is, [3,4,2])>0));
                                
% vnos                            = 
% v3(:)                           = v3(is,    [1,3,4,2])
% 
% vi                              = 
% 
% 
% v4                              = zeros(size(v0,1),         4);
% vvv                             = VOIdef(v0);
% for i=1:1:size(vvv.anm,1);      v4(i,   :)                  = VOIdef(vvv.anm(i, :));                end;
% % v4
% vv0                             = sum(v4(:, [1,3,4]),       2);
% vnos                            = zeros(max(vv0),           1);
% vnos(vv0,   :)                  = 1;
% vv1                             = find(vnos);
% vv2                             = VOIdef(vv1);
% if strcmpi(snmval,'on');
%     for i=1:1:size(vv2.anm,1);  
%         vv2c{i}               	= [deblank(vv2.anm(i, :)),' (',deblank(vv2.snm(i, :)),')'];         end;
%     vv2.anm                     = char(vv2c);                                                       end;
% [vv3, is]                       = sortrows(vv2.anm);
% vv1(:)                          = vv1(is,   :);
% 
% vnos                            = [v2+100, v2+200, v2];
% vsss                            = zeros(size(vnos));
% cwUD.gHs                        = zeros(size(vnos));
% for i=1:1:3;                    vv                          = consolidVOINos(v0,vnos(:,i));
%                                 vsss(:, i)                  = vv(:, 2);                             end;
% cwUD.vss                        = vsss;
% cwUD.vnos                       = vnos;
%                            
nrs                             = 1 + size(v2,1) + nbrval(1);
if ceil(nrs./2)~=floor(nrs./2); nrs                         = nrs + 1;                              end;
[wNo, bpos]                     = bWindow([], ...
                                'nrs',                      nrs,                ...
                                'bwd',                      bwdval,             ...
                                'ttl',                      ttlval,             ...
                                'mnc',                      mncval);

set(wNo,'ToolBar','none',       'MenuBar','none',           'Tag',ttlval,       ...
                                'Resize','off',             'CloseRequestFcn',' ');
%
% first row GUIs ---------------------------------------------------------:
ic                              = 1;
bbb                             = [12,3;1,3];
bHs                             = postJBs(wNo,              'B',bpos(ic,:), bbb);
set(bHs(1),                     'String',                   strval,             ...
                                'Tag',                      'infoB4voiSelector',...
                                'Backgroundcolor',          iv2_bgcs(1));
sst                             = 'LRW';
for i=1:1:3;    
    set(bHs(i+1),               'String',                   sst(i),             ...
                                'userData',                 0,                  ...
                                'Backgroundcolor',          iv2_bgcs(1),        ...
                                'Callback',                 'vL2_selectVOIs(''bycolumn'',[]);');    end;
                                
cwUD.bHs                        = bHs;
set(cwUD.bHs(end),              'userData',                 1);
% VOI GUIs ---------------------------------------------------------------:
for i=1:1:size(v2,1);
    ic                          = ic + 1;
    bHs                         = postJBs(wNo,              'B',bpos(ic,:), bbb);
    cwUD.gHs(i, :)              = bHs(2:end)';
    set(bHs(1),                 'String',                   deblank(vv3(i, :)));                    end;
%
set(cwUD.gHs(cwUD.vss==0),    	'Enable','off');
set(cwUD.gHs,                   'Style','radiobutton');
set(cwUD.gHs(cwUD.vss(:,3)>0, 3),       'Value',1);

% Bottom row GUIs --------------------------------------------------------:
cwUD.qHs                        = zeros(nbrval(1),          sum(brcval(2,:)));
for i=1:1:nbrval(1);
    ic                          = ic + 1;
    cwUD.qHs(i, :)              = postJBs(wNo,              'B',bpos(ic,:), brcval)';               end;
%
cwUD.i2                         = i2;
set(wNo,                        'userData',                 cwUD);
% ic
% size(bpos,1)
if ic==size(bpos,1);                                                                return;         end;
% repeating labels one more time since bpos(end,:) is availble:
bHs                             = postJBs(wNo,              'B',bpos(end,:), bbb);
set(bHs(1),                     'String',                   strval,                 ...
                                'Tag',                      'infoB4voiSelector_2',  ...
                                'Backgroundcolor',          iv2_bgcs(1));
sst                             = 'LRW';
for i=1:1:3;    
    set(bHs(i+1),               'String',                   sst(i),                 ...
                                'userData',                 0,                      ...
                                'Backgroundcolor',          iv2_bgcs(1),            ...
                                'Callback',                 'vL2_selectVOIs(''bycolumn'',[]);');    end;
%
return;
%%

function                        local_bycolumn(fNo,oNo);
%%
ii                              = find(get(g0o,'String')=='LRW');
if length(ii)~=1;                                                                   return;         end;
cwUD                            = get(fNo,                  'userData');
if any(cell2mat(get(cwUD.gHs(:,ii(1)),'Value'))==1);
    set(cwUD.gHs(:,ii(1)),      'Value',0);
else;
    im1                         = umo_cstrs(char(get(cwUD.gHs(:,ii(1)),'Enable')),'on ','im1');
    set(cwUD.gHs(im1,ii(1)),    'Value',1);                                                         end;
return;
%%

function                        local_set4voisets(fNo, gNo);
%% set VOI sets (can be used from other figures as well

h                               = findobj(groot, 'Tag','vL2 VOI selector');
if isempty(h);                                                                      return;         end;
figure(h);
%
set(findobj(gcf,  'Tag','infoB4voiSelector'),   'Value',1,  'Style','popupmenu',    'String',       ...
    {'Ready to set region sets',' - highlight all VOIs to include',                                 ...
    ' - hit L/R/W GUI to select/deselect by column',' * Select this tab when done'},                ...
  	'UserData',[],  'BackgroundColor',iv2_bgcs(11),'CallBack','vL2_selectVOIs(''set_selected_s1'',[]);'); 
return;
%%

function                        local_set_selected_s1(fNo, gNo);
%%
disp('.not implemented');

return;
%%
