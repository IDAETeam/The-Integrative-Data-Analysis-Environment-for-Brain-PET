function    vL2_defACPC(i1); 

% vL2_defACPC:      
%       
%       usage:      open a vL2Land session. Then ...
%                   vL2_defACPC('output.acpc')
%       
% 
% (cL)2010    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

if ischar(i1);                  local_s0(i1);
elseif i1(1)==1;
    fff                         = ['local_s1v',int2str(get(gco, 'Value'))];
    if ~isempty(which(fff));    feval(fff,i1);                                                      end;
elseif i1(1)==0;                
    str                         = get(gco,  'String');
    im1                         = umo_cstrs(['record AC';'show AC  ';'record PC';
                                                                    'show PC  ';'Save     '],str,'im1');
    if im1;                     feval(['local_s1v',int2str(im1+1)], i1);                            end;
else;                           local_s2v0(i1);                                                     end;
return;
%%

function                        local_s0(i1);
%% setting up defACPC GUI window:

clear global g4vL2defACPC;
global g4vL2;
set(g4vL2{double(gcf)}.iHs,   	'ButtonDownFcn',            'vL2_getXYc(1,''vL2_defACPC(2)'');');
local_s1v1(0);

return;
%%

function                        local_s2v0(i2);
%%

global g4vL2defACPC g4vL2;
vNo                             = find(g4vL2{double(gcf)}.aHs==gca);
if ~vNo;                                                                            return;         end;

% whdn AC or PC coordinaets were entered:
xyz                             = zeros(1,      3);
if size(i2,2)==3;               xyz(:)                      = round(i2);
else;                           rxyz                        = [1,2,3;   1,3,2;  2,3,1];
                                p                           = get(gca,  'CurrentPoint');
                                xyz(:,  rxyz(vNo,1:2))      = round(p(1,1:2));
                                xyz(:,  rxyz(vNo,3))        = g4vL2{double(gcf)}.inos(rxyz(vNo, 3));
                                i2                          = xyz;
                                i2(:,   rxyz(vNo,1:2))      = p(1,  1:2);                           end;

g4vL2{double(gcf)}.inos      	= xyz;
vL2_IJs('updatetM',1);
vL2_IJs('updatecM',1);
vL2_IJs('updatesM',1);

need2plot                       = 0;
if isempty(g4vL2defACPC);       need2plot                   = 1;
else;
    if isempty(g4vL2defACPC{double(gcf)});
                                need2plot                   = 1;                            end;    end;
rxy                             = [1,2;     1,3;    2,3];
if need2plot;
    g4vL2defACPC{double(gcf)}.pHs                           = zeros(3,      1);
    for i=1:1:3;                    
        set(gcf,'CurrentAxes',  g4vL2{double(gcf)}.aHs(i));
        g4vL2defACPC{double(gcf)}.pHs(i, :)               	= plot(i2(1,rxy(i,1)),i2(1,rxy(i,2)),'r.');
                                                                                                    end;
else;
    for i=1:1:3;                
        set(g4vL2defACPC{double(gcf)}.pHs(i),'xData',i2(1,rxy(i,1)),'yData',i2(1,rxy(i,2)));      	end;
                                                                                                    end;
return;
%%

function                        local_s1v1(i2);
%% showing info on the infoBoard

bH                              = findobj(gcf,'Tag',        'vL2InfoB');
if isempty(bH);                                                                     return;         end;

set(bH,'String',                [10,'  Defining ACPC points',10,10,  ...
                                '  A red dot appears when point & click LtMBut.',       10, ...
                                '  1.  Bring the red dot at AC (exact)',                10, ...
                                '      Point & click or using the number key pad',      10, ...
                                '      (arrows [Trans-axial] & PgUp/Dn in Z)',          10, ...
                                '  2.  Select ''record AC'' to record AC in memory',    10, ...
                                '  3.  Select ''show AC'' to review AC, if desired',    10, ...
                                '  4.  Repeat the same for PC.',                        10, ...
                                '  5.  Select ''save'' to save ACPC points to the file'],   ...
                                'FontName',                 'Courier New');

return;
%%

function                        local_s1v2(i2);
%% recording AC point
oNo                             = gco;
h                               = findobj(gcf,'String',     'Save');
if isempty(h);                                                                      return;         end;
ud                              = get(h,    'UserData'); 
global g4vL2defACPC;
ud.ac                           = [get(g4vL2defACPC{double(gcf)}.pHs(1),'xData'),   ...
                                    get(g4vL2defACPC{double(gcf)}.pHs(1),'yData'),  ...
                                    get(g4vL2defACPC{double(gcf)}.pHs(2),'yData')];
set(h,'UserData',               ud);
set(oNo,    'BackgroundColor',  iv2_bgcs(6));
return;
%%

function                        local_s1v3(i2);
%% displaying AC Point

h                               = findobj(gcf,'String',     'Save');
if isempty(h);                                                                      return;         end;
ud                              = get(h,    'UserData'); 
if isempty(ud);                                                                     return;         end;
if ~isfield(ud,'ac');                                                               return;         end;
if isempty(ud.ac);                                                                  return;         end;
% 
local_s2v0(ud.ac);

return;
%%

function                        local_s1v4(i2);
%% recording PC point

oNo                             = gco;
h                               = findobj(gcf,'String',     'Save');
if isempty(h);                                                                      return;         end;
ud                              = get(h,    'UserData'); 

global g4vL2defACPC;
ud.pc                           = [get(g4vL2defACPC{double(gcf)}.pHs(1),'xData'),   ...
                                    get(g4vL2defACPC{double(gcf)}.pHs(1),'yData'),  ...
                                    get(g4vL2defACPC{double(gcf)}.pHs(2),'yData')];
set(h,'UserData',               ud);
set(oNo,    'BackgroundColor',  iv2_bgcs(6));
return;
%%

function                        local_s1v5(i2);
%% displaying PC Point

h                               = findobj(gcf,'String',     'Save');
if isempty(h);                                                                      return;         end;
ud                              = get(h,                    'UserData');
if isempty(ud);                                                                     return;         end;
if ~isfield(ud,'pc');                                                               return;         end;
if isempty(ud.pc);                                                                  return;         end;

local_s2v0(ud.pc);

return;
%%

function                        local_s1v6(i2);
%% saving ACPC points to the file:

ud                              = get(gco,                  'UserData');
fnms                            = fieldnames(ud);
fstr                            = str2mat('fname','ac','pc','ifln');
im1                             = umo_cstrs(char(fnms),fstr,'im1');
if any(~im1);                                                                       return;         end;

si                              = struct('h2s',32,'c',mfilename,'p',ud.ifln,'cp','m');
ststus                          = um_save(ud.fname,[ud.ac; ud.pc],si,[],    ...
                                'dataUnit',                 'voxel',            ...
                                'ACPCPoints',               [ud.ac; ud.pc]);

[odx, onm, oex]                 = fileparts(ud.fname);
set(findobj(gcf,'Tag','vL2InfoB'),  'String',   {' ',' Info: ACPC points were saved ..',    ...
    ['  folder: ',odx],['    file: ',onm,oex],' now OK to exit'},   'FontName','Courier New');
clear global g4vL2defACPC
return;
%%