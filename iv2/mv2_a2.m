function                        mv2_a2(i1,i2); 

% To perform tasks of first row GUIs of IDAE Level 2 Window (L2W)
%       
%       usage:      cst     = mv2_a2('fun',[L2W,gco])
%
%   fun     string of any of 1strow GUIs or task class GUIs
%          	[] is valid if both gcf & gco are correct
%   L2W    	Matlab Fig# of IDAE level 2 GUI window
%   gco    	the GUI handle of tha task to do
%
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
if isempty(i1);                 i1                          = get(gco,  'String');                  end;
if nargin==1;                   i2                          = [double(gcf),     double(gco)];       end;
if numel(i1)==1;                local_highlight(i2);    
else;                           feval(['local_',lower(i1)], i2);                                    end;
return;

function                        local_highlight(fos);
%% highligh/de-highlight task class GUIs:

bgcs                            = [iv2_bgcs(0); iv2_bgcs(7); get(fos(2),'BackgroundColor')];
[v, imax]                       = max(sqrt(((bgcs(1:2, :) - bgcs([3;3],  :)).^2)*ones(3,1)));
set(fos(2),'BackgroundColor',   bgcs(imax,  :));

return;
%%

function                        local_update(fos);
%%
global g4iv2;

cwUD                            = get(fos(1),              	'userData');
if isempty(cwUD);                                                                   return;         end;
% updating .fck, .ick, & .ock for one ibase (=cwUD(3)):
mv2_fck(cwUD(1),cwUD(2),        find(g4iv2.ppp(:,4)==cwUD(3)));
mv2_a0(fos(1));

h                               = findobj(fos(1),           'Tag','iv2_L2W_taskGUIs');
if ~isempty(h);                 set(h,'Enable', 'on');                                              end;

return;
%%

function                        local_perform(fos);
%% called from L2W's 'perform' GUI
global g4iv2;

% cwUD  = [L1W figure handle, subject#, iBase#]
cwUD                            = get(fos(1),               'userData');
if isempty(cwUD);                                                                   return;         end;

h                               = findobj(fos(1),           'String','Task descriptions');
if isempty(h);                                                                      return;         end;
coUD                            = get(h,                    'userData');
im1                             = umo_cstrs(char(get(coUD,  'Enable')),'on', 'im1');
if any(im1>0);                  set(coUD(im1),'Enable',     'off');                                 end;
drawnow;

fNo                             = cwUD(1);
sNo                             = cwUD(2);
% performing 'a' and 's' processes in the iBase:
i0                              = double(g4iv2.ppp(:,1)==97) + double(g4iv2.ppp(:,1)==115).*2;
ii                              = find(i0 & g4iv2.ppp(:,4)==cwUD(3));
% 
q1                              = zeros(1,  size(g4iv2.irq,2));
q2                              = zeros(size(q1));
q3                              = [1:1:size(q1,2)-1,1];
% looping over a/s process:
for i=ii(:)';
    mv2_fck(fNo,sNo,            i);
    q1(:)                       = g4iv2.ick{sNo}(i,:)==g4iv2.irq(i,:);
    q2(:)                       = g4iv2.ock{sNo}(i,:)<g4iv2.orq(i,:) & g4iv2.orq(i,:)>0;
    %
    disp(int2str([i,q1(1,:),q2(1,:)]));
    %
    for j=find(q1==1 & q2>0);    
        feval(g4iv2.ifl{g4iv2.ppp(i,3)}, g4iv2.ppp(i,2), [fNo,sNo,q3(j),g4iv2.ppp(i,:)]);         	end;
        % updating g4iv2{fNo}.fck for Subject sNo:
    mv2_fck(fNo,sNo,        i);                                                                     end;
% updating subject x iBase GUIs strings:
mv2_a0(fos(1));                                                                                
%
% updating L2W tasks GUIs:
mv2_fck(fNo,sNo,                find(g4iv2.ppp(:,4)==cwUD(3)));
mv2_a0(fos(1));
if any(im1>0);                  set(coUD(im1),'Enable',     'on');                                  end;

return;
%%

function                        local_exit(fos);
%%
set(fos(1),'Visible','off');

% cwUD = [L1W fing.#, subject#, iBase#]
cwUD                            = get(fos(1),               'userData');
figure(cwUD(1));

global g4iv2;
% updating input/output file status:
mv2_fck(cwUD(1),cwUD(2),find(g4iv2.ppp(:,4)==cwUD(3)));
% updating status flag:
mv2_a0([cwUD(1),cwUD(2),0,cwUD(3)]);

return;
%%

function                        local_next(fos)
%%
h                               = findobj(fos(1),           'Tag','iv2_L2W_subject');
if isempty(h);                                                                      return;         end;
p0                              = get(fos(1),               'Position');
% cwUD = [L1W fing.#, subject#, iBase#]
cwUD                            = get(fos(1),               'userData');
f0UD                            = get(cwUD(1),              'userData');
global g4iv2;
cwUD(1, 2)                      = cwUD(1, 2) + 1;
if size(g4iv2.yyy.snm,1)<cwUD(1, 2);
    postQ({'No more subjects to work on',' '},              []); 
    p1                          = get(gcf,                  'Position');
    set(gcf,    'Position',     [p0(1),round(p0(2)+p0(4)./2),p1(3:4)]);             return;         end;
% more subjects are present (but not on display):
if ~f0UD{3}(cwUD(1,2),2);
    h2                          = findobj(cwUD(1),'Tag',    'L1W_pagedown');
    if ~isempty(h2);            figure(cwUD(1));
                                set(cwUD(1),'CurrentObject',h2(1));
                                mv2_a1(-1);
                                f0UD                        = get(cwUD(1),      'userData');
    else;                       disp('.unexpected problem @local_next@mv2_a2.m');   return;         end;
end;
% Attemting to move L2W to the next subject:              
mv2_p1(cwUD);
%
if strcmpi(g4iv2.yyy.ipk(1,1:4),'prep');
    if strcmpi(get(findobj(gcf, 'Tag','L2W_gUseR0'),'String'),  ...
            'Review / approve outputs of FS for multi-MRIs');
        mv2_run_FS('multi_check',[],[]);                                            return;         end;
    % ezr/*mmm.mat
%     ifl                         = mv2_genfln(fullfile('ezr',[mv2_pmp2code(  ...
%                                 g4iv2.xxx(1).pmp),'.mat']),        [cwUD(1), cwUD(2), 1]); 
%     feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'ig',  ifl,[cwUD(1), cwUD(2), 1],[]);
%     feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'disp',ifl,[cwUD(1), cwUD(2), 1],[]);
    % now adjusted for 
%     feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'ig',  [cwUD(1), cwUD(2), 1],ifl,[]);
%     feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'disp',[cwUD(1), cwUD(2), 1],ifl,[]);
% when iPack is TAC2MPE???:
elseif strcmpi(g4iv2.yyy.ipk(1,1:4),'tac2');
    h                           = findobj('Tag',            'Display MPE results/plots');
    if ~isempty(h);
    % updating subject #:
        ud2                     = get(h(1),                 'userData');
        set(findobj(h(1),       'Tag','dispRes_subjectField'), ...
                                'String',                 deblank(g4iv2.yyy.snm(cwUD(1,2),:)));
        ud2.fbc(1, 2)           = cwUD(1, 2);
        set(h(1),               'userData',                 ud2);                                   end;
end;
%
% updating results/maps display windows, if any:
mv2_p1(cwUD,'update');
return;
%%

function                        local_set4jump2(fos);
%% pewpare for jumping to a desired subject
%
mv2_w4L2Wguis('resetall', gcf);
global g4iv2;
[snm, is]                       = sortrows(g4iv2.yyy.snm);
set(findobj(gcf, 'Tag','L2W_gUseR0'),   'String','Select a subject to jump to:');
set(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'Value',1,  'Style','popupmenu',    'UserData',is,  ...
                                'String',snm,       'CallBack','mv2_a2(''jump2'',[]);');
return;
%%

function                        local_jump2(fos);
%% jump to a selected subject
if ~strcmpi(get(gco,'Tag'),'L2W_gUseR1C1');                                         return;         end;
global g4iv2;
ud                              = get(gcf,      'UserData');
ud2                             = get(gco,      'UserData');
ud(1,   2)                      = ud2(get(gco,  'Value'));
set(gcf,    'UserData', ud);
set(findobj(gcf, 'Tag','iv2_L2W_subject'),      'String',deblank(g4iv2.yyy.snm(ud(2), :)));
set(gcf,    'CurrentObject',findobj(gcf, 'String','Update'));
mv2_a2([]);
return;
%%
