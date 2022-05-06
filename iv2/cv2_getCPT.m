function    cv2_getCPT(i1, i2); 
% To aid to idenfy/convert local plasma/HPLC files to IDAE-ready files
%   When called, this code sets the lower GUIs of L2W to work for plasma
%   and HPLC data
%
%   Callback jobs of the 'plasma' and 'HPLC' GUIs are:       
%     Plasma:     dxetc4xxx('get_plasma', input_2)
%     HPLC:       dxetc4xxx('get_hplc',   input_2)
%     	where input_2.ofl = IDAE-supplied plasma and HPLC files (outputs)
%      	and   input_2.fbc = [figure # of L1W, subject #, pet #]
% 
%   Set your dxetc4xxx.m for the plasma data as follows:
%     Prepare a two-column matrix (any variable name, for example 'cpt')
%       cpt = [sampling times, radioactivity values];
%     Then, save it in 8-digit ASCII format:
%       save(input_2.ofi, 'cpt', '-ascii');
%     Add the following line to update the 'plasma' GUI
%       cv2_getCPT('check',1);
%
%   Set your dxetc4xxx.m for the HPLC data as follows:
%     Prepare a two-column matrix (any variable name, for example 'hplc')
%       hplc = [sampling times, fractions of metabolites in %];
%         2nd column = 100 - parent fraction in %
%     Then, save it in 8-digit ASCII format:
%       save(input_2.ofi, 'hplc', '-ascii');
%     Add the following line to update the 'HPLC' GUI
%       cv2_getCPT('check',2);
%
% (cL)2022    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

feval(['local_',lower(i1)], i2);
return;
%%


function                        local_set(iii);
%% set GIUs of L2W
%
% figure(findobj(groot, 'Tag','iv2L2W'));
mv2_w4L2Wguis('resetall',gcf);
set(findobj(gcf, 'Tag','L2W_gUseR0'),      ...
   	'String',['Select input files (green=done) for PET #',int2str(iii{3})], 'FontWeight','bold');
set(findobj(gcf,'Tag','L2W_gUseR1C1'), 'String','Plasma',   ...
                                'Callback','cv2_getCPT(''get'',1);',    'UserData',iii);
%
set(findobj(gcf,'Tag','L2W_gUseR1C2'), 'String','HPLC',     ...
                                'Callback','cv2_getCPT(''get'',2);',    'UserData',iii);
%
set(findobj(gcf,'Tag','L2W_gUseR1C3'), 'String','Cancel',   ...
                                'Callback','mv2_w4L2Wguis(''resetall'',[]);');
%
set(findobj(gcf,'Tag','L2W_gUseR2C1'), 'String','Restart Plasma',       ...
                                                            'Callback','cv2_getCPT(''sover'',1);');
set(findobj(gcf,'Tag','L2W_gUseR2C2'), 'String','Restart HPLC',       ...
                                                            'Callback','cv2_getCPT(''sover'',2);');

local_check([1,2]);
return;
%%

function                        local_get(iii);
%%
global g4iv2;
job                             = {'get_plasma',    'get_hplc'};
udL2W                           = get(gcf,      'UserData');        
udR1C1                          = get(gco,     	'UserData');
[f1, g1]                        = mv2_genfln(udR1C1{iii(1)}, [udL2W(1, 1:2), udR1C1{3}]);
if g1>0;                        set(gco,    'BackgroundColor',iv2_bgcs(12));        return;         end;
feval(g4iv2.yyy.lds,    job{iii(1)}, struct('ofl',f1,   'fbc',[udL2W(1, 1:2), udR1C1{3}]));
%
return;
%%

function                        local_sover(iii);
%%

udL2W                           = get(gcf,      'UserData');        
udR1C1                          = get(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'UserData');
[f1, g1]                        = mv2_genfln(udR1C1{iii(1)}, [udL2W(1, 1:2), udR1C1{3}]);
if g1>0;                        delete(f1);                                                         end;
set(findobj(gcf, 'Tag',['L2W_gUseR1C',int2str(iii(1))]),    'BackgroundColor',iv2_bgcs(0));
return;
%%

function                        local_check(iii);
%%
figure(findobj(groot, 'Tag','iv2L2W'))
if ~strcmpi(get(findobj(gcf, 'Tag','L2W_gUseR1C1'),  'Callback'),'cv2_getCPT(''get'',1);');
                                                                                    return;         end;
udL2W                           = get(gcf,      'UserData');        
udR1C1                          = get(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'UserData');
for i=iii(:)';
    set(findobj(gcf, 'Tag',['L2W_gUseR1C',int2str(i)]), 'BackgroundColor',iv2_bgcs(0));
    [f1, g1]                    = mv2_genfln(udR1C1{i}, [udL2W(1, 1:2), udR1C1{3}]);
    if g1>0;                    set(findobj(gcf, 'Tag',['L2W_gUseR1C',int2str(i)]),     ...
                                    'BackgroundColor',iv2_bgcs(12));                        end;    end;
return;