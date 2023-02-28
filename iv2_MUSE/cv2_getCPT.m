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
   	'String',['Select input files (green=done) for PET #',int2str(iii{end})], 'FontWeight','bold');
set(findobj(gcf,'Tag','L2W_gUseR1C1'), 'String','Plasma',   ...
                                'Callback','cv2_getCPT(''get'',1);',    'UserData',iii);
%
set(findobj(gcf,'Tag','L2W_gUseR1C2'), 'String','HPLC',     .../.,
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
job                             = { 'cpt',  'hplc' };
udL2W                           = get(gcf,      'UserData');        
udR1C1                          = get(gco,     	'UserData');
fbc                             = [udL2W(1, 1:2), udR1C1{end}];
%
[f1, g1]                        = mv2_genfln(udR1C1{1}, fbc);
pet_is                          = feval(g4iv2.yyy.lds, 'pet_d2i',   f1);
[f2, g2]                      	= mv2_genfln(udR1C1{iii(1)+1}, fbc);
if g2>0;                        set(gco,    'BackgroundColor',iv2_bgcs(12));        return;         end;
my_getCPT(job{iii(1)},  pet_is, f2);
% if exist(pet_is,'dir')~=7;                                                          return;         end;
% if g2>0;                        set(gco,    'BackgroundColor',iv2_bgcs(12));        return;         end;
% [f3, g3]                        = mv2_genfln(fullfile('pet', [job{iii(1)},'_source.m']), fbc);
% if g3>0;                        f3b                         = umo_getptf(f3,0,[]);
% else;                           f3b                         = [];                                   end;
% %
% if isempty(f3b) || ~exist(f3b,'file');
%     [fname, path_x]           	= uigetfile(fullfile(pet_is,'*'));
%     if ~ischar(fname);                                                              return;         end;
%     [fdx, fnm, fex]          	= fileparts(fname);
%     [f3x, f3m]                  = fileparts(f3);
%     f3b                         = fullfile(f3x, [f3m,fex]);
%     disp(['> copying ',job{iii(1)},' source file..']);
%     disp(['  input: ',path_x,fname]);
%     disp([' output: ',f3b]);
%     copyfile(fullfile(path_x, fname),   f3b);
%     write2ptf(f3, f3b);                                                                             end;
% %
% my_getCPT(job{iii(1)},  {f3b}, {f2});
%
if exist(f2,'file');            set(gco,    'BackgroundColor',iv2_bgcs(12));                        end;
return;
%%

function                        local_sover(iii);
%%
global g4iv2;
udL2W                           = get(gcf,      'UserData');        
udR1C1                          = get(findobj(gcf, 'Tag','L2W_gUseR1C1'), 'UserData');
[f1, g1]                        = mv2_genfln(udR1C1{iii(1)+1}, [udL2W(1, 1:2), udR1C1{end}]);
if g1>0;                        delete(f1);                                                         end;
[f2, g2]                        = mv2_genfln(fullfile('pet', [g4iv2.xxx(udR1C1{end}).cpt,'.cpt']),  ...
                                                            [udL2W(1, 1:2), udR1C1{end}]);
if g2>0;                        delete(f2);                                                         end;
[f3, g3]                        = mv2_genfln(fullfile('pet', ...
                                [g4iv2.xxx(udR1C1{end}).cpt,'_ok.txt']),  ...
                                                            [udL2W(1, 1:2), udR1C1{end}]);
if g3>0;                        delete(f3);                                                         end;

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
    [f1, g1]                    = mv2_genfln(udR1C1{i+1}, [udL2W(1, 1:2), udR1C1{end}]);
    if g1>0;                    set(findobj(gcf, 'Tag',['L2W_gUseR1C',int2str(i)]),     ...
                                    'BackgroundColor',iv2_bgcs(12));                        end;    end;
return;