function    d                   = iv2_load_outputs(i1,i2,i3,ref); 

% To load and extract data (subject x vois) for one variable
%
%  To save data in mat file 
%   set up an IDAE session for the project, then
%   >> iv2_load_outputs('full\path\outputs_database.mat')
%
%     need to name the output database file as follows:
%       file = fullfile('whatever',petx,analysisID.mat');
%         where x of petx stands for scan x (PET x) of the project
%         analysis ID is
%         to display existing on (for a TAC2MPE package):
%           Hit 'Existing RTMs/PIMs outputs' GUIs of the display-result module
%
%       usage:      d           = iv2_load_outputs('full\path\outputs_database.mat',var,vnos[,ref])
%       
%   var -   all variables in .ezd files may be used. 
%           addtional valid variables:
%               Ki      K1*k3/(k2+k3)
%               VT      K1/k2 if OTCM (k3 & k4 = 0)
%                       K1/k2/(1 + k3/k4), if TTCM
%               BP      VT/VR - 1, if OTCM or TTCM
%               BP2     k3/k4, if TTCM
%               DVR     returns BPND + 1
%               SUVR    returns 'ratio' if present
%               AIC     Akaike information criteria, if RRS are recorded
% 
% Notes:
%   1. enter [] for 2nd and 3rd inputs for checking presence of output.mat,
%       and create one if not present. called from iv2_save_outputs.m
%   2. return additional variables (by calculation)
%       vnm_added   = iv2_load_outputs([],vnm,mfg{2},ref)
%       	vnm     names of variables recored in a .ezd file (n by m)
%                   vnm = getLseg(gei(ifl, 'orientation'), [0,1]);
%           mfg     mfg{1} = 'PIMs' or 'RTMs'
%                   mfg{2} = s4mpe.minfo{i}{2} (that correspond to 'ifl')
%                   ignored if not OTCM or TTCM*
%           ref     VOIIDNos of reference regions
%       
%
% (cL)2017    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin==4 && isempty(i1);    d                           = local_vnm(i2,i3,ref); return;         end;
if nargin<4;                    ref                         = 0;                                    end;
if nargin==1;                   
    if exist(i1,'file');        d                           = local_check(i1);
                                if d.uptodate>0;                                    return;         end;
    else;                       d.uptodate                  = 0;                                    end;
    if d.uptodate>0;                                                                return;         end;
    local_save(i1);
    if nargin>0;                d                           = load(i1);                             end;
                                                                                    return;         end;
%
if nargin<margin;               helq(mfilename);                                    return;         end;
d                               = [];
if isempty(i2);                                                                     return;         end;
global g4iv2;
x                               = [];
if exist(i1,'file')==2;         
    x                           = local_check(i1);
    if isempty(x);              disp(['.wrong input: ',i1]);                        return;         end;
                                                                                                    end;
if isempty(x);                 	local_save(i1);
                                x                           = load(i1);                             end;
if isempty(x);                  disp('.error! ??? ');                               return;         end;
if isempty(x.mflg);             x.mflg                      = '     ';                              end;
% checking reference region(s):
if ref>0;                       vr                          = consolidVOINos(x.vnos,ref);
    if any(~vr(:,2));           disp('.problem! not all requested VOIs are present (=0)');
                                vvr                         = VOIdef(ref(:));
                                dispCharArrays(1,int2str(vr(:,2)),2,vvr.anm);       return;         end;
else;                           vr                          = [];                                   end;
%
im1                             = umo_cstrs(lower(x.vnm),lower(i2), 'im1');
if isempty(i3);                 i3                          = x.vnos;                               end;
vi                              = consolidVOINos(x.vnos,i3(:));
d                               = nan(size(x.res{1},1),     size(vi,1));
if im1>0;                       
    d(:, vi(:,2)>0)             = x.res{im1}(:,vi(vi(:,2)>0,2));                    return;         end;
%
v2                              = lower(char('Ki','VT','BP','BP2','DVR','SUVR','AIC'));
im2                             = umo_cstrs(v2,[lower(i2),' '], 'im1');
vtc                             = umo_cstrs(lower(x.vnm),'vt ',  'im1');
bpc                             = umo_cstrs(lower(x.vnm),'bp ',  'im1');
%
[idx, inm]                      = fileparts(i1);
if im2==1;
% Ki
    if sum(x.mflg(1:3)=='i')==3;
        x.res{1}                = x.res{1}.*x.res{2}/(x.res{2}+x.res{3});
        d(:, vi(:,2)>0)         = x.res{1}(:,vi(vi(:,2)>0,2));
    else;                       d                           = local_ng(inm,i2);                     end;
    return;
elseif im2==2;
% VT: 
    if ~any(x.mflg(1:4)=='i');
        x.res{1}                = x.res{1}./x.res{2}.*(1+x.res{3}./x.res{4});
        d(:, vi(:,2)>0)         = x.res{1}(:,vi(vi(:,2)>0,2));
    elseif ~any(x.mflg(1:2)=='i');
        x.res{1}                = x.res{1}./x.res{2};
        d(:, vi(:,2)>0)         = x.res{1}(:,   vi(vi(:,2)>0,2));
    else;                       d                           = local_ng(inm,i2);                     end;
    return;
elseif im2==3;
% BP:
    if isempty(vr);             d                           = local_noref(inm,i2);  return;         end;
    if vtc>0;
        vtr                     = nanmean(x.res{vtc}(:,vr(:,2)),2);
        d(:, vi(:,2)>0)         = x.res{vtc}(:,vi(vi(:,2)>0,2))./vtr(:, ones(1,sum(vi(:,2)>0))) - 1;
    % trying OTCM or TTCM:
    elseif ~any(x.mflg(1:4)=='i');
        x.res{1}                = x.res{1}./x.res{2}.*(1+x.res{3}./x.res{4});
        vtr                     = nanmean(x.res{1}(:,vr(:,2)),2);
        d(:, vi(:,2)>0)         = x.res{1}(:,vi(vi(:,2)>0,2))./vtr(:, ones(1,sum(vi(:,2)>0))) - 1;
    elseif ~any(x.mflg(1:2)=='i');
        x.res{1}                = x.res{1}./x.res{2};
        vtr                     = nanmean(x.res{1}(:,vr(:,2)),2);
        d(:, vi(:,2)>0)         = x.res{1}(:,vi(vi(:,2)>0,2))./vtr(:, ones(1,sum(vi(:,2)>0))) - 1;
    else;                       d                           = local_ng(inm,i2);                     end;
    return;      
elseif im2==4;
% BP2:
    if ~any(x.mflg(1:4)=='i');
        x.res{1}                = x.res{3}./x.res{4};
        d(:, vi(:,2)>0)         = x.res{1}(:,   vi(vi(:,2)>0,2));
    else;                       d                           = local_ng(inm,i2);                     end;
    return;
elseif im2==5;
% DVR:
    if bpc>0; 
        d(:, vi(:,2)>0)         = x.res{bpc}(:, vi(vi(:,2)>0,2)) + 1;               return;         end;
    if isempty(vr);             d                           = local_noref(inm,i2);  return;         end;
    if vtc>0;
        vtr                     = nanmean(x.res{vtc}(:,vr(:,2)),2);
        d(:, vi(:,2)>0)         = x.res{vtc}(:,vi(vi(:,2)>0,2))./vtr(:, ones(1,sum(vi(:,2)>0)));
    % trying OTCM or TTCM:
    elseif ~any(x.mflg(1:4)=='i');
        x.res{1}                = x.res{1}./x.res{2}.*(1+x.res{3}./x.res{4});
        vtr                     = nanmean(x.res{1}(:,vr(:,2)),2);
        d(:, vi(:,2)>0)         = x.res{1}(:,vi(vi(:,2)>0,2))./vtr(:, ones(1,sum(vi(:,2)>0)));
    elseif ~any(x.mflg(1:2)=='i');
        x.res{1}                = x.res{1}./x.res{2};
        vtr                     = nanmean(x.res{1}(:,vr(:,2)),2);
        d(:, vi(:,2)>0)         = x.res{1}(:,vi(vi(:,2)>0,2))./vtr(:, ones(1,sum(vi(:,2)>0)));
    else;                       d                           = local_ng(inm,i2);                     end;
    return;
elseif im2==6;
% SUVR:     
    ims                         = umo_cstrs(lower(x.vnm),'ratio',   'im1');
    if ~ims;                    d                           = local_ng(inm,i2); 
    else;                       d(:, vi(:,2)>0)             = x.res{ims}(:, vi(vi(:,2)>0,2));       end;
elseif im2==7;
% AIC:     
    imx                         = umo_cstrs(lower(x.vnm),'rss ',    'im1');
    if ~imx;                    d                           = local_ng(inm,i2); 
    else;                       d(:, vi(:,2)>0)             = x.res{imx}(:, vi(vi(:,2)>0,2));       end;
end;
return;
%%

function    out                 = local_check(ifl);
%%
[idx, inm]                      = fileparts(ifl);
[jdx, jnm]                      = fileparts(idx);
out                             = [];
if ~strncmpi(jnm,'pet',3);                                                          return;         end;
pno                             = str2num(jnm(1,4:end));
if isempty(pno);                disp('.problem! petxxx (x=integers) expected');
                                disp([' input: ',jnm]);                             return;         end;
global g4iv2;
x                               = load(ifl);
x.uptodate                      = 0;
%
im1                             = umo_cstrs(x.sid,  g4iv2.yyy.did,  'im1');
%
[fls, fi]                       = makefarrays('res',[inm,'.ezd'],'fbc',[1,0,pno]);
im2                             = umo_cstrs(char(x.ezd), fls,   'im1');

if sum(im1(:,1)~=[1:1:size(im1,1)]')>0 || sum(im2(:,1)~=[1:1:size(im2,1)]')>0 ||    ...
                                                            max(mv2_get_dnum(fls))>mv2_get_dnum(ifl);
  	disp(['.updating: ',ifl]);
   	local_save(ifl); 
    %
    if max(mv2_get_dnum(fls))>mv2_get_dnum(ifl);
        disp('> update failed (aborting)');                                         return;         end;
    %
    x                           = load(ifl);                                                
    x.uptodate              	= 1;                                                                end;
        
%
disp(['.loading: ',ifl]);
out                             = x;
return;
%%

function                        local_save(ifl);
%%
global g4iv2;
[idx, inm]                      = fileparts(ifl);
[jdx, jnm]                      = fileparts(idx);
[f1, g1]                       = makefarrays('res',[inm,'.ezd'],'fbc',[1,0,str2num(jnm(1,4:end))]);
snm                             = g4iv2.yyy.snm;
sid                             = g4iv2.yyy.did;
if ~any(g1>0);                  disp('.problem! not ready for any subjects (no outputs)');
                                disp([' sought: ',ifl]);                            return;         end;
for i=1:1:size(f1,1);
    if g1(i)>0;                 ezd{i}                      = deblank(f1(i, :));
    else;                       ezd{i}                      = fileparts(f1(i, :));         end;    end;
% 
[x0, mflg]                      = gei(ezd{find(g1>0,1)},    'orientation','modelused');
vnm                             = local_vars(x0);
if ~isempty(mflg);              mflg                        = local_mflg(mflg);                     end;
ic                              = 0;
for j=find(g1'>0);
    ic                          = ic + 1;
    d                           = ged(ezd{j},               1);
    x                           = gei(ezd{j},               'orientation');
    if ic==1;                   x0                          = x;
                                vnm                         = local_vars(x);
                                vc                          = umo_cstrs(lower(vnm),'voiidno',  'im1');
        if ~vc;                 disp('.error! column of VOIIDNo not found');        return;         end;
        for k=1:1:size(vnm,1);  res{k}                      = nan(size(f1,1), size(d,1));          end;
                              	vnos                        = d(:,  vc);
                                vi                          = zeros(size(vnos,1), 2);               end;
    %
    if ~strcmpi(x,x0);          disp('.error! inconsistent data orientation');
                                disp([' ',x0, ' for 1st subject']);
                                disp([' ',x,' for subject #',int2str(j)]);          return;         end;
    vi(:)                       = consolidVOINos(d(:, vc),  vnos);
    for k=1:1:size(vnm,1);      res{k}(j, vi(:,2)>0)        = d(vi(vi(:,2)>0,2), k)';               end;
                                                                                                    end;
%
if ~exist(fileparts(ifl),'dir');                            mkdir(fileparts(ifl));                  end;
save(ifl,                       'vnm', 'res', 'vnos', 'ezd', 'snm', 'sid', 'mflg');
disp('.saved! (MPE output database)');
disp([' file: ',ifl]);
return;
%%

function    out                 = local_vars(x);
%%
x(x=='[' | x==',' | x==']')     = ' ';
out                             = char(line2mat(x));
return;
%%

function    out                 = local_noref(inm,i2);
%%
out                             = [];
disp('.error! enter VOIID# of reference region [=4th input]');
disp([' requested variable: ',i2]);
disp([' output file: ',inm]);
return;
%%

function    out                 = local_ng(inm,i2);
%%
out                             = [];
disp(['.error! ',i2,' may not be available for this output type']);
disp([' output file: ',inm]);
return;
%%

function    out                 = local_uc(i2);
%%
out                             = [];
disp(['.error! under construction']);
disp([' requested variable: ',i2]);
return;
%%

function    out                 = local_mflg(i2);
%% return modeling flag (5 character string) for the input (such as k4v0-model)
% check get_mstr.m perioedically for dstrs and mstrs
dstrs                           = str2mat( ...
                                'k2-model (v0-ignored)','k2v0-model (v0-fixed)','k2v0-model',   ...
                                'k3-model (v0-ignored)','k3v0-model (v0-fixed)','k3v0-model',   ...
                                'k4-model [v0-ignored; Ve(K1/k2)',          ...
                                'k4v0-model [v0-fixed; Ve(K1/k2)',          ...
                                'k4v0-model [Ve(K1/k2) ',                    ...
                                'k4-model (v0-ignored)','k4v0-model (v0-fixed)','k4v0-model');
im1                             = umo_cstrs(dstrs(:,1:size(i2,2)+2),[i2,'  '],   'im1');
if ~im1;                        disp('.problem! a nonregistered modeling flag');
                                disp([' flag: ',i2]);
                                disp(' request hkuwaba1@jhmi.edu to register this flag');
                                out                         = [];               return;             end;
if length(im1)>1;               disp('.problem??? more than one corresponding flags');
                                dispCharArrays(1,dstrs(im1, :));
                                disp(' request hkuwaba1@jhmi.edu to fix this problem');
                                out                         = [];               return;             end;
%
mstrs                           = [ 'eeiiii';   'eeiici';   'eeiiei';
                                    'eeeiii';   'eeeici';   'eeeiei';
                                    'ereeii';   'ereeci';   'ereeei';
                                    'eeeeii';   'eeeeci';   'eeeeei'];
out                             = mstrs(im1,    1:5);
return;
%%

function    out                 = local_vnm(vnm,mfg,ref);
%%
rstr                            = [];
if ~isempty(ref);               vv                          = VOIdef(ref(:));
    for i=1:1:size(vv.anm,1);   rstr{i}                     = deblank(vv.snm(i, :));    	end;    end;
ic                              = 0;
if strcmpi(mfg{1},'pims');
    if mfg{2}(1)=='e';          ic                          = ic + 1;
                                add{ic}                     = 'VT';                                 end;
    if size(mfg{2},2)>=4 && mfg{2}(4)=='e';          
                                ic                          = ic + 1;
                                add{ic}                     = 'BP2 (=k3/k4)';                       end;
    for i=1:1:numel(rstr);      ic                          = ic + 1;
                                add{ic}                     = ['BP (ref: ',rstr{i},')'];
                                ic                          = ic + 1;
                                add{ic}                     = ['DVR (ref: ',rstr{i},')'];           end;
% for RTMs:
else;
    sss                         = ['bp  ';'rati'];
    im1                         = umo_cstrs(lower(vnm),sss,	'im1');
    for i=1:1:size(sss,1);
        if im1(i)>0;            eval([sss(i, :),'           = 1;']);
        else;                   eval([sss(i, :),'           = 0;']);                        end;    end;
    if bp>0;                    ic                          = ic + 1;
                                add{ic}                     = 'DVR (=BP+1)';                        end;
    if rati>0;                  ic                          = ic + 1;
                                add{ic}                     = 'SUVR (=ratio)';                      end;
                                                                                                    end;
%
rss                             = umo_cstrs(lower(vnm),'rss ',  'im1');
if rss>0;                       ic                          = ic + 1;
                                add{ic}                     = 'Akaike IC';                          end;
if ic>0;                        im1                         = umo_cstrs(char(vnm),char(add), 'im1');    
                                out                         = char(vnm,'<calculated>',char(add(~im1)));
else;                           out                         = vnm;                                  end;
return;
%%