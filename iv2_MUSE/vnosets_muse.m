function    [out, out2, out3, out4] = vnosets_muse(i1, varargin); 

% To generate VOIIDNo set vectors, MUSE added by jmei in 2023
%       
%       usage:      vnos            = vnosets('set_flag')
%                   [vnos, vflg]    = vnosets('set_flag')
%
%   >> vnoset   ... will display defined VOIs sets.
%                   Add 'u' to 'set_flag' (e.g., 'p10u' for 'p10') to unite
%                   left and right VOIs.
%   vflg        -   voi string flags in a few ~ several characters.
%
% Options:
%   'dsp','on'  -   to display anatomical names of VOIs
%   'ref',val   -   to remove selected regions (=val by ID.No) from vnos
%                   Use this option to remove the reference region from vnos
%   'arc','on'  -   to archive existing VOI sets into one file
%                   (with left and right VOIs merged)
%                   output: vnosets_archived.mat in which('vnosets');
% New!
%   my_vnosets      Individuals can make their own vnosets
%                   Request hkuwaba1@jhmi.edu to send a template file, if you are interested
%   
% (cL)2004~10   hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);
                                disp(local_vsets);                                  return;         end;

%% Options:
dspval                          = 'off';
refval                          = [];
arcval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
if strcmpi(arcval,'on');        local_archive;                                      return;         end;
if isnumeric(i1);               out                         = i1;                   
    if nargout>1;               out2                        = local_out2(out);                      end;
                                                                                    return;         end;
if size(i1,2)<3;                disp('.error! set_flag = >=3 characters');          return;         end;
if i1(1)=='$';                  out                         = mv2_iv2VOIs([]);      return;         end;
% special VOIs that were defined after 5/4/2010
if ischar(i1) && size(i1,2)>=3;
    i1                          = lower(i1);
    myfile                      = fullfile(fileparts(which('scratch')),'my_vnosets.m');
    if exist(myfile,'file');
        out                     = my_vnosets(i1(1,  1:3));
        % disp('.using my_vnosets');
        if ~isempty(out);
        if size(i1,2)>3;        out                         = local_unite(out);                     end;
        if nargout>1;           out2                        = local_out2(out);                      end;
                                                                                    return;         end;
                                                                                                    end;

    newly_added                 = ['fsl';'f81';'mus';'f45';'s28';'rat';'p18';'p20';'m89';'m37'];
    im0                         = umo_cstrs(newly_added,i1(1,1:3),'im1');
    if im0(1)>0;                out                         = feval(['local_',newly_added(im0,:)]);
        % disp('.using newly_added');
        if ~isempty(out);
        if size(i1,2)>3;        out                         = local_unite(out);                     end;
        if nargout>1;           out2                        = local_out2(out);                      end;
                                                                                    return;         end;
                                                                                            end;    end;
% disp('.using vnosets.m');
out                             = [];
out2                            = [];
out3                            = [];
out4                            = [];

% an existing file is entered:
if ~isnumeric(i1) &&   ~isempty(dir(i1));
    out                         = local_file(i1);       
    if refflg;                  out                         = out(out~=refval(1));                  end;
    if nargout>1;               out2                        = local_out2(out);                      end;
    if nargout>2;               out3                        = local_out3(out, out2);                end;
                                                                                    return;         end;



%% when input is a numeric vector:
if isnumeric(i1);               out                         = i1;
    if nargout>1;               out2                        = local_out2(i1(:));                    end;
    if nargout>2;               out3                        = local_out3(out, out2);                end;
                                                                                    return;         end;

fff                             = which(mfilename);
qqq                             = load([fff,'at']);

uflg                            = 0;
if lower(i1(end))=='u';         uflg                        = 1;
                                i1                          = i1(1,     1:end-1);                   end;

im1                             = umo_cstrs(lower(qqq.vstr),lower(i1),'im1');
if ~im1(1);                     disp(['Unregistered VOI set label ... ',i1]);       return;         end;

out                             = find(qqq.vmat(:,      im1(1)));
[v, is]                         = sort(qqq.vmat(out,    im1(1)));
out(:)                          = qqq.vnos(out(is),     :);

if uflg;                        out                         = local_unite(out);                     end;

if strcmpi(i1,'b17') && uflg;   out(out==59000)             = 66000;                                end;

if refflg;                      out                         = out(out~=refval(1));                  end;
if nargout>1;                   out2                        = local_out2(out);                      end;
if nargout>2;                   out3                        = local_out3(out, out2);                end;

return;


function    out2                = local_out2(vnos);
%% calculating out2, if nargout>1:
vvv                             = VOIdef(vnos);
out2                            = [];
for i=1:1:size(vnos,1);         out2{i}                     = deblank(vvv.snm(i,    :));            end;
return;
%%

function    out3                = local_out3(out, out2);
%% converting out2 to out3:

out30                           = char(zeros(size(out,1),   10) + 32);
for i=1:1:size(out,1);          istr                        = char(out2(i));
                                out30(i, 1:size(istr,2))    = istr;                                 end;
    
out31                           = out30';
out3                            = out31(:)';

return;
%%


function    vsets               = local_vsets;
%% defined vnosets:

disp(['The list is as of 5/4/2010 (no further updates)',10,  ...
    'Because of the increased number of registered VOI sets, labels cannot be tracked.',10,    ...
    'Thus, <<set_vnosets>> checks each user-defined set against current list of sets ',10,    ...
    'and create one (in a xnn format) if the set does not exist yet.']);

vsets                           = [ 's01        striatum                                        '; ...
                                    's02        striatum (left and right)                       '; ...
                                    's03        putamen, caudate, and VS (L+R)                  '; ...
                                    's04        caudate nucleus & putamen per side      (n=4)   '; ...
                                    'c04        s04 + cerebellum                                '; ...
                                    's4t        s04 + thalamus                                  '; ...
                                    's4g        s04 + globus pallidus                           '; ...
                                    'stc        s4t + cerebellum                                '; ...
                                    'sto        s4t + occipital lobe                            '; ...
                                    'str        s4t + cerebellum + occipital lobe               '; ...
                                    'm10        s4t + thalamus + hippocampus + amygdala         '; ...
                                    'obm        s04 + thalamus + hippocampus (opiateBM)         '; ...
                                    'bm0        hippocumpus + amygdala + DLPFC + VS (opiateBM)  '; ...
                                    'cfn        striatum + thalamus + cerebllum (for CFN)       '; ...
                                    'o16        Am, vS, Pu, GP, In, Th, CN, and Cg (lt & Rt)    '; ...
                                    'tpm        ant/port thalamus, pons, medulla                '; ...
                                    's06        caudate nucleus, putamen, & VS per side         '; ...
                                    's10        VS, and/post CDN & PUT, left then right (n=10)  '; ...
                                    'p10        VS, and/post CDN & PUT, side by side    (n=10)  '; ...
                                    'p12        p10 + SM (n=12)                                 '; ...
                                    'p08        VS, and/post PUT & ant CDN, side by side (n-8)  '; ...
                                    'cln        VOIs defined on stdVOIsets(''cln'');              '; ... 
                                    'icb        VOIs defined on stdVOIsets(''icb'');              '; ...
                                    'm56        icb + m10 - striatum                            '; ...
                                    'm57        m56 + dorsal tegmentum nucleus (left + right)   '; ...
                                    'm58        m56 + dorsal tegmentum nucleus (left & right)   '; ...
                                    'stt        struatum + thalamus (left & right)              '; ...
                                    'thh        thalamus + hippocampus                          '; ...
                                    'sth        p10 + thalamus + hippocampus                    '; ...
                                    'sah        p10 + hippocampus + amygdala                    '; ...
                                    'tah        thalamus + amygdala + hippocampus               '; ...
                                    'c31        standard 31 VOIs                                '; ...
                                    'f81        Freeserfer-based 81 VOIs (detailed)             '; ...
                                    'f45        truncated version of f81                        '; ...
                                    'fsl        subcortical regions + Cb given by FSL/FIRST     '; ...
                                    'a29        c31 minus pons and hippocampus (for anil)       '; ...
                                    'c30        standard 31 VOIs minus pons                     '; ...
                                    'c29        cortical 29 VOIs (eliminating VS from ''c31'')    '; ...
                                    't31        cortical 31 VOIs (c29 with ant/postTh + medulla)'; ...
                                    'c28        cortical 28 VOIs (c29 minus pons)               '; ...
                                    's28        c28 + vS (14 structures; no Cb)                 '; ...
                                    'p28        c28 + L/R Cb (PET10288)                         '; ...
                                    'a28        c28 + amygdala - parahippocampus                '; ...
                                    'm28        a28 + ventral striatum (MOR)                    '; ...
                                    'w28        c28 + WM (90000) (for CB1)                      '; ...
                                    'q28        c28 + pons (61000) + Cb (71000) (for CB1)       '; ...
                                    'w30        c28 + WM + pons (for CB1)                       '; ...
                                    'm20        c28 minus Cb/Pu/CN/GP                           '; ...
                                    'g28        c28 + corpus callosum (for Georg)               '; ...
                                    'q38        q28 (30 VOIs) + bm0 (opiateBM)                  '; ...
                                    'q39        q28 (30 VOIs) + bm0 (opiateBM) +WM              '; ...
                                    'm30        c28 + raphe (69002) + hypothalamus              '; ... 
                                    'm32        c30 + raphe (69002) + hypothalamus              '; ... 
                                    'c25        cortical 29 VOIs minus striatal 4               '; ...
                                    'b08        s04 + Fr + Tm   (4 VOIs per side)               '; ...
                                    'b20        20 brain regions for baboons/monkeys            '; ...
                                    'b12        discontinued. use b20u insted                   '; ...
                                    'b17        17 brain regions for baboons/monkeys            '; ...
                                    'b18        b17 + left/right hippocampus                    '; ...
                                    'b19        b17 + left/right hippocampus & amygdala         '; ...
                                    'b21        b19 + VS + hypothalamus                         '; ...
                                    'b16        b17 minus pons                                  '; ...
                                    'b15        b17 minus whole cortex and pons                 '; ...
                                    'mg5        b17 + ant/post cingulate                        '; ...
                                    'cng        ant/mid/post cingulate cortices                 '; ...
                                    'rmg        b16 + amygdala + ant/post cyngulate (roch mGlu5)'; ...
                                    'iif        21 brain regions for pig standard VOI template  '; ...
                                    'all        to report [] (may cause errors in some codes)   '; ...
                                    'rba        for evaluation of RBA methods                   '; ...
                                    'spm        for evaluation of spm RBA errors (f06-f20)      '; ...
                                    '5ht        ctx + choroind plexsus                          '; ...
                                    'fpw        b17 + (left + right) parahippocumpus            '; ...
                                    'n10        spheric VOIs numbered from 1 to 10              '; ...
                                    'n18        spheric VOIs numbered from 1 to 18              '; ...
                                    'n20        spheric VOIs numbered from 1 to 20              '; ...
                                    'n40        spheric VOIs numbered from 1 to 40              '; ...
                                    'soc        stratum, Oc, CB, and WM (L+R)                   '; ...
                                    'pib        c28 + n20                                       '; ...
                                    'lly        for Lilly dosimetry                             '; ...
                                    'phm        for spheric phantom scans                       '; ...
                                    'm02        thalamus + cerebellum for Horti studies         '; ...
                                    'ish        59001 + 59002                                   '; ...
                                    'rat        six regions (left+right) for the rat brain      '; ...
                                    'r09        rat + olfactory bulb + cortex + SN              '; ...
                                    'r10        rat + olfactory bulb + cortex + SN + Lt&Rt Hipp '; ...
                                    'gwf        gray/white/CSF                                  '; ...
                                    'e06        MCA/ACA/PCA                                     '; ...
                                    'e14        MCA/ACA/PCA + Put/SN/Thalamus + C.Semiovale     '; ...
                                    'nb9        striatal + cerebellum listed in mic4nico        '; ...
                                    'm19        19 brain regions for NB (umich)                 '; ...
                                    'v14        p10 + thalamus + GP + IC (for PVC-rac)          '; ...
                                    'v11        p10 + thalamus (for PVC-rac)                    '; ...
                                    '200        200 random spheric VOIs (1-200)                 '; ...
                                    'css        COSS 14 VOIs                                    '; ...
                                    'bmx        for OMAR bone marro (mid=31100/mand=31200)      '; ...
                                    'lrh        left & right hemispheres                        '; ...
                                    'ccs        anterior/mid/posterior CC sub-dividions         '; ...
                                    'dsm        a standard set for dosimetry scans              '; ...
                                    'mkd        organs for MK7916 dosimetry scans               '; ...
                                    'ga5        for GABA-benzodiazepine alpha 5 (ver.Baboon)    '; ...
                                    'ga6        for GABA-benzodiazepine alpha 5 (ver.high)      '; ...
                                    'ctx        hemispheres                                     '; ...
                                    'ics        50 brain VOIs as defined in stdVOIset(''ics'')    ';...
                                    'c10        10 cingulate subdividions (5 per hemisphere)    '; ...
                                    'vc1        VOIs for voi defining certification course 01   '; ...  
                                    'pag        PAG & hypothalamus + PAG(1) (See opioidBM)      '; ...
                                    'qqq        guest string (use this for temporary sets)      '];
return;
%%




function    out                 = local_unite(i1);
%% uniting left and right VOIs
out                             = [];
if isempty(i1);                                                                     return;         end;
i1s                             = int2str(i1);
if size(i1s,2)~=5;                                                                  return;         end;
i1s0                            = i1s;
s0                              = '0123456789';
s1                              = '0000444777';
for i=1:1:length(s0);           i1s(i1s(:,end-2)==s0(i),        end-2)      = s1(i);                end;
cm1                             = umo_cstrs(i1s,[],     'cm1');
out                             = str2num(i1s(cm1(:,2)>=1, :));

return;
%%


function    out                 = local_file(i1);
%% when a file is entered > get VOIID#s from the file:
f0                              = umo_getptf(i1,0,1);
out                             = str2num(f0);
if ~isempty(out);                                                                   return;         end;
fH                              = fopen(i1, 'r');
i0                              = char(fread(fH,    Inf,'uint8')');
fclose(fH);
eval(i0);
if exist('vfl','var');
    out                         = [];
    for i=1:1:length(vfl);      
        if isnumeric(vfl(i).vnos);
            out                 = [out; vfl(i).vnos];
        else;
            if exist(vfl(i).vnos,'file');
                                out                         = [out; local_file(vfl(i).vnos)];
            else;               o00                         = eval(vlf(i).vnos);
                                out                         = [out;     o00];                       end;
                                                                                            end;    end;
elseif exist('tfl','var');
    out                         = [];
    for i=1:1:length(tfl);      out                     = [out; tfl(i).vnos];                       end;
                                                                                                    end;

return;
%%

function    out                 = local_check;
%% make a list of local_* (not accessible?); use get_localFun.m instead
fH                              = fopen(which(mfilename),   'r');
ic                              = 0;
while 1;                        tline                       = fgetl(fH);
    if ~ischar(tline);                                                              break;          end;
    if size(tline,2)>20;        s                           = findstr(tline,'= local_');
        if length(s)==1 && size(tline,2)>=s(1)+11 && tline(1, s(1)+11)==';';
                                ic                          = ic + 1;
                                ff{ic}                      = tline(1,  s(1)+8:s(1)+10);            end;
                                                                                            end;    end;
fclose(fH);
out                             = char(ff);
return;
%%

function    out                 = local_fsl;
%% fls/first VOIs, less Cb and Po (Source: run_fsl_anat.m @local_voi)
out                             = [ 81100;  81200;  82100;  82200;  80180;  80280;  83100;  83200;  
                                    84100;  84200;  66100;  66200;  67100;  67200];
return;
%%

function    out                 = local_f81;
%% FS81 VOIs (modified: 8/31/2017; source: FreeSurferColorLUT.m)
c3                              = umo_getptf(which('FreeSurferColorLUT.m'), 0,3);
cm1                             = umo_cstrs(c3,[],          'cm1');
out                             = str2num(c3(cm1(:,2)>0,    :));
disp('f81');
return;
%%

function    out                 = local_mus;
%% MUSE VOIs (modified: 2/15/2023; source: MUSESurferColorLUT.m)
c3                              = umo_getptf(which('MUSEColorLUT.m'), 0,3);
cm1                             = umo_cstrs(c3,[],          'cm1');
out                             = str2num(c3(cm1(:,2)>0,    :));
disp('mus');
return;


function    out                 = local_f45;
%% Freesurfer 45 (a merged version of f81; source: fs_vnos4f45.m):
v0                              = fs_vnos4f45(1);
cm1                             = umo_cstrs(int2str(v0(:, 2)),[],   'cm1');
out                             = v0(cm1(:,2)>0,        2);
return;
%%

function    out                 = local_s28
%% c28 with vS added (good for 5HT6R studies)

out                             = [
       51100
       51200
       52100
       52200
       53100
       53200
       54100
       54200
       56100
       56200
       58100
       58200
       63100
       63200
       66100
       66200
       78100
       78200
       81100
       81200
       82100
       82200
       80180
       80280
       83100
       83200
       84100
       84200];
return;
%%

function    out                 = local_rat;
%%
out                             = [
       59100
       59200
       60000
       66100
       66200
       69000
       71000
       80100
       80200
       84100
       84200
       85000];
return;
%%

function    out                 = local_p18;
%%
out                             = [
81100                 
81200
82100
82200
81130
81230
82130
82230
81120
81220
82120
82220
81170
81270
82170
82270
80180
80280];
return;
%%

function    out                 = local_p20;
%% striatum subdivision VOIs, including merged VOIs such as whole striatum
out                             = [
       81100
       81200
       81130
       81230
       81120
       81220
       81170
       81270
       82100
       82200
       82130
       82230
       82120
       82220
       82170
       82270
       80180
       80280
       80100
       80200];
return;
%%
function    out                 = local_p10;
%% striatum subdivision VOIs
out                             = [
       80180
       80280
       81120
       81220
       81170
       81270
       82120
       82220
       82170
       82270];
return;
%%

function    out                 = local_m89;
%% VOIs of MC89 (checked: 8/31/2017)
if ~isempty(which('my_MRICloud89_LUT'));
    c1                          = umo_getptf(which('my_MRICloud89_LUT'),0,1);
else;
    c1                          = umo_getptf(which('MRICloud89_LUT'),0,1);                          end;
cm1                             = umo_cstrs(c1,[],          'cm1');
out                             = str2num(c1(cm1(:,2)>0,:));
return;
%%

function    out                 = local_m37;
%% VOIs of MC37 (checked: 8/31/2017)
c1                              = umo_getptf(which('MRICloud37_LUT'),0,1);
cm1                             = umo_cstrs(c1,[],          'cm1');
out                             = str2num(c1(cm1(:,2)>0,:));
return;
%%

function                        local_archive;
%%
ofl                             = fullfile(fileparts(which('vnosets')),'vnosets_archived.mat');
if exist(ofl,'file');
    f1                          = dir(ofl);
    f2                          = dir(which('vnosets'));
    if ~isempty(which('my_vnosets'));
                                f3                          = dir(which('my_vnosets'));
    else;                       f3                          = f2;                                   end;
    if f1.datenum>max([f2.datenum, f3.datenum]);
        disp('.vnosets_archived.mat is up to date');
        disp(' i.e., newer than vnosets.m or my_vnosets.m');                        return;         end;
                                                                                                    end;
%
[c1, c2]                        = umo_getptf(which(mfilename),0, 1);
if isempty(c1);                                                                     return;         end;
im1                             = umo_cstrs(c1,'function ', 'im1');
ss4                             = char(zeros(length(im1), 4)+'_');
for i=1:1:length(im1);
    s1                          = strfind(c2(im1(i),:),'local_');
    if ~isempty(s1);            ss4(i,  :)                  = c2(im1(i),s1(1)+5+[1:4]);     end;    end;
ss4(ss4(:, 4)==' ', 4)          = ';';
ss4s                            = ss4(ss4(:, 4)==';',   :);
ss4s(:, 4)                      = 'u';
vnos                            = zeros(9999,               1);
for i=1:1:size(ss4s,1);         vnos(vnosets(ss4s(i,:)),:)  = 1;                                    end;
y                               = load([which(mfilename),'at']);
for i=1:1:size(y.vstr,1);
    vnos(vnosets([y.vstr(i,:),'u']),    :)                  = 1;                                    end;

if ~isempty(which('my_vnosets'));
    [m1, m2]                    = umo_getptf(which('my_vnosets'),0,1);
    im2                         = umo_cstrs(m1,'function ', 'im1');
    mm4                         = char(zeros(length(im2), 4)+'_');
    for i=1:1:length(im2);
        s1                      = strfind(m2(im2(i),:),'local_');
        if ~isempty(s1);        mm4(i,  :)                  = m2(im2(i),s1(1)+5+[1:4]);     end;    end;
    mm4(mm4(:,4)==' ', 4)       = ';';
    mm4s                        = mm4(mm4(:,4)==';',    :);
    mm4s(:, 4)                  = 'u';
    for i=1:1:size(mm4s,1);     vnos(vnosets(mm4s(i,:)), :) = 1;                                    end;
    vvv                         = char(ss4s,y.vstr,mm4s);
else;                           vvv                         = char(ss4s,y.vstr);                    end;
vvv(:,  4)                      = 'u';
cm1                             = umo_cstrs(vvv,[],         'cm1');
vv2                             = vvv(cm1(:,2)>0,   :);
vno2                            = zeros(9999,               size(vv2,1));
for i=1:1:size(vv2,1);          vno2(vnosets(vv2(i,:)), i)  = 1;                                    end;
vst.vnos                        = find(sum(vno2,2)>0);
vst.mat                         = vno2(sum(vno2,2)>0,   :);
vst.vst                         = vv2;
vst.datenum                     = now;
vst.info                        = '.vnos = VOIID#s (lt/rt merged); .mat = .vno x set; .vst = setnames';
save(ofl,    'vst')
disp('.done! (archived vnosets)');
disp([' output .. ',ofl]);

return;
%%


