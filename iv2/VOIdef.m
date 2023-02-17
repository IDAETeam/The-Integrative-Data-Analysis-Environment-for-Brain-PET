function    out     = VOIdef(i1)

% VOIdef:       The definition file of VOIIDNos for IDAE
%       
%       usage:      out         = VOIdef(i1)
%       
%   i1          -   [] to get 'out' described below
%                   0 to check for duplications in VOIIDNos
%                   other numbsers (a vector) =>    out.anm = anatomical definition
%                                                   out.snm = short anat. defs.
%   Anatomical structires:
%   out.anm     -   anatomical definitions of VOIs
%   out.snm     -   initials of VOIs (in ~3 characters)
%
% Brodmann areas:
%   594xx where xx = BA#. Thus, 59632 = Rt. Brodmann area 32
%   See help messages of vL2_VOIUtility.m to learn how to use BAs.
% 
% (cL)2007/13    hkuwaba1@jhmi.edu 

%   out.vnos    -   VOIIDNos    (format: xx0yz, xx4yz, or xx7yz)
%                   left = +100 / right = + 200
%                   y = the position for descriptive terms
%                       descriptive terms may not be added if y is not 0
%                       (e.g., 51050 = medial frontal lobe)
%                   z = the position for the suffixes (for defining more than 1 VOIs)
%                       for a structure (e.g., 81200 & 81201~81209)
%                       suffixes are not allowed for structures when z is not 0
%                       (e.g., 69002 for raphe's nucleus)
%   out.sanms   -   short descriptions of anatomical names
%   Descriptive terms:
%   out.dnms    -   definitions of descriptive terms (e.g., superior)
%   out.dis     -   initials of descriptive terms
%   out.dnos    -   VOIIDNo to add to VOIIDNos for descriptive terms
%   out.sdnms   -   short descriptions of descriptive terms
% 
% (cL)2007    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

if isnan(i1);                   
    out                         = struct('anm','dammy',     'snm','?');             return;         end;
xxx                             = local_getdef;
if size(i1,2)==4 && i1(1)=='$'; i1                          = i1(1, 2:4);                           end;
if isempty(i1);                 out                         = xxx;                  return;         end;
if ischar(i1);                  out                         = local_str2vno(xxx,i1);
else;                           
    
    if ~i1(1);                  local_checkDupl(xxx);
                                out                         = xxx;                  return;         end;
                                out                         = local_vno2str(xxx,i1);
        if any(i1==68450);      out.anm(i1==68450, 1:20)    = 'median raphe nucleus';               end;
                                                                                                    end;
return;

function        out             = local_getdef;
%%
%
% structture #s are xx000, xx400, or xx700 (e.g.,  51400 = medial frontal lobe);
% left = +100 / right = +200
% 
% Looking for unsed VOIID#s for brain structures? See the end of this code


% anatomical strings
%long name           sym     IDNo       short name
astrs               = [
'frontal lobe               Fr      51000    frontal        b';
'frontal association area   FAA     51005    FAA            b';
'DLPFC                      DLF     51400    DLPFC          b';
'prefrontal cortex          pFC     76400    precrontalC    b';
'orbital gyrus              OG      73000    orbital        b';
'frontal operculum          FO      73400    f.operculum    b';
'orbital operculum          OO      73700    orb.operc      b';
'pars triangularis          PT      74400    pars.triang    b';
'rectal gyrus               RG      76000    rectal         b';
'parietal operculum         POp     77400    p.operculum    b';
'central operculum          COp     98700    MUSE-defined   b';
'cerebral exterior          CEx     98400    MUSE-defined   b';
'basal forebrain            BFB     89700    MUSE-defined   b';
'olfactory bulb             Ob      51700    olfactory      b';
'frontal pole               FP      51020    frontal pol    b';
'caudal middle frontal      cmFr    51050    caudal.m.fr    b';
'rostral middle frontal     rmFr    51015    rostral.m.f    b';
'temporal pole              TP      52020    temporal po    b';
'planum temporale           PlT     77700    planum.Tp      b';
'planum polare              PlP     78400    planum.Polare  b';
'fronto-orbital             FOb     31400    fronto-orbi    b';
'cingulate                  Cg      58000    cingulate      b';
'anterior cingulate         AC      58400    ant.cingul     b';
'isthmus/cingulate          IC      58700    isth.cingul    b';
'precentral gyrus           Pc      75000    preCnt         b';
'insula                     In      78000    insula         b';
'parietal lobe              Pa      53000    parietal       b';
'parietal association area  PAA     53005    PAA            b';
'postcentral gyrus          PS      74000    postCnt        b';
'paracentral                pC      53700    paraCnt        b';
'precuneus                  Pr      79000    precuneus      b';
'angular gyrus              AG      72000    angular        b';
'supramarginal              SM      94000    supmarg        b';
'temporal lobe              Tp      52000    temporal       b';
'temporal association area  TAA     52005    TAA            b';
'banks STS                  bT      52075    banks STS      b';
'transverse temporal        tT      52400    tv.temporal    b';
'tempotal pole              Tt      52700    temporal po    b';
'fusiform gyrus             Fs      56000    fusiform       b';
'occipital lobe             Oc      54000    occipital      b';
'cuneus                     Cu      53400    cuneus         b';
'lingual gyrus              LG      57000    lingual        b';
'calcarine cortex           V1      57400    calcarine      b';
'pericalcarine ctx          pCx     57700    pericalcari    b';
'cortex                     Cx      59000    cortex         b';
'Brodmann area              BA      59400    cortex         b';
'cerebellum                 Cb      71000    cerebellum     b';
'striatum                   S       80000    striatum       b';
'putamen                    Pu      81000    putamen        b';
'caudate-putamen            CP      81400    rodent         b';
'caudate nucleus            CN      82000    caudate        b';
'globus pallidus            GP      83000    glob.pal       b';
'fornix                     Fx      82400    fornix         b';
'thalamus                   Th      84000    thalamus       b';
'geniculate body            GB      84400    geniculate     b';
'hippocampus                Hp      66000    hippocamp      b';
'amygdala                   Am      67000    amygdala       b';
'parahippocampus            PH      63000    paraHipp       b';
'entorhinal area            ER      66400    entorhinal     b';
'subiculum                  Sb      66700    subiculum      b';
'substantia nigra           SN      68000    subs.nigra     b';
'brainstem                  BS      69000    brainstem      b';
'raphe nucleus              RA      68400    raphe          b';
'midbrain                   MB      60000    midbrain       b';
'pons                       PO      61000    pons           b';
'red nucleus                RN      61400    red nucleus    b';
'medulla oblongata          MO      68700    medulla        b';
'hypothalamus               HT      85000    hypothal       b';
'white matter               WM      90000    whiteMat       b';
'centrum semiovale          CS      95000    C.semioval     b';
'internal capsule           IC      95400    int.capsule    b';
'cerebellar WM              CW      90400    Cb.white.m     b';
'cerebral peduncle          CP      90700    crus.Cereb     b';
'corpus callosum            CC      92000    corp.call      b';
'ventricles                 Vn      91000    ventricles     b';
'dura mater                 Dm      93000    meninges       x';
'Sulcal CSF                 CSF     91700    Sulcal CSF     b';
'tegmentum area             TA      62700    tegmentum      b';
'periaqueductal gray        PAG     69400    PAG            b';
'peri-ventr.gray            PVG     69700    PVG            b';
'pituitary grand            PG      93700    pituitgrand    b';
'skull                      SK      99700    skull          b';
'subcallosal area           SCA     99400    subCC.area     b';
'colliculus                 Cc      93400    colliculus     b';
'Braak I                    B1      45010    Braak I        b';
'Braak II                   B2      45020    Braak II       b';
'Braak III                  B3      45030    Braak III      b';
'Braak IV                   B4      45040    Braak IV       b';
'Braak V                    B5      45050    Braak V        b';
'Braak VI                   B6      45060    Braak VI       b';
'bone marrow                BM      31000    bone marrow    x';
'venous sinus               VS      91400    v. sinus       x';
'arterial blood             AB      10000    Ca(t)          x';
'reference region           eR      11000    simulated      x';
'2ndary markings            2m      20000    secondary      x';
'choroid plexus             CPx     44700    choroid ple    x';
'ACA territory              ACA     41000    ACA            b';
'MCA territory              MCA     42000    MCA            b';
'PCA territory              PCA     43000    PCA            b';
'COSS CS                    CCS     41400    CsCS           b';
'COSS PA                    PA      41700    CsPA           b';
'COSS GSM                   GSM     42400    CsGSM          b';
'COSS F2                    F2      42700    CsF2           b';
'COSS F3                    F3      43400    CsF3           b';
'COSS T2                    T2      43700    CsT2           b';
'COSS T1                    T1      44000    CsT1           b';
'unknown                    TBD     44400    unknown        x';
'WM hyperintesity           WMH     92400    WMH            b';
'water shed area            WS      96000    water shed     b';
'SPM Cluster                SPM     99000    spm cluster    b';
'pineal body                PB      92700    pineal         b';
'body                       B       11700    body           x';
'retina                     RT      11400    retina         x';
'nasal cav.                 NC      12000    nasal cav.     x';
'aorta                      AR      18000    aorta          x';
'heart                      HRT     16000    heart          x';
'lung                       Lg      14000    lung           x';
'liver                      Lv      21000    liver          x';
'gall bladder               GB      23000    gall bldder    x';
'spleen                     SP      22000    spleen         x';
'kidney                     KD      32000    kidney-ctx     x';
'bile duct                  BD      23400    bile duct      x';
'small Intst                SI      24000    small Intst    x';
'colon                      CLN     24400    colon          x';
'spine                      SPN     38000    spine          x';
'u. bladder                 UB      35000    u. bladder     x';
'prostate                   PRS     37000    prostate       x';
'abdomen                    ABD     28000    abdomen        x';
'stomach                    STM     28400    stomach        x';
'parotid grand              PGd     12400    parotid g      x';
'reproductive organ         RO      15400    reprod m&f     x';
'submandibular grand        SMG     12700    submand g      x';
'thyloid grand              TYG     13400    tyloid         x';
'diaphragm                  DPH     14400    diaphragm      x';
'tumor                      TM      13700    tumor          x';
'muscle                     Ms      15000    muscle         x';
'lymph nodes                LYN     14700    lymph nodes    x'];
%12345678901234567890123456789012345678901234567890123456789012
%
dstrs               = [ ...
'superior            s       10         sup        ';
'inferior            i       60         inf        ';
'anterior            a       20         ant        ';
'posterior           p       70         post       ';
'dorsal              d       30         dors       ';
'ventral             v       80         vent       ';
'medial              md      40         medial     ';
'lateral             l       90         lateral    ';
'rostral             r       15         rost       ';
'caudal              c       65         caud       ';
'mesial              ms      35         mesial     ';
'middle              m       50         mid        ';
'anterior middle     am      45         ant mid    ';
'posterior middle    pm      55         post mid   ';
'dorsal posterior    dp      75         dors post  ';
'ventral posterior   vp      85         vent post  ';
'pregenual           pg      25         pregenual  ';
'subgenual           sg      95         subgenual  '];

% 'medial frontal lobe mdF     51059      med.frontal    b';
% 'Lateral orbital lobeLOL     51700      Lat.Orbital    b';


%123456789012345678901234567890123456789012345678901234567890
out.anms            = deblank(astrs(:,  1:27));
out.ais             = deblank(astrs(:,  28:35));
out.vnos            = str2num(astrs(:,  36:40)); 

out.sanms           = deblank(astrs(:,  45:49)); 
out.bx              = astrs(:,      end);

out.dnms            = deblank(dstrs(:, 	1:20)); 
out.dis             = deblank(dstrs(:, 	21:28)); 
out.dnos            = str2num(dstrs(:,  29:39));
out.sdnms           = deblank(dstrs(:,  40:end));

return;
%%

function                        local_checkDupl(out);
%% checking duplications in structure names:
cm1                 = umo_cstrs(char(out.anms),[],  'cm1');
if any(cm1(:,2)>1); disp('Duplications in structure names');
                    k                           = find(cm1(:,2)>1);
    for i=1:1:length(k);
                    disp('***');
                    q                           = find(cm1(:,1)==cm1(k(i),1));
                    disp([out.anms(q,:),char(zeros(length(q),3)+32),int2str(out.vnos(q,:))]);       end;
else;               disp('No duplications in structure names');                                     end;
%% checking duplications in VOIIDNos:
d                   = zeros(max(out.vnos),1);
for i=1:1:length(out.vnos);         d(out.vnos(i),:)        = d(out.vnos(i),1) + 1;                 end;
if any(d>1);        disp('Duplications in VOIIDNos');
                    ii                          = find(d>1);
    for i=1:1:length(ii);           
                    j                           = find(out.vnos==ii(i));
        for k=1:1:length(j);        disp([out.anms(j(k),:),' ... ',int2str(out.vnos(j(k)))]);       end;
                                                                                                    end;
else;               disp('No duplications in VOIIDNos');                                            end;
%% return;


function    out                 = local_vno2str(xxx,i2);
%%
out                             = [];
i2s                             = int2str(i2);
if size(i2s,2)~=5;            	disp('.error! ID#s must be 5 digits or less');      return;         end;

i2ss                            = i2s(i2s(:,1)~=' ', :);
i2x                             = zeros(size(i2ss,1),   4);
i2x(:,  2)                      = i2(i2s(:,1)~=' ',     1);
i2x(i2ss(:,3)=='1' | i2ss(:,3)=='5' | i2ss(:,3)=='8', 1)    = 1;
i2x(i2ss(:,3)=='2' | i2ss(:,3)=='6' | i2ss(:,3)=='9', 1)    = 2;
i2x(:,  2)                     	= i2x(:, 2) - i2x(:, 1).*100;
%
i2x(:,  3)                      = floor(str2num(i2ss(:, 4:5))./5).*5;
i2x(:,  4)                      = str2num(i2ss(:, 4:5)) - i2x(:, 3);
i2x(:,  2)                      = i2x(:, 2) - i2x(:, 4);
%
for i=1:1:size(i2x,1);          
    if i2x(i,4)>0;              add{i}                      = ['-',int2str(i2x(i,4))];   
    else;                       add{i}                      = '';                           end;    end;
% 
sid                             = {'', 'left ', 'right '};
sids                            = {'', 'L.', 'R.'};
%
im1                             = umo_cstrs(int2str(xxx.vnos), int2str(i2x(:,2)), 'im1');
im2                             = umo_cstrs(int2str(xxx.vnos), int2str(i2x(:,2)-i2x(:,3)), 'im1');
for i=1:1:size(i2x,1);
    if im1(i)>0;
            anm{i}            	= [sid{i2x(i,1)+1},deblank(xxx.anms(im1(i),:)),add{i}];
            snm{i}              = [sids{i2x(i,1)+1},deblank(xxx.ais(im1(i),:)),add{i}];
    else;
            anm{i}            	= [sid{i2x(i,1)+1},deblank(xxx.dnms(xxx.dnos==i2x(i,3), :)),    ...
                                    ' ',deblank(xxx.anms(im2(i),:)),add{i}];
            snm{i}            	= [sids{i2x(i,1)+1},deblank(xxx.dis(xxx.dnos==i2x(i,3), :)),   ...
                                    deblank(xxx.ais(im2(i),:)),add{i}];                     end;    end;
%
out                             = struct('anm',char(anm),   'snm',char(snm));              
            
    

return;
%xxx.vnos -


if size(i2s,2)<5;
    for i=1:1:size(i2,1);       anm{i}                      = ['region ',int2srt(i2(i))];
                                snm{i}                      = ['R',int2srt(i2(i))];                 end;
    out                         = struct('anm',char(anm),   'snm',char(snm));       return;         end;
%
imq                             = umo_cstrs(lower(xxx.anms), ['spm cluster';'brodmann ar'], 'im1');
%
vxx                             = int2str(xxx.vnos);
% i2s(:,3) indicates whole/Lt/Rt > giving whole since vxx are always whole:
i2s3                            = str2num(i2s(:,3));
i2s3x                          	= zeros(size(i2s3));
i2s3x(i2s3>=4 & i2s3<7,  :)    	= 4;
i2s3x(i2s3>=7,   :)             = 7;
% [i2s3, i2s3x]
% i2s(:,5) can indicate default or copy x > giving the defaults for now: 
i2s5                            = str2num(i2s(:,5));
i2s5x                           = zeros(size(i2s5));
i2s5x(i2s5<5,   :)              = 0;
i2s5x(i2s5>=5,  :)              = 5;
% searching by characters 1-5:
im15                            = umo_cstrs(vxx,[i2s(:,1:2),int2str(i2s3x),i2s(:,4),    ...
                                                            int2str(i2s5x)],    'im1');
% searching by characters 1-3 + '00':
im1300                          = umo_cstrs(vxx, [i2s(:,[1,2]),int2str(i2s3x),          ...
                                                            repmat('0',size(i2s,1),2)], 'im1');
% marking SPM clusters & Brodmann areas, if any:
im13x                           = umo_cstrs(int2str(umo_cstrs(lower(xxx.anms), ...
                                    ['spm cluster';'brodmann ar'], 'im1')), int2str(im1300),    'im1');
% search criterion #1 = by 5 digits with SPM & BA removed:
s1                              = im15.*double(im13x<1);
% search criterion #2 = SPM & BA alone:
s2                              = im1300.*double(im13x>0);
% search criterion #3 = 
s3                              = im1300.*double(s1<1).*double(s2<1);
% descriptive terms:
im45                            = umo_cstrs(int2str(xxx.dnos), [i2s(:,4),int2str(i2s5x)], 'im1');
% some regions are using all 5 digits 
%
% vxr                             = vxx;
% vxr(vxx(:,4)~='0' | vxx(:,5)~='0',  :)                      = 'x';
% [im15, im13(:,1), im13x]                
%
sid                             = {'','left ','right '};
sis                             = {'','Lt','Rt'};
i2s3x(:)                        = i2s3 - i2s3x + 1;
% i2s45                           = str2num(i2s(:,4)).*10 + i2s5x;
i2s5(:)                         = i2s5 - i2s5x;
for i=1:1:size(s1,1);
    if s1(i)>0;
        anm{i}              	= [sid{i2s3x(i)}, deblank(xxx.anms(s1(i),:))];
        snm{i}                  = [sis{i2s3x(i)}, deblank(xxx.ais(s1(i),:))];
        if i2s5(i)>0;           anm{i}                      = [anm{i},'#',int2str(i2s5(i))];
                                snm{i}                      = [snm{i},int2str(i2s5(i))];            end;
    elseif s2(i)>0;
        anm{i}                  = [deblank(xxx.anms(s2(i),:)),' ',i2s(i, 4:5)];
        snm{i}                  = [deblank(xxx.ais(s2(i),:)),i2s(i, 4:5)];
    elseif s3(i)>0;             anm{i}                      = sid{i2s3x(i)};
                                snm{i}                      = sis{i2s3x(i)};
        if im45(i)>0;           anm{i}                      = [anm{i},deblank(xxx.dnms(im45(i),:)),' '];
                                snm{i}                      = [snm{i},deblank(xxx.dis(im45(i),:))];	end;
        %
        anm{i}                  = [anm{i}, deblank(xxx.anms(s3(i),:))];
        snm{i}              	= [snm{i}, deblank(xxx.ais(s3(i),:))];
        if i2s5(i)>0;           anm{i}                      = [anm{i},'#',int2str(i2s5(i))];
                                snm{i}                      = [snm{i},int2str(i2s5(i))];            end;
    else;                       anm{i}                      = ['unknown (',i2s(i, :),')'];
                                snm{i}                      = '?';                          end;    end;
%
out                             = struct('anm',char(anm),   'snm',char(snm));                                
return;
%%

function        out             = local_str2vno(xxx,istr);
%% from a string to VOIIDNo componemts
%
%   out = [VOIIDNo, side#, descriptive#, add#]

istr                            = deblank(istr);
out                             = zeros(1,  4);
%
[s1, s2]                        = getLseg(istr, 1);
out(:,  2)                      = umo_cstrs(['left ';'right'],istr, 'im1').*100;
if out(1,2)>0;                    istr                        = s2;                                   end;

% when istr is Brodmann area:
if strncmpi(istr,'brodmann area',13);
    out(:,  1)                  = 59400;
    out(:,  3)                  = str2num(istr(1,   14:end));                       return;         end;
%
dd                              = zeros(size(xxx.dnms,1),   1);
for i=1:1:size(dd,1);
    if strncmp(istr,xxx.dnms(i,:),sum(xxx.dnms(i,:)~=' ')+1);
                                dd(i,   :)                  = sum(xxx.dnms(i,:)~=' ');      end;    end;
%
if any(dd>0);                   [v, im]                     = max(dd);
                                out(:,  3)                  = xxx.dnos(im);
                                istr(1, 1:v+1)            	= ' ';
                                qq                          = deblank(istr(end:-1:1));
                                istr                        = deblank(qq(end:-1:1));            	end;
% checking if more-than-one-VOI-per-structure flag (=add#) is present:
if any(istr=='(');
    out(:,  4)                  = str2num(istr(1, find(istr=='(',1)+1:find(istr==')',1)-1));
    istr(1, find(istr=='(',1):end)                          = ' ';                                  end;
%
im1                             = umo_cstrs(xxx.anms,istr,  'im1');
out(:,  1)                      = xxx.vnos(im1,     :);
return;
%%

v0 = zeros(150,1); 
ic = 0; 
for i=50:1:99; for j=[0,4,7];   ic                          = ic + 1; 
                                v0(ic, :)                   = i.*1000+j.*100;               end;    end
%
v                               = VOIdef([]);
vi                              = consolidVOINos(v.vnos, v0);
% to display available VOIID#s:
vi(vi(:,2)<1, 1)