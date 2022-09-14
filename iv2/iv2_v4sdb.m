function    sss               = iv2_v4sdb(i1); 

% iv2_v4sdb:    To return/dispaly items for scanDB.m (ver.iv2)
%       
%       usage:      sss         = iv2_v4sdb(flg)
% 
% Inputs:
%   flg 1   -   to return registered scanDB items and their descriptions
%       0   -   to display registered scanDB items and their descriptions  
%       3   -   to check if short and long strings (defined here) are unique 
%       xs  -   to display one-per-scan variables alone
%       xa  -   to display variables of any # of columns (unform across subjects) 
%
% Outputs:
%   sss{:,1}    - short strings
%   sss{:,2}    - long strings
%   sss{:,3}    - three characters showing ...
%               character/numner &  must/reserved/not(=x) & perProject/Subject/Condition/Any
%   sss{:,4}    - descriptions
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

% if ischar(i1) | length(i1)<1;                                                       return;         end;

% sss{:,1}  short strings
% sss{:,2}  long strings
% sss{:,3}  character/numner &  must/reserved/not(=x) & perProject/Subject/Condition/Any
% sss{:,4}  descriptions

% strings that are reserved for the umo file format
[a, b]                      = um_info(4,    1);


for i=1:1:size(a,1);        umo{i, 1}               = deblank(a(i,  :));
                            umo{i, 2}               = umo{i, 1};
                            umo{i, 3}               = 'xrx';
                            umo{i, 4}               = 'See <<um_info>>';                            end;

sss                         = {
'%       ', '(long name)    ','   ',   '... descriptions',
'lds',      'lds4dxetc',      'cmp',   'local directory system file',
'cnd',      'condNames',      'cmp',   'experimental conditions (cndx) names (short)',
'snm',      'subjNames',      'cms',   'study subject IDs',
'g',        'groupNames',     'cms',   'group initials',
'g2',       'gName2nd',       'cxs',   'secondary group initials (see groupNames)',
'did',      'databaseIDs',    'cms',   'subject IDs for PET variable database (one per subject)', 
'pet',      'PETStatus',      'nxs',   '0=Not done; 1=Plasma data -; 2=PD +; 7/8=TACs alone; 9=Not usable',
'%',        '  ',             '   ',   '...',
'cnd0',     'condDescript',   'crp',   'short descriptions of experimental conditions',
'unm',      'IDAEuser',       'crp',   'IDAE username (a restricted item)',
'ipj',      'projectName',    'crp',   'IDAE project name (taken from project folder name)',
'dxs',      'dxs4mri',        'crp',   'reserved for condition-independent directories (=mri)',
'dx4',      'dx4c***',        'crp',   'reserved for condition-dependent directories (=pet)',
'tnm',      'radioligands',   'crp',   'reserved for radioligand names for PET conditions',
'%',        '  ',             '   ',   '...',
'age',      'subjAge',        'nxs',   'age @1st PET',
'sex',      'subjSex',        'cxs',   'Gender initials',
'race',     'subjRace',       'cxs',   'race in 1 character (W/B/A etc)',
'bw',       'subjWeight',     'nxs',   'subject body weight in kg',
'sns',      'scanNotes',      'cxs',   'notes on scans (e.g., troubles)',
'%',        '  ',             '   ',   '...',
'ssf',      'subjSelect',     'nxs',   'subject selection field',
'sid',      'subjIDstr',      'cxs',   'subject identification string (=snm but not shown)',
'ids',      'subjID02',       'cxs',   'subject identification string (in addition to snm)',
'cid',      'subjCTstr',      'cxs',   'subject CT string for microPET studies',
'%',        '  ',             '   ',   '...',
'mass',     'coldmass',       'nxc',   'injected cold mass in micrograms',
'rad',      'RADoseInjected', 'nxc',   'injected dose (in site-specific unit)',
'sra',      'SA(mCi/umol)',   'nxc',   'specific radio-activity in (mCi/micromol)',
'%',        '  ',             '   ',   '...',
'doseA',    'dose4all',       'nxp',   'doses of test drugs for conditions common to all subject',
'doseS',    'dose4eachsubj',  'nxs',   'doses of test drugs for conditions for each subject',
'PK',       'plasmaConc',     'nxc',   'plasma concentration of the test drug',
'blocker',  'blockingAgent',  'cxs',   'name of the blocking agent',
'%',        '  ',             '   ',   '...',
'dtTx',     'durationOfTx',   'nxc',   'time since the start of a treatment',
'dtDx',     'sinceDiagnosis', 'nxc',   'years since diagnosis',
'%',        '  ',             '   ',   '...',
'mri',      'haveMRI',        'nxp',   'have MRI (=1) or not (=0)',
'hrr',      'HRRTorNot',      'nxp',   'PET study done with HRRT (=1) or not (=0)',
'ham',      'subjSpec',       'cxp',   'subject species; h=human; a=non-human; m=micriPET',
'%',        '  ',             '   ',   '...',
'cDA',      'clnDatA',        'cxs',   'generic field for clinical data A (character, any length)',
'cDB',      'clnDatB',        'cxs',   'generic field for clinical data B',
'cDC',      'clnDatC',        'cxs',   'generic field for clinical data C',
'cD0',      'clnDat0',        'nxa',   'generic field for clinical data #0 (numeric, any length)',
'cD1',      'clnDat1',        'nxa',   'generic field for clinical data #1',
'cD2',      'clnDat2',        'nxa',   'generic field for clinical data #2',
'cD3',      'clnDat3',        'nxa',   'generic field for clinical data #3',
'cD4',      'clnDat4',        'nxa',   'generic field for clinical data #4',
'cD5',      'clnDat5',        'nxa',   'generic field for clinical data #5',
'cD6',      'clnDat6',        'nxa',   'generic field for clinical data #6',
'cD7',      'clnDat7',        'nxa',   'generic field for clinical data #7',
'cD8',      'clnDat8',        'nxa',   'generic field for clinical data #8',
'cD9',      'clnDat9',        'nxa',   'generic field for clinical data #9',
'bDA',      'behDatA',        'cxs',   'generic field for behavioral data A (character, any length)',
'bDB',      'behDatB',        'cxs',   'generic field for behavioral data B',
'bDC',      'behDatC',        'cxs',   'generic field for behavioral data C',
'bD0',      'behDat0',        'nxa',   'generic field for behavioral data #0 (numeric, any length)',
'bD1',      'behDat1',        'nxa',   'generic field for behavioral data #1',
'bD2',      'behDat2',        'nxa',   'generic field for behavioral data #2',
'bD3',      'behDat3',        'nxa',   'generic field for behavioral data #3',
'bD4',      'behDat4',        'nxa',   'generic field for behavioral data #4',
'bD5',      'behDat5',        'nxa',   'generic field for behavioral data #5',
'bD6',      'behDat6',        'nxa',   'generic field for behavioral data #6',
'bD7',      'behDat7',        'nxa',   'generic field for behavioral data #7',
'bD8',      'behDat8',        'nxa',   'generic field for behavioral data #8',
'bD9',      'behDat9',        'nxa',   'generic field for behavioral data #9',   
'tDA',      'tecDatA',        'cxs',   'generic field for technical data (character, any length)',
'tDB',      'tecDatB',        'cxs',   'generic field for technical data B',
'tDC',      'tecDatC',        'cxs',   'generic field for technical data C',
'tD0',      'tecDat0',        'nxa',   'generic field for technical data #0 (numeric, any length)',
'tD1',      'tecDat1',        'nxa',   'generic field for technical data #1',
'tD2',      'tecDat2',        'nxa',   'generic field for technical data #2',
'tD3',      'tecDat3',        'nxa',   'generic field for technical data #3',
'tD4',      'tecDat4',        'nxa',   'generic field for technical data #4',
'tD5',      'tecDat5',        'nxa',   'generic field for technical data #5',
'tD6',      'tecDat6',        'nxa',   'generic field for technical data #6',
'tD7',      'tecDat7',        'nxa',   'generic field for technical data #7',
'tD8',      'tecDat8',        'nxa',   'generic field for technical data #8',
'tD9',      'tecDat9',        'nxa',   'generic field for technical data #9',   
'm4p1',     'm4pet_1',        'cxa',   'for MRI directories MRI / PET paires (longitudinal)',
'm4p2',     'm4pet_2',        'cxa',   'for MRI directories MRI / PET paires (longitudinal)',
'm4p3',     'm4pet_3',        'cxa',   'for MRI directories MRI / PET paires (longitudinal)',
'm4p4',     'm4pet_4',        'cxa',   'for MRI directories MRI / PET paires (longitudinal)',
'm4p5',     'm4pet_5',        'cxa',   'for MRI directories MRI / PET paires (longitudinal)',
'm2p',      'mri2pet',        'nxp',   '1 x # of PETs showing MRI #s of m2pi for PETs',
'p2p'       'pet4p2p',        'nxa',   '1 x # of PETs of PET group #s for PET2PET coreg'};


% 'IDAEuser', 'IDAEuser',       'crp',   'reserved for IDAE user',

if ischar(i1);
    if size(i1,2)~=2;           disp(['.wrong input: ',i1,' (must be 1 x 2)']);     return;         end;
    c3                          = char(sss(:,3));
    im1                         = umo_cstrs(c3(:, 2:3),i1(:,1:2),   'im1');
    if ~im1;                    disp(['.wrong input: ',i1]);                        return;         end;
    if i1(2)=='s';              disp('#one-per-scan variables');
    elseif i1(2)=='a';          disp('#variables of any # of columns (unform across subjects');     end;
    disp(' enter 1st column variables to scanDB.m followed by their values');
    ss8                         = char(zeros(1, 8) + 32);
    k                           = im1(find(c3(im1,1)=='c',1));
    if ~isempty(k);             ss8(1, 1:size(sss{k,1},2))  = sss{k,    1};
                                disp([' e.g., ',ss8,'abc']);                                        end;
    k                           = im1(find(c3(im1,1)=='n',1));
    if ~isempty(k);             ss8(:)                      = ' ';
                                ss8(1, 1:size(sss{k,1},2))  = sss{k,    1};
                                disp([' e.g., ',ss8,'3.2,14.12,3.56']);                             end;
    disp(' enter values separated by , without spaces if numbers (e.g., 3.2,14.12,3.56)');
    disp(' 2nd colmn: n=numbers; c=characters');
    dispCharArrays(char(sss(im1,1)),3,c3(im1,1),3,char(sss(im1,4)));                return;         end;
%
if i1(1)==3;
    str                         = {'some short strings are not unieque',    
                                    'some long strings (1-7) are not unieque'};
    dOK                         = ones(2,   1);
    for i=1:1:2;
        cx                      = char(char(umo{1:7,i}),char(sss{1:7,i}));
        if i==1;                n2u                         = find(cx(:,1)=='%');                   end;
        cm1                     = umo_cstrs(cx,[],          'cm1');
        cm1(n2u,    2)          = 1;
        if any(cm1(:,2)>1);     disp(str{i});
                                dOK(i,  :)                  = 0;
                                xxx                         = zeros(max(cm1(cm1(:,2)==0,1)),    1);
                                xxx(cm1(cm1(:,2)==0,1), :)  = 1;
                                ii                          = find(xxx);
            for j=1:1:length(ii);   
                                k                           = find(cm1(:,1)==ii(j));
                                disp(['duplication ',int2str(j),' ... ']);
                                disp(cx(k,:));                                                      end;
                                                                                            end;    end;
    if ~any(~dOK);              disp('short and long strings are unique (=OK)');                    end;
    return;                                                                                         end;

if ~i1(1);                      mmm                         = char(zeros(size(sss,1), 6)+32);
                                mmm(:,  2:4)                = char(sss{:,   3});
                                mmm(:,  1)                  = '(';
                                mmm(:,  5)                  = ')';
                                mm1                         = char(sss{:,   1});
                                mmm(mm1(:,1)=='%',  :)      = ' ';
                                disp(['*** items for scanDB.m',10,'(c1=short name to enter; ',   ...
                                'c3=character/numner|must/reserved/x|perProject/Subject/Any)']);
                                disp([mm1, char(sss{:,2}), mmm, char(sss{:,4})]);                   end;

