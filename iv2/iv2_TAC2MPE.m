function    iv2_TAC2MPE(i1,i2); 

% To generate TAC2MPEs (iPacks for TAC generation + modeling) for IDAE.iv2 
%       
%       usage:      iv2_TAC2MPE(i1,fNo)
%       
%   i1      -   .mat file of information of VOIs for TAC generation 
%               (generated at Step 4 of VOI prep procedure)
%   fNo     -   figure # of IDAE Level 1 Window (L1W)
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
if ~exist(i1,'file');                                                               return;         end;

global g4iv2;
% some ifiles depends of # of MRIs, if a longitudinal study:
if isfield(g4iv2.yyy,'nMRI');  	mx_add                      = ['L',int2str(g4iv2.yyy.nMRI)];
else;                           mx_add                      = '';                                   end;
%
% retieving values of parameters from prepMPxxx:
%   to cope with cases when this code is called while non-prepMPxxx is up:
q23                             = load(fullfile(g4iv2.yyy.idx,'iv2',[g4iv2.xxx(1).pmp,'.mat']));
% items to copy from the stage-1 package to TAC2MPE:
ss0                             = ['mid';'pid';'ifc';'pmp'];
ss1                             = ['mid';'pid';'pio';'pmp'];
for i=1:1:size(ss0,1);
    eval(['belowL{i}         	= [''',ss1(i, :),repmat(' ',1,13),''',q23.prepMP.',ss0(i,:),'];']); end;
%
%
[idx, inm, iex]               	= fileparts(i1);
load(i1);
% 
ic                              = 0;
if strcmpi(Info4TACs.hmc,'hmcMIT');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4HMC        Head Motion Correction (SPM12-MIT=ecc)';
    belowL{end+1}               = 'IDAE4HMC        iv2_hmcMIT';
    belowL{end+1}               = 'ifc             *pio_hmcMIT';
elseif strcmpi(Info4TACs.hmc,'hmcMITnmi');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4HMC        Head Motion Correction (SPM12-MIT=nmi)';
    belowL{end+1}               = 'IDAE4HMC        iv2_hmcMITcfn';
    belowL{end+1}               = 'cfn             nmi';
    belowL{end+1}               = 'ifc             *pio_hmcMIT*cfn';
elseif strcmpi(Info4TACs.hmc,'hmcMIT6mm');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4HMC        Head Motion Correction (SPM12-MIT=nmi6mm)';
    belowL{end+1}               = 'IDAE4HMC        iv2_hmcMITcfn';
    belowL{end+1}               = 'cfn             nmi6mm';
    belowL{end+1}               = 'ifc             *pio_hmcMIT*cfn';
elseif strcmpi(Info4TACs.hmc,'hmcMITncc');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4HMC        Head Motion Correction (SPM12-MIT=ncc)';
    belowL{end+1}               = 'IDAE4HMC        iv2_hmcMITcfn';
    belowL{end+1}               = 'cfn             ncc';
    belowL{end+1}               = 'ifc             *pio_hmcMIT*cfn';
elseif strcmpi(Info4TACs.hmc,'hmcP2M');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4HMC        Head Motion Correction (MIT + frame2MRI)';
    belowL{end+1}               = 'IDAE4HMC        iv2_hmcP2M';
    belowL{end+1}               = 'ifc             *pio_hmcP2M';
elseif strcmpi(Info4TACs.hmc,'hmcSPM');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4HMC        Head Motion Correction (MIT + spm_realign)';
    belowL{end+1}               = 'IDAE4HMC        iv2_hmcSPM';
    belowL{end+1}               = 'ifc             *pio_hmcSPM';
else;
    belowL{end+1}               = 'ifc             *pio';                                           end;
%
belowL{end+1}                   = ['avr             ',Info4TACs.avr];
%
% IDAE4TACs:
if strcmpi(Info4TACs.m2p,'m2p');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4TACs        VOI transfer & TAC generation';
    if strcmpi(g4iv2.xxx(1).pmp,'prepMPbab');
        belowL{end+1}           = 'IDAE4TACs        iv2_bab_genTACs';
    else;
        belowL{end+1}           = ['IDAE4TACs     	 iv2_MPcoreg',mx_add];
        belowL{end+1}           = ['IDAE4TACs        iv2_genTACs_v2',mx_add];                     	end;
%     belowL{iq}                  = 'IDAE4TACs       iv2_mpsVOIs';
    belowL{end+1}               = ['v4t             ',fullfile('iv2', [inm,iex])];
    belowL{end+1}               = 'mpp             pet\*ifc_*avr.ezi';                              end;
% iq                              = iq + 1;
belowL{end+1}                   = ['eza             ',Info4TACs.eza];
%
% inserting SN section, if requested:
if ~strcmpi(Info4TACs.sm,'NoSN');
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4SN         Spacial normalization';
    belowL{end+1}               = ['IDAE4SN         iv2_SN_fMaps_',Info4TACs.sm,'_',mx_add];        end;
%
% the field 'snu' is needed whether to perform SN or not:
belowL{end+1}                   = ['snu',repmat(' ',1,13),Info4TACs.sm];
%
% inserting IDAE4PIMs, if not nocpt:
if ~strcmpi(Info4TACs.cpt,'nocpt');  
    ic                          = ic + 1;
    aboveL{ic}                  = 'IDAE4PIMs       Plasma input methods';
    %
    belowL{end+1}               = ['IDAE4PIMs       ',feval(g4iv2.yyy.lds,'cpt',Info4TACs.cpt)];
    belowL{end+1}               = 'IDAE4PIMs       iv2_doPIMs';                                     end;


% tissue reference methods
ic                              = ic + 1;
aboveL{ic}                      = 'IDAE4RTMs       Reference tissue methods';
belowL{end+1}                   = ['IDAE4RTMs       iv2_doRTMs',mx_add];
% 'cpt' is needed whether to perform PIMs or not:
belowL{end+1}                   = ['cpt             ',Info4TACs.cpt];
% writiung the new iPack:
ofln                            = fullfile(idx,             [inm,'.m']);
fH                              = fopen(ofln,               'w');
if fH<0;                        disp('.problem: unable to open new iPack to write (no output)');
                                disp([' file: ',ofln]);                             return;         end;
%
fwrite(fH,  ['% generated by ',mfilename,' (',datestr(now),')',10],         'char');
for i=1:1:numel(aboveL);        fwrite(fH,  [aboveL{i},10], 'char');                                end;
fwrite(fH,  ['***',10],         'char');
for i=1:1:numel(belowL);       	fwrite(fH,  [belowL{i},10], 'char');                                end;
fclose(fH);
%
disp('.done! (iPack for TAC generation + MPE)');
disp([' output: ',ofln]);
%%
% getting lines of the iPack list file (.ev2):
disp('.revising the package list file');
ffc                             = umo_getptf(g4iv2.yyy.ev2, 1,[]);
c1                              = getLseg(ffc,              1);

% looking for g4iv2.xxx(1).pmp (such as prepMPmcx):
im1                             = umo_cstrs(lower(c1),[lower(g4iv2.xxx(1).pmp),' '],'im1');
if ~im1(1);                     disp(['.error! ',mfilename,'@L147']);               return;         end;
%
% there is only one lins that start with q23.prepMP.pmp
im1                             = umo_cstrs(c1, [q23.prepMP.pmp,' '],   'im1');
fH                              = fopen(g4iv2.yyy.ev2,  'w');
for i=1:1:im1(1);               fwrite(fH,      [deblank(ffc(i, :)),10],    'char');                end;
%
fwrite(fH,                      [inm,'      ',Info4TACs.sstr,10],           'char');
for i=im1(1)+1:1:size(ffc,1);   fwrite(fH,      [deblank(ffc(i, :)),10],    'char');                end;
fclose(fH);
%
disp('.done! (the package list file)');
return;
%%
