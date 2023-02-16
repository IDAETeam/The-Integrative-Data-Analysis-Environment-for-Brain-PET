function    out                 = iv2_tracerList(i1); 

% To return if radiotracers in input1 are registered for IDAE 
%       
%       usage:      x           = iv2_tracerList(listInCharArray)
%       
%   output x is [] if one of input1 is not registered. Otherwise, x will be
%   the same as input1 (except , will be removed)
%   To list current list:
%       >> iv2_tracerList([]);
% 
% Missing tracers?
%   Create a function called my_tracerList.m as follows:
%   Open a matlab editor blank page and copy/paste/edit the following lines there:
%
%       function    out         = my_tracerlist(i1);
%   %
%   out                         = {'[11C]tracerA',  value of k2R,   'target system'
%                                   '[18F]tracerB,  value of k2R,   'target system'
%                                   ...
%                                   '[18F]tracerZ,  nan,            'target system'};
%
%       value of k2R is the value of k2 of reference tissue, if known
%       if not known, enter nan (=not a number).
%
%   save it in the same folder as iv2_tracerList.m
%       (>> which iv2_tracerList will list the full path)
%
% (cL)2016    hkuwaba1@jhmi.edu 

margin                          = 1;
out                             = [];
if nargin<margin;               helq(mfilename);                                    return;         end
ttt                             = {
                                '[11C]rac ',            0.163,          'Dopamine D2/D3'
                                '[11C]raclopride',      0.163,          'Dopamine D2/D3'
                                '[18F]FPEB',            0.043,          'mGluR5'
                                '[11C]DTBZ',            nan,            'VMAT-2'
                                '[18F]RO-948',          [0.05;0.18],    'tau'
                                '[18F]AV1451',          nan,            'tau'
                                '[18F]ASEM',            nan,            'a7AChRs'
                                'MRI',                  -9,             'MRI'
                                '[18F]XTRA',            nan,            'a4b2AChRs'
                                '[18F]AV45',            nan,            'amyloid-beta'
                                '[18F]Flutemetamol',    nan,            'amyloid-beta'
                                '[11C]PiB',             0.149,          'amyloid-beta'
                                '[11C]OMAR',            -9,             'CB1R'
                                '[11C]CFN',             0.104,          'mu-opioid'
                                '[18F]FTP',             nan,            'Dopamine D3?'
                                '[11C]NNC-112',         nan,            'Dopamine D1'
                                '[18F]FDG',             -9,             'cMRglc'
                                '[11C]RO5378',          nan,            'GABAgamma1Rs'
                                '[11C]Flumazenil',      nan,            'GABA-BZPRs'
                                '[15O]gas',             nan,            'cMRoxy'
                                '[15O]water',           nan,            'CBF'
                                '[11C]MeNER',           nan,            'NorEpiNephTs'
                                '[11C]FLB457',          nan,            'Dopamine D2/D3'
                                '[11C]AZ9369',          nan,            '5HT1BR'
                                '[11C]DASB',            0.048,          'SERT'
                                '[18F]AZAN',            -9,             'a4b2AChRs'
                                '[18F]2FA',             -9,             'a4b2AChRs'
                                '[11C]MDL',             -9,             '5HT2A'
                                '[11C]MP',              0.051,          'DAT'};
%
if ~isempty(which('my_tracerList'));
    t2                          = my_tracerList(0);                      
    im1                         = umo_cstrs(char(t2(:,1)), char(ttt(:,1)), 'im1');
    ttt                         = [t2; ttt(im1<1, :)];                                              end

%
if iscell(i1);  
    out                         = ttt{umo_cstrs(char(ttt(:,1)),i1{1},'im1'), 2};    return;         end
%    
[ts, is]                        = sortrows(char(ttt(:,3)));
ttt                             = ttt(is, :);
if isempty(i1);
    if nargout;                 out                         = ttt;                  return;         end
    disp('.registered radioligands:');
    dispCharArrays(char(ttt(:,1)),1,char(ttt(:,3)));
    disp('< end of the list');
    disp('.ask your IDAE administrator to missing tracers to the list');            return;         end
%
i1(i1==',')                     = ' ';
im1                             = umo_cstrs(lower(char(ttt(:,1))),lower(i1),     'im1');
if any(~im1);                   disp('.foreign radiotracer(s) detected (marked by 0)');
                                disp([i1, char(zeros(size(i1,1),1)+32), int2str(im1(:,1))]);
                                iv2_tracerList([]);                                 return;         end
out                             = char(ttt(im1,1)); 
return;