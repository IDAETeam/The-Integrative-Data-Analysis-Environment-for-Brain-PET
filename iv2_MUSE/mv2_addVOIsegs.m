function    out                 = mv2_addVOIsegs(i1,i2, varargin); 

% To ... 
%       
%       usage 1:    out         = mv2_addVOIsegs([],[])
% 
%   To return rules on VOI routines (Need to be synchronized with iv1_prepMP.m)
%   out     -   a cell array of cells that list (element = out{i}{j}):
%           out{1}  -   VOI set strings of VOI routines 
%           out{2}  -   VOI routine strings 
%           out{3}  -   iBase names of VOI routines 
%           out{4}  -   names of sbdivision routines
%           out{5}  -   iBase names of sbdivision routines
% 
%       usage 2:    out         = mv2_addVOIsegs(input1,NM2)
%   
%   To check for required VOIs for sbdivision routines
%   input1  -   any of usage_1_out{4}
%   NM2     -   n by (1 + m) numeric matrix. input2(:,1) = VOIIDNos
%               input2(:,2:end) = required VOI competion status (0/1/2/3)
%   out     -   [] if all required VOIs are marked by 2 or 3 
%               Otherwise, it returns required VOIs in VOIIDNos.
%
%       usage 3:    out         = mv2_addVOIsegs(input1,'full/path/iPack.m')
%
%   input1  -   any of usage_1_out{4}
%   input2  -   iPack with full path
%   out     -   0 if input1 is not found in the iPack (=not activated)
%               >0 if input1 is active in the iPack
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if ~isempty(i1) && ~isempty(which(['local_',lower(i1)]));
    out                         = feval(['local_',lower(i1)],i2);                   return;         end;

% available VOI routines. Always check mv2_info.m
% current non-subdivision VOI sets (={1})& their labels (={2})
out{1}                          = {'fslu',  'f81u', 'f45u', 'pet'};
out{2}                          = {'FSL',   'FS81', 'FS45', 'onPET'};
% current iBases for regular VOI rourines:
out{3}                          = {'IDAE4FSL', 'IDAE4FS',  'IDAE4FS',  'IDAE4pVOIs'};
% current sbdivision routines:
out{4}                          = {'P10',       'C14',      'I08'};
% of these, 
out{5}                          = {'IDAE4vs',   'IDAE4Cng', 'IDAE4Ins'};

return;

function    out                 = local_p10(i2);
%%
% dfine required VOIs for individual subdividion routines:
out                             = [81000; 82000];
if isempty(i2);                                                                     return;         end;
if isnumeric(i2);               vi                          = consolidVOINos(i2(:,1),   out);
    if any(~vi(:,2));                                                               return;         end;
    if ~any(sum(i2(vi(:,2), 2:end)>1,1)==size(vi,1));                               return;         end;
    out                         = []; 
    return;
else;
    sss                         = mv2_addVOIsegs([],[]);
    im1                         = umo_cstrs(char(sss{4}),'P10', 'im1');
    % reading iBase from the iPack (=i2):
    c1                          = umo_getptf(i2,            1,1);
    out                         = umo_cstrs(c1,sss{5}{im1}, 'im1');                                 end;

return;
%%
% Neet generate c14, i08 and so on as they become available

function    out                 = local_onpet(i2);
%%
sss                             = mv2_addVOIsegs([],[]);
im1                             = umo_cstrs(char(sss{2}),'onPET',   'im1');
if isnumeric(i2);               out                         = i2(i2(:,im1+1)>0,     im1+1);
else;
    % reading iBase from the iPack (=i2):
    c1                          = umo_getptf(i2,            1,1);
    out                         = umo_cstrs(c1,sss{5}{im1}, 'im1');                                 end;
return;
