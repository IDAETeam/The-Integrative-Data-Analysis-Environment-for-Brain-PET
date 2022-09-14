function    [out1, out2]        = makefarrays(i1,i2, varargin); 

% makefarrays:      To make arrays of file names for IDAE
%       
%       usage:      [out1, ..., outn, finfo] = makefarrays(dstr,fstr)
%  
%   dstr        -   directory string (pet/res/mri/ezr) to look for files
%   fstr        -   file name string (to look for *fstr)
%
%   outi        -   matrix of full file names for scan condition i when cNo==0.
%                   one output if cNo~=0.
%   finfo       -   No. of subjects by 1. 1 if found, 0 otherwise.
%
% Options:      
%   'fNo',val   -   Figure number of the IDAE overview chart.
%   'bNo',val   -   Subject button No in the chart. 0 for all subjects (=default)
%   'cNo',val   -   Experimental condition No. 0 to include all conditions (=default)
%   'fbc',val   -   To enter fNo, bNo, and cNo in one (over-write them)
%   'all','off' -   Report missing file(s) as empty
%                   default: 'on'  - Not to report if one or more files are missing
%                   in any condition.

margin                          = 2;
if nargin<margin;               help makefarrays;                                   return;         end;
% -----------------------------------------------------------------------------------------------------;

 
fnoval                          = 1;
bnoval                          = 0;
cnoval                          = 0;
fbcval                          = [1,0,1];
allval                          = 'on';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
% -----------------------------------------------------------------------------------------------------;
if fbcflg;                      fnoval                      = fbcval(1);
                                bnoval                      = fbcval(2);
                                cnoval                      = fbcval(3);                       
else;                           fnoval                      = fnoval(1);
                                bnoval                      = bnoval(1);
                                cnoval                      = cnoval(1);                            end;
allflg                          = strncmp(lower(allval),'on',2);
% -----------------------------------------------------------------------------------------------------;

% when called from iv2:
global g4iv2;
if ~isempty(g4iv2) && numel(g4iv2)>=fbcval(1) && ~isempty(g4iv2);     
    [out1, out2]                = local_iv2(i1,i2,          fbcval);                return;         end;
%

fns                             = int2str(fnoval(1));
eval(['global g4b2idae',fns]);
eval(['g4b2idae                 = g4b2idae',fns,';']);
if ~bnoval;                     bnoval                      = [1:1:size(g4b2idae.bHs,1)]';          end;
if ~cnoval;                     cMat                        = gei(g4b2idae.ifl,     'condFlags');
                                cnoval                      = [1:1:size(cMat,1)]';                  end;
% -----------------------------------------------------------------------------------------------------;

cL                              = length(cnoval);
bL                              = length(bnoval);
fLs                             = zeros(bL,     cL);
if bL==1;
% -----------------------------------------------------------------------------------------------------;
    for j=1:1:cL;
        ifl(j).name             = findfls4idae(i1,i2,       [fnoval,bnoval(1),cnoval(j)]);
        fLs(1,  j)              = length(ifl(j).name);                                              end;
    % -------------------------------------------------------------------------------------------------;


    if ~max(fLs(:));            out1                        = [];
    % none was foound:
                                out2                        = fLs;                  return;         end;
    % -------------------------------------------------------------------------------------------------; 

    out1                        = char(zeros(cL,    max(fLs(1,:))) + 32);
    out2                        = fLs(1,    :)';

    for j=1:1:cL;   if fLs(j);  out1(j,  1:fLs(1,j))        = ifl(j).name;                  end;    end;
    if ~allflg;                 out1                        = out1(find(out2),    :);               end;
    out2(find(out2))            = 1;
    return;                                                                                         end;
% -----------------------------------------------------------------------------------------------------;

for j=1:1:cL;                   ifl                         = [];
% -----------------------------------------------------------------------------------------------------;
    for i=1:1:bL;       
    
        ifl(i).name             = findfls4idae(i1,i2,       [fnoval,bnoval(i),cnoval(j)]);
        fLs(i,  j)              = length(ifl(i).name);                                              end;

    out                         = char(zeros(bL,    max(fLs(:,j))) + 32);
    ii                          = find(fLs(:,   j));
    for i=1:1:bL; if fLs(i,j);  out(i,  1:fLs(i,j))         = ifl(i).name;                  end;    end;
    if allflg;                  eval(['out',int2str(j),'    = out;']); 
    else;                       eval(['out',int2str(j),'    = out(ii,   :);']);                     end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
fLs(find(fLs))                  = 1;
eval(['out',int2str(cL+1),'     = fLs;']);

return; 
%%

function    [o1, o2]            = local_iv2(i1,i2,          fbc);
%%
global g4iv2 g4dxs;
if ~isempty(i2);                i1                          = fullfile(i1,  i2);                    end;
%
o2                              = zeros(size(g4dxs.mri,1),  1);
for i=1:1:size(o2,1);
    [ooo{i}, o2(i,:)]           = mv2_genfln(i1,            [fbc(1),i,fbc(3)]);                     end;
%
o1                              = char(ooo);
return;
%%
