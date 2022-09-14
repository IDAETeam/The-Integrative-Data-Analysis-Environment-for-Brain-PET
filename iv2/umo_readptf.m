function    [out, out2]    = umo_getptf(i1,i2,i3, o1,v1,o2,v2,o3,v3,o4,v4,o5,v5); 

% umo_getptf:       To retreive columns from plain text files, using <<getLseg>>
%       
%       usages:     [ci, ci+1-n]    = umo_getptf('ptffile',m,i)
%                   cs              = umo_getptf('ptffile',m,i:j)
%                   cs              = umo_getptf('ptffile',m,[i,j,k, ...])
%
%   ptffile     -   file in a plain text format (ver. UMO)
%                   each line is [c1, c2, ... cn]
%                   c1 may have (but only one) *** as a marker
%                   (Enter 0 to m if no *** in c1)
%                   blank lines and lines starting with % will be removed.
%                   2nd and later *** will be ignored
%   m           -   0 to report all lines, 1/2 for lines above/below ***
%   i, etc      -   column number to extract
%   ci, c1+1-n  -   When a single integer (=i) is given, output will be i-th column
%                   (=ci) and i+1-th to the last columns (=ci+1-n, as one string per row)
%   cs          -   When an integer vector is given, the output will be a structure
%                   array (=cs), cs(i).mat and cs(i).cno (for column number)
%
%
% (cL)2004  hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
% -----------------------------------------------------------------------------------------------------;
i2                              = i2(1);
if i2<0 | i2>2;                 disp('m must be 0/1/2');                            return;         end;

out                             = [];
out2                            = [];

% putting lines of i1 in a matrix, after eliminating comment lines and comments in lines:
% also eliminating initial spaces on each line
cmat                            = local_getLines(i1);

% selecting all/above/below *** 
dmat                            = local_selectLs(cmat,  i2);
if isempty(dmat);                                                                   return;         end;

% putting columns to respective matrices:
%[out, out2]                     = local_getCmat(dmat,i3);
[out, out2]                     = local_test(dmat,i3);
return;
%% ====================================================================================================;

function    [out, out2]         = local_test(dmat,i3);
%%

out                             = [];
out2                            = [];

dmat2                           = zeros(size(dmat,1),   size(dmat,2)+2) + 32;
dmat2(:,    2:end-1)            = dmat;

L                               = max(i3) + 1;
s                               = zeros(size(dmat,1),   max(i3)+1);
e                               = zeros(size(dmat,1),   max(i3)+1);
for i=1:1:size(dmat,1);
    p                           = find(dmat2(i, :)==abs('%'));
    if ~isempty(p);             dmat2(i,    p(1):end)       = 32;                                   end;
    k                           = find(dmat2(i,1:end-1)<=32 & dmat2(i,2:end)>32);
    L(:)                        = min(length(k),max(i3) + 1);
    if L;                       s(i,    1:L)                = k(1:L);                               end;
    k                           = find(dmat2(i,1:end-1)>32 & dmat2(i,2:end)<=32);
    L(:)                        = min(length(k),max(i3) + 1);
    if L;                       e(i,    1:L)                = k(1:L);                               end;
                                                                                                    end;
ii                              = find(s(:, 1));
cc                              = zeros(size(ii));
ss                              = s(ii,     :);
ee                              = e(ii,     :);
for i=1:1:length(i3);           k                           = find(ss(:, i3(i)));
    if isempty(k);              out(i).mat                  = [];
                                out(i).cno                  = i3(i);
    else;                       cc(:)                       = ee(:,i3(i)) - (ss(:,i3(i))+1) +1;
                                L(:)                        = max(cc(k));
                                out(i).mat                  = char(zeros(size(ii,1), L) + 32);
                                out(i).cno                  = i3(i);
        for j=1:1:length(k);   
            out(i).mat(k(j),    1:cc(k(j)))     ...
                                = char(dmat2(ii(k(j)),ss(k(j),i3(i))+1:ee(k(j),i3(i))));            end;
                                                                                            end;    end;
ee(:,   end)                    = size(dmat2,2) -1;
k                               = find(ss(:,    end));
if isempty(k);                  out2                        = [];
else;                           out2                        = char(zeros(size(ii,1), L) + 32);
                                cc(:)                       = ee(:,end) - (ss(:,end)+1) +1;
                                L(:)                        = max(cc(k));
    for j=1:1:length(k);
        out2(k(j),  1:cc(k(j))) = char(dmat2(ii(k(j)),ss(k(j),end)+1:ee(k(j),end)));        end;    end;

if size(out2,2)>1;              out2                        = squeeze_smat(out2);                   end;

return;
%%


function    [out, out2]         = local_getCmat(dmat,i3);
%% putting columns to respective matrices 

out                             = [];
out2                            = [];
% finding string starting points:
[a, b]                          = find(dmat(:,1:end-1)==32 & dmat(:,2:end)>32);
%size(dmat)
sps                             = [[[1:1:size(dmat,1)]',  ones(size(dmat,1),1)]; [a(:), b(:)+1]];
[v, is]                         = sort(sps(:,   1));
sps(:)                          = sps(is,       :);

% finding string end points:
[c, d]                          = find(dmat(:,1:end-1)>32 & dmat(:,2:end)==32);
eps                             = [c(:), d(:)];
[v, is2]                        = sort(c);
eps(:)                          = eps(is2,      :);

% 
if size(sps,1)~=size(eps,1);    disp(['     Error:  uneven string start/end pos']); return;         end;
if any(sps(:,1)~=eps(:,1));     disp(['     Error:  row counting error']);          return;         end;

cc                              = zeros(size(sps,1),    1);
clear c;
for i=1:1:size(sps,1);
    j                           = sps(i,    1);
    cc(j,    :)                 = cc(j, 1) + 1;
    if cc(j)<=max(i3);          c(cc(j)).r(j).d             = char(dmat(j,  sps(i,2):eps(i,2)));
    elseif cc(j)==max(i3)+1;    c(cc(j)).r(j).d             = char(dmat(j,  sps(i,2):end-1));       end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
if length(c)>=max(i3)+1;        out20                       = str2mat(c(max(i3)+1).r.d);
                                out2                        = out20(:, find(sum(abs(out20)-32,1))); end;

if length(i3)==1;               out                         = str2mat(c(i3).r.d);
else;
    for i=1:1:length(i3);       j                           = i3(i);
    if j<=length(c);            out(i).mat                  = str2mat(c(j).r.d);
                                out(i).cno                  = j;                            end;    end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
return;
%% ----------------------------------------------------------------------------------------------------;

function    out                 = local_selectLs(cmat,  i2);
%% selecting all/above/below *** ----------------------------------------------------------------------;

out                             = [];
if ~i2;                         out                         = cmat;                 return;         end;
if size(cmat,2)<3;                                                                  return;         end;
ii                              = find(cmat(:,1)==abs('*') & cmat(:,2)==abs('*') & cmat(:,3)==abs('*'));
% reporting above *** (and all if *** is not found);
if i2==1;
    if isempty(ii);             out                         = cmat;
    else;                       out                         = cmat(1:ii(1)-1,   :);                 end;
else;
    if isempty(ii);             disp(['     Error:  *** not found']);               return;         end;
    if ii(end)==size(cmat,1);   disp(['     Error:  None below *** ']);             return;         end;
    out                         = cmat(ii(end)+1:end,       :);                                     end;

return;
%% ----------------------------------------------------------------------------------------------------;


function    cmat                = local_getLines(i1);
%% putting lines of i1 in a matrix, after eliminating comment lines and comments in lines -------------;

cmat                            = [];
fH                              = fopen(i1,                 'r');
if fH<0;                        disp(['     Error:  Unable to open ',i1]);          return;         end;
a                               = fread(fH,Inf,             'uint8');
fclose(fH);

cs                              = find(a>=32);
% lsps  -   line start positions:
% lsps may miss the 1st line but will not miss the last line
lsps                            = find(a(1:end-1)<32 & a(2:end)>=32) + 1;
if isempty(lsps);               lsps                        = cs(1);                                end;
if lsps(1)~=cs(1);              lsps                        = [cs(1);   lsps];                      end;

% leps  -   line end positions:
% leps will not miss the 1st line but may miss the last line
leps                            = find(a(1:end-1)>=32 & a(2:end)<32);
if isempty(leps);               leps                        = cs(end);                              end;
if leps(end)~=cs(end);          leps                        = [leps;    cs(end)];                   end;

if length(lsps)~=length(leps);  disp(i1);
                                disp(['length(lineStartPos)~=length(lineEndPos)']);
                                disp(int2str([length(lsps),length(leps)]));
                                return;
else;
    if any(leps-lsps<0);        disp(i1);
                                disp('lsps<leps in some lines.');                   return;         end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;

% Eliminating lines starting with % (comment lines):
ii                              = find(a(lsps)~=abs('%'));
cmat                            = zeros(length(ii),     max(leps(ii) - lsps(ii)+2)) + 32;
% Erasing after space + % in each line:
for i=1:1:length(ii);

    % eliminating initial spaces:
    ics                         = find(a(lsps(ii(i)):leps(ii(i)))~=32);
    if isempty(ics);            ics                         = 1;                                    end;
    L                           = leps(ii(i)) - lsps(ii(i)) + 2 - ics(1);
    cmat(i,         1:L)        = a(lsps(ii(i))+ics(1)-1:leps(ii(i)))';
    
    sps                         = find(cmat(i,  1:end-1)==32 & cmat(i,  2:end)==abs('%'));
    if ~isempty(sps);           cmat(i, sps(1):end)         = 32;                                   end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
% Eliminating blank lines (spaces only):
ss                              = sum(abs(cmat - 32),2);
if any(~ss);                    cmat                        = cmat(find(ss),    :);                 end;

return;
%% ----------------------------------------------------------------------------------------------------;
