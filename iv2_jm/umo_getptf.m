function    [out, out2]    = umo_getptf(i1,i2,i3); 

% umo_getptf:       To retreive columns from plain text files, using <<getLseg>>
%       
%       usages:     [ci, ci+1-n]    = umo_getptf('ptffile',m,i)
%                   cs              = umo_getptf('ptffile',m,i:j)
%                   cs              = umo_getptf('ptffile',m,[i,j,k, ...])
%                   all             = umo_getptf('ptffile',0,[]);
%                   lines           = umo_getptf('ptffile',1,[]);
%
%   ptffile     -   file in a plain text format (ver. UMO)
%                   each line is [c1, c2, ... cn]
%                   c1 may have (but only one) *** as a marker
%                   (Enter 0 to m if no *** in c1)
%                   blank lines and lines starting with % will be removed.
%                   2nd and later *** will be ignored
%   m           -   0 to report all lines, 1/2 for lines above/below ***
%                   7 to read .csv files
%   i, etc      -   column number to extract
%                   new -i to report from the back
%   ci, c1+1-n  -   When a single integer (=i) is given, output will be i-th column
%                   (=ci) and i+1-th to the last columns (=ci+1-n, as one string per row)
%   cs          -   When an integer vector is given, the output will be a structure
%                   array (=cs), cs(i).mat and cs(i).cno (for column number)
%   all         -   in a character array (1 x n; including line feeds)
%   lines       -   in a character matrix (n x m)
%
% (cL)2004~18   hkuwaba1@jhmi.edu 

%% when i1 is a matrix (not file name), called from geetLseg
if nargin==2;                   i1(find(abs(i1)<32))        = ' '; 
                                [out, out2]                 = local_getCmat(abs(i1),    i2);
                                return;                                                             end;

margin                          = 3;
if nargin<margin;               helq(mfilename);                                    return;         end;
if isempty(i3);                 
    if i2==1;                   out                         = local_lines(i1); 
    else;                       out                         = local_all(i1);                        end;
                                                                                    return;         end;
if ~any(i3>0);                  [out, out2]                 = local_back(i1,i3);    return;         end;
if i2(1)==7;                    out                         = local_csv(i1,i3);     return;         end;
%
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
[out, out2]                     = local_getCmat(dmat,i3);

return;
%%

function    [out, out2]         = local_getCmat(dmat,i3);
%% putting columns to respective matrices

out                             = [];
out2                            = [];

% adding a blank column in front of dmat:
emat                            = [zeros(size(dmat,1),  1) + 32,    dmat, zeros(size(dmat,1),  1)+32];
% disp(char(emat));

% finding string starting points:
[a, b]                          = find(emat(:,1:end-1)==32 & emat(:,2:end)>32);
%size(emat)
% do not remove (:) to cope with cases where size(emat,1)==1:
sps                             = [a(:), b(:)+1];
[v, iss]                        = sort(sps(:,   1));
sps(:)                          = sps(iss,      :);

% finding string end points:
[c, d]                          = find(emat(:,1:end-1)>32 & emat(:,2:end)==32);
eps                             = [c(:), d(:)];
[v, ise]                        = sort(c);
eps(:)                          = eps(ise,      :);

% 
if size(sps,1)~=size(eps,1);    disp(['     Error:  uneven string start/end pos']);
    disp(char(emat));
    sps
    eps;                                                                            return;         end;
if any(sps(:,1)~=eps(:,1));     disp(['     Error:  row counting error']);          return;         end;

cc                              = zeros(size(sps,1),    1);
clear c;
for i=1:1:size(emat,1);         k                           = find(sps(:,1)==i);
    for j=1:1:length(i3);       m                           = i3(j);
        if length(k)>=m;        c(j).r(i).d                 = char(emat(i,  sps(k(m),2):eps(k(m),2)));
        else;                   c(j).r(i).d                 = '';                           end;    end;
    j                           = length(i3)+1;
    m                           = i3(end) + 1;
    if length(k)>=m;            c(j).r(i).d                 = char(emat(i,  sps(k(m),2):eps(max(k),2)));
    else;                       c(j).r(i).d                 = '';                                   end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
%length(c)
out2                            = str2mat(c(end).r.d);

if length(i3)==1;               out                         = str2mat(c(1).r.d);
else;
    for i=1:1:length(i3);       out(i).mat                  = str2mat(c(i).r.d);
                                out(i).cno                  = i3(i);                                end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;

if isempty(out2);               out2                        = char(zeros(size(emat,1), 3) + 32);    end;
return;
%% 

function    out                 = local_selectLs(cmat,  i2);
%% selecting all/above/below *** 

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
%%

function    cmat                = local_getLines(i1);
%% putting lines of i1 in a matrix, after eliminating comment lines and comments in lines:

cmat                            = [];
fH                              = fopen(i1,                 'r');
if fH<0;                        disp(['.error!  unable to open ',i1]);              return;         end;
a                               = fread(fH,Inf,             'uint8');
fclose(fH);
if ~any(a>32);                  disp(['A empty file ... ',i1]);                     return;         end;
% replacing tabs with spaces:
a(find(a<10))                   = 32;

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
%% 

function        out             = local_all(i1);
%
fH                              = fopen(i1,                 'r');
if fH<0;                        disp(['.error!  unable to open ',i1]);              return;         end;
out                             = char(fread(fH,Inf,        'uint8')');
fclose(fH);
return;
%%

function        out             = local_lines(i1);
%
out                             = [];
fH                              = fopen(i1,                 'r');
if fH<0;                        disp(['.error!  unable to open ',i1]);              return;         end;
ic                              = 1;
while 1;                        c{ic}                       = fgetl(fH);
    if isnumeric(c{ic});                                                            break;          end;
    ic                          = ic + 1;                                                           end;
fclose(fH);
out                             = char(c(1:ic-1));
return;
%%

function        [c, c2]         = local_back(i1,i3);
%% report columns from the back
fH                              = fopen(i1,                 'r');
if fH<0;                        disp(['.error!  unable to open ',i1]);              return;         end;

ic                              = 0;
while 1;                        ic                          = ic + 1;
                                qqq{ic}                     = fgetl(fH);
    if ~ischar(qqq{ic});                                                            break;          end;
    if size(qqq{ic},2)>1;       qqq{ic}                     = qqq{ic}(1,    end:-1:1); 
    else;                       qqq{ic}                     = ' ';                          end;    end;
%
fclose(fH);
qqc                             = char(qqq(1:ic-1));
[c, c2]                         = getLseg(qqc,              abs(i3));
if isstruct(c);
    for i=1:1:numel(c);         c(i).mat(:)                 = c(i).mat(:, end:-1:1);
                                c(i).mat(:)                 = getLseg(c(i).mat, 1);                 end;
else;
    c(:)                        = c(:,      end:-1:1);
    c(:)                        = getLseg(c,    1);                                                 end;
%
c2(:)                           = c2(:,     end:-1:1);
c21                             = char(zeros(size(c2,1), 2));
c21(:,1)                        = '*';
[a, c2(:)]                      = getLseg([c21,c2], 1);
return;
%%

function        c               = local_csv(i1,i3);
%% work on .csv files
fH                              = fopen(i1,                 'r');
if fH<0;                        disp(['.error!  unable to open ',i1]);              return;         end;

ic                              = 0;
while 1;                        qqq                         = fgetl(fH);
    if ~ischar(qqq);                                                                break;          end;
    ic                          = ic + 1;
    qq2                         = [',',qqq,', '];
    ii                          = find(qq2==',');
    jc                          = 0;
    for i=i3(:)';               jc                          = jc + 1;
                                qq3{ic,jc}                  = qq2(ii(i)+1:1:ii(i+1)-1);     end;    end;
%
fclose(fH);
if length(i3)==0;               c                           = char(qq3(:,1));       return;         end;
for i=1:1:size(qq3,2);          c(i).mat                    = char(qq3(:,i));                       end;
return;
%%
    