function   qstr  =  mv2_i2m(i1,i2,i3); 

% To transform IDAE i-files to excutable m-files (ver.iv2)     
%       
%       usage:      mv2_i2m('i-file',[])
%       
%   i-file  -   the name (abc of idx/abc.m) of i-file to work on
%   i2      -   output file name with path and extension
%   output  -   by default (i.e., by entering [] as i2), mv2_abc.m
%
% In IDAE, source i-files are special i-files to reuse in multiple i-files
%   with ### and input parameters (c3). An example line would be:
%       ### iv2_MPcoregLx   m4p=_m2;s2m=m2;pm0=prepMPfsx
%
% Consider this usage when i-files are consisted of source i-file lines alone
%   when one line only:
%     >> mv2_i2m('source i-file','out i-file','name=val combinations');
%
% 
%   e.g., ev2_MPcoregL2.m and ev2_MPcoregL3.m maybe generated without
%   preparing iv2_MPcoregL2.m and iv2_MPcoregL3.m as follows:
%     >> mv2_i2m('iv2_MPcoregLx','iv2_MPcoregL2','m4p=_m2;s2m=m2;pm0=prepMPfsx');
%     >> mv2_i2m('iv2_MPcoregLx','iv2_MPcoregL3','m4p=_m3;s2m=m3;pm0=prepMPfsx');
%   Note that iv2_MPcoregL2.m (abolished) was consisted of above example line
%
% Use cells (1st & 3rd inputs) in multi-lines cases:
%     >> mv2_i2m({'i-file_1',i-file_2', ... }','out i-file',    ...
%           {'name=val_for_i-file_1','name=val_for_i-file_2', ..});
% 
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

%
if nargin==3;                   local_i3(i1,i2,i3);                                 return;         end;

ifl1                            = which(i1);
if isempty(ifl1);               disp(['.problem! not found: ',i1]);                 return;         end;
fH                              = fopen(ifl1,               'r');
ic                              = 0;
while 1;                        inL                         = fgetl(fH);
    if ~ischar(inL);                                                                break;          end;
    if size(inL,2)>=1 && any(inL~=' ') && inL(min(find(inL>=37)))~='%';
                                ic                          = ic + 1;
                                fff{ic}                     = inL;                                  end;
end;
% b                               = fread(fH,Inf,             'char');
fclose(fH);

if isempty(i2);                 
    [idx, inm, iex]             = fileparts(ifl1);
    i2                          = fullfile(fileparts(idx),  'mv2',  ['m',inm(2:end),iex]);          end;
ii                              = 0;
qqq                             = [];
for i=1:1:ic;
    if ~isempty(strfind(fff{i},'###'));
                                mmm                         = local_inc(fff{i});
        for j=1:1:numel(mmm);   ii                          = ii + 1;
                                qqq{ii}                     = mmm{j};                               end;
    else;                       ii                          = ii + 1;
                                qqq{ii}                     = fff{i};                               end;
                                                                                                    end;
clear fff mmm;
charqqq                         = char(qqq);

% checking for o-string:
[ir, ic]                        = find(charqqq(:,   1:end-3)=='*');

if ~isempty(ir);
    rrr                         = char(zeros(length(ir), 3) + 32);
    for i=1:1:length(ir);       rrr(i,  :)                  = charqqq(ir(i),    ic(i)+[1,2,3]);     end;
    % disp(rrr);
    
    % eliminating *xxx that span lines, just in case:
    rrr                         = rrr(sum(rrr>32,2)==3,     :);
    
    cm1                         = umo_cstrs(rrr,[],         'cm1');
    rrr                         = rrr(cm1(:,    2)>0,       :);
    c1                          = umo_getptf(which('iv2_ostrs'),0,1);
    im1                         = umo_cstrs(c1,rrr,         'im1');
    if any(~im1);               disp('Warning');
                                disp(' Un-registered o-stings are noted (marked by 0)');
                                disp([char(rrr),char(zeros(size(rrr))+32),int2str(im1)]);
                                disp(['- if this is an mistake, revise ',i1]);
                                disp('- or add them to iv2_ostrs.m');               return;         end;
else;                           rrr                         = [];                                   end;

if ~isempty(rrr);
    for i=1:1:length(qqq);
        if qqq{i}(1)~='#' && qqq{i}(1)~='$' && any(qqq{i}=='*');
            % replacing ostrings with g4iv2{fbc(1)}.xxx(fbc(3)).
            if qqq{i}(1)=='!';  qqq{i}                      = local_qq2(qqq{i},rrr);
            else;               qqq{i}                      = local_qqq(qqq{i},rrr);        end;    end;
                                                                                            end;    end;                               
qstr                            = char(qqq);
% disp(qstr);
% return;
disp('***')
c1                              = getLseg(qstr,             1);
sss                             = '!#$';
% ss2 = [!/in/out, file#@process, process#, file#@ifile, taskID#]
ss2                             = zeros(size(c1,1),         5);
s0                              = find(qstr(:,1)=='!');
for i=1:1:length(s0);           ss2(s0(i):end,  3)          = i;
                                ss2(s0(i):end,  5)          = abs(c1(s0(i),2));                     end;
for i=1:1:3;                    ss2(c1(:,1)==sss(i),    1)  = i;
    if i>1;                     ss2(c1(:,1)==sss(i),    2)  = str2num(c1(c1(:,1)==sss(i),2:end));   end;
                                                                                                    end;
% ss5 = [in/out, file#@process, process#, file#@ifile, taskID#] for files alone 
ss5                             = ss2(ss2(:,1)>1,           :);
ss5(:,  1)                      = ss5(:,    1) - 1;
c2                              = getLseg(qstr(ss2(:,1)>1,  :),2);
cm1                             = umo_cstrs(c2,[],          'cm1');
ii                              = find(cm1(:,   2)>0);
for i=1:1:length(ii);           ss5(cm1(:,1)==cm1(ii(i),1), 4)  = i;
%                                fff{i}                      = mv2_genfln(deblank(c2(ii(i),:)),[]);  end;
                                fff{i}                      = deblank(c2(ii(i),:));                 end;
% disp(num2str(ss2))
disp(qstr(ss2(:,1)==1,   :));
% disp(char(fff));
% now *xxx are replaced by g4iv2{fnc(1)}.xxx(fbc(3)).xxx in file names: 
% o2(i) is 1 if fff{i} includes *xxx
[ff2, o2]                       = local_fln(fff);
% disp(char(ff2));

disp('all looks ok. now writing in m-file format');
fH                              = local_open(i1,   i2);
if fH<0;                        disp(['.problem! unable to open: ',i2]);            return;         end;
disp(['output: ',i2]);
local_opt(fH,   rrr);
disp('opt ... done');
local_fun(fH,   qstr(ss2(:,1)==1,   :));
disp('fun ... done');
local_fls(fH, ff2,  o2, ss5);
disp('fls ... done');

for i=1:1:max(ss2(:, 3));       local_cls(fH,qstr(ss2(:,1)==0 & ss2(:,3)==i,:), i, ss5);            end;

fclose(fH);
return;
%%

function    out                 = local_qqq(sss,rrr);
%% replacing o=strings with 
str                             = 'g4iv2.xxx(fbc(3)).';
ss2                             = [' ',sss,' '];
s0                              = [find(ss2(1,   1:end-3)=='*'),    size(ss2,2)];
rr2                             = char(zeros(length(s0)-1, 3) + 32);
for i=1:1:length(s0)-1;         rr2(i,  :)                  = ss2(1,    s0(i)+[1:3]);               end;
im1                             = umo_cstrs(rrr, rr2,       'im1');
s0                              = [s0(im1>0),   s0(end)];
ss3                             = ss2(1,    1:s0(1)-1);
for i=1:1:length(s0)-1;         ss3                         = [ss3, str, ss2(s0(i)+1:s0(i+1)-1)];   end;
out                             = ss3(1,    2:end);

return;
%%

function    out                 = local_qq2(sss,rrr);
%% dealing with task declaration lines:
% these lines start with !
str                             = 'g4iv2.xxx(fbc(3)).';
[s1, s2]                        = getLseg(sss,1);
s0                              = find(s2(1,   1:end-3)=='*');
rr2                             = char(zeros(numel(s0), 3) + 32);
for i=1:1:numel(s0);            rr2                         = s2(1, s0(i)+[1:3]);                   end;
im1                             = umo_cstrs(rrr, rr2,       'im1');
if ~any(im1>0);                 out                         = sss;                  return;         end;
s0                              = [s0(im1>0),   size(s2,2)];
if s0(1)==1;                    s3                          = ['[',str];
else;                           s3                          = ['[''',s2(1, 1:s0(1)-1),''',',str];   end;
for i=1:1:numel(s0)-1;   
    s3                          = [s3, s2(1, s0(i)+[1:3]), ',''',s2(1, s0(i)+4:s0(i+1))];           end;
if s0(end-1)+3<s0(end);         out                         = [s1,'  ',s3,''']'];
else;                           out                         = [s1,'  ',s3,']'];                     end;

return;
%%

function                        local_cls(fH,fff,n,ss5);
%% sorting out command lines

str0                            = [ ...
'function                        local_f',intstr(n,3),'(fbc);',    10, ...
'%% command lines for step ',int2str(n),10,'% ',10,];
fwrite(fH,str0,                 'char');

fwrite(fH,['global g4iv2;',10,'% ',10],                     'char');
fwrite(fH,['[iii, ooo]                      = mv2_genfln(fbc,   1);',10],       'char');
fwrite(fH,                      ['% command lines start here: ',10,'% ',10],    'char');
for i=1:1:size(fff,1);          fwrite(fH,[fff(i,:),10],    'char');                                end;

fwrite(fH,  ['return',10,'%% ',10],                         'char');
return;
%%

function                        local_cls_old(fH,fff,n,ff2,ss5);
%% sorting out command lines

str0                            = [ ...
'function                        local_f',intstr(n,3),'(fbc);',    10, ...
'%% command lines for step ',int2str(n),10,'% ',10,];
fwrite(fH,str0,                 'char');

fwrite(fH,['global g4iv2;',10,'% ',10],                     'char');
ii                              = find(ss5(:,1)==1 & ss5(:,3)==n);
if ~isempty(ii);                fwrite(fH,                  ['% input files: ', 10],    'char');    end;
for i=1:1:length(ii);
    fwrite(fH,  ['iii{',int2str(i),'}                          = mv2_genfln(',  ...
                                ff2{ss5(ii(i),4)},',fbc(1,1:4));',10],  'char');                    end;

ii                              = find(ss5(:,1)==2 & ss5(:,3)==n);
if ~isempty(ii);                fwrite(fH,                  ['% output files: ', 10],   'char');    end;
for i=1:1:length(ii);
    fwrite(fH,  ['ooo{',int2str(i),'}                          = mv2_genfln(',  ...
                                ff2{ss5(ii(i),4)},',fbc(1,1:4));',10],  'char');                    end;
fwrite(fH,                      ['% command lines start here: ',10,'% ',10],    'char');
for i=1:1:size(fff,1);          fwrite(fH,[fff(i,:),10],    'char');                                end;

fwrite(fH,  ['return',10,'%% ',10],                         'char');
return;
%%


function                        local_fls(fH, ff2,  o2, ss5);
%% writing file-related info:

str0                            = [ ...
'function        out             = local_fls(fbc);',        10, ...
'%% List of input/output files:',10,'% ',10];
fwrite(fH,str0,                 'char');

if any(o2>0);                   fwrite(fH,                  ['global g4iv2;',10],   'char');        end; 
% writing file names:
fwrite(fH,                      ['out                             = {',10], 'char');
for i=1:1:length(ff2)-1;
    fwrite(fH,                  [ff2{i},    10],            'char');                                end;
fwrite(fH,                      [ff2{end},'};',10,'%',10],  'char');
fwrite(fH,                      ['return;',10,'%% ',10],    'char');

str0                            = [ ...
'function        out             = local_fns(fbc);',        10, ...
'%% [in/out, file#@process, process#, file#@ifile, taskID#]',10,'% ',10];
fwrite(fH,str0,                 'char');
% writing [in/out, file#@process, process#, file#@ifile, taskID#]
fwrite(fH,                      ['out                             = [',10], 'char');
for i=1:1:size(ss5,1)-1;
    fwrite(fH,                  [int2str(ss5(i,:)),10],     'char');                                end;
fwrite(fH,                      [int2str(ss5(end,:)),'];',10],  'char'); 
fwrite(fH,                      ['return;',10,'%% ',10],    'char');
return;
%%

function                        local_fun(fH,fff);
%% copying step category initials and descriptions

[c1, c2]                        = getLseg(fff,              1);
% c2                              = local_fixc2(c2);

str0                            = [ ...
'function        out             = local_fun(fbc);',        10, ...
'%% returning step descriptions',10,'% ',10];
fwrite(fH,str0,                 'char');
c32                             = char(zeros(1,32) + 32);
if any(c2(:,1)=='[');           fwrite(fH,  ['global g4iv2;',10],   'char');                        end;
for i=1:1:size(c2,1);
    c32(:)                      = ' ';
    c32(1,  1:length(int2str(i))+5)                         = ['out{',int2str(i),'}'];
    
    if c2(i,1)=='[';
        fwrite(fH,              [c32,'= ',deblank(c2(i,  :)),';',10],       'char');
    else;
        fwrite(fH,              [c32,'= ''',deblank(c2(i,:)),''';',10],     'char');                end;
end;
fwrite(fH,                      ['return;',10,'%% ',10],    'char');

return;
%%

function    out                 = local_fixc2(c2);
%%
out                             = c2;
% for i=1:1:size(c2,1);
%     q1                          = strfind(c2(i,:),'g4iv2');
%     if numel(q);
%                 


return;
%%

function                        local_opt(fH, rrr);
%% listing option strings of this i-file:

str0                            = [ ...
'function        out             = local_opt(fbc);',        10, ...
'%% listing option strings of this i-file:',10,'% ',10,'% ',10];
fwrite(fH,str0,                 'char');

if isempty(rrr);     
    fwrite(fH,                  ['out                             = [];',10],       'char'); 
                                                                                    return;         end;
fwrite(fH,                      ['out                             = [',10], 'char');
for i=1:1:size(rrr,1)-1;
    fwrite(fH,                  ['''',rrr(i,:),'''',10],    'char');                                end;
fwrite(fH,                      ['''',rrr(end,:),'''];',10],'char');  
fwrite(fH,                      ['return;',10,'%% ',10],    'char');
return;
%%

function    fH                  = local_open(i1,i2);
%%

fH                              = fopen(i2,                 'w');
if fH<0;                        disp(['Unable to open ... ',i2]);                   return;         end;

[idx, inm]                      = fileparts(i2);

str0                            = [ ...
'function        out             = ',inm,'(i1,fbc);',                               10, ...
'%% generated by ',mfilename,' on ',datestr(now),10,'% source file: ',i1,10,'% ',   10, ...
'if isnumeric(i1);               feval([''local_f'',intstr(i1(1),3)],fbc);',        10, ...
'else;                           out                         ',                         ...
'= feval([''local_'',i1],fbc);             end;',10,'return;',10,'%',10];
fwrite(fH,str0,                 'char');

return;
%%

function    [out, o2]           = local_fln(fff);
%%
o2                              = zeros(length(fff),        1);
str                             = 'g4iv2.xxx(fbc(3)).';
for i=1:1:length(fff);          % disp(fff{i});
    if any(fff{i}=='*');
        o2(i,   :)              = 1;
        if size(fff{i},2)==4;   out{i}                      = [str,fff{i}(2:4)];
        else
        vvv                     = [' ',fff{i},' '];
        ss                      = find(vvv=='*');
        for j=1:1:length(ss);   yyy{j}                      = vvv(ss(j)+1:ss(j)+3);
                                vvv(ss(j):ss(j)+3)          = ' ?? ';                               end;

        s0                      = find(vvv(1:end-1)==' ' & vvv(2:end)~=' ')+1;
        s1                      = find(vvv(1:end-1)~=' ' & vvv(2:end)==' ');
        ii                      = 0;
        xxx                     = '';
        for j=1:1:length(s0);
            if vvv(s0(j))=='?'; ii                          = ii + 1;
                                xxx                         = [xxx, str,yyy{ii},','];
            else;               xxx                         = [xxx, '''',vvv(s0(j):s1(j)),''','];   end;
        end;
        xxx(end)                = ']';
        out{i}                  = ['[', xxx];
        end;
    else; 
        out{i}                  = ['''', fff{i},    ''''];                                          end;
end;
return;
%%

function    ggg                 = local_inc(a);
%%
ggg                             = [];
[c1, c2]                        = getLseg(a,2);

c2(c2==';')                     = ' ';
c2x                             = [' ',c2,  ' '];
s1                              = find(c2x(1:end-1)==' ' & c2x(2:end)~=' ')+1;
s2                              = find(c2x(1:end-1)~=' ' & c2x(2:end)==' ');

fH                              = fopen(which(c1),          'r');
if fH<0;                        disp(['Unable to locate ... ',c1]);
                                ggg                         = []; 
                                disp([' ... aborting ',mfilename]);                 return;         end;
ic                              = 0;
while 1;                        rrr                         = fgetl(fH);
    if ~ischar(rrr);                                                                break;          end;
    if size(rrr,2)>1 && rrr(min(find(rrr>=37)))~='%';
                                ic                          = ic + 1;
                                ggg{ic}                     = rrr;                                  end;
end;
fclose(fH);
ppp                             = [];
for i=1:1:ic;
    for j=1:1:length(s1);       ss                          = strfind(ggg{i},['*',c2x(s1(j):s1(j)+2)])+1;
        if ~isempty(ss);        sss                         = [' ',ggg{i},' '];
                                ss                          = [ss, size(sss,2)];
                                s0                          = sss(1:ss(1)-1);
            for k=1:1:length(ss)-1;   
                s0              = [s0,c2x(s1(j)+4:s2(j)),sss(ss(k)+4:ss(k+1)-1)];                   end;
                                ggg{i}                      = s0(1, 2:end);                         end;
    end;
end;
return;
%%

function    out                 = local_repl(a,        rstrs,rvals);
%%
% istrs     -   replacing strings:
%   istrs must be n by 3 matrix
disp('replacing as follows (column 1 by column 2)');
disp([rstrs,char(zeros(size(rstrs))+32),rvals]);
out                             = [];
if length(a)<4;                 out                         = a;                    return;         end;
b                               = zeros(length(a)-3,        4);
b(:)                            = [a(1:end-3), a(2:end-2), a(3:end-1), a(4:end)];
im1                             = umo_cstrs(b,rstrs,        'im1');
if ~any(im1(:));                out                         = a;                    return;         end;

qqq                             = ones(size(im1));
rrr{1}                          = deblank(rvals(1,   :));
for i=2:1:size(im1,1);          qqq(i,  :)                  = i;
                                rrr{i}                      = deblank(rvals(i,  :));                end;
[v, is]                         =  sort(im1(im1>0));
q1                              = qqq(im1>0);
q2                              = q1(is);

% i-files may not start with *xxx:
s1                              = [v-2;     length(a)];
s0                              = [1;   v+3];
for i=1:1:length(s1)-1;         out                         = [out; a(s0(i):s1(i)); rrr{q2(i)}'];   end;
if s1(end)>s0(end);             out                         = [out; a(s0(end):s1(end))];            end;

return;
%%

function                        local_i3(i1,i2,i3);
%%
if iscell(i1);
    L2w                         = '';
    for i=1:1:numel(i1);        [idx, inm]                  = fileparts(which(i1{i}));
                                L2w                         = [L2w,'### ',inm,' ',i3{i},10];     	end;
else;                           [idx, inm]                  = fileparts(which(i1));
                                L2w                         = ['### ',inm,'   ',i3];                end;
%
[jdx, jnm]                      = fileparts(i2);
tfl                             = tmpfln([],    'm');
write2ptf(tfl, L2w);
mv2_i2m(tfl, fullfile(fileparts(idx),  'mv2', ['m',jnm(1,2:end),'.m']));            
delete(tfl);
if exist(which(jnm),'file');
    disp('.suggestion! no longer needed. delete it after confirming');
    disp(['>> delete ',which(jnm)]);                                                                end;
return;
%%
