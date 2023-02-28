function    um_AttEdit(i1,i2); 

% To replace values of attributes in UMO format files
%       
%       usage:      um_AttEdit('full/path/umo.ezx',i2)
%
%   i2{i,1}     -   variable name of i-th attribute
%   i2{i,2}     -   values (in matrix) 
% 
% (cL)2015    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

if ~iscell(i2);                 disp('.second input has to be a cell array');       return;         end;
if size(i2,2)~=2;               disp('.second input must be n by 2');               return;         end;

c                               = zeros(size(i2,1),         1);
for i=1:1:size(i2,1);
    v                           = gei(i1,                   i2{i,1});
    c(i,    :)                  = length(v);                                                        end;
if any(~c);
    disp('.the following attributes (requested) are not present .. ');
    disp([char(zeros(size(c))+32),char(i2(c==0))]);
    q                           = 'x';
    while q(1)~='a';
        q                       = input('.continue (=ret) or abort (=a)? ','s');
        if isempty(q);                                                              break;          end;
        if q(1)=='a';                                                               return;         end;
    end;
end;
%
[idx, inm, iex]                 = fileparts(i1);
tfl                             = tmpfln([],                iex(1,2:end));
disp(['.copying the input file to ',tfl])
copyfile(i1,    tfl);
%
d                               = gei(i1,                   'dataInfo');
di                              = zeros(size(d,1),          1);
df                              = zeros(size(d,1),          10);
%
si                              = struct('h2s',32,'c',mfilename,'p',i1,'cp','a');
ss                              = '[fH, ii]                 = um_save(i1,[],si,[],';
for i=1:1:size(i2,1);
    ss                          = [ss,'i2{',int2str(i),',1},i2{',int2str(i),',2},'];                end;
%
ss                              = [ss(1, 1:end-1),');'];
eval(ss);
% 
if fH<0;                        disp(['.unable to create .. ',i1]);
                                disp('.no change to the file');
                                copyfile(tfl,               i1);                    return;         end;
%
for i=1:1:size(d,1);
    [di(i, :), df(i, :)]        = um_save(fH, ged(tfl,i),   si.h2s, []);                            end;

status                          = um_save(fH, ii, di, df);
disp('.done! (requested attributes revised)');
disp([' output .. ',i1]);
delete(tfl);
return;