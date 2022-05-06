function    mv2_p2(i1); 

% To carry on callback jobs of 'process' GUIs of L2W (IDAE.iv2)    
%       
%       usage:      mv2_p2(1)
%       
%   This code is accessible only from 'process' GUIs of L2W 
%   (i.e., both gcf(=L2W) and gco (=GUI) are used.)
% 
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

% cwUD  = [L1W figure handle, subject#, iBase#]
cwUD                            = get(gcf,                  'userData');
L2W                             = double(gcf);
% coUD  = [fbc, g4iv2{fNo}.ppp], but fbc(2) is a dummy
fbc                             = get(gco,                  'userData');
fbc(1,  2)                      = cwUD(2);
% 'd' (=decesion points) and 'm' (=common to all) are treate as i here);
if fbc(4)==abs('d') | fbc(4)==abs('m');    
                                fbc(4)                      = abs('i');                             end;
%
% ['local_',char(fbc(4))]
if ~isempty(which(['local_',char(fbc(4))]));
                                feval(['local_',char(fbc(4))],fbc,L2W,get(gco,'String'));           end;
return;
%%

function                        local_a(fbc,L2W,g0os);
%% automated processes:
global g4iv2;
% fbc is now with correct subj# (from L2W):
mv2_fck(fbc(1),fbc(2),          fbc(end));
% if fbc(3)==size(g4iv2{fbc(1)}.irq,2);   fbc(3)                  = 1;                                end;
% cc                              = [fbc(3),  size(g4iv2{fbc(1)}.irq,2)];
% is                              = fbc(end);
% process # in the ipack:
i                               = fbc(end);
q1                              = g4iv2.ick{fbc(2)}(i,:)==g4iv2.irq(i,:);
q2                              = g4iv2.ock{fbc(2)}(i,:)<g4iv2.orq(i,:)     ...
                                                            & g4iv2.orq(i,:)>0;
q3                              = [1:1:size(q1,2)-1,1];
fbx                             = fbc;
fbx(1,  3)                      = q3(fbc(3));
% g4iv2.ifl{fbc(6)}
% fbc(5)
if q1(fbc(3)).*q2(fbc(3));      feval(g4iv2.ifl{fbc(6)},fbc(5),fbx);
else;                           local_disp(fbc);                                                    end;
% q1                              = 
% q1                              = g4iv2{fbc(1)}.ick{fbc(2)}(fbc(3), cc)==g4iv2{fbc(1)}.irq(fbc(3), cc);
% q2                              = abs(g4iv2{fbc(1)}.ock{fbc(2)}(is, cc) - g4iv2{fbc(1)}.orq(is, cc))>0;
% %
% % disp(int2str([q1,q2]));
% % return;
% if min(q1)>0 && sum(q2)>0;      feval(g4iv2{fbc(1)}.ifl{fbc(6)},fbc(5),fbc);
% else;                           local_disp(fbc);                                                    end;
% updating input/output file status:
mv2_fck(fbc(1),fbc(2),          fbc(end));
% updating status flag:
mv2_a0(L2W)
return;
%%

function                        local_s(fbc,L2W,g0os);
%% automated processes:
global g4iv2;
% fbc is now with correct subj# (from L2W):
mv2_fck(fbc(1),fbc(2),          fbc(end));
if fbc(3)==size(g4iv2.irq,2);   fbc(3)                  = 1;                                end;
cc                              = [fbc(3),  size(g4iv2.irq,2)];
is                              = fbc(end);
q1                              = abs(g4iv2.ick{fbc(2)}(is, cc) - g4iv2.irq(is, cc))==0;
q2                              = abs(g4iv2.ock{fbc(2)}(is, cc) - g4iv2.orq(is, cc))>0;
if min(q1)>0;                   feval(g4iv2.ifl{fbc(6)},fbc(5),fbc);                        end;
local_disp(fbc);
% updating input/output file status:
mv2_fck(fbc(1),fbc(2),          fbc(end));
% updating status flag:
mv2_a0(L2W)
return;
%%

function                        local_c(fbc,L2W,g0os)
%%
global g4iv2;
% adjusting fbc(3) 
if fbc(3)>size(g4iv2.yyy.cMat,1);
                                fbc(3)                      = 1;                                    end;
if g0os=='-';                   local_disp(fbc);                                    return;         end;
feval(g4iv2.ifl{fbc(6)},fbc(5),fbc);
return;
%%

function                        local_o(fbc,L2W,g0os)
%%
% disp('yes')
global g4iv2;
% adjusting fbc(3) 
if fbc(3)>size(g4iv2.yyy.cMat,1);
                                fbc(3)                      = 1;                                    end;
if g0os~='r';                   local_disp(fbc);                                    return;         end;
feval(g4iv2.ifl{fbc(6)},fbc(5),fbc);
return;
%%

function                        local_i(fbc,L2W,g0os)
%% interactive process:
% fbc
global g4iv2;
% adjusting fbc(3) 
if fbc(3)>size(g4iv2.yyy.cMat,1);
                                fbc(3)                      = 1;                                    end;
if g0os=='-';                   local_disp(fbc);                                    return;         end;
feval(g4iv2.ifl{fbc(6)},fbc(5),fbc);
return;
%%

function                        local_r(fbc,L2W,g0os)
%% reparative process:
global g4iv2;
% adjusting fbc(3) 
if fbc(3)>size(g4iv2.yyy.cMat,1);
                                fbc(3)                      = 1;                                    end;
if g0os~='c';                   local_disp(fbc);                                    return;         end;
feval(g4iv2.ifl{fbc(6)},fbc(5),fbc);
return;
%%

function                        local_disp(fbc);
%% displaying input / output files:
global g4iv2;
if fbc(3)>size(g4iv2.yyy.cMat,1);
                                fbc(3)                      = 1;                                    end;
ppp                             = feval(g4iv2.ifl{fbc(6)},'fun',fbc);
disp(['.leaving process of: ',ppp{fbc(5)}]);
[iii, ooo]                      = mv2_genfln(fbc,           1);
disp('Input files (2=present; 0=not present) ... ');
for i=1:1:numel(iii);           disp([int2str(exist(iii{i},'file')),': ',iii{i}]);                  end;
disp('Output files (2=present; 0=not present) ... ');
for i=1:1:numel(ooo);           disp([int2str(exist(ooo{i},'file')),': ',ooo{i}]);                  end;
return;
%%

