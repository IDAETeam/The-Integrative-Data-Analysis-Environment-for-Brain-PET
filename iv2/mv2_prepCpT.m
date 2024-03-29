function        out             = mv2_prepCpT(i1,fbc);
%% generated by mv2_i2m on 15-Aug-2022 10:30:03
% source file: iv2_prepCpT
% 
if isnumeric(i1);               feval(['local_f',intstr(i1(1),3)],fbc);
else;                           out                         = feval(['local_',i1],fbc);             end;
return;
%
function        out             = local_opt(fbc);
%% listing option strings of this i-file:
% 
% 
out                             = [
'ifc'
'cpt'
'avr'];
return;
%% 
function        out             = local_fun(fbc);
%% returning step descriptions
% 
out{1}                          = 'preparation of plasma & HPLC files';
out{2}                          = 'correct for radio-active metabolites';
out{3}                          = 'plot parent fraction';
out{4}                          = 'plot SUV Ca(t)';
out{5}                          = 'display status of plasma TACs across subjects/scans';
return;
%% 
function        out             = local_fls(fbc);
%% List of input/output files:
% 
global g4iv2;
out                             = {
['pet\',g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).avr,'.ezi']
['pet\cpt_',g4iv2.xxx(fbc(3)).cpt,'.m']
['pet\hplc_',g4iv2.xxx(fbc(3)).cpt,'.m']
['pet\',g4iv2.xxx(fbc(3)).cpt,'.cpt']
['pet\',g4iv2.xxx(fbc(3)).cpt,'_ok.txt']};
%
return;
%% 
function        out             = local_fns(fbc);
%% [in/out, file#@process, process#, file#@ifile, taskID#]
% 
out                             = [
1    1    1    1  105
2    1    1    2  105
2    2    1    3  105
1    1    2    2  105
1    2    2    3  105
2    1    2    4  105
2    2    2    5  105
1    1    3    4  111
1    1    4    4  111
1    1    5    4  111];
return;
%% 
function                        local_f001(fbc);
%% command lines for step 1
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
cv2_getCPT('set',{fullfile('pet',[g4iv2.xxx(fbc(3)).ifc,'_',g4iv2.xxx(fbc(3)).avr,'.ezi']), ...
                                fullfile('pet',['cpt_',g4iv2.xxx(fbc(3)).cpt,'.m']),        ...
                                fullfile('pet',['hplc_',g4iv2.xxx(fbc(3)).cpt,'.m']),fbc(3)}); 
return
%% 
function                        local_f002(fbc);
%% command lines for step 2
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
setCPT(iii{1},ooo{1}, 	'met',struct('fln',iii{2},  'ftype','mat',  'dtype',1), 'apr',ooo{2});  
return
%% 
function                        local_f003(fbc);
%% command lines for step 3
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
plotcpt(iii{1},                 'met','on');                                                   
return
%% 
function                        local_f004(fbc);
%% command lines for step 4
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
plotcpt(iii{1},                 'suv',fbc,'mag','on');                                         
return
%% 
function                        local_f005(fbc);
%% command lines for step 5
% 
global g4iv2;
% 
[iii, ooo]                      = mv2_genfln(fbc,   1);
% command lines start here: 
% 
mv2_getTACstatus(fullfile('pet',[g4iv2.xxx(fbc(3)).cpt,'.cpt']),         fbc);                 
return
%% 
