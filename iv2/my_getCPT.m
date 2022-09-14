
function    my_getCPT(i1,iii,ooo); 

% to read .xls* files of plasma TACs and HPLC      
%       
%       usage:      my_getCPT('full/directory/pet.file','subjID_scanID',cpt)
%
%   cpt     -   plasma file string (xxx of full/path/whatever_xxx.cpt)
%       
%
%   my_getCPT('fun',
%
% (cL)2010    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
feval(['local_',lower(i1)], iii,ooo);
return;
%%

function                        local_cpt(iii,ooo)
%%
disp(['.reading: ',iii{1}]);
strs                            = [ 'Diagnosis                          ';
                                    'Scan Date                          ';
                                    'Scan #                             ';
                                    'Elapsed Time Since Inj.(hh:mm:ss)  ';
                                    'Conversion Factor(cpm/mCi)         ';
                                    'Sample Volume(cc)                  ';
                                    'Tracer/Study                       ';
                                    'Half-Life(minutes)                 ';
                                    'Gamma Counter File                 ';
                                    'Number of plasma samples           '];
sout                            = {'Diagnosis','scanDate','ScanNumber','scan2gamma',    ...
                                    'convFactor','SampleVolume','RadioTracer',          ...
                                    'HalfLife','GammaFile','NoOfSamples'};
sout2                           = {'Diagnosis','scanDate','ScanNumber','scan2gamma(min)',   ...
                                    'convFactor(cpm/mCi)','SampleVolume(ml)','RadioTracer', ...
                                    'HalfLife(min)','GammaFile','NoOfSamples'};

cnos                            = [5,5,5,6,5,5,5,5,5,5];

[a, b, c]                       = xlsread(iii{1},           1);
im1                             = umo_cstrs(char(c{1:13,2}),strs,   'im1');
for i=1:1:size(strs,1);
    if im1(i);                  eval(['sss.',sout{i},'      = c{im1(i),cnos(i)};']);
    else;                       eval(['sss.',sout{i},'      = [];']);                       end;    end;

for i=1:1:size(c,1);
    if isnumeric(c{i,2});       c{i,    2}                  = ' ';                          end;    end;
im2                             = umo_cstrs(char(c{:,2}),'Start',   'im1');


for i=1:1:size(c,2);
    if isnumeric(c{im2(1),i});  c{im2(1),i}                 = ' ';                          end;    end;

im3                             = umo_cstrs(char(c{im2(1),:}),['Time  ';'Activi'],'im1');
cpt                             = zeros(sss.NoOfSamples,    2);
for i=1:1:sss.NoOfSamples;      
    cpt(i,  :)                  = [c{i+im2(1), im3(1)},     c{i+im2(1), im3(2)}];                   end;
% removing negative values:
cpt(cpt<0)                      = 0;
% converting from microCi/ml to nCi/ml:
cpt(:,  2)                      = cpt(:,    2).*1000;

save(ooo{1}, 'cpt', '-ascii');
disp('.done! (plasma TACs [time, nCi/mL] in ASCII)');
disp([' output: ',ooo{1}]);
return;
%%

function                        local_hplc(iii,ooo);
%%
winopen(iii{1});
drawnow;
disp('.manually input HPLC data as follows:');
disp(' first, copy and paste time & decay corrected %parent from the pdf file as follows:');
disp([' hplc = [ ',10,' 0 95',10,'  copy and paste ',10,' ]']);
disp(' then, copy and paste the following 2 lines');
disp([' ofl = ''',ooo{1},''';']);
disp(' my_getCPT(''hplc_2'', {hplc},{ofl});');
return;
%%

function                        local_hplc_2(iii,ooo);
%%
disp('.converting HPLC data from parent % to metabolite % ..');
met                             = iii{1};
met(:,2)                        = 100 - met(:,2);
save(ooo{1}, 'met', '-ascii');
disp('.done! (HPLC data in ascii file)');
disp([' output: ',ooo{1}]);
cv2_getCPT('check', 2);
return;
%%

function                        local_disp(i1,ss,str,i2);
%%
global g4iv2;
w                               = gei(g4iv2.yyy.ifl,        'subjSpec');
if w(1)=='a';                   local_bab(i1,ss,str,i2);                            return;         end;
i1s2                            = [i1(ss(2)+1),'*',i1(ss(2)+2),'*',i1(ss(2)+3:ss(3)-1)];
disp(['.checking .. ',fullfile('Q:', lower(i1(ss(2)+1)), i1s2)]);
idx                             = fullfile('Q:', lower(i1(ss(2)+1)));
pdx                             = dir(fullfile(idx,         i1s2));
for i=1:1:numel(pdx);
    fff                         = dir(fullfile(idx,pdx(i).name,i1(ss(3)+1:ss(5)-1),'plasma','*'));
    % disp(fullfile(idx,pdx(1).name,i1(ss(3)+1:ss(5)-1),'plasma'));
    disp('.potential plasma files .. ');
    for j=1:1:numel(fff);
        if fff(j).name(1)~='.';
            disp(['winopen ''',fullfile(idx,pdx(i).name,i1(ss(3)+1:ss(5)-1),    ...
                                'plasma',fff(j).name),'''']);                               end;    end;
    disp('.check the plasma file and save it as ... ');
    disp(fullfile(i1,           [i2,'_',str,'.xlsx']));                                             end;
return;
%%

function                        local_bab(i1,ss,str,i2);
%%
s1                              = find(i1==filesep);
s2                              = i1(s1(2)+3:s1(3)-1);
s3                              = find(s2~='0',1);
pdx                             = dir(fullfile('R:',['B*_',s2(1, s3(1):end)]));
disp(['.checking .. ',fullfile('R:',['B*_',s2(1, s3(1):end)])]);
for i=1:1:numel(pdx);
    fff                         = dir(fullfile('R:',pdx(i).name,i1(ss(3)+1:ss(5)-1),'plasma','*'));
    disp(fullfile('R:',pdx(1).name,i1(ss(3)+1:ss(5)-1),'plasma'));
    for j=1:1:numel(fff);
        if fff(j).name(1)~='.';
            disp(['winopen ''',fullfile('R:',pdx(i).name,i1(ss(3)+1:ss(5)-1),   ...
                                'plasma',fff(j).name),'''']);                               end;    end;
    disp('.check the plasma file and save it as ... ');
    disp(fullfile(i1,           [i2,'_',str,'.xlsx']));                                             end;
return;
%%    
