function    my_getCPT(i1,i2, i3); 

% To convert plasma files for Wash U 
%       
%       usage:      my_getCPT()
%       
% 
% (cL)2022    hkuwaba1@jhmi.edu 

feval(['local_',lower(i1)], i2, i3);
return;
%%

function                        local_cpt(pet_is,ofl)
%% 
% pet_is = full/path/of/plasma.file
% ofl    = full/path/of/scan_ID_cpt_*cpt.m
if isempty(pet_is);                                                                 return;         end;

[odx, onm]                      = fileparts(ofl);
global g4iv2;
xfls                            = dir(fullfile(odx, [onm(1, 1:end-size(g4iv2.xxx(1).cpt,2)-1),'.xls*']));
if numel(xfls)>1;
    disp(['> problem! more than one ',onm(1, 1:end-size(g4iv2.xxx(1).cpt,2)-1),'.xls*']);
    disp(['  delete excessive one(s)']);
    disp(['  folder: ',odx]);
elseif numel(xfls)<1
    [fname, path_x]             = uigetfile(fullfile(fileparts(pet_is),'plasma','*.xls*'));
    if ~ischar(fname);                                                              return;         end;

    [fdx, fnm, fex]             = fileparts(fname);
    xfln                        = fullfile(odx, [onm(1, 1:end-size(g4iv2.xxx(1).cpt,2)-1),fex]);
    disp(['> copying plasma source file..']);
    disp([' input: ',fullfile(path_x,fname)]);
    disp(['    to: ',xfln]);
    copyfile(fullfile(path_x,fname), xfln);
else;
    xfln                        = fullfile(xfls(1).folder, xfls(1).name);                           end;
%
if exist(xfln,'file')~=2;                                                           return;         end;
disp(['> reading: ',xfln]);
%
qqq                             = readcell(xfln, 'Sheet',1);
ccc                             = zeros(size(qqq));
for i=1:1:size(qqq,1);
    for j=1:1:size(qqq,2);      
        if ~ismissing(qqq(i,j));
            ccc(i, j)           = double(ischar(qqq{i,j})) + double(isnumeric(qqq{i,j})).*2; 
                                                                                    end;    end;    end;
%
strs                            = [    'Start(min)      '
    'Start(sec)      '
    'Stop(min)       '
    'Stop(sec)       '
    'Bcq             '
    'Time            '
    'Activity(mCi/cc)'];
%
for k=find(sum(ccc==1,2)==7)';
    for j=find(ccc(k,:)~=1);    qqq{k, j}                   = ' ';                                  end;
    im1                         = umo_cstrs(char(qqq(k, :)), strs, 'im1');
    if ~any(im1<1);
        ccc(1:k, :)             = 0;
        cpt                     = cell2mat(qqq(sum(ccc(:,im1)==2,2)==7, im1(end-1:end)));   end;    end;
%
if exist('cpt','var')~=1;                                                           return;         end;
% removing negative values:
cpt(cpt<0)                      = 0;
% converting from microCi/ml to nCi/ml:
cpt(:,  2)                      = cpt(:,    2).*1000;

save(ofl, 'cpt', '-ascii');
disp('.done! (plasma TACs [time, nCi/mL] in ASCII)');
disp([' output: ',ofl]);
return;
%%


function                        local_hplc(pet_is,ofl);
%%
if isempty(pet_is);                                                                 return;         end;

[odx, onm]                      = fileparts(ofl);

global g4iv2;
xfls                            = dir(fullfile(odx, [onm(1, 1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf']));
if numel(xfls)>1;
    disp(['> problem! more than one ',onm(1, 1:end-size(g4iv2.xxx(1).cpt,2)-1),'.pdf']);
    disp(['  delete excessive one(s)']);
    disp(['  folder: ',odx]);
elseif numel(xfls)<1
    [fname, path_x]             = uigetfile(fullfile(fileparts(pet_is),'plasma','*.pdf'));
    if ~ischar(fname);                                                              return;         end;

    [fdx, fnm, fex]             = fileparts(fname);
    xfln                        = fullfile(odx, [onm(1, 1:end-size(g4iv2.xxx(1).cpt,2)-1),fex]);
    disp(['> copying HPLC source file..']);
    disp([' input: ',fullfile(path_x,fname)]);
    disp(['    to: ',xfln]);
    copyfile(fullfile(path_x,fname), xfln);
else;
    xfln                        = fullfile(xfls(1).folder, xfls(1).name);                           end;
%
% if exist(xfln,'file')~=2;                                                           return;         end;
disp(['> opening: ',xfln]);
winopen(xfln);
drawnow;
disp('.manually input HPLC data as follows:');
disp(' first, copy and paste time & decay corrected %parent from the pdf file as follows:');
disp([' hplc = [ ',10,' 0 95',10,'  copy and paste ',10,' ]']);
disp(' then, copy and paste the following 2 lines');
disp([' ofl = ''',ofl,''';']);
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
