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

function                        local_cpt(i2,i3)
%% 
xfls                            = dir(fullfile(i2,'*.xlsx'));
ccc                             = zeros(numel(xfls),    1);
for i=1:1:numel(xfls);
    [a, b]                      = xlsfinfo(fullfile(xfls(i).folder,xfls(i).name));
    ccc(i, :)                   = umo_cstrs(char(b), 'Twilite Data', 'im1');                        end;
% 
k                               = find(ccc>0, 1);
if isempty(k);                 
    disp('.problem! no *.xlsx has ''Twilite Data'' sheet');                         return;         end;
%
a                               = readcell(fullfile(xfls(k).folder,xfls(k).name), ...
                                                            'Sheet','Twilite Data', 'Range','A1:D1')
c14                             = readmatrix(fullfile(xfls(k).folder,xfls(k).name), ...
                                                            'Sheet','Twilite Data');
cpt                             = [c14(:,2)./60, c14(:,3).*37];
save(i3, 'cpt', '-ascii');
return;
%% 


function                        local_hplc(i2,i3)
%% 
xfls                            = dir(fullfile(i2,'*.xlsx'));
ccc                             = zeros(numel(xfls),    1);
for i=1:1:numel(xfls);
    [a, b]                      = xlsfinfo(fullfile(xfls(i).folder,xfls(i).name));
    ccc(i, :)                   = umo_cstrs(char(b), 'Data 1', 'im1');                              end;
% 
k                               = find(ccc>0, 1);
if isempty(k);                  disp('.problem! no *.xlsx has ''Data 1'' sheet');   return;         end;
%
c                               = readcell(fullfile(xfls(k).folder,xfls(k).name), 'Sheet','Data 1');
cc2                             = zeros(size(c,1), 2);
for i=1:1:size(c,1);
    if ischar(c{i,1});          
        cc2(i, 1)               = umo_cstrs(...
                                char('Radiometabolites','Time (TOI+t)','hematocrit'), c{i,1},'im1');
    elseif isnumeric(c{i, 1});  cc2(i, 2)                   = 1;                            end;    end;
%
cc2(1:find(cc2(:,1)==2,1),  2)  = 0;
cc2(find(cc2(:,1)==3,1):end, 2) = 0;
%
data                            = cell2mat(c(cc2(:, 2)>0, 1:2));
data(:)                         = [data(:,1),100 - data(:,2).*100];
save(i3, 'data', '-ascii');
return;
%%
