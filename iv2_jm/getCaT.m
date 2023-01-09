function    [cpt, tad]          = getCaT(cptfln, varargin); 

% getCaT:       To recover Ca(t) from files (all versions)
%       
%       usage:      cpt         =  getCaT('cptfln')
%
% Output:
%   cpt     -   [time,parent,total] of plasma TACs by default
%               use clm 
%
% Options:
%   'clm',val   -   To control the column order.    
%       default: [time,parent,total], see above
%       total/measured      -   Ca(t) of total radioactivity  
%       met.corr/met-corr   -   Ca(t) with metabolite correction
%       parent is also valid for this
%       met x               -   Ca(t) for metaboliete x (x=1,2,3, ...)
%       val = {'time','total','met.corr','met 1'}
%       will return [time,total(or measured),met.corr(or met-corr),'met 1'] (n by 4)
%       When one or more requested columns are missing, output cpt will be [];
%
% (cL)2006~15    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

clmval                          = char('time','parent','total');
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
% -----------------------------------------------------------------------------------------------------;

[cpt, tad]                      = local_retCPT(cptfln,clmval);    
return;
%% the rest is ignored for now:
if clmflg;                      
    [cpt, tad]                  = local_retCPT(cptfln,clmval);                      return;         end;

% files generate by <<prepCPT>> have 'm4tad' and 'p4tad' fields but not 'T&DEstimates':
if exist(cptfln,'file')==2;
    [tad, m4tad, p4tad]         = gei(cptfln,               'T&DEstimates','model4tad','param4tad');
    cpt                         = ged(cptfln,               1);                                     end;
% cpt
% when T&D estimate is not completed.
if isempty(tad);                tad                         = [0,0,1,0];
    if isempty(m4tad);          disp('Info: T&D correction not performed.');
    else;                       disp('Info: using T&D corrected Ca(t)');
                                tad(:)                      = [p4tad(6),0,1,0];                     end;
else;                           disp(['T&D correction (TAD only) = ',num2str(tad,3),' (min)']);
                                cpt(:,1)                    = cpt(:,1) + tad(1);                    end;

return;
%%

function    [out, tad]          = local_retCPT(i1,clmval);   
%%
%   output cpt is always [time,total,met.corr] (+ met i etc, if any)
%
out                             = [];
tad                             = [0,0,1,0];
cpt                             = ged(i1,   1);
% if isempty(clmval);                                                                 return;         end;
% if ~iscell(clmval);             disp('Wrong ''val'' for ''clm'' ooption (not a cell array)');
%                                 return;                                                             end;
%%
qqq                             = gei(i1,                   'orientation');
qqq(qqq=='/')                   = ' ';
qq0                             = '[],';
for i=1:1:length(qq0);          qqq(qqq==qq0(i))            = ' ';                                  end;
qq2                             = line2mat(qqq);
ss2                             = '[';
for i=1:1:size(qq2,1);          ss2                         = [ss2,deblank(qq2(i,:)),', '];         end;
disp(['.columns of the cpt file .. ',ss2(1, 1:end-2),']']);
strs                            = str2mat('time','total','measur','met.co','met-co','parent',   ...
                                'met 1','met 2','met 3','met 4','met 5','met 6');
catnos                          = [1,2,2,3,3,3,4,5,6,7,8,9]';
im0                             = umo_cstrs(strs,lower(char(clmval)),   'im1');
clmc                            = char(clmval);
ss0                             = '[';
for i=1:1:size(clmc,1);         ss0                         = [ss0,deblank(clmc(i,:)),', '];        end;
disp(['.requested entities are  .. ',ss0(1, 1:end-2),']']);

if any(~im0);
    disp('.error in the value of ''clm'' option');
    disp(' .fix it using the following terms');
    disp(char(strs));                                                               return;         end;
%
s4                              = zeros(length(im0),        1);
for i=1:1:length(im0); 
    imx                         = umo_cstrs(lower(qq2), strs(catnos==catnos(im0(i)),:), 'im1');
    if any(imx>0);              s4(i,   :)                  = imx(find(imx>0,1));           end;    end;
% 
if any(~s4);
    disp('.following requested entities are missing ..');
    for i=find(s4'==0);         disp(clmc(i,    :));
                                disp(' .which could be one of ..');
                                disp(strs(catnos==catnos(im0(i)),   :));
                                disp('  in the cpt file');                                          end;
    disp('.check the cpt file and retry');                                          return;         end;
%
ss4                             = '[';
for i=1:1:size(s4,1);           ss4                         = [ss4,deblank(qq2(s4(i),:)),', '];     end;
disp(['.reporting columns are   .. ',ss4(1, 1:end-2),']']);
out                             = cpt(:,    s4);
out(out<0)                      = 0;
return;
%%
