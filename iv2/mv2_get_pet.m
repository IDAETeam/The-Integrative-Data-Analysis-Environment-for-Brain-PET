function        out             = mv2_get_pet(lds, i2,i3); 
% To conver source MRI in the image server to IDAE-preferred format
%       
%       usage:      out         = mv2_get_pet('dxetc4xxx',input_2,input_3)
%       
%   dxetc4xxx     local directory adapter file (the name of)
%   input_2/3     date of scan / format, respectively
%                 e.g.. '7/29/2022' & 'mm/dd/yyyy'
%
%   To actually convert:        >> mv2_get_pet([],out);
%       
% 
% (cL)2022    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

out                             = [];
if isempty(lds) && ischar(i2);  local_get_pet_mri(i2);                              return;         end;
if isstruct(i2);                local_convert(i2);                                  return;         end;

x                               = feval(lds, 'get_pet', []);
%
if isdatetime(i2);              dnum                        = datenum(i2);
else;                           dnum                        = datenum(i2,i3);                       end;
%
for i=find(x.mqq>0);
    x.search_q{i}(x.search_q{i}=='$')                       = datestr(dnum, x.search_t{i});         end;
% converting search_q to a single path string:
qdx                             = x.search_q{1};
for i=2:1:size(x.mqq,2);        qdx                         = fullfile(qdx, x.search_q{i});         end;
%
% qdx
qqq                             = dir(fullfile(qdx, '*'));
%                               
if isempty(qqq);                disp(['> PET not ready for: ',qdx]);
                                out                         = [];                   return;         end;
%
cm1                             = umo_cstrs(char(qqq.folder),[], 'cm1');
%
k                               = 1;
if sum(cm1(:,2)>0)>1;           
    disp(['> PET scans done on: ',datestr(dnum)]);
    dispCharArrays(char('ID#',int2str(find(cm1(:,2)>0))),2, ...
                                                        char('Folders:',char(qqq(cm1(:,2)>0).folder)));
    while 1;
        q                       = input('> which to transfer? (by ID#; Ret=No): ','s');
        if isempty(q);          out                         = [];                   return;         end;
        k                       = str2double(q);
        if ~isempty(k) && k<size(cm1,1) && cm1(k,2)>0;                              break;          end;
                                                                                          	end;    end;
%
[fname, path_x]                 = uigetfile(fullfile(qqq(k).folder,'*'));
if ~ischar(fname);                                                                  return;         end;
if path_x(end)==filesep;        path_x                      = path_x(1, 1:end-1);                   end;
%
odx                             = feval(lds, 'pet_i2d', path_x);
odx_c                           = feval(lds, 'char2cell', odx);
%
clear out;
out                             = struct('ifl', fullfile(path_x, fname),    ...
                                    'pnm', odx_c{x.subject_seg_no_ds}, 'odx',odx);
out.done                        = numel(dir(fullfile(out.odx, '*_tra_sum.ezi')));
return;
%%

function                        local_convert(out);
%%
if exist(out.odx,'dir')~=7;     mkdir(out.odx);  
else;                           f                           = dir(fullfile(out.odx, '*_tra_sum.ezi'));
    if numel(f)==1;             disp('> previously done');
                                disp([' present: ',fullfile(f(1).folder,f(1).name)]);
                                return;
    elseif numel(f)>1;          disp('> critial problems! dupplicated *_sum.ezi');
                                disp([' folder: ',f(1).folder]);
                                dispCharArrays(2,'files:',1,char(f(1).name));
                                return;                                                     end;    end;
%
[idx, inm, iex]                 = fileparts(out.ifl);
if strcmpi(iex,'.v');           v2ezm(out.ifl,  'pnm',out.pnm,   'odx',out.odx);    return;         end;
%
f                               = dir(fullfile(fileparts(out.ifl), ['*',iex]));
dcm2umo4pet(f, out.pnm, 'odx',out.odx);
return;
%%
