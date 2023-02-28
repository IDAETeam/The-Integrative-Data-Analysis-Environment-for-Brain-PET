function    out                = mv2_get_mri(lds,i2,i3); 
% To conver source MRI in the image server to IDAE-preferred format
%       
%       usage:      out         = mv2_get_mri('dxetc4you','mri_date','date_format)
%                   
% 
% (cL)2022    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

out                             = [];
if isstruct(i2);                out                         = local_convert(i2);    return;         end;

x                               = feval(lds, 'get_mri', []);
if isdatetime(i2);              dnum                        = datenum(i2);
else;                           dnum                        = datenum(i2,i3);                       end;
for i=find(x.mqq>0);
    x.search_q{i}(x.search_q{i}=='$')                       = datestr(dnum, x.search_t{i});         end;
%
qdx                             = x.search_q{1};
for i=2:1:size(x.mqq,2);     	qdx                      	= fullfile(qdx, x.search_q{i});      	end;
%
qqq                             = dir(fullfile(qdx, '*'));
if isempty(qqq);                
    disp(['> no MRI folders for: ',datestr(datenum(i2,i3))]);                     	return;         end;
%
cm1                             = umo_cstrs(char(qqq.folder),[], 'cm1');
%
k                               = 1;
if sum(cm1(:,2)>0)>1;           
    disp(['> MRI scans done on: ',datestr(dnum)]);
    dispCharArrays(char('ID#',int2str(find(cm1(:,2)>0))),2, ...
                                                        char('Folders:',char(qqq(cm1(:,2)>0).folder)));
    while 1;
        q                       = input('> which to transfer (by ID#; Ret=No): ','s');
        if isempty(q);                                                           	return;     	end;
        k                       = str2double(q);
        if ~isempty(k) && k<size(cm1,1) && cm1(k,2)>0;                              break;          end;
                                                                                          	end;    end;
%    
[fname, path_x]                 = uigetfile(fullfile(qqq(k).folder,'*'));
if ~ischar(fname);                                                                  return;         end;
if ~ischar(fname);              out                         = [];                   return;         end;
if path_x(end)==filesep;        path_x                      = path_x(1, 1:end-1);                   end;
[idx, inm, iex]                 = fileparts(fname);
%      
%
path_c                          = feval(lds, 'char2cell', path_x);
%
if numel(path_c)==x.series_seg_no_is-1;
                                path_c{end+1}               = '*';                  
elseif numel(path_c)<x.series_seg_no_is-1;
    disp(['>fatal error! # of path segments of input file (n=',int2str(numel(path_c)),  ...
                                ') too short than expected (n=',int2str(x.series_seg_no_is-1),')']);
                                                                                    return;         end;
% 
for j=x.series_seg_no_is:1:numel(path_c);
                                path_c{j}                   = '*';                                  end;
%
idx                             = path_c{1};
for j=2:1:numel(path_c);        idx                         = fullfile(idx, path_c{j});             end;
%
mmm                             = dir(fullfile(idx, ['*',iex]));
mm2                             = dir(fullfile(fileparts(idx), ['*',iex]));
%
if numel(mmm)>=numel(mm2);
    if strcmpi(iex,'.par') || strcmpi(iex,'.rec');
                                out                         = local_par(mmm, lds, x, path_c);
    else;                       out                         = local_dcm(mmm, lds, x);               end;
else;
    if strcmpi(iex,'.par') || strcmpi(iex,'.rec');
                                out                         = local_par(mm2, lds, x, path_c);
    else;                       out                         = local_dcm(mmm, lds, x);       end;    end;
return;
%%

function    out                 = local_dcm(mmm,lds,x);
%%
out                             = [];
% checking if inputs (=mmm) are sorted into separate folders:
cm2                             = umo_cstrs(char(mmm.folder),[], 'cm1');
%
nnn                             = zeros(sum(cm2(:,2)>0), 1);
ic                              = 0;
clear dd2 dd1
for j=find(cm2(:,2)'>0);
    ic                          = ic + 1;
    h                           = dicominfo(fullfile(mmm(j).folder,mmm(j).name));
    dd2{ic}                     = [h.SeriesDescription,'; ',h.StudyDescription];
    mdc                         = feval(lds, 'char2cell', mmm(j).folder);
    nnn(ic, :)                	= numel(dir(fullfile(feval(lds,'mri_i2d',mmm(j).folder),'*_tmp.mri')));
    dd1{ic}                     = mdc{x.series_seg_no_is};                                          end;
%
disp('> available MRI series:');
dispCharArrays(char('Series #s',char(dd1)),2,char('Descriptions',char(dd2)),   ...
    2,char('# preset:',int2str(cm2(cm2(:,2)>0, 2))),2,char('done:',int2str(nnn)));
if ic==1;                       q                       	= dd1{1};
else;
    while 1;
        q                       = input('> which to transfer (c/p @column 1; Ret=No): ','s');
        if isempty(q);                                                           	return;     	end;
        if umo_cstrs(char(dd1), q, 'im1')>0;                                        break;          end;
                                                                                  	        end;    end;
%
m2cp                            = umo_cstrs(char(dd1), q, 'im1');
ii                              = find(cm2(:,2)'>0);
clear out;
out.mmm                         = mmm(cm2(:,1)==cm2(ii(m2cp),1));
% 
out.odx                         = feval(lds, 'mri_i2d', mmm(ii(m2cp)).folder);
odx_c                           = feval(lds, 'char2cell', out.odx);
out.snm                         = odx_c{x.subject_seg_no_ds};
out.done                        = nnn(m2cp);
if out.done>0;                  disp(['> previously done for: ',q]);                return;         end;

return;
%%

function    out                 = local_par(par, lds, x, path_c);
%%
qqq                             = [];
nnn                             = zeros(numel(par),     3);
v2r                             = {'Acquisition nr', 'Reconstruction nr', ...
                                                            'Max. number of sl','Protocol name'};
%
for i=1:1:numel(par);
    clear qqq;
    qqq                       	= umo_getptf(fullfile(par(i).folder, par(i).name), 1,[]);
    [c1, c2]                    = getLseg(qqq(qqq(:,1)=='.', :), 1);
    im1                         = umo_cstrs(c2,char(v2r), 'im1');
    nnn(i, 1:2)                 = [str2double(c2(im1(1), find(c2(im1(1),:)==':',1)+1:end)).*100 ...
                                    + str2double(c2(im1(2), find(c2(im1(2),:)==':',1)+1:end)),  ...
                                    str2double(c2(im1(3), find(c2(im1(3),:)==':',1)+1:end))];
    path_c{x.series_seg_no_is}  = int2str(nnn(i,1));
    clear idx;
    idx                         = path_c{1};
    for j=2:1:numel(path_c);    idx                         = fullfile(idx, path_c{j});             end;
    odx_c{i}                    = feval(lds, 'mri_i2d', idx);
    nnn(i, 3)                   = numel(dir(fullfile(odx_c{i}, '*_tmp.mri')));
    ppp{i}                      = deblank(c2(im1(4), find(c2(im1(4),:)==':',1)+1:end));             end;
%
disp('> protocol names of inquired .par files');
disp([' folder: ',par(1).folder]);
dispCharArrays(1,char('#',int2str([1:1:size(nnn,1)]')),2,char('Protocol Names:',char(ppp)),2,   ...
    char('# of slices:',int2str(nnn(:,2))),2,char('Done:',int2str(nnn(:,3))),   ...
                                2,char('files:',char(par.name)));
disp('<end of the list');

while 1;
    res                         = input('Select the one to convert (1st column) [no=Ret]: ','s');
    if isempty(res);                                                                return;         end;
    ic                          = floor(str2double(res));
    if ~isempty(ic) && ic>0 && ic<=size(nnn,1);                                  	
        if nnn(ic, 3)>0;        disp('> previously done');
                                disp([' folder: ',odx_c{ic}]);
                                out                         = [];                   return;         end;
                                                                                    break;          end;
                                                                                                    end;
%
out.par                         = fullfile(par(ic).folder, par(ic).name);
out.odx                         = odx_c{ic};
out.snm                         = path_c{x.subject_seg_no_ds};
return;
%%

function    out                 = local_convert(iis);
%%
out                             = [];
if exist(iis.odx,'dir')~=7;     mkdir(iis.odx);                                                     end;
%
if isfield(iis, 'par');
    out                         = rec2umo(iis.par, '.nii', 'pnm',iis.snm, 'odx',iis.odx);
                                                                                    return;         end;
%
iis.done                        = numel(dir(fullfile(iis.odx, '*_tmp.mri')));
if exist(iis.odx,'dir')~=7;     mkdir(iis.odx);                                                     end;
if iis.done<1
    out                         = dcm2umo4mri(iis.mmm, iis.snm, 'idx',iis.mmm(1).folder,    ...
                                                            'odx',iis.odx,    'nii','on');
else;                           disp(['> previousely done', 10,' for: ',iis.odx]);  return;         end;
% %
% if exist(out,'file')==2;
%     h                           = dicominfo(fullfile(iis.mmm(1).folder,iis.mmm(1).name));
%     [odx, onm]                  = fileparts(out);
%     write2ptf(fullfile(odx, [onm,'.mri']),  h.SeriesDescription);
%     disp('.done! (symbolic MRI file with SeriesDescription)');
%     disp([' output: ',fullfile(odx, [onm,'.mri'])]);                                                end;
return;
%%


