function	[out, out2, out3]   = um_info(i1,i2,i3);

% um_info:  To provide information on the UMO format (not on individual files)
% 
%   usage 1:    To return UMO version statements (current and previous)
%               >> verInfo          = um_info(1,0); 
% 
%       verInfo.strs    -   registered version statements
%       verInfo.prec    -   precision to read version statemtns
% 
%   usage 2:    To convert UMO precision #s to Matlab precision strings
%               >> prec             = um_info(2,prec#);
% 
%       prec#       -   one of precision #s defined for the UMO format
%       prec.str    -   Matlab precision string for prec#
%       prec.len    -   bit length (1/2/4/8 for 8/16/32/64 bits)
% 
% 
%   usage 3:    To convert precision # to equations to write/read
%               >> eqs              = um_info(3,precNo);
%               To display descriptions of precNos
%               >> um_info(3,0);  
% 
%       eqs.e2w     -   equation to convert main/image data to the format to save
%       eqs.e2r     -   equation to convert main/image data in file to usable form
%       eqs.p       -   Matlab precision string to read/write main/image data
% 
%   usage 4:    To return reserved UMO attribute names 
%               >> [rNames, flg]    = um_info(4,1);
%               To display reserved attribute names
%               >> um_info(4,0);    
% 
%       rNames      -   reserved attribute names in a character matrix
%       flg         -   0 for attNames not permitted for users to specify 
%                       1 for attNames that may be copied from a parent file
%                       2/3 for special attNames (history/pFileName)
% 
% (cL)2005  hkuwaba1@jhmi.edu 


%% Notes:
% Some functions were transferred as follows:
%   original usage 3:   To get version No of UMO format files
%                       -> verNo                        = um_open(fH,0);
%   original usage 4:   To retrieve file bit positions and attribute names from an UMO file
%                       -> [foot, attNames7, attNames]  = um_finfo(fH,0)
%   original usage 5:   To retrieve equations for writing/reading image data (main data)
%                       -> replaced by new usage 6

margin                          = 2;
if nargin<margin;               help um_info;                                       return;         end;


if isnumeric(i1);               fnm                         = ['local_',int2str(i1(1))];
else;                           fnm                         = ['local_',i1];                        end;
out                             = [];
out2                            = [];
if ~isempty(which(fnm));        [out, out2]                 = feval(fnm,    i2);                    end;
%%
return;


function    [out, out2]         = local_1(i2);
%% usage 1: inquiring version strings 

out.strs                        = str2mat('ezfFormatVer2.0','ezfFormatVer2.1','umoFormatVer3.0');
out.prec                        = str2mat('uint16','uint32','uint32');
out.L                           = [2,4,4]';
out2                            = [];
%%    
return;


 function   [out, out2]         = local_2(i2);
%% usage 2: precision No to Matlab prec string conversion


out                             = [];
estrs                           = local_pstrs;
c                               = find(estrs.v3==i2(1));
out2                            = [];
if isempty(c);
    if ischar(i2);              disp(['*** Wrong usage. Enter prec#   ... ']); 
    else;                       disp(['*** Wrong prec #   ... ',int2str(i2)]);                      end;
    return;                                                                                         end;
L                               = (estrs.v3(c(1), :) - floor(estrs.v3(c(1), :)./100).*100)./8;

out                             = struct('str',         deblank(estrs.p(c(1),   :)), ...
                                        'pno',          estrs.v3(c(1),    :), ...
                                        'len',          L);
%%    
return;
    
function    [out, out2]         = local_3(i2);
%% usage 3: returning equations to write/read 

out                             = [];
out2                            = [];
estrs                           = local_pstrs;

% displaying precision Nos (of ver3.x) and 
if isnumeric(i2) & ~i2(1);
 	disp('*** precision Nos of UMO format (ver3.x) ****');
    for i=1:1:size(estrs.v3,1); 
    if estrs.v3(i);             disp([int2str(estrs.v3(i)),'   ... ',estrs.d(i,:)]);
                                disp([' eq2read:  ',estrs.r(i,:)]);
                                disp([' eq2write: ',estrs.w(i,:)]);                         end;    end;
    disp('*** end of precision No list');                                           return;         end;

% ver2.x uses precision strings:
if ischar(i2);                  L                       = size(estrs.v2,    2);
                                i2x                     = char(zeros(1, L) + 32);
                                s                       = min([L, size(i2,2)]);
                                i2x(1,  1:s)            = i2(1, 1:s);
                                cc                      = umo_cstrs(estrs.v2, lower(i2x),   'im1');

    % i2 (prec.string) is not in the prec.string list (=estrs.v2):
    if ~cc;                     disp(['*** Wrong prec. string   ... ',i2]);
                                disp('use one of the following strings:');
                                disp([estrs.v2,char(zeros(size(estrs.v2,1),4)+32),estrs.d]);
                                disp('*** end of prec. string list');               return;         end;
        
% ver3.x uses precision numbers:
else;                           cc                      = find(estrs.v3==i2(1));

    if isempty(cc);             disp(['*** Wrong prec. number   ... ',int2str(i2(1))]);
                                um_info(6,  0);                                     return;         end;
                                                                                                    end;

out.e2r                         = deblank(estrs.r(cc(1),:));
out.e2w                         = deblank(estrs.w(cc(1),:));
out.p                           = deblank(estrs.p(cc(1),:));
%%
return;

function    [out, out2]         = local_4(i2);
%% returning reserved attribute names
%
% out2 is 1 to allow copying, 0 not to permit users to specify, 2/3 for history/pFileName:

if i2(1);                       [out, out2]             = local_rNames;
else;                           [out, out2, d]          = local_rNames;
                                disp('*** Reserved UMO attribute names');
    for i=1:1:size(out,1);      disp([int2str(out2(i)),'   ',out(i,:),d(i,:)]);                     end;
                                disp('*** end of attribute name list');
                                disp('0   attNames not permitted for users to specify');
                                disp('1   attNames that may be copied from a parent file');
                                disp('2/3 for special attNames (history/pFileName)');               end;
%%
return;

function    [out, out2]         = local_pstrs(i2)
%% precision strings/numbers and conversion (read/write) equations ------------------------------------;

out2                            = [];

    %        ver2 umo   matlab  ver3    equations to write;     
    %                                   equations to read;
    %                                   descriptions
eqs                             = ...
        [   'scaled8    int8    0       vM(:) = round((vM-minv)./(maxv-minv).*127);          ', ...
                                        'vM(:) = vM./127.*(maxv-minv) + minv;          ',   ...
                                        'scaled, using 0-127 (int8)                    ';
            'scaled8a   int8    2108    vM(:) = round((vM-minv)./(maxv-minv).*255) - 128;    ', ...
                                        'vM(:) = (vM+128)./255.*(maxv-minv) + minv;    ',   ...
                                        'scaled, using -128-127 (int8)                 ';
            'unscaled8  int8    108     vM(:) = round(vM);                                   ', ... 
                                        'vM(:) = vM;                                   ',   ... 
                                        'rounded, using -128-127 (int8)                ';
            'unscaled8a uint8   208     vM(:) = round(vM);                                   ', ...
                                        'vM(:) = vM;                                   ',   ... 
                                        'rounded, using 0 and 255 (uint8)              ';
            'character  uchar   1208    vM(:) = vM;                                          ', ...
                                        'vM = char(vM);                                ',   ...
                                        'characters, using 0 and 255 (uchar)           ';
            'scharacter schar   1108    vM(:) = vM;                                          ', ...
                                        'vM = char(vM);                                ',   ...
                                        'characters, using -128 and 127 (schar)        ';
            'scaled16p  int16   0       vM(:) = round((vM-minv)./(maxv-minv).*32767);        ', ...
                                        'vM(:) = vM./32768.*(maxv-minv) + minv;        ',   ... 
                                        'scaled, using 0-32767 (int16)                 ';
            'scaled16a  int16   2116    vM(:) = round((vM-minv)./(maxv-minv).*65535) - 32768;', ...
                                        'vM(:) = (vM+32768)./65535.*(maxv-minv) + minv;',   ... 
                                        'scaled, using -32768-32767 (int16)            ';
            'scaled16*n int16   0       vM(:) = round(vM.*nv);                               ', ...
                                        'vM(:) = vM./nv;                               ',   ... 
                                        'rounded, after dividing with nv (int16)       ';
            'unscaled16 int16   116     vM(:) = round(vM);                                   ', ...
                                        'vM(:) = vM;                                   ',   ... 
                                        'rounded, using -32768-32767 (int16)           ';
            'unscaled16aint16   216     vM(:) = round(vM);                                   ', ...
                                        'vM(:) = vM;                                   ',   ... 
                                        'rounded, using 0-65535 (uint16)               ';
            'unscaled32sint32   132     vM(:) = round(vM);                                   ', ...
                                        'vM(:) = vM;                                   ',   ... 
                                        'rounded, using -2147483647-2147483648 (int32) ';
            'unscaled32uuint32  232     vM(:) = round(vM);                                   ', ...
                                        'vM(:) = vM;                                   ',   ... 
                                        'rounded, using 0-256.*256.*256.*256-1 (uint32)';
            'unscaled32 float32 32      vM(:) = vM;                                          ', ...
                                        'vM(:) = vM;                                   ',   ... 
                                        'using float32                                 ';
            'unscaled64 float64 64      vM(:) = vM;                                          ', ...
                                        'vM(:) = vM;                                   ',   ... 
                                        'using float64                                 '];


out                             = struct('v2',              eqs(:, 1:11),           ...
                                        'p',                eqs(:, 12:18),          ...
                                        'v3',               str2num(eqs(:, 20:23)), ...
                                        'w',                eqs(:,  28:80),         ...
                                        'r',                eqs(:,  81:126),        ...
                                        'd',                eqs(:,  127:end));
%%
return;

function    [rNames, cflg, estrs]       = local_rNames(i2);
%% 
%
%   cflg:   0   -   reserved attNames which users are not allowed to use
%           1   -   items to copy from pFile name, if info.p is a file
%           2   -   'history' which require a specific treatment
%           3   -   'pFileName' which require a specific treatment

rNames                          = str2mat(  ...
                                    'PatientName',  'studyIDNo',    'imagesize',    'voxelsize',    ...
                                    'orientation',  'imageType',    'dataUnit',     'history',      ...
                                    'pFileName',    'infoIndex',    'dataInfo',     'dataIndex');
cflg                            = ones(size(rNames,1),      1);
cflg(8, :)                      = 2;                        % giving 2 to history
cflg(9, :)                      = 3;                        % marking pFileName with 3
cflg(10:end,    :)              = 0;                        % assigning 0 to not2inherit variables

% definitions of must items:
if nargout==3;
    estrs                       = [ '   ...   Subject name        (characters)                    ';
                                    '   ...   Study ID number     (integer or characters)         ';
                                    '   ...   XYZ image box sizes (1 by 3 integer matrix)         ';
                                    '   ...   XYZ voxel sizes     (1 by 3 in mm)                  ';
                                    '   ...   Depends, tra [1,1,1], if main data sets are images  ';
                                    '   ...   Description of the main data sets                   ';
                                    '   ...   Data unit (of the main data sets)                   ';
                                    '   ...   Code name           (characters, reserved)          ';
                                    '   ...   File name from which attributes were copied from    ';
                                    '   ...   Bit positions in file of attributes (reserved)      ';
                                    '   ...   Number of data sets by 10   (9:10 for min/max)      ';
                                    '   ...   Bit positions in file of data sets  (reserved)      ']; 
                                                                                                    end;
%%
return;
