function	[bits, out2]      = um_write(fH,an,vM,i4);

% um_write:     To write attName/attValue pairs in UMO-format files
%
%   usage 1:    To write an attName/attValue pair
%               >> pos      = um_write(fH,attName,attValue);
%
%       fH          -   The handle of opened umo-format file to write
%                       Bring the file bit position where to start writing
%       attName     -   The name of data attribute      (character string)
%       attValue    -   The value of attName            (string, integer, real)
%       pos         -   The bit position of start writing attName
%
%   usage 2:    To write attName only (special for the version statement)
%               >> pos      = um_write(fH,'eziformatver2.0',[]);
%
%   usage 3:    To write a main/image data set (=vM)
%               >> pos      = um_write(fH,[],vM,pNo)
%
%       pNo         -   UMO precision # to write/read
%                       To display precision #s     ... >> um_info(3,0)
%
% Examples:
%   pos         = um_write(fH,'patientname','John F Kennedy');
%   pos         = um_write(fH,'voxelsize',[1.24,1.43,3.18]);
%
% Notes:
%   1.  Bring file bit position to where to start writing
%   2.  The file bit position will be where it is after writing
%
% (cL)2005  hkuwaba1@jhmi.edu 

margin                          = 3; 
if nargin<margin;               help um_write;                                      return;         end;
if nargin<4;                    i4                          = [];                                   end;

bits                            = ftell(fH);
out2                            = [];

out2                            = ...
                                feval(['local_',int2str(isempty(an)),int2str(isempty(vM))],fH,an,vM,i4);
return;
%%

%%

function    out2                = local_11(fH,an,vM,i4);
%% undefined usage:

out2                            = [];
disp('.error! wrong usage (um_write)');         
return;
%%

function    out2                = local_10(fH,an,vM,i4);
%% When attName is empty - writing data only (usage 3):

out2                            = [];
% This usage requires h2w:
if isempty(i4);                 disp('.problem! specify precision to write');
                                disp('>> bits  = um_write(fH,[],data,p2s);');       return;         end;

% eqs.e2w dictates how to treat vM before saving:
eqs                             = um_info(3,    i4(1));
if isempty(eqs);                                                                    return;         end;

% Matlab can recover NaN and Inf when they are saved by floatxx:
if any(any(isnan(vM))) | any(any(isinf(vM)));
    fwrite(fH,                  [32,size(vM)],          'uint32');
    fwrite(fH,                  vM,                     'float32');                 return;         end;

% eqs.e2w uses minv and maxv for scaling:
minv                            = min(vM(:));
maxv                            = max(vM(:));
out2                            = [minv, maxv];
eval(eqs.e2w);

% eqs.p is the matlab precision string file reading/writing:
fwrite(fH,                      [i4(1),size(vM)],       'uint32');
fwrite(fH,                      vM,                     eqs.p);
return;
%%

function    out2                = local_01(fH,an,vM,i4);
%% When attVal is empty - writing attName only (usage 2):

out2                            = [];
if ~ischar(an);                 disp(['*** attName must be a character string']);   return;         end;

fwrite(fH,                      [1208,size(an)],        'uint32');
fwrite(fH,                      an,                     'uchar');
return;
%%

function    out2                = local_00(fH,an,vM,i4);
%% Both attName and attVal are entered (usage 1):
out2                            = [];
% saving attName:
fwrite(fH,                      [1208,size(an)],        'uint32');
fwrite(fH,                      an,                     'uchar');

% saving attValue - when vM are characters:
if ischar(vM);      
    fwrite(fH,                  [1208,size(vM)],        'uint32');
    fwrite(fH,                  vM,                     'uchar');                   return;         end;

% Matlab can recover NaN and Inf when they are saved by floatxx:
if any(isnan(vM(:))) || any(isinf(vM(:)));
    fwrite(fH,                  [32,size(vM)],          'uint32');
    fwrite(fH,                  vM,                     'float32');                 return;         end;

% when vM is made up of real numbers:
if sum(abs(ceil(vM(:)) - floor(vM(:))));
    if strcmpi(an,'dataInfo');
        % disp([an,': float64']);
        fwrite(fH,          	[64,size(vM)],          'uint32');
        fwrite(fH,            	vM,                     'float64');         
    else;
        fwrite(fH,            	[32,size(vM)],          'uint32');
        fwrite(fH,             	vM,                     'float32');                                 end;
                                                                                    return;         end;
    
% vM is made up of integers - using int8:
if max(abs(vM(:)))<128;
    fwrite(fH,                  [108,size(vM)],         'uint32');
    fwrite(fH,                  vM,                     'int8');                    return;         end;

% vM is made up of integers - using int16:
if max(abs(vM(:)))<128.*128;
    fwrite(fH,                  [116,size(vM)],         'uint32');
    fwrite(fH,                  vM,                     'int16');                   return;         end;

fwrite(fH,                      [132,size(vM)],         'uint32');
fwrite(fH,                      vM,                     'int32'); 
return;
%% 