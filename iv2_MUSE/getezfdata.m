function    [imM, p2r] = getezfdata(efln,dNo);

% getezfdata:   To recover data from ezf files.
%
%       usages:     data = getezfdata(efffln)           % the first data
%                   data = getezfdata(ezffln,dNo);      % to get data No "dNo"
%                   data = getezfdata(ezffln,[]);       % last data
%                   [data, p2r] = ...
%
%Notes:
%   1. datainfo(:,9:10) are preserved for the min and max values.
%   2. Faster for unscaled data:
%   3. ezffln can be the handle of an opened ezffln. In such a use, ...
%       a. Open 'ezffln' using <<ezopen>>.
%       b. The file will be left open. Thus, close it somewhere else.

imM                         = [];
margin                      = 1; 
if nargin<margin;           help getezfdata;                                        return;         end;

if nargin==1;               dNo                         = 1;                                        end;
dNo                         = dNo(1);

if ischar(efln);
% -----------------------------------------------------------------------------------------------------;

    vNo                     = um_finfo(efln,1);
    if vNo==3;              imM                         = um_gdata(efln,    dNo);   return;         end;


    fH                      = ezopen(efln,'r');
    if fH<0;                                                                        return;         end;

    [idx, iname, iext]      = fileparts(efln);
    if strncmp(iext,'.img',4);
    % spm .img? ---------------------------------------------------------------------------------------;

        global defaults;
        if isempty(defaults);   spm_defaults;                                                       end;
        defaults.analyze.flip   = 0;

        V                       = spm_vol(efln);
        vM                      = spm_read_vols(V,0);
        imM                     = reshape(vM,   prod(V.dim(1,1:2)),V.dim(3));

        return;                                                                                     end;
    % -------------------------------------------------------------------------------------------------;


else;

    vNo                     = um_finfo(efln,1);
    if vNo==3;              imM                         = um_gdata(efln,    dNo);   return;         end;
    

    fH                  = efln;                                         end;





[dIndex, dInfo, givs] ...
                    = getezfinfo(fH,'dataindex','datainfo','getimavals');
if isempty(dIndex); disp(['Error: Unable to read ',efln]);  return;     end;

dL                      = size(dIndex,1);
if dNo<1 | dNo>dL;

    disp(['Error: Invalid data No ... ',int2str(dNo),'/',int2str(dL)]);     
    return;                                                             end;


if isempty(dNo);    dNo    = length(dIndex);                            end;
datainfo            = dInfo(dNo,:);

precs               = ['uint16';'uint32'];
v21ORnot            = isv21(fH)+1;
prec                = precs(v21ORnot,:);

fseek(fH,dIndex(dNo,1),'bof');
adm                 = fread(fH,3,prec); 
p2r                 = convprec(adm(1));
imM                 = fread(fH,[adm(2),adm(3)],p2r); 

if ischar(efln);    fclose(fH);                                         end;
if givs(1)=='u';                                            return;     end;

% Re-scaling using <<imagevaltype>>
minv                = dInfo(dNo,9); 
maxv                = dInfo(dNo,10);
[stype, e2w, e2r]   = imagevaltype(givs); 
eval(e2r);



