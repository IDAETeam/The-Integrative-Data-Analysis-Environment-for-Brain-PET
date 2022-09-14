function    [vM, value]  = um_read(fH,fpos,verNo,mmx,giv)

% um_read:      To read attName/Value pairs from UMO files.
% 
%   usage 1:    To read attName/Value pairs from UMO files:
%               >> [vname, value]       = um_read(fH,fpos,ver#); 
% 
%       fH      -   Handle of opend UMO format file (using um_open)
%       fpos    -   To move file position to fpos before reading (optional)
%                   By defaults, um_read Starts reading where it is 
%                   After reading fH will be left where it is after reading
%       ver#    -   UMO version # as given by um_finfo(4,1);
% 
%   usage 2:    To read image/main data or attName only:
%               >> vM                   = um_read(fH,fpos,ver#,mmx);        ... ver 3.0 
%               >> vM                   = um_read(fH,fpos,ver#,mmx,giv);    ... ver 2.0 & 2.1
% 
%       mmx     -   min/max values for scaled data sets
%                   'getimaval' starting with 's' (for scale*) for ver 2.0 * 2.1
%                   adm(1) is 2**** for ver 3.0
%       giv     -   value of getimavals (required if fH is ver 2.0 or 2.1)
%                   >> giv  = um_finfo(fH,'getimaval');
% 
% (cL)2005  hkuwaba1@jhmi.edu 
%%
margin                          = 3;
if nargin<margin;               help um_read;                                       return;         end;
% -----------------------------------------------------------------------------------------------------;
if isempty(fH);                                                                     return;         end;
p                               = ftell(fH);

% getting precision string (=prec) to read administrative parts of UMO files:
verInfo                         = um_info(1,        0);
prec                            = deblank(verInfo.prec(verNo,   :));

% bringing file bit position to the specified position:
fseek(fH,   fpos,   'bof');

% reading the administrative part of the attName:
adm1                            = fread(fH,     3,  prec);
% not a UMO file:
if isempty(adm1);               disp(['.not a UMO file (or collapsed?)! ',10,' file: ',fopen(fH)]);
                                fclose(fH);
                                vM                          = [];
                                value                       = [];                   return;         end;
% adm(1) is the precision to read in number - translating the number to the matlab precision string:
h2r1                            = um_info(2,    adm1(1)); 
vM                              = fread(fH,[adm1(2),adm1(3)],  h2r1.str);

%% reading the values only when requested:
if nargin>=4;                  

    if nargin<5;                giv                     = [];                                       end;
    if isempty(giv);            giv                     = adm1(1);                                  end;
    if isempty(mmx);            mmx                     = [-Inf, Inf];                              end;
                                eqs                     = um_info(3,        giv);
                                minv                    = min(mmx(:));
                                maxv                    = max(mmx(:));
                                eval(eqs.e2r);                                      return;         end;
% -----------------------------------------------------------------------------------------------------;

%% vname is character array by defaults:
vM                              = char(vM);
if nargout==1;                                                                      return;         end;

% reading the administration part of the attValue:
adm2                            = fread(fH,     3,  prec);
h2r2                            = um_info(2,    adm2(1)); 
value                           = fread(fH,[adm2(2),adm2(3)],       h2r2.str);

% converting the valut to characters when requested:
if floor(adm2(1)./1000)==1;     value                   = char(value);                              end;
% -----------------------------------------------------------------------------------------------------;

return;
% -----------------------------------------------------------------------------------------------------;