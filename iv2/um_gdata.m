function    [vM, isz]       = um_gdata(i1,i2,varargin);

% um_gdata:     To retrieve image/main data sets from UMO format files
% 
%   usage 1:    For image matrices:
%               >> [vM, isz]    = um_gdata('umofln',frm#);                  ... in nCi/ml
%               >> [vM, isz]    = um_gdata('umofln',frm#,   'unt','Bq/ml'); ... in Bq/ml
%               >> [vM, isz]    = um_gdata('umofln',frm#,   'dsc',64);      ... for display
% 
%       umofln  -   UMO image file name
%       frm#    -   frame #
%       vM      -   image matrix (isz(1).*isz(2) by isz(3))
%                   UMO format permit only one voxel order scheme:
%                   firstest component          -   left to right
%                   second firstest component   -   inion to nasion (back to front)
%                   slowest component           -   base to cranium (bottom to top)
%       isz     -   image box sizes (1 by 3)
% 
%   usage 2:    For non-image data sets:
%                >> [data, dsz]  = um_gdata('umofln',data#)
% 
% Options:
%   'unt',val   - unit of image data:   'nCi/ml'(=default), 'mcirCi/ml', 'Bq/ml', or 'MBq/ml'
%   'dsc',val   - iM in integers for display. val = length of colormap 
%   'mmx',val   - to use val(1)<=iM<val(2).
% 
% (cL)2005  hkuwaba1@jhmi.edu 

margin                          = 2; 
if nargin<margin;               help um_gdata;                                      return;         end;
% -----------------------------------------------------------------------------------------------------;

% default values for options;
dscval                          = []; 
untval                          = 'nCi/ml';
mmxval                          = [];
opt                             = ['dsc';'unt';'mmx'];
n0                              = nargin;
um_options; 
if ~OptionsOK;                                                                      return;		    end;
% -----------------------------------------------------------------------------------------------------;


%% opening input umo file:
[fH, cflg]                      = um_open(i1,   'r');
if ~fH;                                                                             return;         end;
verNo                           = cflg;
if ~verNo;                      verNo                       = um_finfo(fH,  1);                     end;
% -----------------------------------------------------------------------------------------------------;

%%
isz                             = um_finfo(fH,              'imagesize');
[vM, mmxval]                    = local_getdata(fH,         i2(1),mmxval,verNo);

if untflg;                      vM(:)                       = local_unit(vM,    untval);            end;
if mmxflg | dscflg;             vM(:)                       = local_scale(vM,   mmxval,dscval);     end;

if cflg;                        fclose(fH);                                                         end;

return;
%%

function    [vM, mmx]           = local_getdata(fH,dNo,mmx,vNo);
%%
    giv                         = [];
    if vNo<3;
        [dIndx, dInfo, giv]     = um_finfo(fH,              'dataIndx','dataInfo', 'getimaval');
    else;
        [dIndx, dInfo]          = um_finfo(fH,              'dataIndx','dataInfo');                 end;
    % UMO version #:
    vInfo                       = um_info(1,    0);
    % vInfo.prec = precision for reding administrative poart of data set:

    if isempty(mmx);            mmx                         = dInfo(dNo,    9:10);                  end;

    % reading the main/image data set:
    vM                          = um_read(fH,   dIndx(dNo), vNo,mmx,giv);

return;
%%

function    vM                  = local_unit(vM,unt);
% convertinf from nCi/ml to microCi/ml or Bq/ml:

    ustrs                       = str2mat('nCi/ml','microCi/ml','Bq/ml','MBq/ml');
    ii                          = umo_cstrs(lower(ustrs),   lower(unt),     'im1');
    if ii==2;                   vM(:)                       = vM./1000;
    elseif ii==3;               vM(:)                       = vM.*37000;
    elseif ii==4;               vM(:)                       = vM.*37./1000;                         end;

return;
%%

function    vM                  = local_scale(vM,mmx,dsc);
% rescaling vM for image display ----------------------------------------------------------------------;


    vM(find(vM>max(mmx(:))))    = max(mmx(:));
    vM(find(vM<min(mmx(:))))    = min(mmx(:));

    vM(find(isnan(vM)))         = min(mmx(:));
    
    vM(:)                       = ceil((vM - min(vM(:)))./(max(vM(:)) - min(vM(:))).*dsc(1));
    vM(find(~vM))               = 1;

return;
%%

function    local_spmimg(i1);
if ischar(i1);
% To see if i1 is a SPM.img -----------------------------------------------------------------------;

    [idx, iname, iext]          = fileparts(i1);
    if strncmp(iext,'.img',4);
    % i1 is a SPM.img -----------------------------------------------------------------------------;

        % SPM.img must be a single frame file for now:
        [iM, isz]               = ged(i1,       1);

        if dscflg;              [iM(:), isz]                = rescImM(iM,mmxval,dscval); 
        elseif bcqflg;          iM(:)                       = iM.*37000.*bcqval;                    end;

        return;                                                                             end;    end;
% -----------------------------------------------------------------------------------------------------;

return;
%%
