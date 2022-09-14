function    acpc    = getACPC(i1,i2, o1,v1,o2,v2,o3,v3,o4,v4,o5,v5); 

% getACPC:      To recover ACPC point XYZ coordinates from an UME file
%       
%       usage:      acpc        =   getACPC('file',outflg)
%                   or
%                   acpc(:)     =   getACPC(acpc,outflg,[isz,vsz]);
%       
%   file    -   file containing ACPCpoints as an attribute or as a data set
%   outflg  -   'px' in pixels or 'mm' in mm, eucladien
%
% Options:      


margin                          = 2;
if nargin<margin;               help getACPC;                                       return;         end;
% -----------------------------------------------------------------------------------------------------;

if nargin==3;                   v1                          = o1;
                                o1                          = 'szs';                                end;

szsval                          = [];
dspval                          = [];
opt                             = ['szs';'dsp'];
n0                              = nargin;
options;
if ~OptionsOK;                                                                      return;         end;
% -----------------------------------------------------------------------------------------------------;

acpc                            = [];

oNo                             = whichstr(['px';'mm'],i2);
if ~oNo;                        disp('Error: Invalid ''outflg (getACPC)');          return;         end;
% -----------------------------------------------------------------------------------------------------;

 
if ischar(i1);
% -----------------------------------------------------------------------------------------------------;
    
    [isz, vsz, acpc]            = gei(i1,                   'imagesize','voxelsize','acpcPoints'); 

    if isempty(acpc);           d                           = ged(i1,           1);

        if size(d,1)<4 & size(d,2)==3;                          
                            
                                acpc                        = d;
    
        else;                   disp(['Error: ACPC not found (getACPC)']);          return;         end;
    % -------------------------------------------------------------------------------------------------;
                                                                                                    end;

else;
% -----------------------------------------------------------------------------------------------------;

    acpc                        = i1;
    if isempty(szsval);                                                             return;         end;
    isz                         = szsval(1,:);
    vsz                         = szsval(2,:);
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;

cog                             = (isz + 1)./2;
if abs(acpc(1,1) - cog(1)) > abs(acpc(1,1) - 0);
% acpc is in mm ---------------------------------------------------------------------------------------;

    if oNo==1;                  acpc(:,1:3)                 = mm2pixels(acpc(:,1:3),isz,vsz,'mm');  end;

else;
% acpc is in pixels -----------------------------------------------------------------------------------;

    if oNo==2;                  acpc(:,1:3)                 = mm2pixels(acpc(:,1:3),isz,vsz,'px');  end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
                                                                                                   
if dspflg & exist(dspval,'file');
% -----------------------------------------------------------------------------------------------------;
    
    spm                         = gei(dspval,               'spm_mat');
    if isempty(spm);            dinfo                       = [1,1];
    else;                       dinfo                       = [1,0];                                end;
    p0n                         = unite_dHxSPM(dspval,      dinfo);

    if isempty(spm);            acpc(:, 1:3)                = doDispl2D(acpc(:, 1:3),   p0n.mat);
    else;                       xyzT                        = [acpc(:,1:3)';ones(1,size(acpc,1))];
                                xyzT                        = p0n.mat(:,1:4)\(p0n.mat(:,5:8)*xyzT);
                                acpc(:, 1:3)                = xyzT(1:3, :)';                        end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;
    
