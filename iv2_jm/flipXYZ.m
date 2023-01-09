function    i1     = flipXYZ(i1,i2,i3)

% To flip image volumes in XYZ directions
%
%       usage:      vM(:)   = flipXYZ(vM,isz,flipflg)
%                   flipXYZ(ifl,ofl,    flipflg)
%
%   vM      -   image volume matrix: L(X).*L(Y) by L(Z)
%               file name is also valid (See the second line above)
%                 .nii files will be treated separately. 
%                 .img files are not supported
%   isz     -   XYZ image volume sizes (1 by 3)
%   flipflg -   1 by 3 (=[x,y,z]). Put 1 to flip 0 not to flip
%                   [1,0,1] will flip in X and Z
%
% When vM is [x,y,z] (e.g., SPM vM):
%   Enter [] to isz (=[x,y,z])
%
% (cL)2014~18    hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help flipXYZ;                                       return;         end;

if length(i3)~=3;               disp('Errror: flipflg must be 1 by 3');             return;         end;

if ~i3(1) & ~i3(2) & ~i3(3);                                                        return;         end;

if ischar(i1);                  local_fls(i1,i2,i3);                                return;         end;

if isempty(i2);                 i1                          = local_xyz(i1,i3);     return;         end;

isz                             = i2;
iM                              = zeros(isz(1),             isz(2));    
p                               = [1:1:isz(1).*isz(2)]';
%% flipping in X direction:
if i3(1);                       q                           = reshape(p,    isz(1),isz(2));
                                q(:)                        = q([isz(1):-1:1]',:);
                                p(:)                        = q(:);
    for z=1:1:isz(3);           i1(:,   z)                  = i1(p, z);                             end;
                                                                                                    end;
p(:)                            = [1:1:isz(1).*isz(2)]';
%% flipping in Y direction:
if i3(2);                       q                           = reshape(p,            isz(1),isz(2));
                                q(:)                        = q(:,[isz(2):-1:1]);
                                p(:)                        = q(:);
    for z=1:1:isz(3);           i1(:,   z)                  = i1(p, z);                             end;
                                                                                                    end;
%% flipping in Z direction:    
if i3(3);                       zs                          = [isz(3):-1:1];
                                i1(:)                       = i1(:,zs);                             end;
return;                            
%%

function    vM                  = local_xyz(vM,flg);
%% flip 3D matrices:
if flg(1)>0;
    for z=1:1:size(vM,3);       vM(:,:,z)                   = vM(size(vM,1):-1:1, :,  z);           end;
                                                                                                    end;
if flg(2)>0;
    for z=1:1:size(vM,3);       vM(:,:,z)                   = vM(:, size(vM,2):-1:1,  z);           end;
                                                                                                    end;
if flg(3)>0;                    vM(:,:,:)                   = vM(:,:, size(vM,3):-1:1);             end;
    
return;
%%

function                        local_fls(ifl,ofl,flpflg);
%%
[idx, inm, iex]                 = fileparts(ifl);
if strcmpi(iex,'.nii');         local_nii(ifl,ofl,flpflg);                          return;         end;
if strcmpi(iex,'.img');         disp('.problem! .img no longer supported');         return;         end;

tfl                             = [];
i0                              = ifl;
if strcmpi(ifl,ofl);            disp('.not allowed, input=output');                 return;         end;
if ~exist(ifl,'file');          
    disp('.unable to locate the input file (aborting)');                            return;         end;
%
disp(['.working on .. ',i0]);
if isempty(ofl);
    tfl                         = tmpfln([],                'ezx');
    copyfile(ifl,               tfl);                   
    ofl                         = ifl;
    ifl                         = tfl;
    disp(['.input file copied to .. ',tfl]);                                                        end;
%
[isz, d0]                       = gei(i0,                   'imagesize','dataInfo'); 
vM                              = zeros(isz(1).*isz(2),     isz(3));
si                              = struct('h2s',32,'c',mfilename,'p',i0,'cp','a');
%
di                              = zeros(size(d0,1),         1);
df                              = zeros(size(d0,1),         10);
[fH, h0]                        = um_save(ofl,[],si,[],     ...
                                'flipXYZ',                  flpflg);
for i=1:1:size(d0,1);
    vM(:)                       = flipXYZ(ged(ifl, i),isz,  flpflg);
    [di(i,:), df(i,:)]          = um_save(fH,vM,si.h2s,     []);                                    end;
%
um_save(fH, h0, di, df);
disp('.done! (flipped image volume)');
disp([' output: ',ofl]);
return;
%%

function                        local_nii(ifl,ofl,flpflg);    
%% flip images stored in .nii files

v0                              = spm_vol(ifl);
vM                              = zeros(v0.dim);
vM(:)                           = spm_read_vols(v0);
vM(:)                           = local_xyz(vM,     flpflg);
v1                              = v0;
v1.fname                        = ofl;
for i=1:1:3;                    v1.mat(i, :)                = v1.mat(i, :).*((-1).^(2-flpflg(i)));  end;
v1                              = spm_create_vol(v1);
v1                              = spm_write_vol(v1, vM);
disp('.done! (flipped image volume)');
disp([' output: ',ofl]);
return;
%%


