function            [oM, osz]   =   transViews(iM,isz,viw, varargin); 

% To rearrangea 3-D image matrix as defined below
%
%   Usage 1: to rearrange iM to tM (see below for their definitions)
%
%               [iM, isz]       = transViews(tM,tsz,viw);
%
% Inputs:
%   tM  :   by definition, reshape(tM(:,z),tsz(1),tsz(2)) is a trans-axial image @z 
%   tsz :   [xn,yn,zn] (image dimensions) of tM (=tra-oriented images)
%   viw :   output orientation: 'sag'/'cor' for saggital/coronal orientations 
% Outputs:
%   iM  :   reshape(iM(:,x),isz(1),isz(2)) will be a saggital image @x, if viw = 'sag', or
%           reshape(iM(:,y),isz(1),isz(2)) will be a coronal image @y, if viw = 'cor'
%   isz  -  image dimensions in 'viw' orientation. 
%           [yn,zn,xn] for 'sag', [xn,zn,yn] for 'cor'.
%
%   Usage 2: to rearrange iM to tM (see Usage 1 for their definitions)
%
%               [tM, tsz]       = transViews(iM,isz,viw,'ovw','tra');
%
% Options:
%   'ovw',val   -   output orientation
%
% (cL)2002~18    hkuwaba1@jhmi.edu 

% abolished:
%
%       usage 3:    transViews(ifl,[],viw, 'ovw',val, 'ofl',val)
%
% To transfer image volumes from 'sag' or 'cor' to 'tra' orientations
% Inputs: 
%   viw  -  input orientation (not output orientation)
%           need to use the 'ovw' option
%
%   'ofl',val   -   output file (default = full/path/*_ovw.ext)
%   'flp',val   -   to flip in output orientation 
%                   val = 1 by 3, 1 to flip & 0 not flip
%

margin                          = 3; 
if nargin<margin;               help(mfilename);                                    return;         end;
%
ovwval                          = 'tra';
% oflval                          = [];
% flpval                          = [0,0,0];
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;
if nargin==3;                   ivwval                      = 'tra';
                                ovwval                      = viw;
else;                           ivwval                      = viw;                                  end;
%
if isempty(which(['local_',ivwval,'_2_',ovwval]));
    disp(['.problem! matrix conversion of ',ivwval,' to ',ovwval,' views not supported']);
                                                                                    return;         end;
                                                                                
[oM, osz]                       = feval(['local_',ivwval,'_2_',ovwval],iM,isz);           
return;
%%

function    [iM, isz]           = local_tra_2_sag(tM, tsz);
%% 
% sag = [y, z, x]
isz                             = tsz(1,    [2,3,1]);
iM                              = zeros(isz(1).*isz(2),     isz(3)); 
qM                              = zeros(tsz);
qM(:)                           = reshape(tM(:), tsz(1), tsz(2), tsz(3));
for i=1:1:isz(3);               iM(:, i)                    = reshape(qM(i,:,:),isz(1).*isz(2),1);  end;
return;
%% 

function    [iM, isz]           = local_tra_2_cor(tM, tsz);
%% 
% cor = [x, z, y]
isz                             = tsz(1,    [1,3,2]);
iM                              = zeros(isz(1).*isz(2),     isz(3)); 
qM                              = zeros(tsz);
qM(:)                           = reshape(tM(:), tsz(1), tsz(2), tsz(3));
for i=1:1:isz(3);               iM(:, i)                    = reshape(qM(:,i,:),isz(1).*isz(2),1);  end;
return;
%% 

function    [tM, tsz]           = local_sag_2_tra(iM, isz);
%% 
% sag = [y, z, x]
tsz                             = isz(1,    [3,1,2]);
tM                              = zeros(tsz(1).*tsz(2),     tsz(3)); 
qM                              = zeros(isz);
qM(:)                           = reshape(iM(:), isz(1), isz(2), isz(3));
sM                              = zeros(isz(1), isz(3));
for i=1:1:tsz(3);               sM(:)                       = qM(:,i,:);
                                tM(:, i)                    = reshape(sM',tsz(1).*tsz(2),1);        end;
return;
%% 

function    [tM, tsz]           = local_cor_2_tra(iM, isz);
%% 
% cor = [x, z, y]
tsz                             = isz(1,    [1,3,2]);
tM                              = zeros(tsz(1).*tsz(2),     tsz(3)); 
qM                              = zeros(isz);
qM(:)                           = reshape(iM(:), isz(1), isz(2), isz(3));
for i=1:1:tsz(3);               tM(:, i)                    = reshape(qM(:,i,:),tsz(1).*tsz(2),1); 	end;
return;
%% 

% 
% if ~ovwflg;                     ivwval                      = 'tra';
%                                 ovwval                      = viw;
% else;                           ivwval                      = viw;                                  end;
% % checking input / output orientations (ivwval / ovwval):
% tcs                             = {'tra','cor','sag'};
% vNo                             = umo_cstrs(char(tcs),ivwval,   'im1');
% oNo                             = umo_cstrs(char(tcs),ovwval,   'im1');
% ok                              = 1;
% if vNo<1;                       disp(['.error! unknown input orientation: ',ivwval]);
%                                 disp(' enter one of tra/cor/sag (aborting)');
%                                 ok                          = 0;                                    end;
% if oNo<1;                       disp(['.error! unknown output orientation: ',ovwval]);
%                                 disp(' enter one of tra/cor/sag (aborting)');
%                                 ok                          = 0;                                    end;
% if ~ok;                                                                             return;         end;
% %    
% if ischar(iM);                  local_fls(iM,ivwval,ovwval,oflval,flpval);          return;         end;
% %
% %
% [g1, g2, g3]                    = ndgrid(1:1:isz(1),        1:1:isz(2),     1:1:isz(3));
% vvv                             = [ 'xyz';  'xzy';  'yzx'];
% ixyz                            = vvv(vNo,  :);
% oxyz                            = vvv(oNo,  :);
% 
% osz                             = [isz(ixyz==oxyz(1)),  isz(ixyz==oxyz(2)), isz(ixyz==oxyz(3))];
% xyz0                            = [g1(:),  g2(:),   g3(:)];
% zzz                             = xyz0(:, ixyz==oxyz(3));
% %
% oM                              = zeros(osz(1).*osz(2),     osz(3));
% for rz=1:1:osz(3);              
%     k                           = find(zzz==rz);
%     [v, is]                     = sort(xyz0(k,  ixyz==oxyz(2)));
%     oM(:,   rz)                 = iM(xyz2n(xyz0(k(is), :), isz));                                   end;
% %
% if length(size(iM))>2;          oM                          = reshape(oM, osz);                     end;
% return;
% %%
% 
% function                        local_fls(ifl,iview,oview,  ofl, flpval);
% %%
% [idx, inm, iex]                 = fileparts(ifl);
% if isempty(ofl);                ofl                         = fullfile(idx, [inm,'_',viw,iex]);     end;
% if strcmpi(iex,'.nii');
%     v0                          = spm_vol(ifl);
%     vM                          = zeros(v0.dim);
%     vM(:)                       = spm_read_vols(v0);
%     [oM, osz]                   = transViews(vM,v0.dim,iview,    'ovw',oview);
%     if any(flpval>0);           oM(:)                       = flipXYZ(oM,[],flpval);                end;
%     vsz                         = local_szs(sqrt(sum(v0.mat(1:3,1:3).^2,1)),    iview, oview);
%     v1                          = dimmat(osz,   vsz);
%     v1.fname                    = ofl;
%     v1                          = spm_create_vol(v1);
%     v1                          = spm_write_vol(v1, oM);
% elseif strcmpi(iex,'img');
%     disp('.problem! .img are not supported any more for transViews.m');
%     return;
% else;
%     [isz, vsz, d0]              = gei(ifl,                  'imagesize','voxelsize','dataInfo');
%     vM                         	= zeros(isz(1).*isz(2),     isz(3));
%     %
%     osz                         = local_szs(isz, iview, oview);
%     ovs                         = local_szs(vsz, iview, oview);
%     %
%     si                          = struct('h2s',32,'c',mfilename,'p',ifl,'cp','m');
%     [oH, ii]                    = um_save(ofl,[],si,[],             ....
%                                 'imagesize',                osz,    ...
%                                 'voxelsize',                ovs,    ...
%                                 'orientation',              'tra [1,1,1]');
% 
%     dx                          = zeros(size(d0,1),         1);
%     do                          = zeros(size(d0,1),         1);
%     iM                          = zeros(isz(1).*isz(2),     isz(3));
%     oM                          = zeros(osz(1).*osz(2),     osz(3));
% 
%     for i=1:1:size(d0,1);
%         iM(:)                   = ged(ifl,                  i);
%         oM(:)                   = transViews(iM,isz,iviwe,  'ovw',oview);
%         if any(flpval>0);     	oM(:)                       = flipXYZ(oM,[],flpval);                end;
%         [dx(i,:), do(i,:)]      = um_save(oH,oM,si.h2s,     []);                                    end;
%     %
%     um_save(oH,ii,dx,do);                                                                           end;
% %
% disp('.done! (re-oriented images)');
% disp([' output: ',ofl]);
% return;
% %%
% 
% function    ooo                 = local_szs(iii, iview, oview);
% %% reordering images / voxel sizes:
% tcs                             = {'tra','cor','sag'};
% vNo                             = umo_cstrs(char(tcs),iview,    'im1');
% oNo                             = umo_cstrs(char(tcs),oview,    'im1');
% %
% vvv                             = [ 'xyz';  'xzy';  'yzx'];
% ixyz                            = vvv(vNo,  :);
% oxyz                            = vvv(oNo,  :);
% 
% ooo                             = [iii(ixyz==oxyz(1)),  iii(ixyz==oxyz(2)), iii(ixyz==oxyz(3))];
% return;
