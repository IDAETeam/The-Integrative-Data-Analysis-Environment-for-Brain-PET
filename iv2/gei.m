function    varargout = gei(i1,varargin)

% gei:      To display info-items on screen (for UMO or spm.img/spm.nii)
%       
%       usage:      gei(iflename)
%                   or
%                   [v1,v2, ..., vn]    = gei(iflname,n1,n2, ... , nn);
%
%   n1,n2,  -   info-item names
%   v1,v2,  -   values of corresponding info-items
% 
% (cL)2008    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

%% when the input is given in a structure array ... multiple inputs:
if isstruct(i1);
    if ~isfield(i1,'name');     disp('Use .name filed when entering files in a structure arraty');
                                return;                                                             end;
    for i=1:1:length(i1);
        if exist(i1(i).name,'file');   
                                gei(i1(i).name);
        else;                   disp(['.error! unable to locate: ',i1(i).name]);                  	end;
                                                                                                    end;
    return;                                                                                         end;
% when * is included:
if ~isempty(find(i1=='*'));
    [idx, inm, iex]             = fileparts(i1);
    fls                         = dir(fullfile(idx,     [inm,iex]));
    for i=1:1:length(fls);
        if ~fls(i).isdir;       gei(fls(i).name);                                           end;    end;
    return;                                                                                         end;

if nargin==1 && nargout;         
    disp('Wrong usage of <<gei>>! No outputs if file name alone is given.');        return;         end;
if nargin>1;                    inq                         = char(varargin);
else                            inq                         = [];                                   end;

[idx, inm, iex]                 = fileparts(i1);
if strncmp(iex,'.img',4) | strncmp(iex,'.nii',4);                                   
    varargout                   = local_spm(i1,inq,nargout);                        return;         end;
    
verNo                           = um_finfo(i1,  1);
if isempty(verNo);              varargout{1}                = [];                   return;         end;
if verNo>1 && verNo<4;
    fH                       	= um_open(i1,               'r');
    [foot, inq7]              	= um_finfo(fH,              0);
    no                          = size(inq, 1);
    if isempty(foot);           
        for i=1:1:numel(varargin);
                                varargout{i}                = [];                                   end;
                                                                                    return;         end;
    if isempty(inq);           	inq                         = inq7;
                                no                          = 0;                                    end;
    im1                       	= umo_cstrs(lower(inq7),lower(inq(:,1:7)),  'im1');
    for i=1:1:size(inq,1);
        if im1(i,   1)>0;     	[inq2{i}, varargout{i}]  	= um_read(fH,   foot(im1(i,1),:),3);
        else;                   inq2{i}                     = deblank(inq(i,    :));
                                varargout{i}            	= [];                                   end;
                                                                                                    end;
    fclose(fH);
    if ~no;                   	local_disp(i1,char(inq2),   varargout);
                                varargout{1}              	= [];                        	end;    end;
%                                                                              
return;
%% 


function    out                 = local_spm(i1,inq,no);
%% when input is SPM.img/SPM.nii
global defaults;
if isempty(defaults);           spm_get_defaults;                                                   end;

out                             = [];
strs                            = char('imagesize','voxelsize','description','dataType','spm_mat');
if isempty(inq);                inq                         = strs;
                                no                          = 0;                                    end;
V                               = spm_vol(i1);
if isstruct(V);
    val{1}                      = V(1).dim(1,  1:3);
    val{2}                      = sqrt(sum(V(1).mat(1:3,1:3).^2));
    val{3}                      = V(1).descrip;
    val{4}                      = spm_type(V(1).dt(1));
    val{5}                      = V(1).mat;
else;    
    val{1}                      = V.dim(1,  1:3);
    val{2}                      = sqrt(sum(V.mat(1:3,1:3).^2));
    val{3}                      = V.descrip;
    val{4}                      = spm_type(V.dt(1));
    val{5}                      = V.mat;                                                            end;

if ~no;                         local_disp(i1,inq,val);
                                out{1}                      = [];                   return;         end;
im1                             = umo_cstrs(lower(strs),lower(inq),     'im1');
for i=1:1:size(im1,1);
    if ~im1(i,  1);             out{i}                      = [];
    else;                       out{i}                      = val{im1(i,1)};                        end;
                                                                                                    end;
return;

function                        local_disp(i1,vnm,val);
%% displaying variabl names and thier values:

if size(vnm,2)>13;              vnm                         = vnm(:,    1:13);                      end;
L                               = 16;
Line                            = char(zeros(1,L) + 32);
Line(1,     L-1)                = ':';
disp(['*** ',i1,' ***']);
for i=1:1:size(vnm,1);
    Line(:, 1:L-2)              = ' ';
    ivnm                        = deblank(vnm(i,    :));
    Line(1,     L-3-size(ivnm,2)+1:L-3)                     = ivnm;
    if isnumeric(val{i}); 
        if size(val{i},1)==1 & size(val{i},2)<=5;           add                 = num2str(val{i});
        else;                   s                           = size(val{i});
                                add                         = ['(',int2str(s(1)),' by ',    ...
                                                            int2str(s(2)),'), numeric'];            end;
    else;
        if size(val{i},1)==1;   add                         = val{i};
        else;                   s                           = size(val{i});
                                add                         = ['(',int2str(s(1)),' by ',    ...
                                                            int2str(s(2)),'), characters'];         end;
                                                                                                    end;
    disp([Line,add]);                                                                               end;
return;

