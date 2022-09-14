function    [out1, out2]        = mv2_genfln(i1,i2); 

% To generate file names / check their presence for IDAE (ver.iv2) 
%       
%       usage 1:    fln         = mv2_genfln(finfo,[])
% 
% To sort out theoretical file names for file declaration lines of i-files 
%   finfo   -   abc\ending_str.ext  () are allowed
%               abc\e1.ext(def/e2.ex2) > def/e1_e2.ex2
%               abc\e1.ext(def/#e2.ex2) > def/e1e2.ex2 (to remove _)
%               abc\e1abc.ext(def/#$$e2.ex2) > def/e1ae2.ex2 
%               mri/fsl.nii(#$$$bc.nii) > mri\bc.nii
%   []      -   not used for now
%   fln     -   finfo after removing (and replacing) () and globals.
%
%       usage 2:    [iii, ooo]  = mv2_genfln([fbc,g4iv2{fbc(1).ppp(p#,:)],1)
%       
% To return input and output files (to use from IDAE ver.iv2). 
%   iii     -   input files (=full/path/fileID.ext) in a cell array (e.g., iii{i})
%   ooo     -   output files in a cell array. 
% 
%       usage 3:    [idx, p]    = mv2_genfln('dxs',[fbc,g4iv2{fbc(1).ppp(p#,:)]);
%
%   input_2     1 by 8; input(2) alone is used (=subject #)
%
% To return actual full/path of pet/res/mri/ezr/etc folders
%   p       -   to indicate presence(=1)/absence(=0) of idx
%
%       usage 4:    [fln, p]    = mv2_genfln('dfg/ffg.ext',fbc);
%  
% To return full/path/file.ext (=fln). fbc must be 1 by 3.
%
%       usage 5;    dxs         = mv2_getfln('defdx',[]);
%
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
%%

if isnumeric(i1);               [out1, out2]                = local_fln(i1);        return;         end;
if numel(i2)==8;                [out1, out2]                = local_getdx(i1, i2);  return;
elseif numel(i2)==3;            [out1, out2]                = local_fff(i1,i2);     return;         end;
if ischar(i1) && strcmpi(i1,'defdx');
                                [out1, out2]                = local_defdx([],[]);   return;         end;

out1                            = [];
out2                            = 0;

i1x                             = ['(',i1,')'];
ps0                             = find(i1x=='(');
ps1                             = find(i1x==')');
if length(ps0)~=length(ps1);    disp(['.error! unequal #s of (/): ',i1]);       	return;         end;
while length(ps0)>1;
    j                           = max(find(ps0<ps1(1)));
    [id, in, ix]                = fileparts(i1x(1,   ps0(j-1)+1:ps0(j)-1));
    [jd, jn, jx]                = fileparts(i1x(1,   ps0(j)+1:ps1(1)-1));
    if isempty(id);             id                          = jd;                                   end;
    if isempty(jx);             jx                          = ix;                                   end;
    if ~isempty(jn);            
        if isempty(in);         in                          = jn;
        else;                   in                          = [in,'_',jn];                          end;
                                                                                                    end;
    i1x                         = [i1x(1:ps0(j-1)),fullfile(id, [in, jx]),i1x(ps1(1)+1:end)];       
    ps0                         = find(i1x=='(');
    ps1                         = find(i1x==')');
    if length(ps0)==1;                                                              break;          end;
                                                                                                    end;
ss                              = strfind(i1x,      '_#');
i1x(1,  [1,ss,end])             = '#';
i1                              = i1x(i1x~='#'); 
i1(i1=='/' | i1=='\')           = filesep;
if ~any(i1=='$') && isempty(i2);out1                        = i1;                   return;         end;
%
[odx, onm, oex]                 = fileparts(i1);
n                               = sum(onm=='$');
s2                              = find(onm=='$',1);
if onm(s2-1)=='_';              onm(s2-[sum(onm=='$'):-1:1]-1)  = '$';
else;                           onm(s2-[sum(onm=='$'):-1:1])    = '$';                              end;
%
out1                            = fullfile(odx, [onm(onm~='$'),oex]);
return;
%%

function    [iii, ooo]        	= local_fln(i1);
%% returning input (=iii) and output (=ooo) file names of the specified process (=i1)
global g4iv2 g4dxs;
if isempty(g4iv2) || isempty(g4dxs)                                                 return;         end;
iii                             = [];
ooo                             = [];
% for g4iv2, see iv2_memo:
% i1 = [fNo, bNo, cNo, g4iv2{fNo}.ppp, p#@iPack]
q                               = i1(3);
if q==size(g4iv2.irq,2);        q                           = 1;                                    end;
oflg                            = {'iii', 'ooo'};
for j=1:1:2;
    xxx                         = [];
    for i=find(g4iv2.fff(:,7)==i1(1,end) & g4iv2.fff(:,1)==j)';
        k                       = g4iv2.fff(i,   2);
        fNo                     = g4iv2.fff(i,   6);
        [idx, inm, iex]         = fileparts(deblank(g4iv2.fls{q}(g4iv2.fff(i,6),  :)));
        % pet/res:
        if g4iv2.fpc(fNo,end)==5;
            xxx{k}              = fullfile(local_getdx(idx,i1(1,1:3)), [inm,iex]);
            if ~any(xxx{k}=='*') && ~exist(fileparts(xxx{k}),'dir');  
                                mkdir(fileparts(xxx{k}));                                           end;
        else;
            % multi-mri:
            if g4iv2.fpc(fNo,end)>5;
                eval(['jdx     	= deblank(g4dxs.',idx,'(i1(2), :));']); 
                eval(['xxx{k}   = fullfile(jdx, [deblank(g4dxs.mid',idx(2:end),'(i1(2), :)),inm,iex]);']);
                if ~any(jdx=='*') && ~exist(jdx,'dir');     mkdir(jdx);                             end;
            % mri/ezr
            elseif g4iv2.fpc(fNo,end)>2;
                eval(['jdx    	= deblank(g4dxs.',idx,'(i1(2), :));']); 
                xxx{k}          = fullfile(jdx, [deblank(g4dxs.msid(i1(2),:)),inm,iex]);
                if ~any(jdx=='*') && ~exist(jdx,'dir');     mkdir(jdx);                             end;
            % pet/res
            else;
                eval(['xxx{k} 	= fullfile( deblank(g4dxs.',idx,'{i1(3)}(i1(2), :)), ',     ...
                                    '[deblank(g4dxs.psid{i1(3)}(i1(2),:)),inm,iex]);']);    
            if ~any(xxx{k}=='*') && ~exist(fileparts(xxx{k}),'dir');  
                                mkdir(fileparts(xxx{k}));                                           end;
                                                                                    end;    end;    end;
    %
    if ~isempty(xxx);             eval([oflg{j},'             = xxx;']);                	end;    end;
return;
%%

function    [o1, o2]          	= local_defdx(idx, fbc);
%%
o1                              = ['pet ';'res ';'mri ';'ezr '; 
                                                            'm01 ';'m02 ';'m03 '; 'e01 ';'e02 ';'e03 '];
o2                              = 'directory strings for IDAE (iv2)';
return;
%%

function    [o1, o2]            = local_getdx(idx,fbc);
%%
% attention! copied to mv2_fck.m > copy this section when revised
global g4iv2 g4dxs;
if isempty(g4iv2) || isempty(g4dxs)                                                 return;         end;
idx(:)                          = lower(idx);
im1                             = umo_cstrs(local_defdx([],[]),[idx,' '],    'im1');
if im1==0;
    o1                          = fullfile(g4iv2.yyy.idx,idx);
elseif im1>3;                   
    eval(['o1                   = deblank(g4dxs.',idx,'(fbc(2),:));']); 
else;                           
    eval(['o1                   = deblank(g4dxs.',idx,'{fbc(3)}(fbc(2),:));']);                     end;
%
if ~any(o1=='*') & ~exist(o1,'dir');                        mkdir(o1);            	               	end;
o2                              = 7;
return;
%%

function    [o1, o2]            = local_fff(fff,fbc);
%%
global g4iv2 g4dxs;
if isempty(g4iv2) || isempty(g4dxs)                                                 return;         end;
o1                              = [];
o2                              = 0;
[idx, inm, iex]                 = fileparts(fff);
%
if isempty(idx) | isempty(inm) | isempty(iex);                                      return;         end;
%
im1                             = umo_cstrs(local_defdx([],[]),[idx,' '],    'im1');
k                               = str2num(idx(1, 2:end));
if im1==0;
    o1                          = fullfile(g4iv2.yyy.idx,idx,   [inm,iex]);
elseif im1>4;
    eval(['o1                   = fullfile(deblank(g4dxs.',idx,'(fbc(2),:)),',                  ...
                                    '[deblank(g4dxs.mid',idx(1, 2:end),'(fbc(2),:)),inm,iex]);']);
elseif im1>2;
    eval(['o1                   = fullfile(deblank(g4dxs.',idx,'(fbc(2),:)),',                  ...
                                    '[deblank(g4dxs.msid(fbc(2),:)),inm,iex]);']);
elseif im1>0;
 	eval(['o1                   = fullfile(deblank(g4dxs.',idx,'{fbc(3)}(fbc(2),:)),',          ...
                                    '[deblank(g4dxs.psid{fbc(3)}(fbc(2),:)),inm,iex]);']);          end;
%
o2                              = double(exist(o1,'file')==2);
return;
%%