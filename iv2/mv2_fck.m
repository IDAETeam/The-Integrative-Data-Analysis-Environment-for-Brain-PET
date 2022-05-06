function    mv2_fck(fNo,sNos,pNos); 

% To update statua of input/output files for IDAE (ver.iv2)     
%       
%       usage:      mv2_fck(fNo,subj#s,p#s)
%       
%   fNo     -   figure # fo IDAE L1W
%   subj#s  -   subject #s to update 
%   p#s     -   process #s to update
%               e.g., find(g4iv2{fNo}.ppp(:,4)==i) to update processes of i-th ibase
% 
% (cL)2013    hkuwaba1@jhmi.edu 


margin                          = 3;
if nargin<margin;               help(mfilename);                                    return;         end;
%% 
%
% pNos
global g4iv2 g4dxs;
fNo                             = fNo(1);
% g4iv2.fck{1}
if isempty(g4iv2);                                                                  return;         end;
% selecting files to work (=f2w):
ccc                             = zeros(max(g4iv2.fff(:,7)),   1);
if ~pNos(1);                    pNos                        = g4iv2.ppp(:, end);                    end;
for i=pNos(:)';          
    ccc(g4iv2.fff(g4iv2.fff(:,7)==i,   6),  1)              = 1;                                	end;

% looping over subjects (=i);
for i=sNos(:)';
    % all others:
    for j=find(g4iv2.fpc(:, end)==5)';
        [idx, inm, iex]        	= fileparts(deblank(g4iv2.fls{1}(j, :)));
        g4iv2.fck{i}(j, end) 	= double( exist(fullfile(g4iv2.yyy.idx,idx, [inm,iex]),'file')==2);	end;
    % mxx & exx:
    for j=find(g4iv2.fpc(:, end)>5)';
        [idx, inm, iex]     	= fileparts(deblank(g4iv2.fls{1}(j, :)));
        eval(['g4iv2.fck{i}(j,end)  = double(exist( fullfile( deblank(g4dxs.',idx,'(i,:)),',    ...
            '[deblank(g4dxs.mid',idx(1,2:end),'(i, :)), inm, iex]), ''file'')==2);']);              end;
    % mri:
    for j=find(g4iv2.fpc(:, end)==3)';
        [idx, inm, iex]     	= fileparts(deblank(g4iv2.fls{1}(j, :)));
        g4iv2.fck{i}(j, end) 	= double( exist( fullfile( deblank(g4dxs.mri(i, :)),            ...
                                    [deblank(g4dxs.msid(i, :)), inm, iex]), 'file')==2);            end;
    % ezr:
    for j=find(g4iv2.fpc(:, end)==4)';
        [idx, inm, iex]     	= fileparts(deblank(g4iv2.fls{1}(j, :)));
        g4iv2.fck{i}(j, end) 	= double( exist( fullfile( deblank(g4dxs.ezr(i, :)),            ...
                                    [deblank(g4dxs.msid(i, :)), inm, iex]), 'file')==2);            end;
    % pet/res:
    for j=find(g4iv2.fpc(:, end)<1)';
        % pet:
        for k=find(g4iv2.fpc(j, :)==1);
            [idx, inm, iex]    	= fileparts(deblank(g4iv2.fls{k}(j, :)));
            g4iv2.fck{i}(j, k) 	= double( exist( fullfile( deblank(g4dxs.pet{k}(i, :)),       	...
                                    [deblank(g4dxs.psid{k}(i, :)), inm, iex]), 'file')==2);         end;
        for k=find(g4iv2.fpc(j, :)==2)
            [idx, inm, iex]    	= fileparts(deblank(g4iv2.fls{k}(j, :)));
            g4iv2.fck{i}(j, k) 	= double( exist( fullfile( deblank(g4dxs.res{k}(i, :)),       	...
                                    [deblank(g4dxs.psid{k}(i, :)), inm, iex]), 'file')==2); end;    end;
                                                                                                    end;
% 
% f2w                             = find(ccc>0);
% [ii, jj]                        = find(g4iv2.fpc(ccc>0,:));
% [ii, is]                        = sort(ii);
% i2                              = f2w(ii);
% jj(:)                           = jj(is);
% j2                              = jj;
% j2(jj==size(g4iv2.fpc,2))  = 1;
% 
% for i=1:1:numel(ii);
%     [idx, inm, iex]             = fileparts(deblank(g4iv2.fls{j2(i)}(i2(i), :)));
%     % to cope with mxx/exx:
%   	qNo                         = str2num(idx(2:end));
%     if g4iv2.fpc(i2(i),end)==5;
%         for j=sNos(:)';
%             g4iv2.fck{j}(i2(i),    jj(i))              ...
%                                 = exist(fullfile(g4iv2.yyy.idx,idx, [inm,iex]),'file')==2;          end;
%     % mri/ezr, including mxx/exx:
%     elseif g4iv2.fpc(i2(i),end)>2;
%         for j=sNos(:)';
%         	jdx                	= mv2_genfln(idx,[0,j,zeros(1,6)]);
%             if isempty(qNo);
%                 g4iv2.fck{j}(i2(i), jj(i))                  = double( exist(        ...
%                                 fullfile(jdx, [deblank(g4dxs.msid(j, :)),inm,iex]),'file')==2);
%             else;
%                 g4iv2.fck{j}(i2(i), jj(i))                  = double( exist(        ...
%                                 fullfile(jdx, [deblank(g4dxs.mid4pet{qNo}(j,:)),inm,iex]),'file')==2);  
%                                                                                             end;    end;
%     % pet/res:
%     else;
%         for j=sNos(:)';
%         	jdx                	= mv2_genfln(idx,[0,j,zeros(1,6)]);
%             g4iv2.fck{j}(i2(i), jj(i))                      = double( exist(        ...
%                                 fullfile(jdx, [deblank(g4dxs.psid{jj(i)}(j, :)),inm,iex]),'file')==2);
%                                                                                     end;    end;    end;
% revising .ick and .ock:
fi                              = zeros(size(g4iv2.fls{1},1),     1);
fo                              = zeros(size(fi));
for i=pNos(:)';
    fi(:)                       = 0;
    fi(g4iv2.fff(g4iv2.fff(:, 7)==i & g4iv2.fff(:, 1)==1,   6),  :)      = 1;
    fo(:)                       = 0;
    fo(g4iv2.fff(g4iv2.fff(:, 7)==i & g4iv2.fff(:, 1)==2,   6),  :)      = 1;
    for j=sNos(:)';
        % disp(int2str([i,j]));
        g4iv2.ick{j}(i, :)      = sum(g4iv2.fck{j}(fi>0,   :),1);
        g4iv2.ock{j}(i, :)      = sum(g4iv2.fck{j}(fo>0,   :),1);                         	end;    end;
return;
%%
