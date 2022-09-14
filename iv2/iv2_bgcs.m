function    out                 = iv2_bgcs(i1); 

% To return background colors of GUIs for IDAE (iv2)
%       
%       usage:      out         = iv2_bgcs(0);          % to return the default
%                   out         = iv2_bgcs(i);          % to extract i-th set (i<=9)
%                   out         = iv2_bgcs(-9);         % to return all sets (n by 3)
%                   iv2_bgcs(-1);                       % to display colors
%   
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               helq(mfilename);                                    return;         end;

if ~isnumeric(i1);              i1                          = -1;                                   end;
if isempty(i1);                 i1                          = -1;                                   end;
i1                              = round(i1);
if ~i1;                         
    out                         = get(0,'DefaultUicontrolBackgroundColor');         return;         end;

bgcs                            = [0.9,0.9,0.6;             0.8,0.8,0;          0.4,0.7,0.9;
                                    0.8,1,1;                0.4,0.9,0.9;        0.7,0.9,0.4;
                                    0.9,0.7,0.4;            0.4,0.4,0.9;        0.7,0.4,0.9;
                                    0.8,0,0;                1,0.8,0.9;        0.4,0.8,0;
         0    0.4470    0.7410
    0.6350    0.0780    0.1840
    0.3010    0.7450    0.9330
    0.4660    0.6740    0.1880
    0.4940    0.1840    0.5560
    0.9290    0.6940    0.1250
    0.8500    0.3250    0.0980];
                                    
% bgcs(:)                         = bgcs.*0.9;
if length(i1)>1;                out                         = bgcs(i1(:),   :);     return;         end;
if i1>=1 && i1<=size(bgcs,1);   out                         = bgcs(i1,  :);         return;         end;
if i1(1)==-9;                   out                         = bgcs;                 return;         end;
% 
[fNo, bpos]                     = bWindow([], ...
                                'nrs',                      size(bgcs,  1),     ...
                                'bwd',                      200,                ...
                                'ttl',                      'iv2 BGCs');
set(fNo,                        'Toolbar',                  'none',             ...
                                'Menubar',                  'none'); 
for i=1:1:size(bgcs, 1);
    bHs                         = postJBs(fNo,              'B',bpos(i,:),[1;1]);
    set(bHs(1),                 'String',                   int2str(i),         ...
                                'CallBack', 'disp(num2str(get(gco,''BackgroundColor'')));',         ...
                                'BackgroundColor',          bgcs(i, :));                            end;
return;