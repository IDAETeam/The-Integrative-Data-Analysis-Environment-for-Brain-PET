function    mv2_do_guGUIs(i1); 

% To ... 
%       
%       usage:      mv2_do_guGUIs()
%       
% Options:      
% 
% (cL)2016    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

f2                              = findobj(0,    'Tag','iv2L2W');
if isempty(f2);                                                                     return;         end;
%
if ~isempty(which(['local_',lower(i1)])); 
                                feval(['local_',lower(i1)],f2);
else;                           disp(['.unknown usage .. ',i1]);                                    end;
return;

function                        local_s0(f2);
%% canceling L2W_gUseRiCj GUIs:
for i=1:1:3;
    for j=1:1:3;
        set(findobj(f2, 'Tag',['L2W_gUseR',int2str(i),'C',int2str(j)]),             ...
                                'Value',1,                  'Style','pushbutton',   ...
                                'String',' ',               'CallBack',' ',         ...
                                'BackgroundColor',iv2_bgcs(0));                             end;    end;
return;
%%