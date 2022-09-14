function                        mv2_a0(fbci); 

% To calculate/display data analysis completion status in % on subject x iBase GUI matrix 
%       
%       usage:                  mv2_a0(fbci)
%   
%  To update % completion stata on L1W
%   fbci    -   [L1W figure #, subj#, 0, iBase#]
%               Set subj# to 0 to update for all subjects
%               Set iBase# to 0 to update for all iBase
%  To update DA status of individual steps for one iBase
%   fbci    -   L2W figure# (1 by 1)
% 
% (cL)2013    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;

if size(fbci,2)==1;             local_L2W(fbci(1));                                 return;         end;

%% updating completion status in L1W GUIs:
global g4iv2;
cwUD                            = get(fbci(1),              'userData');
%
if fbci(3);                     s2s                         = fbci(2);
else;                           s2s                         = find(cwUD{3}(:,2)>0);                 end;
if any(~cwUD{3}(s2s,1));                                                            return;         end;
if fbci(4);                     b2b                         = fbci(4);
else;                           b2b                         = [1:1:max(g4iv2.ppp(:,4))];   end;
d                               = 0;
for i=b2b(:)';                      
    d(:)                        = sum(sum(g4iv2.orq(g4iv2.ppp(:, 4)==i, :)));
    ic                          = 0;
    for j=s2s(:)';
        ic                      = ic + 1;
        set(cwUD{4}(ic,i),      'String',                   [int2str(round(     ...
            sum(sum(g4iv2.ock{j}(g4iv2.ppp(:, 4)==i, :)))./d.*100)),'%']);        end;
                                                                                                    end;
return;
%%

function                        local_L2W(i1);
%% updating strings of L2W task GUIs:

cwUD                            = get(i1(1),                'userData');
% cwUD = [L1W figure#, subjct#, iBase#]:
h                               = findobj(i1(1),'String',   'Task descriptions');
if isempty(h);                  disp(['error @local_L2W@',mfilename]);              return;         end;
bb2                             = get(h,                    'userData');
global g4iv2;
ii                              = find(g4iv2.ppp(:, 4)==cwUD(3));
% ick is 1 if all input files are present, 0 otherwise:
ick                             = g4iv2.ick{cwUD(2)}(ii,:)==g4iv2.irq(ii,:);
ick(:, 1:end-1)                 = ick(:, 1:end-1).*ick(:,end+zeros(size(ick,2)-1,1));
% o2p is 1 if any output files is present, 0 otherwise:
o2p                             = g4iv2.ock{cwUD(2)}(ii,:)>0;
% ock is 1 if all output files are present, 0 otherwise:
ock                             = g4iv2.ock{cwUD(2)}(ii,:)==g4iv2.orq(ii,:) &     ...
                                                            g4iv2.orq(ii,:)>0;
%
opt                             = g4iv2.ppp(ii,    ones(1,size(ick,2)));
%
qq2                             = ones(size(ick));
% classifying status into not-ready, ready, pending, & completed:
qq2(ick>0)                      = 2;
qq2(qq2>=2 & o2p>0)             = 3;
qq2(qq2>=3 & ock>0)             = 4;
qq2(g4iv2.irq(ii,:)==0 & ock>0)                    = 4;
qq2(opt~=111 & g4iv2.orq(ii,:)==0)                 = 1;
% qq2(:)                          = qq2.*((g4iv2{cwUD(1)}.irq(ii,:) + g4iv2{cwUD(1)}.orq(ii,:))>0);
% qq2
bHs                             = bb2(1:1:length(ii),       3:end);
smk                             = '-rpcx';
for i=1:1:length(smk);          set(bHs(qq2==i),            'String',smk(i));                       end;
if length(ii)<size(bb2,1);      set(bb2(length(ii)+1,2),    'String','*** No more tasks ***');      end;
return;
%%
