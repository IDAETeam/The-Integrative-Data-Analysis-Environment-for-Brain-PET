function    out                 = mv2_iv2VOIs(i1,i2); 

% To return VOIs that are marked in the 'vL2 VOI selector' window
%       
%       usage:      vnos        = mv2_iv2VOIs(fH)
%       
% 
% (cL)2015    hkuwaba1@jhmi.edu 
margin                          = 1;
if nargin<margin;               help(mfilename);                                    return;         end;
out                             = [];
if isempty(i1);
    h                           = findbyn(0,    'Tag','vL2 VOI selector');
    if isempty(h);              h                           = findbyn(0,'Name','vL2 VOI selector'); end;
    if isempty(h);              disp('.check if ''vL2 VOI selector'' is up');       return;         end;
    if length(h)>1;             delete(h(2:end));                                                   end;
    i1                          = h(1);
else;
    if ~strcmpi(get(i1(1),'Name'),'vL2 VOI selector') && ~strcmpi(get(i1(1),'Tag'),'vL2 VOI selector');
                                disp('.wrong window handle (aborting)');            return;         end;
                                                                                                    end;
%
cwUD                            = get(i1,                   'userData');
if ~isfield(cwUD,'gHs');        disp('.wrong window handle (aborting)');            return;         end;
a                               = reshape(cell2mat(get(cwUD.gHs, 'Value')), size(cwUD.gHs));
b                               = a'.*cwUD.vnos';
c                               = b(:);
out                             = c(c>0);
if nargin==2;                   x                           = out==i2(1);
                                out(:)                      = [out(x==0); out(x>0)];                end;  
return;
%%
