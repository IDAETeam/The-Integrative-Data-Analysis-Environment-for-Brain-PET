function    out = umo_cstrs(m1,m2,w2r); 

% umo_cstrs:    To check duplications or appearance in character matrices
%       
% Usages:   1. To return row Nos with duplications (first appearance only)
%               out     = umo_cstrs(m1,[],'idp');
%           2. out(i,1) indicates the row No where m1(i,:) first appeared
%              out(i,2) No. of duplications (first appearances only), or 0.
%               out     = umo_cstrs(m1,[],'cm1');
%           3. similar to cm1 but out(i,1) are successful integers 
%               out     = umo_cstrs(m1,m2,'cm2');
%           4. out(i,:) lists row Nos of m1 that agreed with m2(i,:)
%               out     = umo_cstrs(m1,m2,'im1');
%           5. To return outputs of 'cm2' and 'im1' in a structure array
%              (i.e., out.cm2 and out.im1)
%               out     = umo_cstrs(m1,m2,'ist');
%
% Notes:    1. Eliminate unwanted rows (such as starting with % or spaces)
%           2. When m1 and m2 differ in widths, the shoeter will be used.
%
% New!  m1/m2 could be a cell array (converted to a character array)
%       enpty cells will be replaced by spaces. 
%
% (cL)2004  hkuwaba1@jhmi.edu 

margin                          = 3;
if nargin<margin;               help umo_cstrs;                                     return;         end;
% -----------------------------------------------------------------------------------------------------;

if isempty(m1);                 out                         = zeros(size(m2,1), 1); return;         end;
if iscell(m1);
    m1c                         = m1;
    clear m1;
    for i=1:1:numel(m1c);
        if isempty(m1c{i})      m1c{i}                      = ' ';                          end;    end;
    m1                          = char(m1c);                                                        end;
if iscell(m2);
    m2c                         = m2;
    clear m2;
    for i=1:1:numel(m2c);
        if isempty(m2c{i})      m2c{i}                      = ' ';                          end;    end;
    m2                          = char(m2c);                                                        end;
if isempty(m2);                 
% -----------------------------------------------------------------------------------------------------;

    cm1                         = umoLocal_checkduplic(m1);
    % cm1 is m1L by 2:
    % cm1(:,    2) >=1 if duplications, 0 otherwise.
    % Thus, k = find(cm1(:,2)); d = find(cm1(:,1)==k(i)); will find i-th duplications

    if strncmp(lower(w2r(1,1:3)),'idp',3);                  out             = find(cm1(:,   2)>1);
    elseif strncmp(lower(w2r(1,1:3)),'cm1',3);              out             = cm1;
    elseif strcmpi(w2r,'cm2');
        out                     = cm1;
        ic                      = 0;
        for i=find(cm1(:,2)');  ic                          = ic + 1;
                                out(cm1(:,1)==cm1(i,1), 1)  = ic;                                   end;
    else;                       disp([mfilename,': Wrong ''w2o''']);                                end;


    return;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;

[cm2, m2inm1]                   = umoLocal_checkm2(m1,m2);
% cm2 is size(m2,1) by 1.       Lists Nos of m2(i,:) appeared in m1
% m2inm1 is size(m2,1) by nn.   Lists row Nos of m1 that agreed with m2(i,:).


if strncmp(lower(w2r(1,1:3)),'ist',3);
                                out.cm2                     = cm2;
                                out.im1                     = m2inm1;
elseif strncmp(lower(w2r(1,1:3)),'cm2',3);
                                out                         = cm2;
elseif strncmp(lower(w2r(1,1:3)),'im1',3);
                                out                         = m2inm1;
else;                           disp([mfilename,': Wrong ''w2o''']);                                end;


return;
% -----------------------------------------------------------------------------------------------------;


function    [o1, o2]            = umoLocal_checkm2(m1,m2);

    n                           = min([size(m1,2),          size(m2,2)]); 

%     o1                          = zeros(size(m2,1),         1);
%     o2                          = zeros(size(m2,1),         size(m1,1));
%     ddd                         = zeros(size(m1,1),         n);
%     for i=1:1:size(m2,1);       ddd(:)                      = abs(m2(i+zeros(size(m1,1),1), 1:n));
%                                 ddd(:,  1)                  = abs(abs(m1(:,1:n))-ddd)*ones(n,1);
%                                 o1(i,   :)                  = sum(ddd(:,1)==0);
        

    L                           = size(m1,  1);
    m1a                         = abs(m1(:, 1:n));
    m1d                         = zeros(size(m1a));
    m2L                         = size(m2,  1);
    m2a                         = abs(m2(:, 1:n));

    o1                          = zeros(m2L,    1);
    o2                          = zeros(m2L,    L);

    for i=1:1:m2L;              m1d(:)                      = m2a(i+zeros(L,1), :);
                                m1d(:,  1)                  = abs( m1a - m1d )*ones(n,1);
                                k                           = find(~m1d(:,1));
        if ~isempty(k);         o1(i,   :)                  = length(k);
                                o2(i,   1:length(k))        = k(:)';                        end;    end;
    % -------------------------------------------------------------------------------------------------;
    
    if max(o1);                 o2                          = o2(:,     1:max(o1));
    else;                       o2                          = o2(:,     1);                         end;


return;
% -----------------------------------------------------------------------------------------------------;



function    out                 = umoLocal_checkduplic(m1);
%%

[L, n]                          = size(m1);
m1a                             = abs(m1);
m1d                             = zeros(size(m1a));
out                             = zeros(L,  2);

for i=1:1:L;    if ~out(i,1);   m1d(:)                      = m1a(i+zeros(L,1), :);
                                m1d(:,  1)                  = abs( m1a - m1d )*ones(n,1);
                                k                           = find(m1d(:,1)==0);
                                out(k,  1)                  = i;
                                out(i,  2)                  = length(k);                    end;    end;
    
return;
%%
