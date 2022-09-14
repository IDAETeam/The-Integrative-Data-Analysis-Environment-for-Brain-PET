function    [o1, o2] = um_save(i1,vM,i3,i4, varargin)
% To save data to an UMO format file.
%
%   usage 1:    To save one single data/image volume (e.g., one frame PET volume):
%
%       >>  status  = um_save(ofln,vM,si,[]     [,infoName/valuePairs] );
%       where   ...
%           si      = struct('h2s',prec#,'c',mfilename,'p',file.name,'cp','m/a');
%           prec#   -   integers for precision and scaling      (must, see below)
%           cName   -   code name (use mfilename as default)    (must)
%           pFile   -   parent file to copy info-items          (optional)
%           m/a     -   to copy all / must items only ('a'/'m') (optional)
%
%   usage 2:    To save mutiple data sets (excluding dynamic PET, VOIs) - 3 steps
%
%       step 1. To open output file (ofln) and save data info-items
%       >>  fH      = um_save(ofln,[],si,[]     [,infoName/valuePairs] );
%           See usage 1 for 'si'
%           fH(1)   -   handle of opened output UMO file
%           fH(2)   -   session ID (used for global of this session)
%
%       step 2. To save individual data sets (i-th data). 
%       >>  [dIndx(i,:), dInfo(i,:)]        = um_save(fH,vM,si.h2s,[]);
%		or
%       >>  [dIndx(i,:), dInfo(i,:)]        = um_save(fH,vM,si.h2s,[], attName/valuePairs)
%           Note that dInfo(:,9:10) are reserved for min/max values, and
%                     dInfo(:,1:5) are reserved in VOI files (See um_save)
%
%       step 3. To add info-items before/after each data set (optional)
%       >>  um_save(fH,[],[],[],      infoName/valuePairs)
%           aIndx will be n by 1 where n is the number of infoName/valuePairs
%
%       step 4. To close the file.
%       >>  status  =   um_save(fH,1,dIndex,dInfo  [,infoName/valuePairs] );
%
%   usage 3:    To save dynamic PET (.ezm) - 6 steps
%       
%       step 1. Open the output file (ofln) and save data info-items
%       >>  fH      = um_save(ofln,[],si,[]     [,infoName/valuePairs] );
%
%       step 2. Save individual frames as .mat
%      	>>  [odx, onm]      = fileparts(ofln);
%       >>  mkdir(fullfile(odx, [onm,'_ezm']));
%           for i=1:1:#ofFrames;
%               save(fullfile(odx, [onm,'_ezm'], [onm,'_frm',int2str(i),'.mat]), 'vM');
%           end;
%           'vM' is the name of the image volume matrix (any name, such as mM) 
%
%       step 4. Record the names of the image volume matrix
%           for i=1:1:#ofFrames; um_save(['!',int2str(i),'!vM!'],208,[]);    end;
%
%       step 5. To add info-items before/after each data set (optional)
%       >>  um_save(fH,[],[],[],      infoName/valuePairs)
%           aIndx will be n by 1 where n is the number of infoName/valuePairs
%
%       step 6. To close the file.
%       >>  status  =   um_save(fH,1,dIndex,dInfo  [,infoName/valuePairs] );
%
%   usage 4:    To save VOI data    >   use save2ezr.m
%
% Variables:
%   h2s:        -   integers defined for precision and scaling or not
%                   >> um_info(3,0) will display current UMO precision numbers.
%   codeName    -   m-file name where um_save is called
%                   mfilename is valid  (i.e., si.c  = mfilename)
%   pfln        -   The file from which to copy data info-items
%                   Omit 'p' (i.e., info.p) to create a file from scratch
%                   Enter 'm' to info.cp to copy 'must' items only
%                   Enter 'a' to info.cp to copy all info-items 
%   infoName    -   the name of a data info-item (a data info item)
%   infoValue   -   the value of infoName. infoName and infoValue must be given in pairs
%
% (cL)2005~16    hkuwaba1@jhmi.edu 

%%
margin                          = 4; 
if nargin<margin;               helq (mfilename);                                   return;         end;

o1                              = -1; 
o2                              = 0;

if ~isempty(vM) && isstruct(i3);
    fH                          = um_save(i1,[],i3,[],      varargin{:});
    if fH(1)<0;                 disp(['.unable to open .. ',i1]);                   return;         end;
    [di, df]                    = um_save(fH,vM,i3.h2s,     []);
    um_save(fH, 1, di, df);                                                         return;         end;

% i1 is output file name (usages 1 or 2-1):
if ischar(i1);                  

    sno                         = round(datenum(clock).*24.*3600.*100);
    sid                         = int2str(sno);
    attin                       = [];
    if ~isempty(varargin);      attin                       = local_getatt(varargin);
        if isempty(attin);                                                          return;         end;
                                                                                                    end;
    att                         = local_cpatt(i3,attin);
    if isempty(att);                                                                return;         end;
    
    % opening output umo-file:
    fH                          = um_open(i1,               'w');
    if fH<0;                                                                        return;         end;

    eval(['global umsave',sid,';']);
    eval(['umsave',sid,'        =                           att;']);
    
    o1                          = [fH; sno];
    o2                          = sno;

    % the main data is entered (save it and close the file):
    if ~isempty(vM);            [dx, di]                    = um_save(o1,vM,i3.h2s,[]);
                                o1                          = um_save(o1,1,dx,di);                  end;

    return;

% i1 is opened file handle:
else;

    jstr                        = [int2str(isempty(vM)),int2str(isempty(i3)),int2str(isempty(i4))];
    if ~isempty(which(['local_',jstr]));    
        [o1, o2]                = feval(['local_',jstr],i1,vM,i3,i4,varargin);
    else;                       o1                          = 0;
                                o2                          = zeros(1,  10);                        end;
                                                                                                    end;
return;
%% 

function    att                 = local_getatt(c);
%% sorting out info-item names and their values from input c (a cell array from varargin):
%
%   att.n   = info-item names
%   att.v   = info-item values
%
% inputs are checked against 'restricted' info-items give by: um_info(4,    1) 

    att                         = [];
    
    attOK                       = 1;
    n                           = size(c,2)./2;
    if ceil(n)~=floor(n);       disp(['attName/attValue NOT in paires']);           return;         end;
    att                         = struct('n','',            'v',[]);

    % copying info-item name/value pairs to att
    for i=1:1:n;                att(i).n                    = c{i.*2-1};
                                att(i).v                    = c{i.*2};                              end;

    attn                        = char(att.n);
    if any(attn(:,7)==' ');     disp(['Some info-item names are too short (<7)']);
                                disp(attn);
                                attOK                       = 0;                                    end;

    cm1                         = umo_cstrs(attn,[],        'cm1');
    if any(cm1(:,2)>1);         disp('Duplications in info-item names');
                                disp(attn(find(cm1(:,2)>1),:));
                                attOK                       = 0;                                    end;

    % checking if reserved attNames are entered:
    [rNames, cflg]              = um_info(4,    1);
    notuse                      = umo_cstrs(lower(attn(:,1:7)),lower(rNames(find(cflg~=1),1:7)),'im1');

    if any(notuse);             disp(['Restricted info-item entered:']);
                                cc                          = find(cflg~=1);
                                disp(rNames(cc(find(notuse)),  :));
                                att                         = [];                   return;         end;

return;
%%

function    att                 = local_cpatt(si,att);
%% copying info-items from the parent file
%
% When si.p (parent file to copy info-items) is given ...
%   1. retrieve 'patt' (info-items) of si.p
%   2. compare 'patt' against 'att' (input). 'att' overwrites 'patt'.
%   3. remove non-must items from 'patt' if ci.cp is 'm' (must items only)
%   4. add 'patt' (after steps 2 and 3) after 'att' to generate new 'att' (output)
%   5. add hitory (si.c) and pFileName ('Created') to output 'att'
%
% When si.p is not given or si.p is empty
%   1. check if all must itesm are in input 'att'
%   2. add hitory (si.c) and pFileName ('Created') to output 'att'

out                             = [];
if ~isfield(si,'h2s') | ~isfield(si,'c');
    att                         = [];
    disp(['.h2s and/or .c are missing from input ''si''']);                         return;         end;

if ~isempty(att);               attn                        = lower(str2mat(att.n));                end;
[rNames, cflg]                  = um_info(4,    1);


if isfield(si,'p') && ~isempty(si.p);
% si.cp (parent file) found:
if ~isfield(si,'cp');
    att                         = [];
    disp(['Need to specify to copy all or must items only (a/m)']);                 return;         end;
                                

if exist(si.p,'file')==2;
%% si.p is present

[idx, inm, iex]                 = fileparts(si.p);
if strcmpi(iex,'.nii'); 
    n                           = length(att);
    vx                          = spm_vol(si.p);
    v                           = vx(1);
    q77                         = char(rNames(cflg>0,:),    'spm0mat');
    if n;                       im1                         = umo_cstrs(char(att.n),q77,'im1');     end;
    att00(1).v                  = v.fname;
    att00(2).v                  = v.fname;
    att00(3).v                  = v.dim(1,  1:3);
    att00(4).v                  = sqrt(sum(v.mat(1:3,1:3).^2));
    att00(5).v                  = 'tra [1,1,1]';
    att00(6).v                  = v.descrip;
    att00(7).v                  = 'unknown';
    att00(8).v                  = [si.c,'/',v.descrip];
    att00(9).v                  = v.fname;
    att00(10).v                 = v.mat;
    if n;                       ii                          = find(im1(:,1)==0);
    else;                       ii                          = [1:1:10]';                            end;
    for i=1:1:length(ii);       att(n+i).n                  = deblank(q77(ii(i),:)); 
                                att(n+i).v                  = att00(ii(i)).v;                       end;
    
else;
%% creating the output by copying must/all info items from an existing file:

    [fH, verNo]                 = um_open(si.p,             'r');
    [apos, patt]                = um_finfo(fH,              0);

    hist0                       = um_finfo(fH,              'history');

    % removing restricted info-items from copying:
    k                           = find(cflg~=1);
    not2cp                      = umo_cstrs(lower(rNames(find(cflg~=1),1:7)),lower(patt),'im1');
    patt(find(not2cp),  1)      = '%';
    
    % removing info-items give in input 'att':
    if ~isempty(att);       
        inatt                   = umo_cstrs(attn(:,1:7),lower(patt),    'im1');
        patt(find(inatt),   1)  = '%';                                                              end;

    % removing those not permitted to copy (flag=0 in um_info(4,1)):
    if ~isempty(att);           
        n2cp                    = umo_cstrs(lower(rNames(find(~cflg),1:7)),lower(patt),'im1');
        patt(find(n2cp),    1)  = '%';                                                              end;

    % removing non-must info items when si.cp when si.cp is 'm' 
    if si.cp(1)=='m';           
        musts                   = umo_cstrs(lower(rNames(find(cflg==1),1:7)),lower(patt),'im1');
        patt(find(~musts),  1)  = '%';                                                              end;

    % copying the remaing info-items
    ii                          = find(patt(:,1)~='%');
    n                           = length(att);
    for i=1:1:length(ii);       [att(n+i).n, att(n+i).v]    = um_read(fH,   apos(ii(i)),verNo);     end;
    fclose(fH);

    % adding history and parent filename:
    n                           = length(att);
    att(n+1).n                  = 'history';
    att(n+1).v                  = [si.c,'/',hist0];
    att(n+2).n                  = 'pFileName';
    att(n+2).v                  = si.p;                                                             end;
    
    return;    
                                
% unable to locate si.cp (parent file):
else;                           disp(['Unable to open   ... ',si.p]);               return;         end;
                                                                                                    end;
%% generating the file from scratch (parent file not given):

if isempty(att);                disp(['Need to enter must info-items']);
                                disp(rNames(find(cflg==1),:));                      return;         end;
musts                           = umo_cstrs(attn(:,1:7),lower(rNames(find(cflg==1),1:7)),'im1');
if any(~musts);                 disp(['Missing ''must'' info-items']);
                                disp(rNames(find(~musts),   :));
                                att                         = [];                   return;         end;
n                               = length(att);
att(n+1).n                      = 'history';
att(n+1).v                      = si.c;
att(n+2).n                      = 'pFileName';
att(n+2).v                      = 'Created';

return;
%% 
    

function    [o1, o2]            = local_111(fH,vM,i3,i4,vs);
%% when only info-items are entered (usage 2-3):

o1                              = 1;
o2                              = 0;
sid                             = int2str(fH(2));

oHs                             = fopen('all');
if length(find(oHs==fH(1)))~=1;                                                     return;         end;
if isempty(vs);                                                                     return;         end;

attin                           = local_getatt(vs);
if isempty(attin);              o1                          = 0;
                                ifln                        = fopen(fH(1));
                                disp(['Error. Closing/deleting ... ',ifln]);
                                eval(['clear global umsave',sid,';']);
                                fclose(fH(1));
                                delete(ifln);                                       return;         end;

eval(['global umsave',sid,';']);
eval(['n                        = length(umsave',sid,');']);
for i=1:1:length(attin);        eval(['umsave',sid,'(n+i)   = attin(i);']);                         end;

return;
%% 

function    [o1, o2]            = local_001(fH,vM,i3,i4,vs);
%% Saving individual data sets (i-th data) (usage 2-2)

o1                              = 0;
o2                              = zeros(1,  10);

oHs                             = fopen('all');
if length(find(oHs==fH(1)))~=1;                                                     return;         end;

% infoItem/itsValue pairs are added when saving individual data:
if ~isempty(vs);                status                      = local_111(fH,0,0,0,vs);
    if ~status;                                                                     return;         end;
                                                                                                    end;
o1                              = um_write(fH(1),           [],vM,i3(1));
p                               = find(~isnan(vM(:)) & ~isinf(vM(:))); 
o2(1,   9:10)                   = [min(vM(p)),  max(vM(p))];

return;
%%


function    [o1, o2]            = local_000(fH,notused,dIndx,dInfo,vs);
%% closing the UMO file (usage 2-4)

o1                              = 0;
o2                              = 0;

oHs                             = fopen('all');
if length(find(oHs==fH(1)))~=1;                                                     return;         end;

% infoItem/itsValue pairs are added when closing the file:
if ~isempty(vs);                [addIndx, addInfo]          = local_111(fH,0,0,0,vs);
    if ~addIndx;                                                                    return;         end;
                                                                                                    end;
[rNames, cflg]                  = um_info(4,    1);
sid                             = int2str(fH(2));
eval(['global umsave',sid,';']);
eval(['n                        = length(umsave',sid,');']);
eval(['attn                     = lower(char(umsave',sid,'.n));']);
cattn                           = ones(size(attn,1),        1);

% removing those att-names whose values are empty:
for i=1:1:n;
    if eval(['isempty(umsave',sid,'(i).v)']);               cattn(i,    :)          = 0;    end;    end;

% checking for duplications of info-item names:
cm1                             = umo_cstrs(lower(attn(end:-1:1,   1:7)),[],   'cm1');
cm1(:)                          = cm1(end:-1:1,     :);

% the last one will be taken when duplications are found;
if any(~cm1(:,2));              disp(['Warning: info-items are given more than once']);
                                disp(attn(find(cm1(:,2)>1), :));
                                disp(['Keeping the last ones']);
                                cattn(find(~cm1(:,2)),  :)  = 0;                                    end;

% Replacing input info-item names by registered ones:
iIndx                           = zeros(n,      1);
rattn                           = rNames(find(cflg),    :);
im1                             = umo_cstrs(attn(:,1:7),lower(rattn(:,1:7)),    'im1');

% just in case must + reserved info-items are not in att: 
if any(~im1);                   
    disp(['Some must/reserved info-items are missing (marked by 0)']);
    for i=1:1:size(im1,1);      disp([rattn(i,  :),' ... ',int2str(im1(i))]);                       end;
    ofln                        = fopen(fH(1));
    disp(['Aborting (no output). Closing/deleting ... ',ofln]);
    eval(['clear global umsave',sid,';']);
    fclose(fH(1));
    delete(ofln);                                                                   return;         end;

% saving 'must' items first:
for i=1:1:size(im1,1);          
    eval(['iIndx(i, :)          = um_write(fH(1),deblank(rattn(i,:)),umsave',sid,'(im1(i)).v);']);  end;
cattn(im1(:,1), :)              = 0;

% disp([attn,int2str(cattn)]);
% saving non-must items:
jj                              = find(cattn);
s                               = size(im1, 1);
for i=1:1:length(jj);
    eval(['iIndx(i+s, :)        = um_write(fH(1),umsave',sid,'(jj(i)).n,umsave',sid,'(jj(i)).v);']);end;

iIndx                           = iIndx(find(iIndx),    :);
eval(['clear global umsave',sid,';']);

% adding data index (dIndx) and dataInfo (=dInfo):
pdx                             = um_write(fH(1),           'dataIndx',dIndx);
pdi                             = um_write(fH(1),           'dataInfo',dInfo);

% saving infoIndx:
fpos                            = um_write(fH(1),           'infoIndex',[iIndx;pdx;pdi]);

% putting the shoes on:
fwrite(fH(1),                   fpos,   'uint32');

o2                              = ftell(fH(1));
o1                              = fclose(fH(1));

return;
%%