function    varargout         = um_finfo(i1,  varargin);

% To retrieve file attributes (info items) from UMO files
% 
%   usage 1:    To display all info items on the Matlab command window
%               >> um_finfo('fileName');
% 
%   usage 2:    To retrieve values of specified attributes
%               >> [v1, .. , vn] = um_finfo('fileName','att1', .. ,'attn');
% 
%       atti    -   name of i-th attribute
%                   Try usage1 to find what attributes are available
%       vi      -   values of atti (i-th attributes)
% 
%   usage 3:    To retrieve file bit positions and names of all attributes
%               >> [foot, attName7]     = um_finfo('fileName',0);
%   
%   usage 4:    To retrieve UMO version #
%               >> verNo                = um_finfo('fileName',1);
% 
% Notes:
%   1.  'fileName' may be replaced by fH (the handle of opened file) if opened
%   2.  The first 7 characters (after converting to lower cases) must be unique.
%       e.g., abcdef7hijk and abCDEf7xyz are treated as same attribute names.
% 
% (cL)2005  hkuwaba1@jhmi.edu   

% adopted from <<getezfinfo>>

%disp(int2str(nargin));
margin                          = 1;
if nargin<margin;               help um_finfo;                                      return;         end;
varargout{1}                    = [];

%% opening the input file:
[fH, cflg]                      = um_open(i1,   'r');
verNo                           = cflg;
if ~fH;                                                                             return;         end;
if ~verNo;                      verNo                       = local_getver(fH);                     end;
p0                              = ftell(fH);

%% for usages 3 and 4 
if nargin==2 & isnumeric(varargin{1});

    flg                         = varargin{1};

    % obtaining foot (file positions of attributes = a1), attNames7 (=a2), and attName (=a3):
    if ~flg(1);                 [a1, a2, a3]                = local_foot(fH,    varargin{1},verNo);
        for i=1:1:3;            eval(['varargout{i}         = a',int2str(i),';']);                  end;
    % obtaining version #:
    else;                       varargout{1}                = verNo;                                end;

    if isempty(fopen(fH));                                                          return;         end;
    if cflg;                    fclose(fH);
    else;                       fseek(fH,                   p0,'bof');                              end;
    return;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;


%% displaying all attribute on screen 
if nargin==1;                   local_display(fH,   verNo);

%% retrieving values of specified attributes:
else;

    if nargout;
    if nargin-1~=nargout;       disp(['#out = ',int2str(nargout),'; #inq = ',int2str(nargin-1)]);
                                disp('n(v''s) must be = n(items).');                return;         end;

    vals                        = local_getatt(fH,          varargin,verNo);
    for i=1:1:nargout;          varargout(i)                = vals{i};                              end;

    else;
    vals                        = local_getatt(fH,          varargin,verNo);
    for i=1:1:nargin-1;         disp([varargin{i},': ']);
        if isnumeric(vals{i});  disp(num2str(vals{i}));
        else;                   disp(vals{i});                                              end;    end;
                                                                                                    end;
% -----------------------------------------------------------------------------------------------------;
                                                                                                    end;

if cflg;                        fclose(fH);
else;                           fseek(fH,                   p0,'bof');                              end;

return;
%% ====================================================================================================;


function                        local_display(fH,vNo);
%% displaying all attributes --------------------------------------------------------------------------;

    foot                        = local_foot(fH,    0); 
    if ~foot(1);                                                                    return;         end;

    nL                          = 10;
    line                        = char(zeros(1,nL+5) + 32);
	for i=1:1:size(foot,1);     line(:)                     = char(zeros(size(line)) + 32);
    % looping over attributes -------------------------------------------------------------------------;

                                [aName, attVal]             = um_read(fH,   foot(i),vNo);
                                L                           = min([nL, size(aName,2)]);
                                line(1, 1:L)                = aName(1,  1:L);
                                vsz                         = size(attVal);
                                vstr                        = [];

        % attVal is numerical -------------------------------------------------------------------------;
        if isnumeric(attVal);

            if vsz(1)==1;       vstr                        = num2str(attVal);
            elseif vsz(2)==1;   vstr                        = num2str(attVal');                     end;

            if length(vstr)>80; vstr                        = [];                                   end;

        % attVal is string (matrix) -------------------------------------------------------------------;
        else;

            if vsz(1)==1;       L                           = min([80, vsz(2)]);
                                vstr                        = attVal(1,     1:L);                   end;
        % ---------------------------------------------------------------------------------------------;
                                                                                                    end;

        % when attVal is too big (in terms of its size) to display:
        if isempty(vstr);       
            vstr                = ['(',int2str(vsz(1)),' by ',int2str(vsz(2)),')'];                 end;
        % ---------------------------------------------------------------------------------------------;
        
        disp([line,': ',vstr]);
    % -------------------------------------------------------------------------------------------------;
                                                                                                    end;
return;
%% ----------------------------------------------------------------------------------------------------;


function    out                 = local_getatt(fH,att,vNo);
%% to retrieve specified info items -------------------------------------------------------------------;
%length(att)
%size(att)
    
    out                         = [];
    % info7 list all attNames stored in fH:
    [foot, info7]               = local_foot(fH,    0);
    info7(:)                    = lower(info7);

    % iname lists all attName to inquire:
    iname                       = str2mat(att{:});

    % when all attNames are too short (<7):
    if size(iname,2)<7;         disp('*** To short inputs (<7)');
                                disp(iname);
                                disp('*** end of requested attName list');          return;         end;

    % finding iname in info7:
    c                           = umo_cstrs(info7,lower(iname), 'im1');

    % iname are identified uniquely:
    if size(c,2)==1;
        for i=1:1:size(c,1);    
            if c(i);            [item, ival]                = um_read(fH,   foot(c(i)),vNo);
                                out{i}                      = {ival};
            else;               disp(['*** Not found ... ',att{i}]);
                                out{i}                      = {''};                         end;    end;

    % more than 1 info7 matched for a iname(s):
    else;

        k                       = find(prod(c,   2));
        disp(['*** duplications (or too short input(s)) found']);
        for i=1:1:length(k);    disp(att{k(i)});                                                    end;
        disp('** end of the list');                                                                 end;
    % -------------------------------------------------------------------------------------------------;

return;
%% ----------------------------------------------------------------------------------------------------;


function    [out, out2, out3]   = local_foot(fH,iNo,vNo)
%% retrieving foot and info item names of UMO format files --------------------------------------------;

     p                           = ftell(fH);
   
    vNo                         = um_finfo(fH,  1);
    vIs                         = um_info(1,    0);
    prec                        = deblank(vIs.prec(vNo,   :));
    L                           = str2num(prec(find(abs(prec)>=48 & abs(prec)<=57)))./8;

    
    % The last integer (whose length is given by prec.len) (==shoes) indicates where foot starts:
    fseek(fH,                   -L,'eof'); 
    shoes                       = fread(fH,1,           prec);

    % foot (=='infoIndex') lists bit positions of attributes (thus, functions as an index):
    [footName, foot]            = um_read(fH,           shoes,vNo);
    % the file will be closed if anything went wrong:
    if isempty(fopen(fH));
       out                      = [];
       out2                     = [];
       out3                     = [];                                               return;         end;
    %
    
    n                           = size(foot,            1);

    % retrieving all atturibute names:
    % out   =   file positions of all attributes (thus indexed) ---------------------------------------;
    if ~iNo;                    
        out                     = foot;

        c                       = struct('name',[]);
        for i=1:1:n;            c(i).name               = um_read(fH,   foot(i),vNo);               end;

        out3                    = str2mat(c.name);
        out2                    = out3(:,               1:7);

    % retieving one simgle info item ------------------------------------------------------------------;
    elseif iNo>0 & iNo<=n;      [out, out2]             = um_read(fH,   foot(iNo),vNo);             end;


return;

function        verNo           = local_getver(fH)
%% obtaining version # of the opened file (=fH) -------------------------------------------------------;

    p                           = ftell(fH);

    % version statements:
    verinfo                     = um_info(1,0);
    % reading the opening statment until one of the registered version statments is met:
    c                           = zeros(size(verinfo.prec,1),   1);
    for i=1:1:size(c,1);        fseek(fH,               0,'bof');
                                adm                     = fread(fH,     3,deblank(verinfo.prec(i,:)));
        if adm(2)==1 & adm(3)==size(deblank(verinfo.strs(i,:)),2);
        % comparing the version statement of the file and registered version statements: 
                                h2r                     = um_info(2,    adm(1));
                                vstr                    = fread(fH,     adm(3),h2r.str)';
            if strncmp(lower(deblank(verinfo.strs(i,:))),lower(char(vstr)),size(vstr,2));
                                c(i,    :)              = 1;
                                break;                                                      end;    end;
    % -------------------------------------------------------------------------------------------------;
                                                                                                    end;
    % one of c must be 1 if any registered version statment is encountered:
    % c
    out                         = find(c);
    if isempty(out);            verNo                   = 0;
    else;                       verNo                   = out(end);                                 end;

return;
%%

