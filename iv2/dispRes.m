function    dispRes(i1, varargin); 

% dispRes:      To resplay results (stored in .ezd)
%       
%       usage:      dispRes(ezdfln)
%       
% Options:      
%   'dno',val   -   data # to display           default: 1
%   'prc',val   -   precision (# of digits)     default: 5
%   'rpl',val   -   to replace columns with whatever given by columns
%                   val = 'equation strings'
%       e.g., val = 'd(:,6) = d(:,1)./d(:,2); d(:,7) = d(:,6).*(1+d(:,3)./d(:,4));'
%       will replace columns 6 and 7 as in the equation string (use d).
%   'hst','on'  to add histogram                default: 'off'
%
% (cL)2004~14  hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               help dispRes;                                       return;         end;
% -----------------------------------------------------------------------------------------------------;

dnoval                          = 1;
prcval                          = 5;
rplval                          = '';
vnoval                          = [];
hstval                          = 'off';
n0                              = nargin;
um_options;
if ~OptionsOK;                                                                      return;         end;

%%
dno                             = dnoval(1);
im1                             = [];
[di, dx, vvv]                   = gei(i1,               'orientation','dataIndx','orientation');
if dno<0 | dno>size(dx,1);      
    disp(['.Not a valid data #: ',int2str(dno),'/',int2str(size(dx,1))]);
else;
    d                           = ged(i1,                   dno);

    if ~isempty(vnoval);
        vvv(vvv==',' | vvv=='[' | vvv==']')                 = ' ';
        s1                      = find(vvv(1,1:end-1)==' ' & vvv(1,2:end)~=' ') + 1;
        s2                      = find(vvv(1,1:end-1)~=' ' & vvv(1,2:end)==' ');
        for i=1:1:length(s1);   qqq{i}                      = vvv(1,    s1(i):s2(i));               end;
        im1                     = umo_cstrs(lower(char(qqq)),'voiidno ','im1');
        if im1;                 vi                          = consolidVOINos(d(:, im1),vnoval(:));
                                d                           = d(vi(vi(:,2)>0,2),    :);
        else;                   disp('''vno'' option ignored (reason: no VOIIDNo column)');         end;
                                                                                                    end;
    
    disp(i1);
    disp(['.columns are:  ',di]);
    if ~isempty(rplval);        eval(rplval);
                                disp(['column(s) replaced: ',rplval]);                              end;
    if isempty(im1);            disp(num2str(d,prcval(1))); 
    else;                       vv                      = VOIdef(d(:, im1));
                                ii                      = [1:1:size(d,2)]~=im1;
                                dispCharArrays(vv.anm,1,num2str(d(:, ii>0),prcval(1)));           	end;
    disp('.Means and SDs: ');
    disp(num2str(nanmean(d,1),prcval(1)));
    disp(num2str(nanstd(d,0,1),prcval(1)));
    disp('.Min & max values: ');
    disp(num2str(nanmin(d,[],1),prcval(1)));
    disp(num2str(nanmax(d,[],1),prcval(1)));                                                        end;
%
if ~strcmpi(hstval,'on');                                                           return;         end;
% adding histogram
vvv(vvv==']' | vvv==']' | vvv==',')                         = ' ';
vnms                            = getLseg(vvv, [0,1]);
im1                             = umo_cstrs(vnms,['VT ';'BP '], 'im1');
if ~any(im1>0);
    disp('.problem! histogram option works for VT or BP alone (not found variable list)');
                                                                                    return;         end;
ii                              = find(im1>0, 1);  
vstr                            = {'V_T (mL/mL)',   'BP_{DN} (unitless)'};
for i=1:1:size(vv.anm,1);       glb.n{i}                    = deblank(vv.snm(i, :));                end;
figure;
HistSD(d(:, im1(ii))',[],    	'glb',glb);
ylabel(vstr{ii});
set(gca,'XTickLabelRotation',90);
disp('> to chnage x-tick angle: set(gca,''XTickLabelRotation'',x);');
return;
%%