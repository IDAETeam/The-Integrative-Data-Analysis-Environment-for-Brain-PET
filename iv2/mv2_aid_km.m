function    mv2_aid_km(i1,i2, i3); 

% To aid to adopt codes from km (Keisuke Matsubara) to iv2
%       
%       usage:      mv2_aid_km('fun',iii,ooo)
%       
% 
% 
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

feval(['local_',lower(i1)],i2,  i3);
return;
%%


function                        local_edit_ev2(iii,ooo);
%%  edit descriptions of iv2 packages while startIDAE_L2W is up
set(gco,    'String','Edit - Done',     'CallBack','mv2_aid_km(''edit_ev2_done'',[],[]);');
set(findobj(gcf, 'Tag','start_iv2_L2W_edits'), 'String',    ['Edit descriptions in cells ',   ...
    'right to packages (green GUIs)']);
set(findobj(gcf,  'Tag','iPacksDescrip'),   'Style','Edit');
return;
%%

function                        local_edit_ev2_done(iii,ooo);
%%  edit descriptions of iv2 packages while startIDAE_L2W is up
h1                              = findobj(gcf,  'Tag','back2iPacks');
h2                              = findobj(gcf,  'Tag','iPacksDescrip');
ud                              = get(gco,  'UserData')

for i=1:1:numel(h1);            h1i{i}                      = get(h1(i), 'UserData');              
                                h2i{i}                      = get(h2(i), 'UserData');
                                ss{i}                       = get(h2(i), 'String');                 end;
im1                             = umo_cstrs(char(h2i),char(h1i),  'im1');

t1                              = umo_getptf(ud,    1,[])
c1                              = getLseg(t1,   1);
c1p2                            = [char(zeros(size(c1,1), 2)+32), c1];
ggg                             = ' ';
for i=1:1:size(c1,1);
    if c1(i,1)=='#';            ggg(:)                      = c1(i, 2);
    else;                       c1p2(i, 1)                  = ggg;                          end;    end;
%
im2                             = umo_cstrs(char(h1i),c1p2,    'im1');
[idx, inm]                      = fileparts(ud);
ofl                             = fullfile(idx,inm,'iv2',[inm,'_ev2_',datestr(now,'mmddyyyy_HHMMSS'),'.m']);
copyfile(ud,    ofl);
disp('.saving the package list file, just in case');
disp([' file: ',ofl]);
fH                              = fopen(ud, 'w');
disp('.wrting the following lines to the list file:');
for i=1:1:size(t1,1);
    if im2(i)<1;                fwrite(fH,  [deblank(t1(i, :)),10], 'char');
                                disp([' ',t1(i,:)]);
    else;                       disp([c1(i, :),' ',ss{im1(im2(i))}]);
                                fwrite(fH,  [c1(i, :),' ',ss{im1(im2(i))},10],  'char');    end;    end;
fclose(fH);
%      
disp('.done! (file of the analysis package list)');
disp([' output: ',ud]);
return;
%%