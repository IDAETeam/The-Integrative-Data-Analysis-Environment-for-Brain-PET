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

h2                              = findobj(gcf,  'Tag','iPacksDescrip');
% iproj.ev2
ud                              = get(gco,  'UserData');
h0                              = gco;
set(h0, 'Enable','off', 'String','Edit');
set(h2, 'Style','pushbutton');
drawnow

for i=1:1:numel(h2);            h2i{i}                      = get(h2(i), 'UserData');               end
%
% reding package list file:
t1                              = umo_getptf(ud,    1,[])
c1                              = getLseg(t1,   1);
c1p2                            = [char(zeros(size(c1,1), 2)+32), c1];
ggg                             = ' ';
for i=1:1:size(c1,1);
    if c1(i,1)=='#';            ggg(:)                      = c1(i, 2);
                                [ss1{i}, ss2{i}]            = getLseg(t1(i, :), 1);
    elseif c1(i,1)=='%';        [ss1{i}, ss2{i}]            = getLseg(t1(i, :), 2);
                                ss1{i}                      = ['% ',ss1{i}];
    else;                       c1p2(i, 1)                  = ggg;
                                [ss1{i}, ss2{i}]            = getLseg(t1(i, :), 1);         end;    end;
%
im2                             = umo_cstrs(char(h2i),c1p2,    'im1');
for i=find(im2'>0)
    istr                        = get(h2(im2(i)), 'String');
    if istr(find(istr>32,1))=='%'
        if ss1{i}(find(ss1{i}>32,1))~='%';
                                ss1{i}                      = ['% ',ss1{i}];                        end
        ss2{i}                  = deblank(istr(1, find(istr~=' ' & istr~='%',1):end));
    else;
        ss2{i}                  = deblank(istr(1, find(istr>32,1):end));                    end;    end
%
ss1char                         = char(ss1);

tfl                             = tmpfln([], 'm');
copyfile(ud,    tfl);
disp('> the package list file copied to a tempral file:');
disp([' file: ',tfl]);
%
fH                              = fopen(ud, 'w');
disp('> making changes to the list file:');
for i=1:1:size(t1,1);           fwrite(fH,  [ss1char(i,:),'   ',ss2{i},10], 'char');
                                disp([ss1char(i,:),'   ',ss2{i}]);                                  end
%          
fclose(fH);
%      
disp('.done!');
disp([' output: ',ud]);
%
delete(tfl);
set(h0, 'Enable','on');

return;
%%