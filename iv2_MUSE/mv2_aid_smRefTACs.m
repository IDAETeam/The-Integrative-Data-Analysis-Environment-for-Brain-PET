function    mv2_aid_smRefTACs(i1,iii,ooo,fbc); 

% To perform RTMs with smoothed reference region TACs 
%       
%       usage:      mv2_aid_smRefTACs()
%       
% Options:      
% 
% (cL)2020    hkuwaba1@jhmi.edu 

margin                          = 4;
if nargin<margin;               help(mfilename);                                    return;         end;

feval(['local_',lower(i1)],     iii,ooo,fbc);
return;
%%

function                        local_run(iii,ooo,fbc);
%%
% #1   res/*ifc_*eza.eza
% #2   res/*ifc_*eza_ok.txt
% #3   mpe/r4RTMs_*eza.mat
% $1   mpe/r4RTMs_*eza_smRefTACs.mat
% $2   res/*ifc_*eza_RTMs_smRefTACs_done.m
global g4iv2
% disp(char(iii))
dti                             = mv2_get_dnum(iii);
dto                             = mv2_get_dnum(ooo);
if dto(1)>dti(3);               mv2_performRTMs(iii{1},ooo{1},ooo(2),fbc);          return;         end;
load(iii{3});
for i=1:1:numel(s4mpe.ext_flg);
    s1                          = strfind(s4mpe.ffg(i,:),s4mpe.ext_flg{i}{3});
    ffg{i}                      = [s4mpe.ffg(i,1:s1(end)-1),'sm',s4mpe.ffg(i,s1(end):end)];
    % adding sm to the reference region string:
    s4mpe.ext_str{i}(s4mpe.ext_str{i}=='/' | s4mpe.ext_str{i}=='\')     = ' ';
    sss                         = getLseg(s4mpe.ext_str{i},[0,2]);
    sss{3}                      = ['sm',sss{3}];
    s4mpe.ext_str{i}            = sss{1};
    for j=2:1:numel(sss);       s4mpe.ext_str{i}            = [s4mpe.ext_str{i},' / ',sss{j}];      end;
    s4mpe.ext_res{i}{3}         = fullfile('res', [g4iv2.xxx(1).ifc,'_',g4iv2.xxx(1).eza,   ...
                                                            '_fit',s4mpe.ext_flg{i}{3},'.mat']);    end;
% 
s4mpe.ffg                       = deblank(char(ffg));
save(ooo{1},    's4mpe');
disp('.done! (adjusted MPE parameters)');
disp([' output: ',ooo{1}]);

mv2_performRTMs(iii{1},ooo{1},ooo(2),fbc);
return;
%%

function                        local_optk2r(iii,ooo,fbc);
%%
if numel(iii)>1;                disp('.problem! need to revise local_optk2r');      return;         end;
h                               = findobj(groot, 'Tag','iv2L2W');
if isempty(h);                                                                      return;         end;
mv2_w4L2Wguis('resetall',h(1));
set(findobj(h(1), 'Tag','L2W_gUseR0'),  'String','Select RTGAs to perform & MRTM2 for k2R');

[idx, inm]                      = fileparts(iii{1});
if exist(fullfile(idx, [inm,'_smRefTACs.mat']),'file');
    iii{2}                      = fullfile(idx, [inm,'_smRefTACs.mat']);                            end;
add                             = {'', '- smRefTACs'};
%
x                               = [];
s1                              = [];
s2                              = [];
for j=1:1:numel(iii);
    clear x s1 s2;
    x                        	= load(iii{j});
    rr                          = [j,0; zeros(size(x.s4mpe.ffg,1),  2)];
    mm                          = rr;
    rtga                       	= umo_cstrs(char(x.s4mpe.ext_mfg),'RTGA ',  'im1');
    if ~rtga(1);
        set(findobj(h(1), 'Tag','L2W_gUseR0'),  'String','Not ready',   'BackgroundColor',iv2_bgcs(11));
        pause(0.5);
        mv2_w4L2Wguis('resetall',h(1));                                         	return;         end;
    %
    s1{1}                    	= ['RTGAs ',add{j},' (selected = *)'];
    ic                         	= 1;
    for i=rtga;                    
        if ~isempty(x.s4mpe.ext_res{i}{4})
                                ic                          = ic + 1;
                                rr(i+1, :)                  = [j, i];
                                s1{ic}                      = [' ',x.s4mpe.ext_str{i}];   	end;    end;
    %
    set(findobj(h(1), 'Tag',['L2W_gUseR1C',int2str(j)]),    'Value',1,  'Style','popupmenu', 	...
        'String',s1,    'CallBack',['s = get(gco,''String''); v = get(gco,''Value''); ',       	...
        'if s{v}(1)=='' ''; s{v}(1)=''*''; set(gco,''String'',s); elseif s{v}(1)==''*''; ',    	...
        's{v}(1)='' ''; set(gco,''String'',s); end'],   'UserData',rr(rr(:,1)>0,:));
    %
    s2{1}                    	= ['MRTM2 ',add{j},' (selected = *)'];
    mrtm2                     	= umo_cstrs(char(x.s4mpe.ext_mfg),'MRTM2 ',  'im1');
    if ~mrtm2(1);
        set(findobj(h(1), 'Tag','L2W_gUseR0'),  'String','Not ready',   'BackgroundColor',iv2_bgcs(11));
        pause(0.5);
        mv2_w4L2Wguis('resetall',h(1));                                         	return;         end;
    
    ic                       	= 1;
    for i=mrtm2;
        if strcmpi(x.s4mpe.ext_flg{i}{4},'optk2R');
                                ic                          = ic + 1;
                                mm(i+1, :)                  = [j, i];
                                s2{ic}                      = [' ',x.s4mpe.ext_str{i}];   	end;    end;
    % 
    set(findobj(h(1), 'Tag',['L2W_gUseR2C',int2str(j)]),    'Value',1,  'Style','popupmenu', 	...
        'String',s2,    'CallBack',['s = get(gco,''String''); v = get(gco,''Value''); ',       	...
        'if s{v}(1)=='' ''; s{v}(1)=''*''; set(gco,''String'',s); elseif s{v}(1)==''*''; ',    	...
        's{v}(1)='' ''; set(gco,''String'',s); end'],   'UserData',mm(mm(:,1)>0, :));               end;
%
set(findobj(h(1), 'Tag','L2W_gUseR1C3'),    'String','Cancel',  ...
                                'CallBack','mv2_w4L2Wguis(''resetall'',[]);');
set(findobj(h(1), 'Tag','L2W_gUseR2C3'),    'String','Save',   	...
    'UserData',{iii,ooo},       'CallBack','mv2_aid_smRefTACs(''save_optk2R'',[],[],[]);');
return;
%%

function                        local_save_optk2r(iii,ooo,fbc);
%%
ud                              = get(gco,  'UserData');
rr                              = [];
mm                              = [];
x                               = [];
for i=1:1:numel(ud{1});
    h1                          = findobj(gcf, 'Tag',['L2W_gUseR1C',int2str(i)]);
    sR1Cx                    	= char(get(h1,	'String'));
    udR1Cx                      = get(h1,       'UserData');
    rr                          = [rr; udR1Cx(sR1Cx(:,1)=='*',:)];
    h2                          = findobj(gcf, 'Tag',['L2W_gUseR2C',int2str(i)]);
    sR2Cx                    	= char(get(h2,	'String'));
    udR2Cx                      = get(h2,       'UserData');
    mm                          = [mm; udR2Cx(sR2Cx(:,1)=='*', :)];
    x{i}                        = load(ud{1}{i});                                                   end;
%
% using structure of s4mpe from ud{1}{1}
load(ud{1}{1});
f2n                             = {'ffg', 'ext_flg', 'ext_res', 'ext_str', 'ext_var', 'ext_mfg'};
for i=1:1:numel(f2n);           eval(['s4mpe.',f2n{i},'     = [];']);                           	end;
ic                              = 0;
for i=1:1:size(rr,1);
    for j=1:1:size(mm,1);
        ic                      = ic + 1;
        mrtm2{ic}               = x{mm(j,1)}.s4mpe.ffg(mm(j,2), :);
        s2                      = find(mrtm2{ic}=='_' | mrtm2{ic}=='.');
        s1                      = find(x{rr(i,1)}.s4mpe.ffg(rr(i,2), :)=='_');
        rtga{ic}                = [x{rr(i,1)}.s4mpe.ffg(rr(i,2), 1:s1(end)),    ...
                                                            'k2R',mrtm2{ic}(1, s2(end-3):end)];
        ext_flg{ic}             = x{rr(i,1)}.s4mpe.ext_flg{rr(i,2)};
        ext_res{ic}             = x{rr(i,1)}.s4mpe.ext_res{rr(i,2)};
        s3                      = find(x{rr(i,1)}.s4mpe.ext_str{rr(i,2)}=='/' |         ...
                                                            x{rr(i,1)}.s4mpe.ext_str{rr(i,2)}=='\');
        ext_str{ic}             = [x{rr(i,1)}.s4mpe.ext_str{rr(i,2)}(1, 1:s3(3)),' ', 	...
                                    'k2R',mrtm2{ic}(1, s2(end-3):s2(end)-1),' ',        ...
                                     	x{rr(i,1)}.s4mpe.ext_str{rr(i,2)}(1, s3(end):end)];
        ext_var{ic}             = x{rr(i,1)}.s4mpe.ext_var{rr(i,2)};  
        ext_mfg{ic}             = x{rr(i,1)}.s4mpe.ext_mfg{rr(i,2)};                        end;    end;
%
s4mpe.ffg                       = deblank(char(rtga));
s4mpe.mrtm2                     = deblank(char(mrtm2));
for i=2:1:numel(f2n);           eval(['s4mpe.',f2n{i},'     = ',f2n{i},';']);                     	end;
% s4mpe
save(ud{2}{1},    's4mpe');
disp('.done! (parameters for individualized k2R for RTGA)');
disp([' output: ',ud{2}{1}]);
return;
%%

function                        local_run_optk2r(iii,ooo,fbc);
%%
% #1   res/*ifc_*eza.eza
% #2   res/*ifc_*eza_ok.txt
% #3   mpe/r4RTMs_*eza_RTGA_optk2R.mat
% $1   res/*ifc_*eza_RTMs_RTGA_optk2R_done.m
global g4iv2;
disp(['.working on Scan #',int2str(fbc(3)),'; Subject: ',g4iv2.yyy.snm(fbc(2),:)]);
dti                             = mv2_get_dnum(iii);
load(iii{3});
%
for i=1:1:size(s4mpe.ffg,1);    
    [f0, g0]                    = mv2_genfln(fullfile('res',deblank(s4mpe.mrtm2(i, :))),    fbc);
    f1                          = mv2_genfln(fullfile('res',deblank(s4mpe.ffg(i, :))),      fbc);
    if mv2_get_dnum({f1})>dti(1);
     	disp(['.previously done: ',s4mpe.ext_str{i}]);
    else;
        if g0>0;              	
            d0                  = ged(f0,   1);
            iok                 = 1;
            if ischar(s4mpe.ext_res{i}{3});
                                [ref, iok]                  = mv2_genfln(s4mpe.ext_res{i}{3},   fbc);
            else;               ref                         = s4mpe.ext_res{i}{3};                  end;
            if iok>0;
                RTGA4BP(iii{1}, s4mpe.ext_flg{i}{2}, f1, 'ref',ref,  'k2r',mean(d0(:, 3))); 
            else;               disp('.not ready? unable to locate file of smoothed reference region');
                                disp([' sought: ',ref]);                                            end;
        % MRTM2 not done yet:
        else;                   disp('.not ready? unable to locate input MRTM2 file');
                                disp([' sought: ',f0]);                         	end;    end;    end;
%
write2ptf(ooo{1},   'done - RTGA by optk2R of MRTM2');
return;
%%
