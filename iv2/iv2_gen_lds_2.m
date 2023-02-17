function    iv2_gen_lds_2(lds);
% To construct a local adapter of IDAE using information gatherd with iv2_gen_lds.m
% This code will be called when the final 'Done' GUI is hit in iv2_gen_lds.m
%
% It is also possible to start from 'saved' information file 
% Let's assume the user started generation of the local adaptor file as follows: 
%       >> iv2_gen_lds('set','yourLocalAdapter')
%   And completed all steps and saved the information (hit 'Save' GUI)
%   The information file woule be yourLocalAdapter.mat
%
%       >> x = load('yourLocalAdapter.mat');
%       >> iv2_gen_lds_2(x.ud);
%
% This usage is for cases when IDAE set a new output to local adaptor files
% So, make sure to keep  yourLocalAdapter.mat
%
% (cL)2022    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               helq(mfilename);                                    return;         end;

%% idae-related
local_idx{1}                    = ['out                             = ',    ....
                                    '[''',replace(lds.idae.symbolic,'$User_ID$',''',i2,'''),'''];'];
%
%
pet_ds_c                      	= lds.pet_ds.symbolic_c;
for i=umo_cstrs(char(lds.pet_ds.symbolic_c),'$','im1');
                                pet_ds_c{i}                 = '*';                                  end;
pdx                             = pet_ds_c{1};
for i=2:1:numel(pet_ds_c);      pdx                         = fullfile(pdx, pet_ds_c{i});           end;
local_str{1}                    = ['out.pet                         = ''',pdx,''';'];
%
mri_ds_c                      	= lds.mri_ds.symbolic_c;
for i=umo_cstrs(char(lds.mri_ds.symbolic_c),'$','im1');
                                mri_ds_c{i}                 = '*';                                  end;
mdx                             = mri_ds_c{1};
for i=2:1:numel(mri_ds_c);      mdx                         = fullfile(mdx, mri_ds_c{i});           end;
local_str{2}                    = ['out.mri                         = ''',mdx,''';'];
                                                        

%% local_pet_i2d

pet_is_sc                       = lds.pet_is.symbolic_c;
for i=umo_cstrs(char(lds.pet_is.symbolic_c), '*', 'im1');
    pet_is_sc{i}                = ['$segment_',int2str(i),'$'];                                     end;
%    
pet_ds_sc                       = lds.pet_ds.symbolic_c;
pet_ds_sq                       = lds.pet_ds.symbolic_c;
ppp                             = zeros(numel(pet_ds_sc),   3);
%
im1                             = umo_cstrs(char(lds.pet_ds.symbolic_c), ['$seg';'$Sub';'$PET'], 'im1');
% 
% setting ppp(i,2) @1 if any - in mri_ds_sc{i} (i.e., -modify) & 
%   ppp(i,3)=1 if the last character of mri_ds_sc{i} is not $ (i.e., to append): 
for i=im1(im1(:)>0)';
    pet_ds_sq{i}                = '*';
    ppp(i, 2:3)               	= [any(pet_ds_sc{i}=='-'), pet_ds_sc{i}(end)~='$'];
    if ppp(i,3)>0;
       	pet_ds_sc{i}            = pet_ds_sc{i}(1, 1:find(pet_ds_sc{i}=='$',1,'last'));      
        pet_ds_sq{i}            = ['*',lds.pet_ds.symbolic_c{i}(size(pet_ds_sc{i},2)+1:end)];
    elseif ppp(i,2)>0;
        pet_ds_sc{i}            = [pet_ds_sc{i}(1, 1:find(pet_ds_sc{i}=='-',1)-1),'$'];     end;    end;
%
ppp(:,1)                        = umo_cstrs(char(pet_is_sc), char(pet_ds_sc),   'im1');
%
                                    
local_pet_i2d{1}                = ['ppp                             = [  ',int2str(ppp(1,:))];
for i=2:1:size(ppp,1);
    local_pet_i2d{i}            = [repmat(' ',1,37),int2str(ppp(i,:))];                             end;
local_pet_i2d{end}              = [local_pet_i2d{end},'];'];
%
ss                              = 'pdc         = {';
for i=1:1:size(ppp,1);          ss                          = [ss, '''',pet_ds_sq{i},''', '];       end;
ss(1, end-1:end)              	= '};';
%
local_pet_i2d{end+1}            = ss;
%
% 
if sum(ppp(:,1))>0;
    local_pet_i2d_2{1}          = 'out                             = local_pet_i2d_2(ppp,pdc,i2c);';
else;
    local_pet_i2d_2{1}          = 'out                             = i2c{1};';
    local_pet_i2d_2{2}          = 'for i=2:1:size(ppp,1);';
    local_pet_i2d_2{3}          = '    out                         = fullfile(out, i2c{i});';
    local_pet_i2d_2{4}          = 'end;';                                                           end;
%
% lines to add to local_pet_i2d_3
local_pet_i2d_3{1}              = '% not applicable to this local adapter';
ic                              = 0;
for i=find(ppp(:,1)>0 & ppp(:,2)>0)';
    ic                          = ic + 1;
    local_pet_i2d_3{ic}         = ['i                               = ',int2str(i),';'];
    ic                          = ic + 1;
    local_pet_i2d_3{ic}         = ['disp(''> enter user-defined rule for: ',    ...
                                                            lds.pet_ds.symbolic_c{i},''');']; 
    ic                          = ic + 1;
    local_pet_i2d_3{ic}         = '% pet_ds_c{i}                     = enter here';                 end;
%
%% local_pet_d2i
is_scc                          = char(lds.pet_is.symbolic_c);
im1                             = umo_cstrs(char(lds.pet_ds.symbolic_c), is_scc, 'im1');
for i=1:1:size(is_scc, 1);      is_out{i}                   = '*';                                  end;
for i=find(im1<1 & is_scc(:,1)~='*' & is_scc(:,1)~='$')';
                                is_out{i}                   = lds.pet_is.symbolic_c{i};             end;
%
local_pet_d2i{1}                = ['rbd                             = [',int2str(im1'),'];'];
rrr                             = 'is_c                            = {';
for i=1:1:numel(is_out);        rrr                         = [rrr,'''',is_out{i},''', '];          end;
rrr(1, end-1:end)               = '};';
local_pet_d2i{2}                = rrr;

%% Aiding PET source file search in the image server:

mqq                             = ones(1,   numel(lds.pet_is.search_c));
mqq(:, umo_cstrs(char(lds.pet_is.format_c(1:size(mqq,2))), '*', 'im1'))    	= 0;
%
search_q                        = lds.pet_is.search_c;
search_t                        = lds.pet_is.search_c;
for i=find(mqq>0);              
    cqq                         = repmat('*',1,size(lds.pet_is.format_c{i},2));
    dqq                         = abs(lds.pet_is.format_c{i} - lds.pet_is.real_c{i});
    cqq(dqq<1)                  = lds.pet_is.format_c{i}(dqq<1);
    % to avoid cases where -_ are missspelled in the directory:
    cqq(cqq=='-' | cqq=='_')    = '*';
    %
    tqq                         = lds.pet_is.format_c{i}(dqq>0);
    ii                          = find(dqq>0);
    cqq(ii(umo_cstrs(('ymdHMS')', tqq', 'im1')>0))          = '$';
    search_q{i}                 = cqq;
    search_t{i}                 = lds.pet_is.format_c{i}(cqq=='$');                               	end;
%
pet_seg_no_is                   = umo_cstrs(char(lds.pet_is.symbolic_c      ...
                                                            (1:numel(search_q))), '$PET', 'im1');
subject_seg_no_ds               = umo_cstrs(char(lds.pet_ds.symbolic_c), '$Sub', 'im1');
%
%
local_get_pet{1}                = ['out.mqqsymbolic_c                         = [',int2str(mqq),'];'];
s0                              = 'out.search_q                    = {';
s1                              = 'out.search_t                    = {';
for i=1:1:numel(search_q);      s0                          = [s0,'''',search_q{i},''', ']; 
                                s1                          = [s1,'''',search_t{i},''', '];         end;
s0(1,   end-1:end)              = '};';
s1(1,   end-1:end)              = '};';
local_get_pet{2}                = s0;
local_get_pet{3}                = s1;
local_get_pet{4}                = ['out.pet_seg_no_is               = ',int2str(pet_seg_no_is),';'];
local_get_pet{5}                = ['out.subject_seg_no_ds           = ',int2str(subject_seg_no_ds),';'];

%% local_mri_i2d
mri_is_sc                       = lds.mri_is.symbolic_c;
imx_mri_is_sc                   = umo_cstrs(char(lds.mri_is.symbolic_c), '*', 'im1');
for i=imx_mri_is_sc(imx_mri_is_sc(:)>0)';
    mri_is_sc{i}                = ['$segment_',int2str(i),'$'];                                     end;
%    
mri_ds_sc                       = lds.mri_ds.symbolic_c;
mri_ds_sq                       = lds.mri_ds.symbolic_c;
mmm                             = zeros(numel(mri_ds_sc),   3);
%
im1                             = umo_cstrs(char(lds.mri_ds.symbolic_c),    ...
                                                            ['$seg';'$Sub';'$MRI';'$Ser'], 'im1');
% setting mmm(i,2) @1 if any - in mri_ds_sc{i} (i.e., -modify) & 
%   mmm(i,3)=1 if the last character of mri_ds_sc{i} is not $ (i.e., to append): 
for i=im1(im1(:)>0)';
    mri_ds_sq{i}                = '*';
    mmm(i, 2:3)               	= [any(mri_ds_sc{i}=='-'), mri_ds_sc{i}(end)~='$'];
    if mmm(i,3)>0;
        mri_ds_sc{i}            = mri_ds_sc{i}(1, 1:find(mri_ds_sc{i}=='$',1,'last'));
        mri_ds_sq{i}            = ['*',lds.mri_ds.symbolic_c{i}(size(mri_ds_sc{i},2)+1:end)];
    elseif mmm(i,2)>0;
        mri_ds_sc{i}            = [mri_ds_sc{i}(1, 1:find(mri_ds_sc{i}=='-',1)-1),'$'];     end;    end;
%
mmm(:,1)                        = umo_cstrs(char(mri_is_sc), char(mri_ds_sc),   'im1');
%
%
local_mri_i2d{1}                = ['mmm                             = [  ',int2str(mmm(1,:))];
for i=2:1:size(mmm,1);
    local_mri_i2d{i}            = [repmat(' ',1,37),int2str(mmm(i,:))];                             end;
local_mri_i2d{end}              = [local_mri_i2d{end},'];'];
%
ss                              = 'mdc         = {';
for i=1:1:size(mmm,1);          ss                          = [ss, '''',mri_ds_sq{i},''', '];       end;
%
ss(1, end-1:end)              	= '};';
local_mri_i2d{end+1}            = ss;
%
%
if sum(mmm(:,1))>0;
    local_mri_i2d_2{1}          = 'out                             = local_mri_i2d_2(mmm,mdc,i2c);';
else;
    local_mri_i2d_2{1}          = 'out                             = i2c{1};';
    local_mri_i2d_2{2}          = 'for i=2:1:size(mmm,1);';
    local_mri_i2d_2{3}          = '    out                         = fullfile(out, i2c{i});';
    local_mri_i2d_2{4}          = 'end;';                                                           end;
%
%
% lines to add to local_mri_i2d_3
local_mri_i2d_3{1}              = '% not applicable to this local adapter';
ic                              = 0;
for i=find(mmm(:,1)>0 & mmm(:,2)>0)';
    ic                          = ic + 1;
    local_mri_i2d_3{ic}         = ['i                               = ',int2str(i),';'];
    ic                          = ic + 1;
    local_mri_i2d_3{ic}         = ['disp(''> enter user-defined rule for: ',    ...
                                                        lds.mri_ds.symbolic_c{i},''');'];
    ic                          = ic + 1;
    local_mri_i2d_3{ic}         = '% mri_ds_c{i}                     = enter here';                 end;

%% Aiding MRI source file search:

mqq                             = ones(1,   numel(lds.mri_is.search_c));
mqq(:, umo_cstrs(char(lds.mri_is.format_c(1:size(mqq,2))), '*', 'im1'))    	= 0;
%
search_q                        = lds.mri_is.search_c;
search_t                        = lds.mri_is.search_c;
for i=find(mqq>0);              
    cqq                         = repmat('*',1,size(lds.mri_is.format_c{i},2));
    dqq                         = abs(lds.mri_is.format_c{i} - lds.mri_is.real_c{i});
    cqq(dqq<1)                  = lds.mri_is.format_c{i}(dqq<1);
    % to avoid cases where -_ are missspelled in the directory:
    cqq(cqq=='-' | cqq=='_')    = '*';
    %
    tqq                         = lds.mri_is.format_c{i}(dqq>0);
    ii                          = find(dqq>0);
    cqq(ii(umo_cstrs(('ymdHMS')', tqq', 'im1')>0))          = '$';
    search_q{i}                 = cqq;
    search_t{i}                 = lds.mri_is.format_c{i}(cqq=='$');                               	end;
%
series_seg_no_is            	= umo_cstrs(char(lds.mri_is.symbolic_c), '$Ser', 'im1');
subject_seg_no_ds               = umo_cstrs(char(lds.mri_ds.symbolic_c), '$Sub', 'im1');
%
%
local_get_mri{1}                = ['out.mqq                         = [',int2str(mqq),'];'];
s0                              = 'out.search_q                    = {';
s1                              = 'out.search_t                    = {';
for i=1:1:numel(search_q);      s0                          = [s0,'''',search_q{i},''', ']; 
                                s1                          = [s1,'''',search_t{i},''', '];         end;
s0(1,   end-1:end)              = '};';
s1(1,   end-1:end)              = '};';
local_get_mri{2}                = s0;
local_get_mri{3}                = s1;
local_get_mri{4}                = ['out.series_seg_no_is            = ',int2str(series_seg_no_is),';'];
local_get_mri{5}                = ['out.subject_seg_no_ds           = ',int2str(subject_seg_no_ds),';'];

%% local_w2d
im2                             = umo_cstrs(char(lds.etc.fs_script_ws_c),   ...
                                                            char(lds.etc.fs_script_ds_c), 'im1');
%
if any(im2<1);
    sw                       = lds.etc.fs_script_ws_c{1};
    for i=2:1:min(im2(im2>0))-1;
        sw                   = [sw,'\',lds.etc.fs_script_ws_c{i}];                                  end;
    sd                       = lds.etc.fs_script_ds_c{1};
    for i=2:1:find(im2>0,1)-1;  
        sd                  	= [sd,'/',lds.etc.fs_script_ds_c{i}];                               end;
    local_w2d{1}                = ['sw                              = ''',sw,''';'];
    local_w2d{2}                = ['sd                              = ''',sd,''';'];
    local_w2d{3}                = 'if strncmpi(i2,sw,size(sw,2));';
    local_w2d{4}                = '    out                         = [sd,i2(1, size(sw,2)+1:end)];';
    local_w2d{5}                = '    out(out==''\'')               = ''/'';';
    local_w2d{6}                = 'else;';
    local_w2d{7}                = '    out                         = [sw,i2(1, size(sd,2)+1:end)];';
    local_w2d{8}                = '    out(out==''/'')               = ''\'';';
    local_w2d{9}                = 'end;';
else;
    local_w2d{1}                = 'out                             = i2;';                        	end;

%% fs
local_fsd{1}                    = ['out.fs.home                     = ',    ....
                                    '[''',replace(lds.etc.fs_script_ws,'$User_ID$',''',i2,'''),'''];'];
local_fsd{2}                    = ['out.fs.linux                    = ',    ....
                                    '[''',replace(lds.etc.fs_script_ds,'$User_ID$',''',i2,'''),'''];'];   
local_fsd{3}                    = ['out.fs.subj                     = ',    ...
                                    '[''',replace(lds.etc.fs_subject_ds,'$User_ID$',''',i2,'''),'''];']; 
%% local_unit
local_unit{1}                   = ['out                             = ''',lds.etc.unit,''';'];


%% generation of local adapter:
y                               = umo_getptf(which('dxetc4xxx_tmp.m'), 1, []);
y1                              = getLseg(y,    1);
y_c                             = umo_cstrs(y1, ['fun';'ret'],  'im1');
ofl                             = fullfile(fileparts(which(mfilename)), [lds.lds,'.m']);
fH                              = fopen(ofl, 'w');
for i=y_c(1, :);
    if i==1;
        fwrite(fH,  [replace(deblank(y(i, :)),'$lds$',lds.lds),10],       'char');                  end;
    %
    y(i+1, 2)                   = ' ';
    fwrite(fH,  [deblank(y(i+1, :)),10],    'char');
    if i>1;                     y(i+1, 2)                   = '%';                                  end;
    if contains(y(i+2, :),'$lds$');
        fwrite(fH,  [replace(deblank(y(i+2, :)),'$lds$',lds.lds),10],     'char');
    else;
        fwrite(fH,  [deblank(y(i+2, :)),10],    'char');                                            end;
    if i>1;                     fwrite(fH,  ['% ',10],  'char');                            end;    end;
%
fwrite(fH,  ['% (cL)',datestr(clock,'yyyy'),'    hkuwaba1@jhmi.edu ',10],   'char');
fwrite(fH,  ['    ',10],    'char');
%
for j=y_c(1,1)+2:y_c(2,1);      fwrite(fH,  [deblank(y(j, :)),10],    'char');                      end;
fwrite(fH,  ['    ',10],    'char');
ic                              = 1;
for i=y_c(1, 2:end);
    ic                          = ic + 1;
    for j=i:1:y_c(2,ic);
        if y1(j,1)=='#';        
            clear j_str;
            eval(['j_str     	= ',getLseg(y(j, :), 2),';']);
            for k=1:1:numel(j_str);
                                fwrite(fH,      [deblank(j_str{k}),10], 'char');                 	end;
        else;                   fwrite(fH,      [deblank(y(j, :)),10],  'char');            end;    end;
    fwrite(fH,  ['%% ',10,'    ',10],   'char');                                                    end;
%
fclose(fH);
disp('.done! (local adapter for youe computer environment)');
disp([' output: ',ofl]);
