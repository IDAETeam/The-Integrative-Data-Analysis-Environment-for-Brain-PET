function    fs_f81Tf45(f81,f45); 

% fs_f81Tf45:   To convert Freesurfer VOIs from f81 to f45 VOIs
%       
%       usage:      fs_f81Tf45('freesurfer.ezr','out.ezr')
%       
% Options:      
% 
% (cL)2011~20    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

disp(['.entering ',mfilename,' ..']);
% look-up tabld for conversion from FS81 to FS45: 
vv                              = fs_vnos4f45(1);
[isz, d, mri]                   = gei(f81,                 	'imagesize','dataInfo', 'mri4vois'); 
mM                              = zeros(isz(1).*isz(2),     isz(3)); 

vi                              = consolidVOINos(d(:,2),    vv(:,1));
vstr                            = VOIdef(vi(:,1));
if any(~vi(:,2));               disp('.info: missing VOI(s)');
                                dispCharArrays(1, vstr.anm(~vi(:,2),:));
                                disp([' file: ',f81]);                                               end;
%
cm1                             = umo_cstrs(int2str(vv(:,2)),[],'cm1');
%
[idx, inm]                      = fileparts(f81);
[odx, onm]                      = fileparts(f45);
if ~exist(fullfile(odx,onm,'vois'),'dir');
                                mkdir(fullfile(odx,onm,'vois'));                                    end;
% saving individual VOIs:
vpw                             = [];
for i=find(cm1(:,2)>0)';        mM(:)                       = zeros(size(mM));
    for j=find(cm1(:,1)==cm1(i,1) & vi(:,2)>0)';
                                clear vpw;
                                load(fullfile(idx,inm,'vois', ['v_',int2str(vv(j)),'.mat']));
                                mM(vpw(:, 1))               = 1;                                    end;
    vpw                         = find(mM(:)>0);
    save(fullfile(odx,onm,'vois', ['v_',int2str(vv(i,2)),'.mat']),  'vpw');                         end;
% saving shared VOI file:
save2ezr(f45,f81,   'mri',mri);
disp('.done! (shared Freesurfer VOI file, v.FS45)');
disp([' output: ',f45]);  
return;