function    eAT                 = interpl_ezm(i1,msiir); 

% To interpolate dynamic PET as specified 
%       
%       usage:      interpl_ezm(input_ezm,msiir)
%       
%   input_ezm   'full/path/input.ezm'
%   msiir       [mid-frame times, mean cortical TAC, frames one before /
%               one after multi-frame inputrpolation (0/1/2), one-frame
%               interpolation frames (0/1 or up), the last frame to keep (0/1 or up).
%               
% Notes:
%
% (cL)2020/21    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
%
disp(['.entering: ',mfilename]);
if sum(sum(msiir(:, 3:5)))<1;                                                       return;         end;

[isz, tim]                     	= gei(i1,                   'imagesize','PETtimes');
if size(tim,1)~= size(msiir,1) || max(abs(tim(:,1)-msiir(:,1)))>10.^-3;
    disp('.critical problems: incomparable frame time vectors');                    return;         end;

[pdx, pnm]                     	= fileparts(i1);
if ~exist(fullfile(pdx, [pnm,'_original.ezm']),'file');
    copyfile(i1, fullfile(pdx, [pnm,'_original.ezm']));                                              
    disp('> original .ezm saved for record (onece alone)');
    disp(['  to: ',fullfile(pdx, [pnm,'_original.ezm'])]);                                        	end;
%
%
iM                              = zeros(isz(1).*isz(2),     isz(3)); 
jM                              = zeros(isz(1).*isz(2),     isz(3)); 
vM                              = zeros(isz(1).*isz(2),     isz(3)); 
%
qqq                             = ones(size(msiir,1),  1);
eAT                             = msiir(:, 2);
%
% multi-frame interpolation:
ips                             = find(msiir(:,3)==1)';
ipe                             = find(msiir(:,3)==2)';
for i=1:1:length(ips);
    disp(['> interpolating frames: ',int2str(ips(i)+1),'-',int2str(ipe(i)-1)]);
    qqq(ips(i)+1:ipe(i)-1,  :)  = 0;
  	iM(:)                       = ged(i1,  ips(i));
  	jM(:)                       = ged(i1,  ipe(i));
  	w                           = interp1(msiir([ips(i),ipe(i)],1),[1;0], msiir(ips(i):ipe(i), 1));
    jc                          = ips(i):1:ipe(i);
    for j=2:1:length(w)-1;
        vM(:)                   = iM.*w(j) + jM.*(1-w(j));
        save(fullfile(pdx, [pnm,'_ezm'], [pnm,'_frm',int2str(jc(j)),'.mat']), 'vM');      	end;    end;
%
qqq(msiir(:,4)>0, :)            = 0;
eAT(qqq<1, :)                   = interp1(msiir(qqq>0,1), msiir(qqq>0,2), msiir(qqq<1,1));
v                               = [];
%
for i=find(msiir(:,4)>0)';
    disp(['> single-frame interpolation, frame: ',int2str(i)]);
    if ~exist(fullfile(pdx, [pnm,'_ezm'], [pnm,'_frm',int2str(i),'_original.mat']),'file');
        copyfile(fullfile(pdx, [pnm,'_ezm'], [pnm,'_frm',int2str(i),'.mat']),   ...
            fullfile(pdx, [pnm,'_ezm'], [pnm,'_frm',int2str(i),'_original.mat']));                  end;
    %
    v                           = load(fullfile(pdx, [pnm,'_ezm'],  ...
                                                            [pnm,'_frm',int2str(i),'_original.mat']));
    vM(:)                       = v.vM./msiir(i,2).*eAT(i);
    clear v;
    save(fullfile(pdx, [pnm,'_ezm'], [pnm,'_frm',int2str(i),'.mat']), 'vM');                        end;
%
tfl                             = tmpfln([],    'ezm');
copyfile(i1, tfl);
disp('> duplicating input .ezm, just in case');
disp(['  file: ',tfl]);
n                               = size(msiir, 1);
di                              = zeros(n,      1);
df                              = zeros(n,      10);
si                              = struct('h2s',208, 'c',mfilename, 'p',i1, 'cp','a');
fH                              = um_save(i1,[],si,[],   'sss.msiir',msiir,                     ...
                                'c1_msiir','mid-frame times',   'c2_msiir','mean cortical TAC', ...
                    'c3_msiir','frames one before / after multi-frame inputrpolation (0/1/2)',  ...
                                'c4_msiir','one-frame interpolation frames (0/1 or up)',        ...
                                'c5_msiir','the last frame to keep (0/1 or up)');         
%
if fH<0;                        disp('.error! unable to create output file');
                                disp([' output: ',i1]);                             return;         end;
%
for i=1:1:n;
  	[di(i,:), df(i,:)]          = um_save(fH,['!',int2str(i),'!vM!'],208,[]);                   	end;
um_save(fH, 1, di, df);
%
%
disp('.done! (interpolation / rescaling of dynamic PET frames)');
disp([' output: ',i1]);
delete(tfl);
return;
%%