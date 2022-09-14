function    vL2_plots(i1,i2); 

% To generate plots for VOILand:
%       
%       usage:      vL2_plots(fNo,input_2)
%       
%   input_2     a .xyz file (data #1 = [x,y,z] to ploe
%               a .mat file (variable name: xyz =[x,t,z]) 
% 
% (cL)2012    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               helq(mfilename);                                    return;         end;

if ~ischar(i2);                 local_xyz(gcf,gco);                                 return;         end;
if ~strcmpi(get(i1,'Name'),'VOILand');                                              return;         end;
% i1 is the figure handle:
local_file(double(i1),i2);

return;

function                        local_file(fNo,i2);
%%
%
[idx, inm, iex]                 = fileparts(i2);
if strcmpi(iex,'.mat');         load(i2); 
else;                           [isz, vsz]                  = gei(i2,       'imagesize','voxelsize');
                                xyz                         = ged(i2,       1);                     end;
if size(xyz,2)~=3;                                                                  return;         end;

rxy                             = [ 1,2,3;  1,3,2;  2,3,1];
global g4vL2;
% recalculating xyz in terms of isz/vsz of the image:
xyz(:)                          = mm2pixels(mm2pixels(xyz, isz,vsz,  'px'), ...
                                                            g4vL2{fNo}.isz,g4vL2{fNo}.vsz,'mm');
for i=1:1:3;
    set(fNo,'CurrentAxes',g4vL2{fNo}.aHs(i));
    ii                          = find(round(xyz(:,rxy(i,3)))==g4vL2{fNo}.inos(rxy(i,3)));
    if isempty(ii);             
        pH                      = plot(1,1,'r.');
        set(pH,'Tag','vL2plot', 'Visible','off');
    else;
        pH                      = plot(xyz(ii,rxy(i,1)),xyz(ii,rxy(i,2)),'r.');
        set(pH,'Tag','vL2plot');                                                                    end;
                                                                                                    end;
g4vL2{fNo}.xyz                  = xyz;
if exist('label','var');        
    g4vL2{fNo}.label            = label;
    vvv                         = 'tcs';
    for i=1:1:3;                h                           = findobj(fNo,'Tag',    [vvv(i),'Mleft']);
        if ~isempty(h);
            set(h,              'CallBack',                 'vL2_plots(0,1)');              end;    end;
                                                                                                    end;
return;
%%

function                        local_xyz(fNo,oNo);
%%
global g4vL2;
if ~isfield(g4vL2{fNo},'label');                                                    return;         end;
% deleting existing labels;
h                               = findobj(gcf,              'Tag','xyz_labels');
if ~isempty(h);                 delete(h);                                                          end;

tflg                            = get(oNo,                  'Tag');
vNo                             = find(tflg(1)=='tcs');
sss                             = [3, 2, 1];
set(fNo,                        'CurrentAxes',              g4vL2{fNo}.aHs(vNo));
q                               = find(g4vL2{fNo}.xyz(:, sss(vNo))==g4vL2{fNo}.inos(sss(vNo)));

xLim                            = get(g4vL2{fNo}.aHs(vNo),  'xLim');
yLim                            = get(g4vL2{fNo}.aHs(vNo),  'yLim');
for i=1:1:length(q);
    qH                          = text(xLim(1)+(xLim*[-1;1])./20,                   ...
                                yLim(1)+(yLim*[-1;1])./30 + (yLim*[-1;1])./14.*i,   ... 
                                g4vL2{fNo}.label{q(i)});                                        
    set(qH,                     'Tag',                      'xyz_labels',           ...
                                'FontSize',                 14,                     ...
                                'FontWeight',               'bold',                 ...
                                'Color',                    [1,1,1]);                               end;

                                
return;
