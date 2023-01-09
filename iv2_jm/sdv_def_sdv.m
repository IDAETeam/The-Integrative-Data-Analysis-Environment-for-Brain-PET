function 	out                 = sdv_def_sdv(i1,sdv);
% definition of critical variable for subdivision routines 
%
%       usage:      out         = sdv_def_sdv(sdv,[])
%                                 sdv_def_sdv([],sdv)
%
%   
%
%
% (cL)2018    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;
%
out                             = [];
sdvRs                           = {'p10','p06','c10'};
% VOIIDNos of subdivision VOIs
vnos                            = {vnosets('p20u'), vnosets('p20u'), vnosets('c10u')};
% VOIIDNos of source VOIs
svois                           = {[81000;82000], [81000;82000;80080], 58000};
sdvflg                          = {'Striatum (p10)','Striatum (p06)','Cingulate (c10)'};
if isempty(i1);
    if isempty(sdv);            out                         = {sdvRs,sdvflg};       return;         end;
    im1                         = umo_cstrs(char(sdvRs),sdv,    'im1');
    if im1<1;                                                                       return;         end;
    out                         = struct('sdvRs',sdvRs{im1}, 'vnos',vnos{im1},  'svois',svois{im1});
                                                                                    return;         end;
% lines to add to the ipack:
%  one line for above *** and othes for below ***
im1                             = umo_cstrs(sdvRs,[lower(i1),' '],  'im1');
if ~im1(1);                     disp('.error! unknown subdivision flag: ',i1);      return;         end;
out                             = {vnos{im1(1)},svois{im1(1)}};
out{3}                          = char('IDAE4sdv        Generation of subdivision VOIs',    ...
                                    'IDAE4sdv        iv2_genURmri',                       	...
                                    ['IDAE4sdv        iv2_sdvVR4',lower(i1)]);
return;
%%
