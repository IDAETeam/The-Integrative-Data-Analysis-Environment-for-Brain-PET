function    ok                  = mv2_check_segetc(vmo,fbc); 

% To check if MRI-coreg and GM segmentation are approved 
%       
%       usage:  ok              = mv2_check_segetc(vmo,fbc)
%       
% Inputs:
%   vmo : fff(i) of the following:     
%       fff                     = feval(mv2_pmp2code(g4iv2.xxx(1).pmp),'vmo',[],[],[]);
%       i                       = umo_cstrs(upper(char(fff.voi_flag)),'vfg ', 'im1');
%   fbc : [1, subject#, scan#] of an IDAE session.
% Outputs:
%   ok  : 2 if both are approved (or not requred to approve)
%
% (cL)2019    hkuwaba1@jhmi.edu 

margin                          = 2;
if nargin<margin;               help(mfilename);                                    return;         end;

ok                              = 0;
if ~any(vmo.seg_ok~=' ');       ok                          = 1;
else;                           [f1, g1]                   	= mv2_genfln(vmo.seg_ok, fbc);
    if g1>0;                    c1                          = umo_getptf(f1, 1,[]);
                                ok                          = any(c1(1)=='at');             
    else;                       disp('..problem! GM segmentation not approved yet');     	end;    end;
%
if ~any(vmo.m2m_ok~=' ');       ok                          = ok + 1;
else;                           [f1, g1]                   	= mv2_genfln(vmo.m2m_ok, fbc);
    if g1>0;                    c1                          = umo_getptf(f1, 1,[]);
                                ok                          = ok + double(any(c1(1)=='at'));             
    else;                       disp('..problem! MRI-coreg not approved yet');              end;    end;
return;
%%        