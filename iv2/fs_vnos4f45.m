function    out                 = fs_vnos4f45(i1); 

% fs_vnos4f45:  Look-up table for full Freesurfer VOIs to selected 45
%       
%       usage:      out         = fs_vnos4f45(1)
%       
%   out(:,  1)  -   VOIIDNos of full Freesurfer VOIs (excluding ventricles)
%   out(:,  2)  -   VOIIDNos of the selected 45 relative to out(:,1) 
%                   (Thus, duplications)
% 
% (cL)2011    hkuwaba1@jhmi.edu 

margin                          = 1;
if nargin<margin;               helq(mfilename);                                    return;         end;

out                             = [
81100  81100
81200  81200
82100  82100
82200  82200
80180  80180
80280  80280
83100  83100
83200  83200
84100  84100
84200  84200
66100  66100
66200  66200
67100  67100
67200  67200
58565  58100
58665  58200
51150  51100
51250  51200
53500  54100
53600  54200
66500  66500
66600  66600
56100  56100
56200  56200
53160  53100
53260  53200
52160  52100
52260  52200
58800  58100
58900  58200
54190  54100
54290  54200
73190  51100
73290  51200
57100  57100
57200  57200
73140  51100
73240  51200
52150  52100
52250  52200
63100  63100
63200  63200
53800  53800
53900  53900
73500  51100
73600  51200
73800  51100
73900  51200
74500  51100
74600  51200
57800  54100
57900  54200
74100  74100
74200  74200
58170  58100
58270  58200
75100  75100
75200  75200
79100  79100
79200  79200
58515  58100
58615  58200
51115  51100
51215  51200
51110  51100
51210  51200
53110  53100
53210  53200
52110  52100
52210  52200
52175  52100
52275  52200
94100  53100
94200  53200
51120  51100
51220  51200
52120  52100
52220  52200
52500  52100
52600  52200
78100  78100
78200  78200
71000  71000
92000  92000
69000  69000];