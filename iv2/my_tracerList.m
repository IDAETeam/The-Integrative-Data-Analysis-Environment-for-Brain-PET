function        ttt             = my_tracerList(i1);
% To return a n x 3 cell array of tracers missing from iv2_tracerList.m
%   mylist(i, :) = {'[isotope]tracerName', k2R, 'Target system};
%       use [11C], [15O], or [F18] for [isotope]
%
%       usage:  my_list         = my_tracerList([]); 
%

ttt                             = {
                                '[18F]D6-FP-DTBZ',      nan,            'VMAT_2'
                                '[18F]RO6958948',       [0.05;0.18],    'tau'
                                '[11C]RO6931643',       nan,            'tau'
                                '[11C]RO6924963',       nan,            'tau'
                                '[11C]RO5011232',       nan,            'mGluR5'
                                '[18F]D3-P16-129',      nan,            'amyloid-beta'
                                '[18F]AV45/[11C]PiB',   nan,            'amyloid-beta'
                                '[11C]GE179',           nan,            'BACE-1'
                                '[18F]JHUxxxxx',        nan             'omnibus'
                                '[18F]GE179/180',       nan,            'BACE-1'
                                '[11C]JHU11602',        nan,            'vasopressin'
                                '[11C]LuAE60157',       nan,            '5HT6Rs'
                                '[11C]RO0154513',       nan,            'GABA-alpha5Rs'};
