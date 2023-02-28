function ret = IDAEHeaderColor(order)

% IDAEHeaderColor:	Return color matrix for IDAE GUI header
%
%		usage:	<Color matrix> = IDAEHeaderColor(<order>);
%
%
    % Color No. for iv2_bgcs
    colorNos = [1, 2, 3, 7, 9, 11];
	ret = iv2_bgcs(colorNos(order));
end
