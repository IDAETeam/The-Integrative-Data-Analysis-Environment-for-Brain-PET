function    err = equalsizes(p);

global mM pL

err                             = norm(sum(mM(:)>p(1)) - pL);

