% complex-b.mp -- do not edit, generated automatically by NIR.tex
input feynmp
require_RCS_revision "1.30";
beginchar(1, 200*1.0pt#, 200*1.0pt#, 0);
"feynmf: 1";
LaTeX_unitlength:=1.0pt;
subgraph (0, 0, w, h);
vinit;
pickup pencircle scaled thin;
straight_galleries;
vleft (__i0, __i1, __i2, __id1, __id2, __i3, __i4, __i5);
vright (__o0, __o1, __o2, __od1, __od2, __o3, __o4, __o5);
vconnect ("fermion,label=$d$", __i1, __o1);
vconnect ("fermion,tension=1.5,label=$\overline {b}$", __v2, __i4);
vconnect ("fermion,label=$\overline {c}$", __o4, __v2);
vfreeze;
vconnect ("fermion", __o2, __v3, __o3);
vconnect ("fermion,label=$\overline {s}$", __o2, __v3);
vconnect ("fermion,label=$c$", __v3, __o3);
vconnect ("photon, tension=2,label=$W^{+}$", __v2, __v3);
vconnect ("phantom,tension=1.5", __i1, __v3);
vvertex ("lab=$V_{cb}^{\ast }$", __v2);
vvertex ("lab=$V_{cs}$,lab.dist=-.1w", __v3);
vfreeze;
vdraw;
endsubgraph;
endchar;
% the end.
end.
endinput;
