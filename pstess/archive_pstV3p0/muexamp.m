G0=[87.8 -86.4;108.2 -109.6];
dyn=nd2sys(1,[75 1]);
Dyn = daug(dyn,dyn);
G=mmult(Dyn,G0);
dynk=nd2sys([75 1],[ 1 1e-5],0.7);
Dynk=daug(dynk,dynk);Kinv=mmult(Dynk,minv(G0));
wp=nd2sys([1 0.2],[0.05 1]);Wp=daug(wp,wp);
wi=nd2sys([1 0.2],[90.5 1]);Wi=daug(wi,wi);
systemnames='G Wp Wi';
inputvar = '[ydel(2);w(2);u(2)]';
outputvar = '[Wi;Wp;-G-w]';
input_to_G='[u+ydel]';
input_to_Wp='[G+w]';input_to_Wi='[u]';
sysoutname='P';
cleanupsysic='yes';
sysic;

N=starp(P,Kinv);omega = logspace(-3,3,61);
Nf = frsp(N,omega);

blk=[1 1;1 1;2 2];
[mubnds,rowd,sens,rowp,rowg]=mu(Nf,blk,'c');
