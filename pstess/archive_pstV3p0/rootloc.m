function rl = rootloc(a,b,c,nstep,sstep)
for k = 0:nstep-1
   an = a;
   an = an - sstep*k*b*c;
   rl(:,k+1)=sort(eig(an));
end