function rl = rltf(a,b,c,tfnt,tfdt,nstep,stsize)
%syntax rl = rltf(a,b,c,tfnt,tfdt,nstep,stsize)

%form state space model of transfer function
% fill out zero vector
nz=length(tfnt);
np=length(tfdt);
if nz==1&tfnt==0;nz=0;end % no zeros
ttz = zeros(1,np);
if nz~=0
   ttz(1:nz)=tfnt;
end
trat = ttz./tfdt;
at = diag(-1./tfdt);
bt = (ones(np,1)-trat')./tfdt';
at(2,1)=(1-trat(2))/tfdt(2);
bt(2)=bt(2)*trat(1);
for k=3:np
   at(k,k-1)=(1-trat(k))/tfdt(k);
   for kk = k-1:-1:1
      bt(k)=bt(k)*trat(kk);
   end
   for kk=k-2:-1:1
      at(k,kk)=at(k,kk+1)*trat(kk+1);
   end
end
ct(1,np)=1;
for k= np-1:-1:1
   ct(1,k)=ct(1,k+1)*trat(k+1);
end
% interconnect at,bt ct with a,b,c
na=length(a(:,1));
NumStates=na+np;
atot = zeros(NumStates);
atot(1:na,1:na)=a;
atot(na+1:NumStates,na+1:NumStates)=at;
atot(na+1:NumStates,1:na)=bt*c;
for n = 0:1:nstep-1
   atot(1:na,na+1:NumStates)=-n*stsize*b*ct;
   rl(:,n+1)=sort(eig(atot));
end
   
