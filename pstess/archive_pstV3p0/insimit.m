function [y,L] = insimit(a,s,g)
%simultaneous iteration on inverse a
%y is the matrix of the s dominant eigenvalues of a
%L is the matrix of s eigenvectors of a
%g is the number of guard vectors
m=s+g;
[n,dummy]=size(a);
U=eye(n,m);
V=zeros(n,m);
UH=U';
G=eye(m);
H=zeros(m);
B=zeros(m);
lambda=zeros(m);
lambda1=zeros(m);
P=eye(m);
[AL,AU]=lu(a);
for i = 0 : s-1
	for k = 1: 30
		V=AL\U;
		V=AU\V;
		G=UH*U;
		H=UH*V;
		B=G\H;
		[P,lambda1]=eig(B);
		lambda2=diag(inv(lambda1));
		[lambda2,ind]=sort(lambda2);
		lambda2=1./lambda2;
		lambda1=diag(lambda2);
		P=P(:,ind);
		U=V*P;
		for j=1:m
			U(:,j)=U(:,j)/norm(U(:,j));
		end
	UH=conj(U');
%check convergence of lambda
		lermax=0.;
		for j=i+1 : s
			ler=abs(lambda1(j,j)-lambda(j,j));
			if lambda1(j,j)>1E-6; 	
				ler=ler/abs(lambda1(j,j));
			end
			if ler>lermax;lermax=ler;end;
		end
		if lermax<1E-6 
			break;
		else
			lambda=lambda1;
		end
	end
end
y=lambda(1:s,1:s);
L=U(:,1:s);





