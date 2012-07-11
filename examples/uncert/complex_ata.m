% Created by Octave 3.2.4, Fri Jul 06 11:39:13 2012 CEST <rmartin@kreide>

n=3;m=5;

%S=ones(n,m)+i

for k=1:n;for l=1:m; S(k,l)=k+i*l;end;end

S

% transpose of complex number is already hermitesch, means complex conjugate

St=S'

% 
A=St*S


