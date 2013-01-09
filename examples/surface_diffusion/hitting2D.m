function fem = hitting2D(fem)

Ncells = size(fem.mesh.p,2); 

Mspecies   = 1;
Mreactions = 1;

A = 1;

fem.urdme.N = sparse([-1]);
fem.urdme.G = sparse(ones(Mreactions, Mspecies+Mreactions));

% Find point closest to origo.
dofs = xmeshinfo(fem,'out','dofs');
p = dofs.coords(:,1:Mspecies:end);

fem.urdme.tspan = 0:0.25:50;
fem.urdme.data = p;

mem = find(fem.urdme.sd==1);
fem.urdme.u0 = zeros(Mspecies,Ncells);
nA = 10000;

% Sample uniformly on the subvolumes. 
vol = fem.urdme.vol(mem);
vol = vol./norm(vol,1);

sv = discreternd(vol,nA);
for i=1:nA
    fem.urdme.u0(A,mem(sv(i))) = fem.urdme.u0(A,mem(sv(i)))+1; 
end

% Suppress output. 
fem.urdme.report = 0;

end

function x = discreternd(p,n)
if nargin == 1
    n = 1;
end
k = length(p);
cmp = cumsum(p)/sum(p);
r = rand(1,n);
x = zeros(1,n);
for i=1:n
    x(i)=find(cmp >= r(i),1);
end
end