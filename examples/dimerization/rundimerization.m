% Run script for the dimerization model

function [X,X2] = rundimerization(N)

% Run the model on a convoluted surface
model = dimerization();
[Ns,Ncells]=size(model.urdme.u0);
Ndofs = Ns*Ncells;
X = zeros(Ndofs,numel(model.urdme.tspan));
for i=1:N
  model = urdme(model,[],{'Propensities','dimerization','Report',0});
  X = X+model.urdme.U;
end
X = X/N;

%Run the same model on a flat, square surface with the same area
area = sum(model.urdme.vol);
L = sqrt(area);
np = 120;

model2 = abc(L,np);
h = L/np;
[Ns,Ncells]=size(model2.urdme.u0);
Ndofs = Ns*Ncells;
X2 = zeros(Ndofs,numel(model2.urdme.tspan));
for i=1:N
  model2 = urdme(model2,[],{'Propensities','dimerization','Report',0});
  X2 = X2+model2.urdme.U;
end
X2 = X2/N;

plot(model.urdme.tspan,sum(X(1:3:end,:)),'-ob',model.urdme.tspan,sum(X2(1:3:end,:)),'-+r');
legend('Folded','Flat');

disp(['Using surface area: ',num2str(area)]);
disp(['Using voxel size: ',num2str(h)]);
