% Model file for reversible dimerization. 
% 
% This is an experimental model file with the aim to describe non-comsol models for urdme2.0. 
%
% 

function model  = dimerization()

% Model parameters
Mspecies   = 3;
Mreactions = 2;

% Define chemical species objects and set their properties. What properties do they need to have??
% The minimal set of properties for compatibility with the micro solver is the diffusion constant (1e-14),
% below, and the reaction radius (1e-9) below. A description of the species as a class at the Matlab level
% can/should be accompanied by the corresponding class or struct in the C code. 
%
% A critical consideration if we go with this solution is how to make it general enough for
% to support spatially variying diffusion. 
%
A = species('A',1e-14,1e-9);
B = species('B',1e-14,1e-9);
C = species('C',1e-14,1e-9);

% The URDME model struct contains a collection of species objects. 
model.species = {A,B,C};

% We can use the same style to define chemical reactions (as classes, not implemented yet)
%R1 = create_reaction('Id','R1','Reactants',{A,B},'Products',{C},'Rate',1e8,'Unit','pMps');
%R2 = create_reaction('Id','R1','Reactants',{C},'Products',{A B},'Rate',1e8,'Unit','pMps');
%model.reactions = {R1,R2};
kr = 1e-14;
kd = 0.05;  % Dissociation rate C-> A+B
model.urdme.parameters = [kr kd]';

% Stoichiometry matrix (Can be embedded in the reaction object). 
model.urdme.N=sparse([ -1   1; 
                       -1   1;
                        1  -1]);    
           
% With an extended reaction desctioption, the information in G
% can be embedded in the reaction object. 
model.urdme.G = sparse(ones(Mreactions,Mspecies+Mreactions));

% Import a Gmsh mesh
model.mesh = importmesh('meshes/surfacef4.msh');

%Scale domain (This is specific to this example and how the geometry was 
%              generated)
model.mesh.p = 1e-6*model.mesh.p;

Ncells = size(model.mesh.p,2);

% Assemble (we need better code for this step...)
[model.urdme.vol,model.urdme.D] = assemble2d(model);
                  
% Set initial populations (we should write wrappers/utility functions 
% for this!)   
Atot = 100;
Btot = 100;
Ctot = 0;

model.urdme.u0 = zeros(Mspecies,Ncells);
ind1  = discreternd(model.urdme.vol,Atot);
for i=1:Atot
  model.urdme.u0(1,ind1(i)) = model.urdme.u0(1,ind1(i))+1;   
end

ind2  = discreternd(model.urdme.vol,Btot);
for i=1:Btot
  model.urdme.u0(2,ind2(i)) = model.urdme.u0(2,ind2(i))+1;   
end

ind3  = discreternd(model.urdme.vol,Ctot);
for i=1:Ctot
   model.urdme.u0(3,ind3(i)) = model.urdme.u0(3,ind3(i))+1;   
end

% Data should be given a default value! 
model.urdme.data  = zeros(0,Ncells);
model.urdme.tspan = 0:100;
model.urdme.report = 0;
model.urdme.sd = ones(1,Ncells);

end

function x = discreternd(p,n)

if nargin == 1
    n = 1;
end
r = rand(1,n);
k = length(p);
cmp = cumsum(p)/sum(p);
x = zeros(1,n);

for i=1:n
    x(i)=find(cmp >= r(i),1);
end

end
