clear all
close all
clc

addpath('./utils/');
addpath('./elements/');
addpath('./meshes/');
addpath('./problems/');
addpath('./kernels/');
addpath('./bcs/');
addpath('./materials/');
addpath('./constraints/');

L = 1;
nx = 100;

Gc = 2.7;
l = 0.025;
c0 = 8/3;
M = Gc/c0/l;
psic = 14.88;

p = FEProblem();

p.setMesh(MeshBAR2(0, L, nx));

p.addMaterial('1 p', ConstantMaterial(p, -1));
p.addMaterial('1 M', ConstantMaterial(p, M));
p.addMaterial('1 psic', ConstantMaterial(p, psic));
p.addMaterial('1 d', DamageProfile(p, L/2, l));
p.addMaterial('1 grad_d', DamageProfileGradient(p, L/2, l));
p.addMaterial('2 g', Degradation(p, '1 d', '1 M', '1 psic', 1));

p.addVariable('u');
p.addVariable('p');

p.setup();

p.addKernel(MatDiffusion(p, 'u', '2 g'));
% p.addKernel(Contact(p, 'u', '2 g', '1 grad_d'));
p.addKernel(CoupledPressure(p, 'u', 'p', '1 grad_d'));
p.addKernel(PenaltyContact(p, 'u', 'p', '2 g', 1));

p.addNodalBC(ConstantBC(p, 'left', 'u', 0));
p.addNodalBC(ConstantBC(p, 'right', 'u', 1));
p.addNodalBC(ConstantBC(p, 'left', 'p', 0));
p.addNodalBC(ConstantBC(p, 'right', 'p', 0));

% p.setUpperBound('p', 0);
 
x = p.solve();

p.plot('u');
p.plot('p');