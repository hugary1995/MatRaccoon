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

p = FEProblem();

p.setMesh(MeshBAR2(0, 1, 100));

p.addMaterial('k', ConstantMaterial(p, 5));

p.addVariable('u');
p.addVariable('v');

p.addKernel(Diffusion(p, 'u'));
p.addKernel(Source(p, 'u', @(x, y) -x^2));
p.addKernel(Diffusion(p, 'v'));
p.addKernel(Source(p, 'v', @(x, y) -x^2+2));
p.addKernel(CoupledReaction(p, 'u', 'v'));

p.addNodalBC(ConstantBC(p, 'left', 'u', 0));
p.addNodalBC(ConstantBC(p, 'right', 'u', 1));
p.addNodalBC(ConstantBC(p, 'left', 'v', 2));
p.addNodalBC(ConstantBC(p, 'right', 'v', 3));

p.setup();
 
x = p.solve();
