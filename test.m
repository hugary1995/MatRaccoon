clear all
close all
clc

addpath('./utils/');
addpath('./elements/');
addpath('./meshes/');
addpath('./problems/');
addpath('./kernels/');
addpath('./bcs/');

p = FEProblem();

p.setMesh(MeshBAR2(0, 1, 1000));

p.addVariable('u');
% p.addVariable('v');

p.addKernel(Diffusion(p, 'u'));
% p.addKernel(CoupledReaction(p, 'u', 'v'));
% p.addKernel(Source(p, 'u', @(x, y) -x^2));
% p.addKernel(Diffusion(p, 'v'));
% p.addKernel(CoupledReaction(p, 'v', 'u'));
% p.addKernel(Source(p, 'v', @(x, y) -x^2+2));

p.addNodalBC(ConstantBC(p, 'left', 'u', 0));
p.addNodalBC(ConstantBC(p, 'right', 'u', 1));
% p.addNodalBC(ConstantBC(p, 'left', 'v', 2));
% p.addNodalBC(ConstantBC(p, 'right', 'v', 3));

p.setup();
 
f = @(x) p.computeObjective(x);
H = @(x, lambda) p.computeHessian(x);
nonlcon = @(x) p.computeConstraint(x);

options = optimoptions(@fmincon,...
  'Algorithm', 'interior-point',...
  'SpecifyObjectiveGradient', true,...
  'SpecifyConstraintGradient', true,...
  'HessianFcn', H,...
  'Display', 'iter');

x = fmincon(f, p.solution, [], [], [], [], [], [], nonlcon, options);