classdef FEProblem < handle
  
  properties
    mesh;
    variables;
    materials;
    kernels;
    nodal_bcs;
    integrated_bcs;
    
    solution;
    objective;
    gradient;
    hessian;
    
    % subject to Aeq * x = Beq
    Aeq;
    Beq;
  end
  
  methods
    
    function this = FEProblem()
      this.variables = {};
      this.materials = containers.Map();
      this.kernels = {};
      this.integrated_bcs = {};
    end
    
    function x = solve(this)
      f = @(x) this.computeObjective(x);
      H = @(x, lambda) this.computeHessian(x);
      
      options = optimoptions(@fmincon,...
        'Algorithm', 'interior-point',...
        'SpecifyObjectiveGradient', true,...
        'SpecifyConstraintGradient', true,...
        'HessianFcn', H,...
        'Display', 'iter');
      
      this.computeLinearConstraints();
      
      x = fmincon(f, this.solution, [], [], this.Aeq, this.Beq, [], [], [], options);
    end
    
    function setMesh(this, mesh)
      this.mesh = mesh;
    end
    
    function addVariable(this, var)
      this.variables{end+1} = var;
    end
    
    function addMaterial(this, mat_name, mat_obj)
      this.materials(mat_name) = mat_obj;
    end
    
    function addKernel(this, kernel)
      this.kernels{end+1} = kernel;
    end
    
    function addNodalBC(this, bc)
      this.nodal_bcs{end+1} = bc;
    end
    
    function addIntegratedBC(this, bc)
      this.integrated_bcs{end+1} = bc;
    end
    
    function setup(this)
      num_nodes = length(this.mesh.nodes);
      num_vars = length(this.variables);
      this.solution = zeros(num_nodes*num_vars, 1);
      this.gradient = zeros(num_nodes*num_vars, 1);
      this.hessian = zeros(num_nodes*num_vars);
      this.Aeq = zeros(num_nodes*num_vars);
      this.Beq = zeros(num_nodes*num_vars, 1);
    end
    
    function id = getVariableId(this, variable)
      for id = 1:length(this.variables)
        if strcmp(this.variables{id}, variable)
          return
        end
      end
      error(['cannot find variable ', variable]);
    end
    
    function dof = globalDoF(this, node, var_id)
      dof = (var_id-1)*length(this.mesh.nodes)+node.id;
    end
    
    function computeMaterials(this, x)
      if norm(x-this.solution) > 0
        this.solution = x;
        for e = this.mesh.elems
          for m = this.materials.values()
            m{1}.reinitElem(e);
            m{1}.computeMaterial();
          end
        end
      end
    end
    
    function [f, g, H] = computeObjective(this, x)
      this.computeMaterials(x);
      this.objective = 0;
      this.gradient = this.gradient*0;
      this.hessian = this.hessian*0;
      
      for e = this.mesh.elems
        for i = 1:length(this.kernels)
          this.kernels{i}.reinitElem(e);
          this.kernels{i}.computeObjective();
          if nargout > 1
            this.kernels{i}.computeGradient();
            if nargout > 2
              this.kernels{i}.computeHessian();
            end
          end
        end
      end
      
      % integrated bcs over side sets
      for boundary = this.mesh.side_sets.keys()
        for b = 1:length(boundary)
          for e = this.mesh.side_sets(boundary{b})
            for i = 1:length(this.integrated_bcs)
              if ~strcmp(boundary{b}, this.integrated_bcs{i}.boundary)
                continue
              end
              this.integrated_bcs{i}.reinitElem(e);
              this.integrated_bcs{i}.computeObjective();
              if nargout > 1
                this.integrated_bcs{i}.computeGradient();
                if nargout > 2
                  this.integrated_bcs{i}.computeHessian();
                end
              end
            end
          end
        end
      end
      
      f = this.objective;
      if nargout > 1
        g = sparse(this.gradient);
      end
      if nargout > 2
        H = sparse(this.hessian);
      end
    end
    
    function H = computeHessian(this, x)
      this.computeMaterials(x);
      this.hessian = this.hessian*0;
      
      for e = this.mesh.elems
        for i = 1:length(this.kernels)
          this.kernels{i}.reinitElem(e);
          this.kernels{i}.computeHessian();
        end
      end
      
      for boundary = this.mesh.side_sets.keys()
        for b = 1:length(boundary)
          for e = this.mesh.side_sets(boundary{b})
            for i = 1:length(this.integrated_bcs)
              if ~strcmp(boundary{b}, this.integrated_bcs{i}.boundary)
                continue
              end
              this.integrated_bcs{i}.reinitElem(e);
              this.integrated_bcs{i}.computeHessian();
            end
          end
        end
      end
      
      H = sparse(this.hessian);
    end
    
    function computeLinearConstraints(this)
      this.Aeq = this.Aeq*0;
      this.Beq = this.Beq*0;
      
      % nodal bcs over node sets
      for boundary = this.mesh.node_sets.keys()
        for b = 1:length(boundary)
          for n = this.mesh.node_sets(boundary{b})
            for i = 1:length(this.nodal_bcs)
              if ~strcmp(boundary{b}, this.nodal_bcs{i}.boundary)
                continue
              end
              this.nodal_bcs{i}.computeConstraint();
            end
          end
        end
      end
    end
    
  end
  
end