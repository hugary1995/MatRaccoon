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
    
    constraint;
    constraint_gradient;
  end
  
  methods
    
    function this = FEProblem()
      this.variables = {};
      this.materials = containers.Map();
      this.kernels = {};
      this.integrated_bcs = {};
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
      this.constraint_gradient = zeros(num_nodes*num_vars);
      this.constraint = zeros(num_nodes*num_vars, 1);
      for m = this.materials.values()
        m.setup(this);
      end
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
    
    function [C, Ceq, GC, GCeq] = computeConstraint(this, x)
      this.computeMaterials(x);
      C = [];
      GC = [];
      this.constraint = this.constraint*0;
      this.constraint_gradient = this.constraint_gradient*0;
      
      % nodal bcs over node sets
      for boundary = this.mesh.node_sets.keys()
        for b = 1:length(boundary)
          for n = this.mesh.node_sets(boundary{b})
            for i = 1:length(this.nodal_bcs)
              if ~strcmp(boundary{b}, this.nodal_bcs{i}.boundary)
                continue
              end
              this.nodal_bcs{i}.computeConstraint();
              this.nodal_bcs{i}.computeConstraintGradient();
            end
          end
        end
      end
      
      Ceq = sparse(this.constraint);
      GCeq = sparse(this.constraint_gradient);
    end
    
  end
  
end