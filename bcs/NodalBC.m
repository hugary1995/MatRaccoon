classdef (Abstract) NodalBC < handle
  
  properties
    problem;
    var_id;
    
    boundary;
    node_set;
  end
  
  methods
    
    function this = NodalBC(problem, boundary, variable)
      this.problem = problem;
      this.var_id = problem.getVariableId(variable);
      this.boundary = boundary;
      this.node_set = problem.mesh.node_sets(boundary);
    end
    
    function computeConstraint(this)
      for i = 1:length(this.node_set)
        dof = this.problem.globalDoF(this.node_set(i), this.var_id);
        this.problem.Beq(dof) = this.computeNodalConstraint(this.node_set(i));
        this.problem.Aeq(dof, dof) = 1;
      end
    end
    
    function modifyGradient(this)
      for i = 1:length(this.node_set)
        dof = this.problem.globalDoF(this.node_set(i), this.var_id);
        this.problem.gradient(dof) = -this.computeNodalConstraint(this.node_set(i));
      end
    end
    
    function modifyHessian(this)
      for i = 1:length(this.node_set)
        dof = this.problem.globalDoF(this.node_set(i), this.var_id);
        this.problem.hessian(dof, :) = 0;
        this.problem.hessian(dof, dof) = 1;
      end
    end
    
  end
  
  methods (Abstract)
  
    computeNodalConstraint(this, n)
    
  end
  
end