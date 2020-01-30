classdef (Abstract) NodalBC < handle
  
  properties
    problem;
    var_id;
    
    boundary;
    node_set;
    
    u;
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
        this.u = this.problem.solution(dof);
        this.problem.constraint(dof) = this.computeNodalConstraint(this.node_set(i));
      end
    end
    
    function computeConstraintGradient(this)
      for i = 1:length(this.node_set)
        dof = this.problem.globalDoF(this.node_set(i), this.var_id);
        this.u = this.problem.solution(dof);
        this.problem.constraint_gradient(dof, dof) = this.computeNodalConstraintGradient(this.node_set(i));
      end
    end
    
  end
  
  methods (Abstract)
  
    computeNodalConstraint(this, n)
  
    computeNodalConstraintGradient(this, n)
    
  end
  
end