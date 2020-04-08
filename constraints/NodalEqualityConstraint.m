classdef (Abstract) NodalEqualityConstraint < handle
  
  properties
    problem;
    boundary;
    
    var_ids;
    values;
    
    node;
  end
  
  methods
    
    function this = NodalEqualityConstraint(problem, boundary)
      this.problem = problem;
      this.boundary = boundary;
      
      this.values = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    end
    
    function id = coupledValue(this, variable)
      id = this.problem.getVariableId(variable);
      this.values(id) = [];
      this.var_ids = unique([this.var_ids, id]);
    end
    
    function reinitNode(this, n)
      this.node = n;
      for jvar = this.values.keys()
        this.values(jvar{1}) = this.reinitCoupledValue(jvar{1});
      end
    end
    
    function v = reinitCoupledValue(this, var_id)
      dof = this.problem.globalDoF(this.node, var_id);
      v = this.problem.solution(dof);
    end
    
    function computeConstraint(this)
      this.problem.Ceq = this.problem.Ceq+this.elem.JxW(qp_)*this.computeNodalConstraint();
    end
    
    function computeConstraintGradient(this)
      for ivar = this.var_ids
        dof = this.problem.globalDoF(this.node, ivar);
        this.problem.GCeq(dof) = this.problem.GCeq(dof)+this.computeNodalConstraintGradient(ivar);
      end
    end
    
    function computeConstraintHessian(this)
      for ivar = this.var_ids
        dof_i = this.problem.globalDoF(this.node, ivar);
        for jvar = this.var_ids
          dof_j = this.problem.globalDoF(this.node, jvar);
          this.problem.hessian(dof_i, dof_j) = this.problem.hessian(dof_i, dof_j)+this.problem.lambda*this.computeNodalConstraintHessian(ivar, jvar);
        end
      end
    end
    
  end
  
  methods (Abstract)
    
    computeNodalConstraint(this)
    
    computeNodalConstraintGradient(this, ivar)
    
    computeNodalConstraintHessian(this, ivar, jvar)
    
  end
  
end