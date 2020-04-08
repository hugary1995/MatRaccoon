classdef (Abstract) EqualityConstraint < handle
  
  properties
    problem;
    
    var_ids;
    values;
    gradients;
    
    i;
    j;
    qp;
    
    elem;
  end
  
  methods
    
    function this = EqualityConstraint(problem)
      this.problem = problem;
      
      this.values = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
      this.gradients = containers.Map('KeyType', 'uint32', 'ValueType', 'any');
    end
    
    function id = coupledValue(this, variable)
      id = this.problem.getVariableId(variable);
      this.values(id) = [];
      this.var_ids = unique([this.var_ids, id]);
    end
    
    function id = coupledGradient(this, variable)
      id = this.problem.getVariableId(variable);
      this.gradients(id) = [];
      this.var_ids = unique([this.var_ids, id]);
    end
    
    function reinitElem(this, e)
      this.elem = e;
      
      for jvar = this.values.keys()
        this.values(jvar{1}) = this.reinitCoupledValue(jvar{1});
      end
      
      for jvar = this.gradients.keys()
        this.gradients(jvar{1}) = this.reinitCoupledGradient(jvar{1});
      end
    end
    
    function v = reinitCoupledValue(this, var_id)
      v = zeros(1,length(this.elem.q_points));
      for i_ = 1:length(this.elem.nodes)
        dof = this.problem.globalDoF(this.elem.nodes(i_), var_id);
        dof_value = this.problem.solution(dof);
        for qp_ = 1:length(this.elem.q_points)
          v(qp_) = v(qp_)+this.elem.test(i_, qp_)*dof_value;
        end
      end
    end
    
    function v = reinitCoupledGradient(this, var_id)
      v(1:length(this.elem.q_points)) = Vector(0, 0);
      for i_ = 1:length(this.elem.nodes)
        dof = this.problem.globalDoF(this.elem.nodes(i_), var_id);
        dof_value = this.problem.solution(dof);
        for qp_ = 1:length(this.elem.q_points)
          v(qp_) = v(qp_)+this.elem.grad_test(i_, qp_)*dof_value;
        end
      end
    end
    
    function computeConstraint(this)
      for qp_ = 1:length(this.elem.q_points)
        this.qp = qp_;
        this.problem.Ceq = this.problem.Ceq+this.elem.JxW(qp_)*this.computeQpConstraint();
      end
    end
    
    function computeConstraintGradient(this)
      for qp_ = 1:length(this.elem.q_points)
        this.qp = qp_;
        for i_ = 1:length(this.elem.nodes)
          this.i = i_;
          for ivar = this.var_ids
            dof = this.problem.globalDoF(this.elem.nodes(i_), ivar);
            this.problem.GCeq(dof) = this.problem.GCeq(dof)+this.elem.JxW(qp_)*this.computeQpConstraintGradient(ivar);
          end
        end
      end
    end
    
    function computeConstraintHessian(this)
      for qp_ = 1:length(this.elem.q_points)
        this.qp = qp_;
        for i_ = 1:length(this.elem.nodes)
          this.i = i_;
          for ivar = this.var_ids
            dof_i = this.problem.globalDoF(this.elem.nodes(i_), ivar);
            for j_ = 1:length(this.elem.nodes)
              this.j = j_;
              for jvar = this.var_ids
                dof_j = this.problem.globalDoF(this.elem.nodes(j_), jvar);
                this.problem.hessian(dof_i, dof_j) = this.problem.hessian(dof_i, dof_j)+this.problem.lambda*this.elem.JxW(qp_)*this.computeQpConstraintHessian(ivar, jvar);
              end
            end
          end
        end
      end
    end
    
  end
  
  methods (Abstract)
    
    computeQpConstraint(this)
    
    computeQpConstraintGradient(this, ivar)
    
    computeQpConstraintHessian(this, ivar, jvar)
    
  end
  
end