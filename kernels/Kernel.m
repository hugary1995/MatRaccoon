classdef (Abstract) Kernel < handle
  
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
    
    function this = Kernel(problem)
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
    
    function computeObjective(this)
      for qp_ = 1:length(this.elem.q_points)
        this.qp = qp_;
        this.problem.objective = this.problem.objective+this.elem.JxW(qp_)*this.computeQpObjective();
      end
    end
    
    function computeGradient(this)
      for qp_ = 1:length(this.elem.q_points)
        this.qp = qp_;
        for i_ = 1:length(this.elem.nodes)
          this.i = i_;
          for ivar = this.var_ids
            dof = this.problem.globalDoF(this.elem.nodes(i_), ivar);
            this.problem.gradient(dof) = this.problem.gradient(dof)+this.elem.JxW(qp_)*this.computeQpGradient(ivar);
          end
        end
      end
    end
    
    function computeHessian(this)
      num_dofs = length(this.elem.nodes)*length(this.var_ids);
      local_hessian = zeros(num_dofs, num_dofs);
      for qp_ = 1:length(this.elem.q_points)
        this.qp = qp_;
        for i_ = 1:length(this.elem.nodes)
          this.i = i_;
          for j_ = 1:length(this.elem.nodes)
            this.j = j_;
            for ivar = this.var_ids
              dof_i = this.problem.globalDoF(this.elem.nodes(i_), ivar);
              for jvar = this.var_ids
                dof_j = this.problem.globalDoF(this.elem.nodes(j_), jvar);
                local_hessian(dof_i, dof_j) = local_hessian(dof_i, dof_j)+this.elem.JxW(qp_)*this.computeQpHessian(ivar, jvar);
              end
            end
          end
        end
      end
      this.problem.hessian = this.problem.hessian+local_hessian;
    end
    
  end
  
  methods (Abstract)
    
    computeQpGradient(this, ivar)
    
    computeQpHessian(this, ivar, jvar)
    
  end
  
end