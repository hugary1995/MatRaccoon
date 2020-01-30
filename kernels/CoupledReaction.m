classdef CoupledReaction < Kernel
  
  properties
    coupled_var;
  end
  
  methods
    
    function this = CoupledReaction(problem, variable, coupled_variable)
      this = this@Kernel(problem, variable);
      this.coupled_var = this.problem.getVariableId(coupled_variable);
      this.coupledValue(coupled_variable);
    end
    
    function f = computeQpObjective(this)
      v = this.coupled_values(this.coupled_var);
      f = this.u(this.qp)*v(this.qp);
    end
    
    function g = computeQpGradient(this)
      v = this.coupled_values(this.coupled_var);
      g = this.elem.test(this.i, this.qp)*v(this.qp);
    end
    
    function h = computeQpHessian(this)
      h = 0;
    end
    
    function h = computeQpOffDiagHessian(this, jvar)
      h = 0;
      if jvar == this.coupled_var
        h = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      end
    end
    
  end
  
end