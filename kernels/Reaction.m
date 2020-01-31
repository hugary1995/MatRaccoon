classdef Reaction < Kernel
  
  properties
    var_id;
  end
  
  methods
    
    function this = Reaction(problem, variable)
      this = this@Kernel(problem);
      this.var_id = this.coupledValue(variable);
    end
    
    function f = computeQpObjective(this)
      u = this.values(this.var_id);
      f = 0.5*u(this.qp)*u(this.qp);
    end
    
    function g = computeQpGradient(this, ivar)
      u = this.values(this.var_id);
      g = this.elem.test(this.i, this.qp)*u(this.qp);
    end
    
    function h = computeQpHessian(this, ivar, jvar)
      h = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
    end
    
  end
  
end