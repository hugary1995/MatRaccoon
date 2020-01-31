classdef Source < Kernel
  
  properties
    var_id;
    fun;
  end
  
  methods
    
    function this = Source(problem, variable, fun)
      this = this@Kernel(problem);
      this.var_id = this.coupledValue(variable);
      this.fun = fun;
    end
    
    function f = computeQpObjective(this)
      u = this.values(this.var_id);
      p = this.elem.mapped_q_points(this.qp);
      f = u(this.qp)*this.fun(p.x, p.y);
    end
    
    function g = computeQpGradient(this, ivar)
      p = this.elem.mapped_q_points(this.qp);
      g = this.elem.test(this.i, this.qp)*this.fun(p.x, p.y);
    end
    
    function h = computeQpHessian(this, ivar, jvar)
      h = 0;
    end
    
  end
  
end