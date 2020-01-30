classdef Source < Kernel
  
  properties
    fun;
  end
  
  methods
    
    function this = Source(problem, variable, fun)
      this = this@Kernel(problem, variable);
      this.fun = fun;
    end
    
    function f = computeQpObjective(this)
      p = this.elem.mapped_q_points(this.qp);
      f = this.u(this.qp)*this.fun(p.x, p.y);
    end
    
    function g = computeQpGradient(this)
      p = this.elem.mapped_q_points(this.qp);
      g = this.elem.test(this.i, this.qp)*this.fun(p.x, p.y);
    end
    
    function h = computeQpHessian(this)
      h = 0;
    end
    
  end
  
end