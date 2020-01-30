classdef Reaction < Kernel
  
  methods
    
    function this = Reaction(problem, variable)
      this = this@Kernel(problem, variable);
    end
    
    function f = computeQpObjective(this)
      f = 0.5*this.u(this.qp)*this.u(this.qp);
    end
    
    function g = computeQpGradient(this)
      g = this.elem.test(this.i, this.qp)*this.u(this.qp);
    end
    
    function h = computeQpHessian(this)
      h = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
    end
    
  end
  
end