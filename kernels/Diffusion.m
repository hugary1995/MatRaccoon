classdef Diffusion < Kernel
  
  properties
    u_var_id;
  end
  
  methods
    
    function this = Diffusion(problem, variable)
      this = this@Kernel(problem);
      this.u_var_id = this.coupledGradient(variable);
    end
    
    function f = computeQpObjective(this)
      grad_u = this.gradients(this.u_var_id);
      f = 0.5*grad_u(this.qp)*grad_u(this.qp);
    end
    
    function g = computeQpGradient(this, ivar)
      grad_u = this.gradients(this.u_var_id);
      g = this.elem.grad_test(this.i, this.qp)*grad_u(this.qp);
    end
    
    function h = computeQpHessian(this, ivar, jvar)
      h = this.elem.grad_test(this.i, this.qp)*this.elem.grad_test(this.j, this.qp);
    end
    
  end
  
end