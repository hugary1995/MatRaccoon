classdef MatDiffusion < Kernel
  
  properties
    u_var_id;
    material;
  end
  
  methods
    
    function this = MatDiffusion(problem, variable, material)
      this = this@Kernel(problem);
      this.u_var_id = this.coupledGradient(variable);
      this.material = material;
    end
    
    function f = computeQpObjective(this)
      mat = this.problem.materials(this.material).data;
      coef = mat{this.elem.id}{this.qp};
      grad_u = this.gradients(this.u_var_id);
      f = 0.5*coef*grad_u(this.qp)*grad_u(this.qp);
    end
    
    function g = computeQpGradient(this, ivar)
      mat = this.problem.materials(this.material).data;
      coef = mat{this.elem.id}{this.qp};
      grad_u = this.gradients(this.u_var_id);
      g = coef*this.elem.grad_test(this.i, this.qp)*grad_u(this.qp);
    end
    
    function h = computeQpHessian(this, ivar, jvar)
      mat = this.problem.materials(this.material).data;
      coef = mat{this.elem.id}{this.qp};
      h = coef*this.elem.grad_test(this.i, this.qp)*this.elem.grad_test(this.j, this.qp);
    end
    
  end
  
end