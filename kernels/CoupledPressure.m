classdef CoupledPressure < Kernel
  
  properties
    var_id;
    p_var_id;
    grad_damage;
  end
  
  methods
    
    function this = CoupledPressure(problem, u, pressure, grad_damage)
      this = this@Kernel(problem);
      this.var_id = this.coupledValue(u);
      this.p_var_id = this.coupledValue(pressure);
      this.grad_damage = grad_damage;
    end
    
    function f = computeQpObjective(this)
      p = this.values(this.p_var_id);
      grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
      u = this.values(this.var_id);
      f = -grad_d*p(this.qp)*u(this.qp);
    end
    
    function value = computeQpGradient(this, ivar)
      p = this.values(this.p_var_id);
      grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
      u = this.values(this.var_id);
      test = this.elem.test(this.i, this.qp);
      
      if ivar == this.var_id
        value = -grad_d*test*p(this.qp);
      elseif ivar == this.p_var_id
        value = -grad_d*test*u(this.qp);
      end
    end
    
    function value = computeQpHessian(this, ivar, jvar)
      if ivar == jvar
        value = 0;
      else
        grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
        value = -grad_d*this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      end
    end
  end
  
end