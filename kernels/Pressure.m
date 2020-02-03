classdef Pressure < Kernel
  
  properties
    var_id;
    pressure;
    grad_damage;
  end
  
  methods
    
    function this = Pressure(problem, u, pressure, grad_damage)
      this = this@Kernel(problem);
      this.var_id = this.coupledValue(u);
      this.pressure = pressure;
      this.grad_damage = grad_damage;
    end
    
    function f = computeQpObjective(this)
      p = this.problem.materials(this.pressure).data{this.elem.id}{this.qp};
      grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
      u = this.values(this.var_id);
      f = -p*u(this.qp)*grad_d;
    end
    
    function value = computeQpGradient(this, ivar)
      p = this.problem.materials(this.pressure).data{this.elem.id}{this.qp};
      grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
      test = this.elem.test(this.i, this.qp);
      value = -p*test*grad_d;
    end
    
    function value = computeQpHessian(this, ivar, jvar)
      value = 0;
    end
  end
  
end