classdef Contact < Kernel
  
  properties
    var_id;
    degradation;
    grad_damage;
  end
  
  methods
    
    function this = Contact(problem, u, degradation, grad_damage)
      this = this@Kernel(problem);
      this.var_id = this.coupledValue(u);
      this.coupledGradient(u);
      this.degradation = degradation;
      this.grad_damage = grad_damage;
    end
    
    function f = computeQpObjective(this)
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
      u = this.values(this.var_id);
      grad_u = this.gradients(this.var_id);
      f = -(1-g)*grad_u(this.qp).x*u(this.qp)*grad_d;
    end
    
    function value = computeQpGradient(this, ivar)
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
      u = this.values(this.var_id);
      grad_u = this.gradients(this.var_id);
      
      test = this.elem.test(this.i, this.qp);
      grad_test = this.elem.grad_test(this.i, this.qp);
      
%       value = -(1-g)*grad_d*test*(grad_u(this.qp).x)-(1-g)*grad_d*(grad_test.x)*u(this.qp);
      value = -(1-g)*grad_d*test*(grad_u(this.qp).x);
    end
    
    function value = computeQpHessian(this, ivar, jvar)
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      grad_d = this.problem.materials(this.grad_damage).data{this.elem.id}{this.qp};
      
      test_i = this.elem.test(this.i, this.qp);
      test_j = this.elem.test(this.j, this.qp);
      grad_test_i = this.elem.grad_test(this.i, this.qp);
      grad_test_j = this.elem.grad_test(this.j, this.qp);
      
%       value = -(1-g)*grad_d*test_i*(grad_test_j.x)-(1-g)*grad_d*(grad_test_i.x)*test_j;
      value = -(1-g)*grad_d*test_i*(grad_test_j.x);
    end
  end
  
end