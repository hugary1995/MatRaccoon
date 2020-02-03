classdef PenaltyContact < Kernel
  
  properties
    u_var_id;
    p_var_id;
    degradation;
    penalty;
  end
  
  methods
    
    function this = PenaltyContact(problem, u, pressure, degradation, penalty)
      this = this@Kernel(problem);
      this.u_var_id = this.coupledGradient(u);
      this.p_var_id = this.coupledValue(pressure);
      this.degradation = degradation;
      this.penalty = penalty;
    end
    
    function value = computeQpObjective(this)
      p = this.values(this.p_var_id);
      grad_u = this.gradients(this.u_var_id);
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      k = this.penalty;
      
      value = 0.5*k*(p(this.qp)-(1-g)*grad_u(this.qp).x)^2;
    end
    
    function value = computeQpGradient(this, ivar)
      p = this.values(this.p_var_id);
      grad_u = this.gradients(this.u_var_id);
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      k = this.penalty;
      test = this.elem.test(this.i, this.qp);
      grad_test = this.elem.grad_test(this.i, this.qp);
      
      if ivar == this.p_var_id
        value = test*k*(p(this.qp)-(1-g)*grad_u(this.qp).x);
      elseif ivar == this.u_var_id
        value = -grad_test.x*k*(1-g)*(p(this.qp)-(1-g)*grad_u(this.qp).x);
      end
    end
    
    function value = computeQpHessian(this, ivar, jvar)
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      k = this.penalty;
      test_i = this.elem.test(this.i, this.qp);
      test_j = this.elem.test(this.j, this.qp);
      grad_test_i = this.elem.grad_test(this.i, this.qp);
      grad_test_j = this.elem.grad_test(this.j, this.qp);
      
      if ivar == this.p_var_id && jvar == this.p_var_id
        value = test_i*k*test_j;
      elseif ivar == this.p_var_id && jvar == this.u_var_id
        value = -test_i*k*(1-g)*grad_test_j.x;
      elseif ivar == this.u_var_id && jvar == this.p_var_id
        value = -grad_test_i.x*k*(1-g)*test_j;
      elseif ivar == this.u_var_id && jvar == this.u_var_id
        value = grad_test_i.x*(1-g)^2*grad_test_j.x;
      end
    end
  end
  
end