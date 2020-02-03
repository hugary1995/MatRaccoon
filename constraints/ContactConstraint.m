classdef ContactConstraint < EqualityConstraint
  
  properties
    p_var_id;
    u_var_id;
    degradation;
  end
  
  methods
    
    function this = ContactConstraint(problem, pressure, displacement, degradation)
      this = this@EqualityConstraint(problem);
      this.p_var_id = this.coupledValue(pressure);
      this.u_var_id = this.coupledGradient(displacement);
      this.degradation = degradation;
    end
    
    function value = computeQpConstraint(this)
      p = this.values(this.p_var_id);
      grad_u = this.gradients(this.u_var_id);
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      
      value = 0.5*(p(this.qp)-(1-g)*grad_u(this.qp).x)^2;
    end
    
    function value = computeQpConstraintGradient(this, ivar)
      p = this.values(this.p_var_id);
      grad_u = this.gradients(this.u_var_id);
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      test = this.elem.test(this.i, this.qp);
      grad_test = this.elem.grad_test(this.i, this.qp);
      
      if ivar == this.p_var_id
        value = test*(p(this.qp)-(1-g)*grad_u(this.qp).x);
      elseif ivar == this.u_var_id
        value = -grad_test.x*(1-g)*(p(this.qp)-(1-g)*grad_u(this.qp).x);
      end
    end
    
    function value = computeQpConstraintHessian(this, ivar, jvar)
      g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
      test_i = this.elem.test(this.i, this.qp);
      test_j = this.elem.test(this.j, this.qp);
      grad_test_i = this.elem.grad_test(this.i, this.qp);
      grad_test_j = this.elem.grad_test(this.j, this.qp);
      
      if ivar == this.p_var_id && jvar == this.p_var_id
        value = test_i*test_j;
      elseif ivar == this.p_var_id && jvar == this.u_var_id
        value = -test_i*(1-g)*grad_test_j.x;
      elseif ivar == this.u_var_id && jvar == this.p_var_id
        value = -grad_test_i.x*(1-g)*test_j;
      elseif ivar == this.u_var_id && jvar == this.u_var_id
        value = grad_test_i.x*(1-g)^2*grad_test_j.x;
      end
    end
    
  end
  
end