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
    
    function value = computeQpConstraint(this, ivar)
      if ivar == this.p_var_id
        p = this.values(this.p_var_id);
        grad_u = this.gradients(this.u_var_id);
        g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
        value = this.elem.test(this.i, this.qp)*(p(this.qp)-(1-g)*grad_u(this.qp).x);
      else
        value = 0;
      end
    end
    
    function value = computeQpConstraintGradient(this, ivar, jvar)
      if ivar == this.p_var_id && jvar == this.p_var_id
        value = this.elem.test(this.i, this.qp)*this.elem.test(this.j, this.qp);
      elseif ivar == this.p_var_id && jvar == this.u_var_id
        g = this.problem.materials(this.degradation).data{this.elem.id}{this.qp};
        value = -(1-g)*this.elem.test(this.i, this.qp)*this.elem.grad_test(this.j, this.qp).x;
      else
        value = 0;
      end
    end
    
  end
  
end