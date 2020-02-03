classdef DamageProfileGradient < Material
  
  properties
    center;
    l;
  end
  
  methods
    
    function this = DamageProfileGradient(problem, center, l)
      this = this@Material(problem);
      this.center = center;
      this.l = l;
    end
    
    function c = computeQpMaterial(this)
      p = this.elem.mapped_q_points(this.qp);
      tau = p.x-this.center;
      if abs(tau) < 2*this.l
        c = -sign(tau)/this.l*(1-abs(tau)/2/this.l);
      else
        c = 0;
      end
    end
    
  end
  
end