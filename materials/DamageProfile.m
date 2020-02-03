classdef DamageProfile < Material
  
  properties
    center;
    l;
  end
  
  methods
    
    function this = DamageProfile(problem, center, l)
      this = this@Material(problem);
      this.center = center;
      this.l = l;
    end
    
    function c = computeQpMaterial(this)
      p = this.elem.mapped_q_points(this.qp);
      tau = abs(p.x-this.center);
      if tau < 2*this.l
        c = (1-tau/2/this.l)^2;
      else
        c = 0;
      end
    end
    
  end
  
end