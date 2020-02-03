classdef Degradation < Material
  
  properties
    damage;
    mobility;
    critical_fracture_energy;
    p;
  end
  
  methods
    
    function this = Degradation(problem, damage, M, psic, p)
      this = this@Material(problem);
      this.damage = damage;
      this.mobility = M;
      this.critical_fracture_energy = psic;
      this.p = p;
    end
    
    function c = computeQpMaterial(this)
      d = this.problem.materials(this.damage).data{this.elem.id}{this.qp};
      M = this.problem.materials(this.mobility).data{this.elem.id}{this.qp};
      psic = this.problem.materials(this.critical_fracture_energy).data{this.elem.id}{this.qp};
      c = (1-d)^2/((1-d)^2+M/psic*d*(1+this.p*d));
    end
    
  end
  
end