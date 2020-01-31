classdef ConstantMaterial < Material
  
  properties
    value;
  end
  
  methods
    
    function this = ConstantMaterial(problem, value)
      this = this@Material(problem);
      this.value = value;
    end
    
    function c = computeQpMaterial(this)
      c = this.value;
    end
    
  end
  
end