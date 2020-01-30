classdef Vector
  
  properties
    x;
    y;
  end
  
  methods
    
    function this = Vector(x,y)
      this.x = x;
      this.y = y;
    end
    
    function r = plus(a,b)
      if isa(a,'Vector') && isa(b, 'Vector')
        r = Vector(a.x+b.x,a.y+b.y);
      else
        error(['cannot add vector and ',class(b)]);
      end
    end
    
    function r = minus(a,b)
      if isa(a,'Vector') && isa(b, 'Vector')
        r = Vector(a.x+-b.x,a.y-b.y);
      else
        error(['cannot substract vector and ',class(b)]);
      end
    end
    
    function r = mtimes(a,b)
      if isa(a,'Vector') && isa(b, 'Vector')
        r = a.x*b.x+a.y*b.y;
      elseif isa(a,'Vector') && isa(b, 'double')
        r = Vector(a.x*b,a.y*b);
      elseif isa(a,'double') && isa(b, 'Vector')
        r = Vector(a*b.x,a*b.y);
      else
        error(['cannot dot vector and ',class(b)]);
      end
    end
    
  end
  
end