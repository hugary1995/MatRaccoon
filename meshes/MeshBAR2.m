classdef MeshBAR2 < Mesh

methods

  function this = MeshBAR2(xmin, xmax, nx)
    s = (xmax-xmin)/nx;
    for x = xmin:s:xmax
      this.nodes = [this.nodes, Node(x,0)];
    end
    
    for i = 1:nx
      this.elems = [this.elems, BAR2([this.nodes(i), this.nodes(i+1)])];
    end
    
    this.node_sets = containers.Map();
    this.node_sets('left') = this.nodes(1);
    this.node_sets('right') = this.nodes(end);
    
    this.side_sets = containers.Map();
    this.side_sets('left') = POINT1(this.nodes(1));
    this.side_sets('right') = POINT1(this.nodes(end));
  end
  
end

end