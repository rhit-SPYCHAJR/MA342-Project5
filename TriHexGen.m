function H = TriHexGen(num_layers)

%Mark Worden - 05/10/2024
%
%H = TriHexGen(num_layers) generates a graph object H depicting a hexagon
%made up of triangles. num_layers dictates how many total layers of
%triangles are generated onto the geometry, with THREE being the minimum
%possible

if(num_layers<3)%enforce input req
    num_layers = 3;
end

H = graph;%init
H = addnode(H,1);%make core node
currnodes = H.numnodes;%number of current nodes
for(i=1:6)
    lower = 1;%set which node to connect outer layer to
    next = mod(i,6)+currnodes+1;%ID of next node to be connected
    H = addedge(H,lower,i+currnodes,1);%connects outer layer to inner layer
    H = addedge(H,next,i+currnodes,1);%connects outer layer together
end

prevnodes = 0;
prevnodes = prevnodes + currnodes;
currnodes = H.numnodes;%set up new counting values for logic

branch_count = 0;
numbranch = 3;
for(i=1:currnodes-prevnodes)%connect new nodes to lower layer

for(j=1:numbranch)
    H = addedge(H,i+prevnodes,branch_count+currnodes+prevnodes,1);
    branch_count = branch_count + 1;
    if(branch_count+currnodes+prevnodes>19)
    branch_count = 0;
    end
end
branch_count = branch_count-1;

end

for(i=1:12)%connect outer layer to eachother
    
    next = mod(i,12)+currnodes+1;
    H = addedge(H,next,i+currnodes,1);

end
corner_start = 1;
for(layer=3:num_layers)

node_layer = 6*layer;

corner_count = corner_start;
prevnodes = currnodes;
currnodes = H.numnodes;%set up new counting values for logic
branch_count = layer-1;
numbranch = 3;
for(i=1:node_layer-6)%connect new nodes to lower layer

for(j=1:2)
    if(branch_count+currnodes>currnodes+node_layer-1)
        branch_count = 0;
    end
    if(j==2&&i==node_layer-6)
        H = addedge(H,prevnodes+node_layer-6,layer+currnodes,1);
    else
        H = addedge(H,i+prevnodes,branch_count+currnodes+1,1);
    end
    branch_count = branch_count + 1;
    if(corner_count == layer-1)
        H = addedge(H,i+prevnodes,branch_count+currnodes+1,1);
        branch_count = branch_count + 1;
        corner_count = 0;
    end

end
corner_count = corner_count + 1;
branch_count = branch_count-1;

end
for(i=1:node_layer)%connect outer layer to eachother
    
    next = mod(i,node_layer)+currnodes+1;
    if(layer>=3)
    H = addedge(H,next,i+currnodes,1);
    end

end
corner_start = corner_start + 1;
end
end