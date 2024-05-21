function A = ConnectOppositeEdges(A)

n = A.numnodes;%how many nodes there are
count = 1;

for(i=1:n)

    if(length(neighbors(A,i))<length(neighbors(A,1)))%assumes 1 is center point and will always have the desired amount of connections
        needsNeighbor(count,1) = length(neighbors(A,i));%how many neighbors this point has
        needsNeighbor(count,2) = i;%indicies of nodes that need an extra neighbor
        count = count+1;
    end

end

halfway = (length(needsNeighbor)/2);%assume length is always even
needsNeighbor_ordered = sortrows(needsNeighbor);%generate array of needs neighbor points in descending order of current neighbors

for(i=needsNeighbor(1,2):needsNeighbor(1,2)+halfway-1)
    A = addedge(A,i,i+halfway,1);
end

for(j=1:length(needsNeighbor))
    i = needsNeighbor_ordered(j,2);
    if(i+halfway-1<needsNeighbor(1,2))
        next = mod(i+halfway-1-needsNeighbor(1,2),length(needsNeighbor))+n;%determine next point using wraparound  
    else
        next = mod(i+halfway-1-needsNeighbor(1,2),length(needsNeighbor))+needsNeighbor(1,2);%determine next point using wraparound        
    end
    if(length(neighbors(A,next))<length(neighbors(A,1))&&length(neighbors(A,i))<length(neighbors(A,1)))
        A = addedge(A,i,next,1);
    end
    next = mod(i+halfway+1-needsNeighbor(1,2),length(needsNeighbor))+needsNeighbor(1,2);%determine next point using wraparound
    
    if(length(neighbors(A,next))<length(neighbors(A,1))&&length(neighbors(A,i))<length(neighbors(A,1)))
        A = addedge(A,i,next,1);
    end

end


end