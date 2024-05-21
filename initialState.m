function sigma = initialState(A)

%A is graph object

n = A.numnodes;%how many points need to be itterated through

sigma = ones(1,n);

for(i=1:n)

    if(rand<0.5)
        sigma(i)=1;
    else
        sigma(i)=-1;
    end

end