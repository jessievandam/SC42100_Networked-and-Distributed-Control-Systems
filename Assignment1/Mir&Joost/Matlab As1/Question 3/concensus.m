function y = concensus(P, stub, u, iter)
    %Extract B from W
    B = [P(:,stub(1)) P(:,stub(2))];
    
    %Remove stubborn nodes from W and B
    for i = 1 : length(stub)
        P(stub(i),:) = [];
        P(:,stub(i)) = [];
        B(stub(i),:) = [];
    end
    Q = P; %store in Q
    
    %initialize y
    y = zeros(size(Q,1),iter);
    y(:,1) = 0.5*ones(size(Q,1),1);

    for i = 1 : iter
        y(:,i+1) = Q*y(:,i) + B*u;    
    end
end

