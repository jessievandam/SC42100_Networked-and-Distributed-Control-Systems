function [stack, pred] = addnode(stack, next, visited, pred)

newnode = setdiff(next{stack(end)}, pred);
possible = arrayfun(@(x) sprintf('%d,', [stack x]), newnode, 'uni', 0);

isnew = ~ismember(possible, visited);

if any(isnew)
    idx = find(isnew, 1);
    stack = str2num(possible{idx});
    pred = stack(end-1);
else
    [stack, pred] = popnode(stack);
end