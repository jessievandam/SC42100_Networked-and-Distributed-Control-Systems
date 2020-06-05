function [stack, pred] = popnode(stack)

stack = stack(1:end-1);
if length(stack) > 1
    pred = stack(end-1);
else
    pred = [];
end

