function c = split_string(str)
% c = split string

n = 1;
not_done = true;
while not_done
    [a b] = strtok(str);
    if (isempty(a)==1)
        not_done = false;
    else
        c(n) = cellstr(a);
        str = b;
        n = n+1;
    end
end

