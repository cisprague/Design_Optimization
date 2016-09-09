function flux = f(a1, a2, a3)
    h = @(x) HE.Ly(a1, a2, a3, x);
    h(1)
    
end