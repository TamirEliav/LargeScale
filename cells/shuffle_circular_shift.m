function t = shuffle_circular_shift(t,ti,shift)
    t = ti(1) + mod(t+shift,diff(ti));
end

