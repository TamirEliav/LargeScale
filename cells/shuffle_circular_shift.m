function t = shuffle_circular_shift(t,ti,shift)
    t = ti(1) + mod(t-ti(1) + shift,diff(ti));
end

