function x = cramer_Jeffrey_HartantoKwee(A, b)
    [r, c] = size(A);
    if (r ~= c)
        error("A is not a square matrix.")
    end
    [rb, cb] = size(b);
    if (r ~= rb)
        error("The number of rows in b does not match A.")
    end
    if (cb ~= 1)
        error("The number of columns in b is not 1.")
    end
    
    det_A = det(A);
    if abs(det_A) < 1e-9
        det_A = 0;
    end
    det_Ai = zeros(1, r);
    for i = 1:1:r
        temp = A;
        temp(:, i) = b;
        det_Ai(i) = det(temp);
        if (abs(det_Ai(i)) < 1e-9)
            det_Ai(i) = 0;
        end
    end

    if det_A == 0
        if det_Ai == zeros(1, r)
            error("A has infinitely many solutions.")
        else
            error("A has no solutions.")
        end
    end
    x = (det_Ai / det_A).';

end