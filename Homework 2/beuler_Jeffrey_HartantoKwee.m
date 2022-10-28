function [tt, yy] = beuler_Jeffrey_HartantoKwee(fun, y0, tf, h)
   tt = [0];
   yy = [y0];
   while (true)
        ti = tt(end);
        yi = yy(end);
        te = ti + h;
        impl = @(y)yi + h * fun(te, y) - y;
        ye = fzero(impl, yi);
        tt(end + 1) = te;
        yy(end + 1) = ye;
        if (te > tf)
            return
        end
   end
end