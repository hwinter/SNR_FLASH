; 
; shocktube -- evolve a shock tube problem 
;    an IDL version of a fortran code by B. Fryxell
;
;  inputs:  left state (rhol, pl, ul)
;           right state (rhor, pr, ur)
;           initial interface posn (xi)
;           time of final solution requested (t)
;           grid where solution is required (xpts) -- assumed to be 
;                       a uniform 1d mesh of points 
;  outputs:  thermodynamic quantities on the mesh at time t:
;              psoln, rsoln, vsoln

;
      function f,p4, p1, p5, rho1, rho5, gamma
;
;     shock tube equation
;
      z = (p4 / p5 - 1.)
      c1 = sqrt (gamma * p1 / rho1)
      c5 = sqrt (gamma * p5 / rho5)
;
      gm1 = gamma - 1.
      gp1 = gamma + 1.
      g2  = 2. * gamma
;
      fact = gm1 / g2 * (c5 / c1) * z / sqrt (1. + gp1 / g2 * z)
      fact = (1. - fact) ^ (g2 / gm1)
;
      f = p1 * fact - p4
;
      return, f
      end





pro shocktube, rhol, pl, ul, rhor, pr, ur, gamma, xi, t, xpts, psoln, rsoln, vsoln

	  subsample = 50                ; # pts/grid point to average over
                                    ; this should be even!
      npts      = n_elements(xpts)
      nsspts    = npts*subsample+1L

      psoln = xpts*0.
      rsoln = xpts*0.
      vsoln = xpts*0.
  
	  rho = findgen(nsspts)*0.
	  u   = findgen(nsspts)*0.
      p   = findgen(nsspts)*0.

      dx = xpts[2] - xpts[1]       ; assuming a uniform mesh!!!
      xl = xpts[0]-.5*dx
      xr = xpts[npts-1]+.5*dx
      xsspts = findgen(nsspts)*(dx/(1.*subsample)) + xl
      
;
      if ((ul ne 0.) or (ur ne 0.)) then begin
         print, 'must have ul = ur = 0.'
      end
;
;--------interval over which to compute solution
;
     if (xr lt xl) then begin
         print, 'xr must be greater than xl'
      end
;
;
;-----begin solution
;
;
      if (pl gt pr) then begin
         rho1 = rhol
         p1   = pl
         u1   = ul
         rho5 = rhor
         p5   = pr
         u5   = ur
      end else begin
         rho1 = rhor
         p1   = pr
         u1   = ur
         rho5 = rhol
         p5   = pl
         u5   = ul
      end
;
;-----solve for post-shock pressure by secant method
;
;--------initial guesses
;
         p40 = p1
         p41 = p5
         f0  = f(p40, p1, p5, rho1, rho5, gamma)
;
;--------maximum number of iterations
;
         itmax = 20
;
;--------maxium allowable relative error
;
         eps = 1.e-8
;
         iter = 1
         error = 100.
         f1 = 314159

         while (iter lt itmax and error gt eps) do begin
            f1 = f(p41, p1, p5, rho1, rho5, gamma)
;
            p4 = p41 - (p41 - p40) * f1 / (f1 - f0)
;
            error = abs (p4 - p41) / p41
;
            p40 = p41
            p41 = p4
            f0  = f1
         end
;
       	 if (error gt eps) then begin
              print, 'did not converge'
         end
;
;--------compute post-shock density and velocity
;
         z  = (p4 / p5 - 1.)
         c5 = sqrt (gamma * p5 / rho5)

         gm1 = gamma - 1.
         gp1 = gamma + 1.
         gmfac1 = 0.5 * gm1 / gamma
         gmfac2 = 0.5 * gp1 / gamma
;
         fact = sqrt (1. + gmfac2 * z)
;
         u4 = c5 * z / (gamma * fact)
         rho4 = rho5 * (1. + gmfac2 * z) / (1. + gmfac1 * z)
;
;--------shock speed
;
         w = c5 * fact
;
;--------compute values at foot of rarefaction
;
         p3 = p4
         u3 = u4
         rho3 = rho1 * (p3 / p1)^(1. /gamma)
;
;--------compute positions of waves
;
      if (pl gt pr) then begin
         c1 = sqrt (gamma * p1 / rho1)
         c3 = sqrt (gamma * p3 / rho3)
;
         xsh = xi + w * t
         xcd = xi + u3 * t
         xft = xi + (u3 - c3) * t
         xhd = xi - c1 * t
;
;
;--------compute solution as a function of position
;
;
;
		print, nsspts
	    i = ulong(1)
         for i = ulong(0), ulong(nsspts-1) do begin
            if (xsspts(i) lt xhd) then begin
               rho(i) = rho1
               p  (i) = p1
               u  (i) = u1
            end else if (xsspts(i) lt xft) then begin
               u(i)   = 2. / gp1 * (c1 + (xsspts(i) - xi) / t)
               fact   = 1. - 0.5 * gm1 * u(i) / c1
               rho(i) = rho1 * fact ^ (2. / gm1)
               p(i)   = p1 * fact ^ (2. * gamma / gm1)
            end else if (xsspts(i) lt xcd) then begin
               rho(i) = rho3
               p  (i) = p3
               u  (i) = u3
            end else if (xsspts(i) lt xsh) then begin
               rho(i) = rho4
               p  (i) = p4
               u  (i) = u4
            end else begin
               rho(i) = rho5
               p  (i) = p5
               u  (i) = u5
            end
         end
      end
;
;-----if pr > pl, reverse solution
;
      if (pr gt pl) then begin
;
         c1 = sqrt (gamma * p1 / rho1)
         c3 = sqrt (gamma * p3 / rho3)
;
         xsh = xi - w * t
         xcd = xi - u3 * t
         xft = xi - (u3 - c3) * t
         xhd = xi + c1 * t
;
         for i = 0, nsspts-1 do begin
            if (xsspts(i) lt xsh) then begin
               rho(i) = rho5
               p  (i) = p5
               u  (i) = -u5
            end else if (xsspts(i) lt xcd) then begin
               rho(i) = rho4
               p  (i) = p4
               u  (i) = -u4
            end else if (xsspts(i) lt xft) then begin
               rho(i) = rho3
               p  (i) = p3
               u  (i) = -u3
            end else if (xsspts(i) lt xhd) then begin
               u(i)   = -2. / gp1 * (c1 + (xi - xsspts(i)) / t)
               fact   = 1. + 0.5 * gm1 * u(i) / c1
               rho(i) = rho1 * fact ^ (2. / gm1)
               p(i)   = p1 * fact ^ (2. * gamma / gm1)
            end else begin
               rho(i) = rho1
               p  (i) = p1
               u  (i) = -u1
            end
         end
      end
;

      for i=ulong(0),ulong(npts-1) do begin
	      rsoln(i) = mean(rho[subsample*i:subsample*(i+1)])
	      psoln(i) = mean(p[subsample*i:subsample*(i+1)])
	      vsoln(i) = mean(u[subsample*i:subsample*(i+1)])
      end

end

