;
; soderror
;  
;   pro soderror, filename, RHOL=rhol, PL=pl, UL=pl, 
;                 RHOR=rhor, PR=pr, UR=ur, GAMMA=gamma, XI=xi, 
;                 pres, comppres, dp
;
;       Input
;         filename -- name of the FLASH file to read in
;         left  state (rhol, pl, ul)
;         right state (rhor, pr, ur)
;         gamma of fluid (gamma)
;         initial posn (xi)
;       Output
;         pres     -- pressure field from semi-analytic soln
;         comppres -- computed field from FLASH output
;         dp       -- difference of the two
;         error    -- integral over area of error


pro soderror, filename, RHOL=rhol, PL=pl, UL=ul, RHOR=rhor, PR=pr, UR=ur, GAMMA=gamma, XI=xi, EXCISESHOCK=exciseshock, EXCISEREGION=exciseregion, pres, comppres, dp, error

    if (not(keyword_set(rhol))) then rhol = 1.
    if (not(keyword_set(rhor))) then rhor = .125
    if (not(keyword_set(pl)))   then pl = 1.
    if (not(keyword_set(pr)))   then pr = .1
    if (not(keyword_set(ul)))   then ul = 0.
    if (not(keyword_set(ur)))   then ur = 0.
    if (not(keyword_set(exciseshock)))   then exciseshock = 0.
    if (not(keyword_set(exciseregion)))   then exciseregion = .1

    if (not(keyword_set(xi)))    then xi = .5
    if (not(keyword_set(gamma))) then gamma=1.4

	loadfile, filename, comppres, /PRES, XCOORDS=x, YCOORDS=y, FILETIME=t
    sizes = size(comppres)
    nx = sizes[1]
    ny = sizes[2]

    shocktube, rhol, pl, ul, rhor, pr, ur, gamma, xi, t, x, psoln, rsoln, vsoln

    pres = findgen(nx,ny)

    if (exciseshock > 0.) then begin
		dp = abs((psoln[1:nx-1]-psoln[0:nx-2])/(psoln[1:nx-1]))
        shockindx = min(where(dp gt .1))+1
        shocklocn = x[shockindx]
        exminind = min(where(x ge x[shockindx]-.5*exciseregion))
        exmaxind = max(where(x le x[shockindx]+.5*exciseregion))
    end

    for j=0,ny-1 do begin
       pres[*,j] = psoln
    end

    dp = comppres - pres
    frac = 1.
    if (exciseshock > 0.) then begin
        print, 'Excising from x=',x[exminind],' to ',x[exmaxind]
		dp[exminind:exmaxind,*] = 0.
	    frac = (1.*nx-1.*(exmaxind-exminind+1.))/(1.*nx)
    end

    dv = (x[2]-x[1])*(y[2]-y[1])
	error = total(abs(dp))*dv/frac

end
