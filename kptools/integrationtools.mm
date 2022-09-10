
#int_x0^x1f(x)=int_-1^1 f(x(s))dx/ds ds

# uses double exponental quadrature to integrate f(x(t)) from one endpoint x0=x(t=-1)=gamma(-1) 

# to the next x1=gamma(1) in the complex plane via the straight line
# gamma(s)=x0(1-s)/2+x1(1+s)
lineint_de:= proc(x0, x1, f, x, DI := Digits)
	local val_lvl, max_crt := -1,	tol := 1^(-2*DI), tol_crt := 1, 
	    lvl_crt := max( base_lvl - 1,0 ), val_tot := 0;
	Digits := DI;

	
    while tol_crt > tol or lvl_crt < global_lvl do 
		#increase lvl of edge
		lvl_crt := lvl_crt + 1;
		val_lvl := lvlk_lineint_de(x0, x1, f, x, lvl_crt);
		max_crt := max( max_crt, abs(val_lvl) ); #too restrictive
		val_lvl := 1/2^(lvl_crt-1)*val_lvl; 
		tol_crt := evalf( abs( val_lvl - val_tot ) );		
		val_tol := evalf( (val_tol + val_lvl)/2 ); # bc 2^blvl valatLvl
	od;
	return val_tol;
end proc;

lvlk_lineint_de := proc(x0,x1, f, x, lvl)
	local tol_crt := 1, tol := 1^(-Digits), val_cum := 0,
        der := diff(x0*(1-s_coord)/2-x1*(1+s_coord)/2,s_coord), 
        k := 0, j, x_smaller, x_larger, w, 
        h_lvl := 2^(-lvl), val_max := -1;
	
    
    if global_lvl < lvl then #should throw error for now      
		global_lvl := lvl;
	fi;
	
	if lvl = base_lvl then #include origin
		w := intData["w"][0];	
		x_larger := subs(s_coord = intData["s"][0], );#fix this
		Digits := Digits + EXTRA2;
		val_cum := val_cum +
			subs(x=x_larger, subs(s_coord=intData["s"][0], der*wt*f));
		Digits := Digits-EXTRA2;
	fi;
		
	while tol_crt > tol do		
		k := k+1;
		j := ifelse( lvl = base_lvl, k, 2*k-1);
		t_j := normal(j*h_lvl);
		if assigned(intData["s"][t_j]) then
			sVal :=intData["s"][t_j];
			w := intData["w"][t_j];
		else
			intData["w"][t_j]:=subs( t =t_j,intData["G_t"]);
			intData["s"][t_j]:=subs( t =t_j,intData["abs"]);
			sVal := intData["s"][t_j];
			wt := intData["w"][t_j];
			global_abscissas_count:=global_abscissas_count+1;				
		fi;	
		if evalf(sVal)=1 then #check that we are avoiding endpoints
			break;
		fi;			
		xPos:=subs(sVar=sVal,p);	
		xNeg:=subs(sVar=-sVal,p);	
		val_Pos:=subs(x=xPos,
				subs(sVar=sVal,
	    	 	der*wt*f));
		val_Neg:=subs(x=xNeg,
				subs(sVar=-sVal,
	    	 	der*wt*f));
		tol_crt :=max(evalf(abs(mPos)),evalf(abs(mNeg)));	
		maxVal:=max(eps,maxVal);
		val_cum:=val_cum+mPos+mNeg;					
	od;
		return val_cum;			
end proc;