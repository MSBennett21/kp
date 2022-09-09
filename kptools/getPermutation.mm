 
 # This procedure gives the permutation that transforms 
# lijst1 to lijst2. Because of numerical roundoff, equalities
# between elements are never used.

  getPermutation:=proc(lijst1::list,lijst2::list)
        local outlijst, m,mm,ii, jj,r,S,d;
            
        outlijst:=NULL;
        S:={seq(ii,ii=1..nops(lijst1))};
        for ii in lijst2 do
            m:=infinity;
            for jj in S minus {outlijst} do
                mm:=abs(ii-lijst1[jj]);
                if evalf(mm-m)<0 then
                    d:=m;
                    m:=mm;
                    r:=jj;
                elif evalf(mm-d) <0 then
                    d:=mm;
                fi;
            od;
            outlijst:=outlijst,r;
            if fnormal(evalf(m*10-d)) >0 then
                print(m*10-d);
                error "Computation inaccurate, use larger value for Digits"
            fi;
        od;
        #print(outlijst);
        [outlijst]
    end;