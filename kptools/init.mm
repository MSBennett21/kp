

#initalize the int scheme at startup so that if integration is to be used abiscusis will have already be produced, reducing computational cost later
#since 1) level-based int scheme is used, 2) x(-t)=x(t) x'(-t)=x'(t) so only half of the N need to be computed
# uses 2p percision so that we have persicion of abiscusses at least p and so that
# we can have more abiscusses at lower levels 
# uses yj=1-xj so that cancelation error can be avodied <-- come back to this   
#will need to add a set of codes to re poplaute when digits are increased
initialize1 := proc()
	local tol:=1^(-ui_digits), tol_crt:=1, t_j;
    Digits := ui_digits_high+EXTRA2; #may need to update
    # idea: int_-1^1 f(x)dx= int_(-inf)^inf f'(x(t))*x'(t)dt
    #thought x(t)x'(t) is odd*even=odd= der x(t)^2/2
    #could replace t->g(x(t)) with oddeven comp? g(x(-t))g'(x(-t))x'(-t)
    intData := table();
    intData["x"] := table();
    intData["dx/dt"] := table(); 
    intData["x(t)"] := tanh(Pi/2 * sinh(t_cord));
    intData["x'(t)"] := Pi/2 * cosh(t_cord) * sech(Pi/2 * sinh(t_cord))^2;     
    intData["y=1-x"] := 1/( exp(Pi/2 * sinh(t_cord)) * cosh(Pi/2 * sinh(t_cord)));
    #sinh=(e(x)-e(-x))/2 -->exp growth
    #cosh=(e(x)+e(-x))/2 --> exp growth
    #sech= 2/(e(x)+e(-x))-->1/exp growth
    #cosh* sech^2(sinh) --> exp growth * 1/exp^2(exp) growth so exponential decay
    # also monotonically decreasing <-prove

    while tol_crt > tol do
        t_j := normal(global_abscissas_count*2^(-global_lvl)); #to ensure uniqueness of tvalue
        intData["x"][t_j] := evalf(subs(t=t_j, intData["y=1-x"])); 
        intData["dx/dt"][t_j] := evalf(subs(t=t_j, intData["x'(t)"]));

        global_abscissas_count:=global_abscissas_count+1;    
        tol_crt := intData["w"][t_j];
    od;

    print(cat("scheme has been preloaded with ", convert(global_abscissas_count,string), " points"));
end proc;