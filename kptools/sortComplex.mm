sortComplex:=proc(ll,tt)
    local s,t,func,b;

	#if not type(ll[1],'float') then		
	#	procname(map(a-> convert(a,float),ll) ,args[2..-1]);
	#else
	if tt=`normal` then #take left most and closest to realaxis 
		func := proc(s::complexcons,t::complexcons) ::boolean;#just lex order
			return evalb( s=t or fnormal(evalf(Re(s) - Re(t) ))< 0) or evalb( fnormal(evalf( Re(s) - Re(t) ) ) = 0. and fnormal( evalf(Im(s)-Im(t)))<0 );
		end proc;
		sort(ll, func )
    elif tt='Arg' then 
		b:=args[3]; #need the b that is the focal point of this ordering
		func := proc(s,t)::boolean; #if s=b or t cannot be be and has larger argument, or same argument but t is farther way
			return evalb( s=b or (t <> b and ( evalf(argument( s- b) - argument( t - b ) ) < 0 or ( evalf(argument( s-b ) - argument( t-b ))= 0 and evalf( abs( s- b)-abs(t-b) ) <= 0))) ) ;
		end proc;   
		sort(ll,func )	
	elif tt='absNorm' then
		b:=args[3]; 
		func := proc(s,t)::boolean;		#closest to real line	
		return evalb( s=t or abs(Im(s-b))<abs(Im(t-b)) or (abs(Im(s-b))=abs(Im(t-b)) and  Re(s-b)< Re(t-b)));
		end proc;   
		sort(ll,func )
	fi;
	#fi;

	
end:

#SortFunction := proc( z,w ) :: boolean; 
#      local Z, W;
#         Z := Complex( ( z - canonicalVertex )[1], ( z - canonicalVertex )[2] );
#         W := Complex( ( w - canonicalVertex )[1], ( w - canonicalVertex )[2] );
         

        
#         if evalf(argument( Z ) - argument( W ))<=0 then  
#            return true;
#         elif evalf(argument( Z ) - argument( W ))<=0 and evalf(abs( Z ) - abs( W )) <=0 then
#            return true;
#         else  
#            return false;
#         fi;
#end proc;  