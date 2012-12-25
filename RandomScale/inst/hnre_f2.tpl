// Half normal detection function with random effects
DATA_SECTION
   init_int n;                        // number of distances
   init_int width;                    // truncation width
   init_int debug;                    // flag for debugging
   init_vector xs(1,n);               // distances
   // Weights applied to the likelhood to obtain normalizing probablity
   // Note that (by mistake) weights run 2,4,6,....
   vector w(1,2*n+2)
   !! w=1.0;
   !! w(2*n+2)=-n;
	
PARAMETER_SECTION 
   init_bounded_number beta(-3,2,1);    // beta parameter for log-sigma;
   init_bounded_number sigeps(-10,1,2); // log(sigma_epsilon) for random effect;             
   random_effects_vector u(1,n+1,2);    // random effect for scale
   !!set_multinomial_weights(w);        // weights to substract denominator
   objective_function_value f;        	// negative log-likelihood

GLOBALS_SECTION

  // Modification of files needed to make set_multinomial_weights work
  #include "minfil.cpp"		// Laplace approx (too inaccurate here); f1b2fnl3.cpp in ADMB code base
  #include "df1b2gh.cpp"	// Gauss-Hermite integration
  #include "xmodelm5.cpp"	// Gauss-Hermite integration

PROCEDURE_SECTION
   int j;
   f=0;
   if(debug>0)
   {
      cout << "beta = " << beta << endl;
      cout << "sigeps = " << sigeps << endl;
   }
   for (j=1;j<=n;j++)                                           // loop over each observation computing numerator
     ll_j(xs(j),beta,sigeps,u(j));                              // which is the average g(x) integrated over epsilon
   if(debug>0)cout << "f = " << f << endl;
   denom(beta,sigeps,u(n+1));                                   // compute constant denominator which is the average_mu 
                                                                // integrated over epsilon and weighted by n
   if(debug>0)cout << "f = " << f << endl;

SEPARABLE_FUNCTION void ll_j(const double x, const dvariable& beta,const dvariable& sigeps,const dvariable& u)
   dvariable eps=u*mfexp(sigeps);                               // random scale component - N(0,exp(sigeps))
   dvariable sigma=mfexp(beta+eps);                             // random scale
   f -= -0.5*square(u)-log(sqrt(2*PI));                         // log of std normal density for epsilon
   f -= -log(sqrt(2*PI))-log(sigma)-0.5*square(x/sigma);        // log of f(x) for half-normal

SEPARABLE_FUNCTION void denom(const dvariable& beta,const dvariable& sigeps,const dvariable& u)
   dvariable eps=u*mfexp(sigeps);                                 // random scale component - N(0,exp(sigeps))
   dvariable sigma=mfexp(beta+eps);                               // random scale	
   f -= -0.5*square(u)-log(sqrt(2*PI));                           // log of std normal density for epsilon
   f -= log(1e-10+cumd_norm(width/sigma)-0.5);	                  // log of integral of f(x)
	
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(250502); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(800);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=5000000;

