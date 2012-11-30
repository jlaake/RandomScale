// Truncated half normal detection function with random effects
DATA_SECTION
   init_int n;                        // number of distances
   init_int width;                    // truncation width
   init_vector xs(1,n);               // distances
   init_int m;                        // number of columns in design matrix
   init_matrix dm(1,n,1,m);           // design matrix for fixed effects
   int i;
   // Weights applied to the likelhood to obtain normalizing probablity
   // Note that (by mistake) weights run 2,4,6,....
   vector w(1,4*n)
   !! w=1.0;
   !! for(i=1;i<=n;i++) w(2*n+2*i)=-1;
	
PARAMETER_SECTION 
   init_vector beta(1,m);             	// beta parameter for log-sigma;
   init_number sigeps; 		            // log(sigma) for random effect;                
   random_effects_vector u(1,2*n);      // random effect for scale; first n for
                                        // numerator and second for integral
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
   for (j=1;j<=n;j++)                                           // loop over each observation computing numerator
     ll_j(xs(j),beta,sigeps,u(j),dm(j));                        // which is the average g(x) integrated over epsilon
   for (j=1;j<=n;j++)                                           // loop over each observation computing numerator
     denom(beta,sigeps,u(n+j),dm(j));                           // compute denominator which is the average_mu jth obs
                                                                // integrated over epsilon and weighted by -1

SEPARABLE_FUNCTION void ll_j(const double x, const dvar_vector& beta,const dvariable& sigeps,const dvariable& u, const dvector& dm)
   dvariable eps=u*exp(sigeps);                                 // random scale component - N(0,exp(sigeps))
   dvariable sigma=exp(dm*beta+eps);                            // detection function scale
   f -= -0.5*square(u)-log(sqrt(2*PI));                         // log of std normal density for epsilon
   f -= -log(sqrt(2*PI))-log(sigma)-0.5*square(x/sigma);        // log of f(x) for half-normal

SEPARABLE_FUNCTION void denom(const dvar_vector& beta,const dvariable& sigeps,const dvariable& u, const dvector& dm)
   dvariable eps=u*exp(sigeps);                                   // random scale component - N(0,exp(sigeps))
   dvariable sigma=exp(dm*beta+eps);                              // detection function scale	
   f -= -0.5*square(u)-log(sqrt(2*PI));                           // log of std normal density for epsilon
   f -= log(1e-10+cumd_norm(width/sigma)-0.5);	                  // cum half std normal
	
TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(250502); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(800);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=5000000;

