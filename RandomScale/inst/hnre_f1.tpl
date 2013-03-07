// Half normal detection function with random effects
DATA_SECTION
   init_int n;                          // number of distances
   init_number width;                   // truncation width
   init_int debug;                      // flag for debugging
   init_vector xs(1,n);                 // distances
PARAMETER_SECTION
   init_bounded_number beta(-3,2,1);    // beta parameter for log-sigma;
   init_bounded_number sigeps(-10,1,2); // log(sigma_epsilon) for random effect;             
   random_effects_vector u(1,n,2);      // random effect for scale
   objective_function_value f;          // negative log-likelihood

PROCEDURE_SECTION
   int j;
// loop over each observation computing sum of log-likelihood values
   f=0;
   if(debug>0)
   {
      cout << "beta = " << beta << endl;
      cout << "sigeps = " << sigeps << endl;
   }
   for (j=1;j<=n;j++)                          // loop over each observation computing log-likelihood
      ll_j(xs(j),beta,sigeps,u(j));   
	  
   if(debug>0) 
      cout << "f = " << f << endl;

  

SEPARABLE_FUNCTION void ll_j(const double x, const dvariable& beta,const dvariable& sigeps,const dvariable& u)
   dvariable eps=u*mfexp(sigeps);                                                   // random scale component - N(0,exp(sigeps))
   dvariable sigma=mfexp(beta+eps);                                                 // detection function scale
   f -= -0.5*square(u)-log(sqrt(2*M_PI));                                             // log of std normal density for epsilon
   f -= -log(sigma*sqrt(2*M_PI)*(cumd_norm(width/sigma)-.5)) - 0.5*square(x/sigma);   // log of f(x) for half-normal

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(250502); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(800);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=5000000;



