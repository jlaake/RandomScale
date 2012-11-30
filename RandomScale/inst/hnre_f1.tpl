DATA_SECTION
   init_int n;                        // number of distances
   init_number width;                 // truncation width
   init_vector xs(1,n);               // distances
PARAMETER_SECTION
   init_number beta;                  // beta parameter for log-sigma;
   init_number sigeps;                // log sigma for random effect;             
   random_effects_vector u(1,n);      // random effect for scale
   objective_function_value f;        // negative log-likelihood

PROCEDURE_SECTION
   int j;
// loop over each observation computing sum of log-likelihood values
   f=0;
   for (j=1;j<=n;j++)
   {
      ll_j(j,beta,sigeps,u(j));
   }  

SEPARABLE_FUNCTION void ll_j(const int j, const dvariable& beta,const dvariable& sigeps,const dvariable& u)
   dvariable eps=u*exp(sigeps);
   dvariable sigma=exp(beta+eps);
   f -= -0.5*square(u)-log(sqrt(2*PI));  
   f -= -log(sigma*sqrt(2*PI)*(cumd_norm(width/sigma)-.5)) - 0.5*square(xs(j)/sigma);  

TOP_OF_MAIN_SECTION
  gradient_structure::set_MAX_NVAR_OFFSET(250502); 
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(800);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(100000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  arrmblsize=5000000;



