//////////////////////////////////////////////////////
// glyphis Wenlock CKMR code /////////////////////////
// includes male and female skip spawning ////////////
// T. Patterson, R. Hillary CSIRO 2017 ///////////////
//////////////////////////////////////////////////////

#include <TMB.hpp>

template<class Type> Type inv_logit(Type x) {return(Type(1)/(Type(1)+exp(-x)));}

template<class Type>
Type objective_function<Type>::operator() ()
{
	PARAMETER(lambda);		// growth rate
	PARAMETER(lN0);					// 	 - log N0
	PARAMETER(iphi); 		// survival probability.
	PARAMETER(lnu); 	       		// litter effect
  PARAMETER(izeta);        		// sex-ratio parameter
  PARAMETER(itheta);              // multiple paternity parameter (relates to females)
  PARAMETER(lgamma);              // multiple female breeding partners per year (relates to males)
  PARAMETER(ipsi);          // female skip-spawning parameter

	DATA_INTEGER(nsex);
	DATA_INTEGER(tmax);  			// number of cohorts covered including tzero
  DATA_INTEGER(tzero); 			// cohort corresponding to time 0
	DATA_SCALAR(p_hsp); 			// critical PLOD (0.92).

	DATA_IVECTOR(c1); 				// cohort 1.
	DATA_IVECTOR(c2); 				// cohort 2
	DATA_IVECTOR(kcode); 			// ktype (0 UP, 1 HSP, 2 FSP)
	DATA_IVECTOR(k); 				// number of kin in that comparison grouping
	DATA_VECTOR(plod); 				// PLOD for mtDNA stuff later on
  DATA_IVECTOR(h1);         // haplotype of fish 1
  DATA_IVECTOR(h2);         // haplotype of fish 2
  DATA_VECTOR(ph1);         // haplotype frequency of fish 1's haplo
  DATA_VECTOR(ph2);         // haplotype frequency of fish 2's haplo

	int nobs = c1.size();

	matrix<Type> N(nsex, tmax);
  matrix<Type> Ntilde(nsex,tmax);
	N.setZero();
  Ntilde.setZero();

	Type theta = inv_logit(itheta);
  Type gamma = exp(lgamma);

	Type nu = exp(lnu);
	Type zeta = inv_logit(izeta);
  Type phi = inv_logit(iphi);
	vector<Type> sexrat(nsex);
  vector<Type> psi(nsex);
  sexrat(0) = zeta;
  sexrat(1) = Type(1)-zeta;
  psi(0) = inv_logit(ipsi);

	// indices.
	// s = 1,2 = fem, male.

	//  POP DYNAM /
    // set up the sex-specific parameters

  Type N0 = exp(lN0);
	for(int t=0; t < tmax; t++) { 		// time
		for(int s=0; s < nsex; s++) {	// sex
			   N(s,t) = N0 * sexrat(s) * exp(lambda * t);
         Ntilde(s,t) = N(s,t) * (Type(2)-psi(s))/Type(2);
		}
	}

  Type Nf = N0*zeta*(Type(1)-psi(0)/Type(2));

  //  CK LIKE CALCS
  Type logl = Type(0);
  Type logl_ckmr = Type(0);
  Type logl_mtdna = Type(0);
  Type pkin;
  Type pfsp;
  Type pmhsp;
  Type pphsp;
  Type psum;
  vector<Type> praw(3);
  vector<Type> pmtdna(3);
  int cmin,cmax;

  for(int i=0; i < nobs; i++) {           // loop over obs.
    // within-cohort cases
    if(c1(i)==c2(i)) {// within cohort comparisons

      cmin = c1(i);
      cmin -= tzero;

      // FSPs

       pfsp = (nu * (Type(1)-theta))/ Ntilde(0, cmin);


      // MHSPs

      pmhsp = (p_hsp * nu * theta)/ Ntilde(0, cmin);

      // PHSPs

     pphsp = (p_hsp * gamma)/ Ntilde(1, cmin);

      pkin = pfsp + pmhsp + pphsp;

      if(pkin > Type(0)) {

        if(k(i) == 0) logl_ckmr += log(Type(1)-pkin); // Female
        if(k(i) == 1) {

          logl_ckmr += log(pkin);

          // mtDNA bit needs calculating for an identified kin pair (PLOD > eta)

          praw(0) = pfsp;
          praw(1) = pmhsp + pphsp * ph2(i);
          praw(2) = pphsp * (Type(1)-ph2(i));
          psum = praw.sum();
          for(int m=0; m<3; m++) pmtdna(m) = praw(m)/psum;

          // is it an FSP?
          if(kcode(i) == 2) logl_mtdna += log(pmtdna(0));
          // is it an HSP which shares a haplo?
          if(kcode(i) == 1 & h2(i) == h1(i)) logl_mtdna += log(pmtdna(1));
          // is it an HSP which doesn't share a haplo?
          if(kcode(i) == 1 & h2(i) != h1(i)) logl_mtdna += log(pmtdna(2));

        }

      }
    }

    // across-cohort case

	  if( c1(i)!=c2(i) ) {

      cmin = c1(i) < c2(i) ? c1(i) : c2(i);
      cmax = c1(i) < c2(i) ? c2(i) : c1(i);
      cmin -= tzero;
      cmax -= tzero;


      // FSPs

      pfsp = Type(0);

      // MHSPs

      if((cmax-cmin) % 2 == 0) {
   	  
        pmhsp = pow(phi,cmax-cmin)/Ntilde(0,cmax);
        pmhsp *= p_hsp;

      } else {

        pmhsp = (Type(1)-psi(0))*pow(phi,cmax-cmin) /(Ntilde(0,cmax));
        pmhsp *= p_hsp; 

      }

      // PHSPs
      
      if((cmax-cmin) % 2 == 0) {

        pphsp = pow(phi,cmax-cmin) /(Ntilde(1,cmax));
        pphsp *= p_hsp;

      } else {

        pphsp = (Type(1)-psi(1))*pow(phi,cmax-cmin) /(Ntilde(1,cmax));
        pphsp *= p_hsp; 

      }

      pkin = pfsp + pmhsp + pphsp;

      if(pkin > Type(0)) {

        if(k(i) == 0) logl_ckmr += log(Type(1)-pkin);
        if(k(i) == 1) {

          logl_ckmr += log(pkin);

          // mtDNA bit needs calculating for an identified kin pair (PLOD > eta)

          praw(0) = pfsp;
          praw(1) = pmhsp + pphsp * ph2(i);
          praw(2) = pphsp * (Type(1)-ph2(i));
          psum = praw.sum();
          for(int m=0; m<3; m++) pmtdna(m) = praw(m)/psum;

          // is it an FSP?
          if(kcode(i) == 2) logl_mtdna += log(pmtdna(0));
          // is it an HSP which shares a haplo?
          if(kcode(i) == 1 & h2(i) == h1(i)) logl_mtdna += log(pmtdna(1));
          // is it an HSP which doesn't share a haplo?
          if(kcode(i) == 1 & h2(i) != h1(i)) logl_mtdna += log(pmtdna(2));
        }
      }

    }
  }

  // overall log-likelihood

  logl = logl_ckmr + logl_mtdna;

  // reporting on natural scale
  ADREPORT(N0);
  ADREPORT(nu);
  ADREPORT(zeta);
  ADREPORT(theta);
  ADREPORT(gamma);
  ADREPORT(zeta*N0);
  ADREPORT(psi(0));
  ADREPORT(Nf);

  return(-logl);
}


