functions{
	real coalescent_log(real []x, real pop){
		real p;
		real n;
		real l;
		int i;
		real temp;
		n = size(x);
		p = 0;
		l = n + 1.0;
		i = 1;

		while(i <= floor(n)){
			temp = ((l*(l-1.0))/2.0)/pop;
			p = p + log(temp) - (temp * x[i]);
			l = l - 1.0;
			i = i + 1;
		}
		return p;
	}
	
	real oneOnX_log(real x){
		return -log(x);
	}
}

data {
    int <lower=0> L;               // alignment length
    int <lower=0> S;               // number of tips
    real<lower=0,upper=1> tipdata[S,L,4];
    int <lower=0,upper=2*S> peel[S-1,3];   // list of nodes for peeling
    int map[S-1,S];
}

transformed data {
    int bcount; // number of branches
    int tauCount;
    bcount = 2*S-2;
    tauCount = S-1;
}

parameters {
	real <lower=0,upper=100> taus[tauCount]; // intercoalescent times
	real <lower=0,upper=100> kappa;
	real <lower=0,upper=100>theta;
}


model {
	vector[4] partials[2*S,L];   // partial probabilities for the S tips and S-1 internal nodes
    matrix[4,4] fttm[bcount]; // finite-time transition matrices for each branch
	real k;
	real r;
	real exp1;
	real exp22;
	vector [bcount] blens; // branch lengths
	
    // set some branch length priors
    theta ~ oneOnX();
    
    taus ~ coalescent(theta);
    
    kappa ~ uniform(0,100);
    //kappa ~ lognormal(1.0,1.25);
	
	for( b in 1:bcount ){
   		blens[b] = 0;
    }
    
    for( i in 1:tauCount ){
    	for( j in 1:S ){
    		if(map[i,j] > 0){
    			blens[map[i,j]] = blens[map[i,j]] + taus[i];
    		}
    	}
    }
    
	k = (kappa + 1.0)/2.0;
	r = 4.0 / (kappa + 2.0);

    for( b in 1:bcount ){
	
		exp1 = exp(-blens[b]*r);
        exp22 = exp(-k * blens[b]*r);
	
		//A
		fttm[b][1,1] = 0.25 * (1. + exp1) + 0.5 * exp22; //A
		fttm[b][1,2] = 0.25 * (1. - exp1);              //C
		fttm[b][1,3] = 0.25 * (1. + exp1) - 0.5 * exp22; //G
		fttm[b][1,4] = 0.25 * (1. - exp1);              //T, U
	
		//C
		fttm[b][2,1] = 0.25 * (1. - exp1);              //A
		fttm[b][2,2] = 0.25 * (1. + exp1) + 0.5 * exp22; //C
		fttm[b][2,3] = 0.25 * (1. - exp1);              //G
		fttm[b][2,4] = 0.25 * (1. + exp1) - 0.5 * exp22; //T, U
	
		//G
		fttm[b][3,1] = 0.25 * (1. + exp1) - 0.5 * exp22; //A
		fttm[b][3,2] = 0.25 * (1. - exp1);              //C
		fttm[b][3,3] = 0.25 * (1. + exp1) + 0.5 * exp22; //G
		fttm[b][3,4] = 0.25 * (1. - exp1);              //T, U
	
		//T, U
		fttm[b][4,1] = 0.25 * (1. - exp1);              //A
		fttm[b][4,2] = 0.25 * (1. + exp1) - 0.5 * exp22; //C
		fttm[b][4,3] = 0.25 * (1. - exp1);              //G
		fttm[b][4,4] = 0.25 * (1. + exp1) + 0.5 * exp22; //T, U
    }
    
    
    // copy tip data into node probability vector
    for( n in 1:S ) {
        for( i in 1:L ) {
            for( a in 1:4 ) {
                partials[n,i][a] = tipdata[n,i,a];
            }
        }
    }
    
    // calculate tree likelihood for the topology encoded in peel
	for( i in 1:L ) {
		for( n in 1:(S-1) ) {
            partials[peel[n,3],i] = (fttm[peel[n,1]]*partials[peel[n,1],i]) .* (fttm[peel[n,2]]*partials[peel[n,2],i]);
        }
        
        // multiply by background nt freqs (assuming uniform here)
        for(j in 1:4){
        	partials[2*S,i][j] = partials[peel[S-1,3],i][j] * 0.25;
        }
        // add the site log likelihood
        target += log(sum(partials[2*S,i]));
    }
}
