#include "BootstrapRandom.h"
#include "MathSupport.h"

// --------------------------------------------------------------------
double BootstrapRandom::bootstrap(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	alocateMemory();
	
	double likelihoodValue = -1000000;
	
	BootstrapType btype = EVOLUTION_STRATEGY;
	//BootstrapType btype = RANDOM_TRIES_SEPARATE_VARS;
	switch(btype)
	{
	case ONLY_RANDOM_TRIES:
		bootstrapRandomly(&likelihoodValue, &aVars[0], mN);
		break;
		
	case RANDOM_TRIES_SEPARATE_VARS:
		bootstrapEachDirectionRandomly(&likelihoodValue, &aVars[0], mN);
		bootstrapEachDirectionRandomly(&likelihoodValue, &aVars[0], 0);
		//bootstrapEachDirectionRandomly(&likelihoodValue, &aVars[0], 0);
		break;
		
	case EVOLUTION_STRATEGY:
		//int numGenerations = static_cast<int> (log(static_cast<double>(mN))*1.5);
		int numGenerations = 5;
		bootstrapEvolutionStrategy(&likelihoodValue, &aVars[0], numGenerations);
		break;
	};
	return -likelihoodValue;
}


// --------------------------------------------------------------------
double BootstrapRandom::evaluateLikelihood(double *x)
{
	memcpy(&x_[0], x, size_vect);
	return evaluateLikelihood(x_);
}


// --------------------------------------------------------------------
double BootstrapRandom::evaluateLikelihood(const std::vector<double> &x)
{
	return mModel->computeLikelihood(x, mTrace);
}


// --------------------------------------------------------------------
void BootstrapRandom::alocateMemory(void)
{
	size_vect = mN*sizeof(double);
	
	x_.resize(mN);
	mSpace.resize(mN);
	
#ifdef BOOTSTRAP_ES
	
	// choose a population size
#if 1
	// we take here a population of:
	// - 15 for small problems
	// - 70 for medium to large problems
	// - linear in between 
	if(mN < 10)
		mPopSize = 15;
	else if(mN < 50)
		mPopSize = static_cast<int>(20 + float(mN-10)*1.25);
	else
		mPopSize = 70;
#else
	// linear model: the population size is linearly dependant 
	// on the number of variables
	 mPopSize = 10 + static_cast<int>( 1.15*static_cast<double>(mN) );
#endif	
	mGASpace.resize( mPopSize*(mN+1) );
	
	mPopFitness = &mGASpace[0];
	mPopPos 	= mPopFitness+mPopSize;
	
#endif // BOOTSTRAP_ES
}


// --------------------------------------------------------------------
double BootstrapRandom::generateRandom(unsigned int i)
{
	double randNumber, low, high;
	
	high = mUpperBound[i];
    low  = mLowerBound[i];
    randNumber = high+1.0;
    
	// branch lengths
	if (i < mNumTimes)
	{
		while(randNumber >= high || randNumber <= low)
        	randNumber = gamma_dist_T( rng );
	}
	else
	{
		switch(i-mNumTimes)
		{
#ifdef USE_ORIGINAL_PROPORTIONS
		case 0:
        case 1:
        	double v0(2.0),v1(2.0);
        	low = 0.0;
        	high = 1.0;
        	while(v0 >= high || v0 <= low)
        		v0 = 1.0 - exp_dist_v0( rng );
			while(v1 >= high || v1 <= low)
        		v1 = 1.0 - gamma_dist_v1( rng );
        	
        	// convert to original proportion variables
        	double a,b, w0, w1;
        	
        	a = v0*v1;
        	b = v0*(1.0-v1);
        	a = a/(1.0-a);
        	b = b/(1.0-b);
        	
        	w0 = a * (1.0+b)/(1.0-a*b);
        	w1 = b * (1.0+a)/(1.0-a*b);
        	
        	if ( (i-mNumTimes) == 0 )	// v0
        	{
        		randNumber = log(w0);
        		std::cout << "v0:"<<randNumber << std::endl;
        	}
        	else						// v1
        	{
        		randNumber = log(w1);
        		std::cout << "v1:"<<randNumber << std::endl;
        	}
        	if(randNumber < mLowerBound[i])
        		randNumber = mLowerBound[i];
        	if(randNumber > mUpperBound[i])
        		randNumber = mUpperBound[i];
        	break;
#else
		case 0: // v0
			while(randNumber >= high || randNumber <= low)
        		randNumber = 1.0 - exp_dist_v0( rng );
        	break;
        	
        case 1: // v1
			while(randNumber >= high || randNumber <= low)
        		randNumber = 1.0 - gamma_dist_v1( rng );
        	break;
#endif
        case 2: // w0
			randNumber = beta_dist_w0( rng ); // always between 0 and 1
        	break;
        
        case 3: // kappa
			while(randNumber >= high || randNumber <= low)
        		randNumber = gamma_dist_k( rng );
        	break;
        	
        case 4: // w2
			while(randNumber >= high || randNumber <= low)
        		randNumber = 1.0 + gamma_dist_w2( rng );
        	break;
        
        default:
        	std::cout << "Error: wrong index of i: " << i << std::endl;
        	break;			
		}
	}
	
	return randNumber;
}

// --------------------------------------------------------------------
void BootstrapRandom::bootstrapRandomly(double *f, double *x, unsigned int numTries)
{
	
	*f = evaluateLikelihood(x);
	double *rand_x = &mSpace[0];
	
	//bool optimize_w2 = (mN-mNumTimes) == 5;
    
    // local variables
    size_t step;
    
	for(step = 0; step < numTries; ++step)
	{	
		// generate the random variables
		
		for(int i=0; i<mN; ++i) rand_x[i] = generateRandom(i);
		  	
        
        // verify if it is a better choice
        double rand_f = evaluateLikelihood(rand_x);
        
        if( rand_f > *f )
        {
        	*f = rand_f;
        	memcpy(x, rand_x, size_vect);
        }
	}
}


// --------------------------------------------------------------------
void BootstrapRandom::bootstrapEachDirectionRandomly(double *f, double *x, int numGlobal)
{
	if (numGlobal > 0)
		bootstrapRandomly(f, x, numGlobal);
	
	double rand_f;
	double *rand_x = &mSpace[0];
	memcpy(rand_x, x, size_vect);
	
	size_t numTriesPerVar = static_cast<size_t>(ceil(log(static_cast<double>(mN))));
	
	// look in each direction the "best" variable
	for(int i=0; i<mN; ++i)
	{
		for(unsigned int j=0; j<numTriesPerVar; ++j)
		{
			rand_x[i] = generateRandom(i);
			memcpy(&x_[0], rand_x, size_vect);
			rand_f = mModel->computeLikelihoodForGradient(x_, false, i);
			
			if( rand_f > *f )
			{
				*f = rand_f;
				x[i] = rand_x[i];
			}
			else
				rand_x[i] = x[i];
		}	
	}
}


// --------------------------------------------------------------------
#ifdef BOOTSTRAP_ES
void BootstrapRandom::bootstrapEvolutionStrategy(double *f, double *x, int maxNumGenerations)
{
	// initialize the state of each individual randomly according to the distribution of the variables
	int best_individual = 0;
	
	double fmin, fmax;
	fmin = fmax = -1e16;
	
	for (int individual_id( 0 ); individual_id < mPopSize; ++individual_id)
	{
		double *individual_pos = &mPopPos[individual_id*mN];
		
		// generate the initial position of individuals
		for (int i(0); i<mN; ++i)
			individual_pos[i] = generateRandom(i);
			
		// compute its likelihood
		double ftmp = evaluateLikelihood(individual_pos);
		mPopFitness[individual_id] = ftmp;
		
		// save the min/max
		if ( fmax<mPopFitness[individual_id] )
		{
			fmax = ftmp;
			best_individual = static_cast<int>(individual_id);
		}
		fmin = (ftmp < fmin) ? ftmp : fmin;
	}
	*f = fmax;
	memcpy(x, &mPopPos[best_individual*mN], size_vect);
	
	if ( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		std::cout << "Finished initialization of the population, evolving generations..." << std::endl;
	
	
	int lambda = mPopSize>>2;
	// selection temp storage
	std::vector<int> selected( lambda<<1 );
	// children temp storage
	std::vector<double> children( lambda*mN );
	
	
	// loop evolving the generations	
	for (int generation=0; generation<maxNumGenerations; ++generation)
	{	
		if ( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
			std::cout << "Starting generation " << generation << std::endl;
			
		// --children generation from random parents
		
		/*
		for(int i=0; i<lambda<<1; ++i)
			selected[i] = int( randFrom0to1()*(mPopSize-1) );			
		*/
		for (int i=0; i<lambda; ++i)
		{
			selected[i<<1] = i<<2; //static_cast<int>( randFrom0to1()*(mPopSize-1) );	
			if (randFrom0to1() < 0.9)
				selected[(i<<1)+1] = best_individual;
			else
				selected[(i<<1)+1] = i<<2 + 1;
		}
		
		for (int childid(0); childid<lambda; ++childid)
		{
			const double proportion_parent_1 = randFrom0to1();
			const double proportion_parent_2 = 1.0 - proportion_parent_1;
			
			int parid = selected[childid<<1];
			double* child = &children[childid*mN];
			double* parent = &mPopPos[parid*mN];
			memcpy(child, parent, size_vect);
			dscal_(&mN, &proportion_parent_1, child, &I1);
			
			parid = selected[(childid<<1) + 1];
		
			parent = &mPopPos[parid*mN];
			daxpy_(&mN, &proportion_parent_2, parent, &I1, child, &I1);
		}
		
		// --mutations
		
		for (int i(0); i<lambda*mN; ++i)
		{	
			// apply a "little" mutation
			const double prop = 0.8 + 0.15*randFrom0to1();
			children[i] = prop*children[i] + (1.0-prop)*generateRandom(i%mN);
		}
		
		// --selection
		
		for (int childid(0); childid<lambda; ++childid)
		{
			double *childPos = &children[childid*mN];
			double childFitness = evaluateLikelihood(childPos);
			
			const int selection_shift = static_cast<const int> (randFrom0to1()*static_cast<double>(mPopSize-1));
			for (int selectTry(0); selectTry<mN; ++selectTry)
			{
				const int selected_oponent = selectTry + selection_shift;
				if (childFitness > mPopFitness[selected_oponent])
				{
					memcpy(&mPopPos[selected_oponent*mN], childPos, size_vect);
					mPopFitness[selected_oponent] = childFitness;
					
					if (fmax < childFitness)
					{
						fmax = childFitness;
						best_individual = static_cast<int>(selected_oponent);
					}
					break;
				}
			}
		}
	}
	
	if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		std::cout << "Finished bootstrap!" << std::endl;
		
	// save the best result
	*f = fmax;
	memcpy(x, &mPopPos[best_individual*mN], size_vect);
}
#endif // BOOTSTRAP_ES



