#include "BootstrapRandom.h"
#include "MathSupport.h"

// --------------------------------------------------------------------
double BootstrapRandom::bootstrap(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	allocateMemory();
	
	double likelihood_value = -1000000;
	
	BootstrapType btype = EVOLUTION_STRATEGY;
	//BootstrapType btype = RANDOM_TRIES_SEPARATE_VARS;
	switch(btype)
	{
	case ONLY_RANDOM_TRIES:
		bootstrapRandomly(&likelihood_value, &aVars[0], mN);
		break;
		
	case RANDOM_TRIES_SEPARATE_VARS:
		bootstrapEachDirectionRandomly(&likelihood_value, &aVars[0], mN);
		bootstrapEachDirectionRandomly(&likelihood_value, &aVars[0], 0);
		break;
		
	case EVOLUTION_STRATEGY:
		//int numGenerations = (mN < 20) ? 0 : ((mN>60) ? 100 : 3);
		int num_generations = static_cast<int> ( static_cast<double>(mN) / 7.0 - 4.0 );
		num_generations = num_generations > 0 ? num_generations : 0;
		bootstrapEvolutionStrategy(&likelihood_value, &aVars[0], num_generations);
		break;
	};
	return -likelihood_value;
}


// --------------------------------------------------------------------
double BootstrapRandom::evaluateLikelihood(const double *aX)
{
	memcpy(&mVarsCopy[0], aX, mSizeVect);
	return evaluateLikelihood(mVarsCopy);
}


// --------------------------------------------------------------------
double BootstrapRandom::evaluateLikelihood(const std::vector<double> &aX)
{
	return mModel->computeLikelihood(aX, mTrace);
}


// --------------------------------------------------------------------
void BootstrapRandom::allocateMemory(void)
{
	mSizeVect = mN*sizeof(double);
	
	mVarsCopy.resize(mN);
	mSpace.resize(mN);
	
#ifdef BOOTSTRAP_ES
	
	// choose a population size
#if 1
	// we take here a population of:
	// - 15 for small problems
	// - 70 for medium to large problems
	// - linear in between 
	if (mN < 10)
		mPopSize = 15;
	else if (mN < 50)
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
double BootstrapRandom::generateRandom(unsigned int aIndexVar)
{
	double randNumber, low, high;
	
	high = mUpperBound[aIndexVar];
    low  = mLowerBound[aIndexVar];
    randNumber = high+1.0;
    
	// branch lengths
	if (aIndexVar < mNumTimes)
	{
		while(randNumber >= high || randNumber <= low)
        	randNumber = mGammaDistT( mUnifRandNumGenerator );
	}
	else
	{
		switch(aIndexVar-mNumTimes)
		{
#ifdef USE_ORIGINAL_PROPORTIONS
		case 0:
        case 1:
        	double v0(2.0),v1(2.0);
        	low = 0.0;
        	high = 1.0;
        	while(v0 >= high || v0 <= low)
        		v0 = 1.0 - mExpDistV0( mUnifRandNumGenerator );
			while(v1 >= high || v1 <= low)
        		v1 = 1.0 - mGammaDistV1( mUnifRandNumGenerator );
        	
        	// convert to original proportion variables
        	double a,b, w0, w1;
        	
        	a = v0*v1;
        	b = v0*(1.0-v1);
        	a = a/(1.0-a);
        	b = b/(1.0-b);
        	
        	w0 = a * (1.0+b)/(1.0-a*b);
        	w1 = b * (1.0+a)/(1.0-a*b);
        	
        	if ( (aIndexVar-mNumTimes) == 0 )	// v0
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
        		randNumber = 1.0 - mExpDistV0( mUnifRandNumGenerator );
        	break;
        	
        case 1: // v1
			while(randNumber >= high || randNumber <= low)
        		randNumber = 1.0 - mGammaDistV1( mUnifRandNumGenerator );
        	break;
#endif
        case 2: // w0
			randNumber = mBetaDistW0( mUnifRandNumGenerator ); // always between 0 and 1
        	break;
        
        case 3: // kappa
			while(randNumber >= high || randNumber <= low)
        		randNumber = mGammaDistK( mUnifRandNumGenerator );
        	break;
        	
        case 4: // w2
			while(randNumber >= high || randNumber <= low)
        		randNumber = 1.0 + mGammaDistW2( mUnifRandNumGenerator );
        	break;
        
        default:
        	std::cout << "Error: wrong index of variable aIndexVar: " << aIndexVar << std::endl;
        	break;			
		}
	}
	
	return randNumber;
}

// --------------------------------------------------------------------
void BootstrapRandom::bootstrapRandomly(double *aF, double *aX, unsigned int aNumTries)
{
	
	*aF = evaluateLikelihood(aX);
	double *rand_x = &mSpace[0];
    
	for (unsigned int step = 0; step < aNumTries; ++step)
	{	
		// generate the random variables
		for (int i=0; i<mN; ++i) rand_x[i] = generateRandom(i);
        
        // verify if it is a better choice
        double rand_f = evaluateLikelihood(rand_x);
        
        if( rand_f > *aF )
        {
        	*aF = rand_f;
        	memcpy(aX, rand_x, mSizeVect);
        }
	}
}


// --------------------------------------------------------------------
void BootstrapRandom::bootstrapEachDirectionRandomly(double *aF, double *aX, int aNumGlobal)
{
	if (aNumGlobal > 0)
		bootstrapRandomly(aF, aX, aNumGlobal);
	
	double rand_f;
	double *rand_x = &mSpace[0];
	memcpy(rand_x, aX, mSizeVect);
	
	const unsigned int num_tries_per_var = static_cast<unsigned int>(ceil(log(static_cast<double>(mN))));
	
	// look in each direction the "best" variable
	for (int i=0; i<mN; ++i)
	{
		for (unsigned int j=0; j<num_tries_per_var; ++j)
		{
			rand_x[i] = generateRandom(i);
			memcpy(&mVarsCopy[0], rand_x, mSizeVect);
			rand_f = mModel->computeLikelihoodForGradient(mVarsCopy, false, i);
			
			if( rand_f > *aF )
			{
				*aF = rand_f;
				aX[i] = rand_x[i];
			}
			else
				rand_x[i] = aX[i];
		}	
	}
}


// --------------------------------------------------------------------
#ifdef BOOTSTRAP_ES
void BootstrapRandom::bootstrapEvolutionStrategy(double *aF, double *aX, int aMaxNumGenerations)
{
	if (aMaxNumGenerations == 0)
		return;
		
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
	*aF = fmax;
	memcpy(aX, &mPopPos[best_individual*mN], mSizeVect);
	
	if ( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		std::cout << "Finished initialization of the population, evolving generations..." << std::endl;
	
	
	int lambda = mPopSize>>2;
	// selection temp storage
	std::vector<int> selected( lambda<<1 );
	// children temp storage
	std::vector<double> children( lambda*mN );
	
	
	// loop evolving the generations	
	for (int generation=0; generation<aMaxNumGenerations; ++generation)
	{	
		if ( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
			std::cout << "Starting generation " << generation << std::endl;
			
		// --children generation from distinct parents / dominant
		for (int i=0; i<lambda; ++i)
		{
			// selected[i<<1] = static_cast<int>( randFrom0to1()*static_cast<double>(mPopSize-1) ); // i<<2; //
			
			if (randFrom0to1() < 0.5)
				selected[i<<1] = i<<2;
			else
				selected[i<<1] = best_individual;
				
			if (selected[i<<1] != best_individual)
			{
				if (randFrom0to1() < 0.9)
					selected[(i<<1)+1] = best_individual;
				else
					selected[(i<<1)+1] = (i<<2) + 1;
			}
		}
		
		for (int childid(0); childid<lambda; ++childid)
		{
			const int parid1 = selected[childid<<1];
			const int parid2 = selected[(childid<<1)+1];
			
			double proportion_parent_1 = mPopFitness[parid2] - fmax;
			double proportion_parent_2 = mPopFitness[parid1] - fmax;
			double tot = proportion_parent_1 + proportion_parent_2;
			proportion_parent_1 /= tot;
			proportion_parent_2 /= tot;
			
			double* child = &children[childid*mN];
			double* parent = &mPopPos[parid1*mN];
			memcpy(child, parent, mSizeVect);
			dscal_(&mN, &proportion_parent_1, child, &I1);
			
			parent = &mPopPos[parid2*mN];
			daxpy_(&mN, &proportion_parent_2, parent, &I1, child, &I1);
		}
		
		// --mutations
		
		for (int i(0); i<lambda*mN; ++i)
		{	
			// apply a "little" mutation
			const double prop = 0.9 + 0.1*randFrom0to1();
			children[i] = prop*children[i] + (1.0-prop)*generateRandom(i%mN);
		}
		
		// --selection
		
		for (int childid(0); childid<lambda; ++childid)
		{
			double *child_pos = &children[childid*mN];
			double child_fitness = evaluateLikelihood(child_pos);
			
			const int selection_shift = static_cast<const int> (randFrom0to1()*static_cast<double>(mPopSize-1));
			for (int select_try(0); select_try<mN; ++select_try)
			{
				const int selected_oponent = select_try + selection_shift;
				if (child_fitness > mPopFitness[selected_oponent])
				{
					memcpy(&mPopPos[selected_oponent*mN], child_pos, mSizeVect);
					mPopFitness[selected_oponent] = child_fitness;
					
					if (fmax < child_fitness)
					{
						fmax = child_fitness;
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
	*aF = fmax;
	memcpy(aX, &mPopPos[best_individual*mN], mSizeVect);
}
#endif // BOOTSTRAP_ES



