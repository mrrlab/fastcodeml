#include "BootstrapRandom.h"
#include "MathSupport.h"
//#include <iomanip>

// --------------------------------------------------------------------
double BootstrapRandom::bootstrap(std::vector<double>& aVars)
{
	mN = static_cast<int>(aVars.size());
	
	if ((mInitStatus & BranchSiteModel::INIT_TIMES_FROM_FILE) == BranchSiteModel::INIT_TIMES_FROM_FILE) // if times initialized from file
		mIndexBegin = mNumTimes;
	else
		mIndexBegin = 0;
	if ((mInitStatus & BranchSiteModel::INIT_PARAMS) == BranchSiteModel::INIT_PARAMS) // if parameters (not times) initialized from command line
		mIndexEnd = mNumTimes;
	else
		mIndexEnd = mN;
	
	if (mIndexBegin == mIndexEnd) // nothing to change, don't perform the bootstrap
		return -evaluateLikelihood(aVars);
	
	allocateMemory();
	
	double likelihood_value = -1000000;
#ifdef BOOTSTRAP_GA
	int num_generations = static_cast<int> ( static_cast<double>(mN) / 7.0 - 4.0 );
	num_generations = num_generations > 0 ? num_generations : 0;
	if (mIndexEnd - mIndexBegin < 6)
	{
		num_generations = (mN > 30) ? 1:0;
	}
	bootstrapGeneticAlgorithm(&likelihood_value, &aVars[0], num_generations);
#else
	int num_generations = static_cast<int> ( static_cast<double>(mN) / 7.0 - 4.0 );
	num_generations = num_generations > 0 ? num_generations : 0;
	if (mIndexEnd - mIndexBegin < 6)
	{
		num_generations = (mN > 30) ? 1:0;
	}
	bootstrapParticlSwarm(&likelihood_value, &aVars[0], num_generations);
#endif // BOOTSTRAP_GA

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
	
#ifdef BOOTSTRAP_GA
	
	// choose a population size
	// we take here a population of:
	// - 15 for small problems
	// - 70 for medium to large problems
	// - linear in between 
	if (mN < 10)
		mPopSize = 15;
	else if (mN < 50)
		mPopSize = static_cast<int>(20 + static_cast<float>(mN-10)*1.25);
	else
		mPopSize = 70;	
	mGASpace.resize( mPopSize*(mN+1) );
	
	mPopFitness = &mGASpace[0];
	mPopPos 	= mPopFitness+mPopSize;
	
#endif // BOOTSTRAP_GA

#ifdef BOOTSTRAP_PSO
	
	// choose a population size
	// we take here a population of:
	// - 15 for small problems
	// - 70 for medium to large problems
	// - linear in between 
	if (mN < 10)
		mPopSize = 15;
	else if (mN < 50)
		mPopSize = static_cast<int>(20 + static_cast<float>(mN-10)*0.75);
	else
		mPopSize = 60;	
		
	mPSOSpace.resize( 3*mPopSize*mN + 3*mN );
	
	mWorkSpace		= &mPSOSpace[0];
	mPositions 		= mWorkSpace + mN;
	mBestPositions	= mPositions + mPopSize*mN;
	mFitnesses		= mBestPositions + mPopSize*mN;
	mBestFitnesses	= mFitnesses + mN;
	mVelocities		= mBestFitnesses + mN;
	
#endif // BOOTSTRAP_PSO
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
        	else								// v1
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
#ifdef BOOTSTRAP_GA
void BootstrapRandom::bootstrapGeneticAlgorithm(double *aF, double *aX, int aMaxNumGenerations)
{
	if (aMaxNumGenerations == 0)
		return;
		
	// initialize the state of each individual randomly according to the distribution of the variables
	int best_individual = 0;
	
	double fmin, fmax;
	fmin = fmax = -1e16;
	
	for (int individual_id(0); individual_id < mPopSize; ++individual_id)
	{
		double *individual_pos = &mPopPos[individual_id*mN];
		memcpy(individual_pos, aX, mSizeVect);
		
		// generate the initial position of individuals
		// add mutations to the tree file
		const double prop_file_init = 0.8;
		for (unsigned int i(0); i<mIndexBegin; ++i)
			individual_pos[i] = prop_file_init*generateRandom(i) + (1.-prop_file_init)*aX[i];
		for (unsigned int i(mIndexBegin); i<mIndexEnd; ++i)
			individual_pos[i] = generateRandom(i);
		for (int i(mIndexEnd); i<mN; ++i)
			individual_pos[i] = prop_file_init*generateRandom(i) + (1.-prop_file_init)*aX[i];

		// compute its likelihood
		double ftmp = evaluateLikelihood(individual_pos);
		mPopFitness[individual_id] = ftmp;
		
		// save the min/max
		if ( fmax<mPopFitness[individual_id] )
		{
			fmax = ftmp;
			best_individual = individual_id;
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
			std::cout << "\tStarting generation " << generation << std::endl;
			
		// --children generation from distinct parents / dominant
		for (int i=0; i<lambda; ++i)
		{
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
			else
			{
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
		// allow mutations on every variable, even if it has been initialized from the files/command line
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
		if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
			std::cout << "\t\tBest individual: f = " << fmax << std::endl;
	}
	
	if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		std::cout << "Finished bootstrap!" << std::endl;
		
	// save the best result
	*aF = fmax;
	memcpy(aX, &mPopPos[best_individual*mN], mSizeVect);
}
#endif // BOOTSTRAP_GA


// --------------------------------------------------------------------
#ifdef BOOTSTRAP_PSO
void BootstrapRandom::bootstrapParticlSwarm(double *aF, double *aX, int aMaxNumGenerations)
{
	// don't proceed to bootstrap if aMaxNumGenerations is 0
	if (aMaxNumGenerations == 0)
		return;
	
	if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		std::cout << "\tInitialize the swarm..." << std::endl;
	
	// initialize the population
	for (int individual_id(0); individual_id<mPopSize; ++individual_id)
	{
		double *individual_pos = mPositions + individual_id*mN;
		double *individual_best_pos = mBestPositions + individual_id*mN;
		double *individual_velocity = mVelocities + individual_id*mN;
		// generate the initial position of individuals
		for (int i(mIndexBegin); i<mIndexEnd; ++i)
			individual_pos[i] = generateRandom(i);
		memcpy(individual_best_pos, individual_pos, mN*sizeof(double));
		// generate velocities
		for (int i(0); i<mN; ++i)
		{
			individual_velocity[i] = generateRandom(i) - individual_pos[i];
		}
		// compute log-likelihood
		double f = evaluateLikelihood(individual_pos);
		mFitnesses[individual_id]		= f;
		mBestFitnesses[individual_id]	= f;
		// save the best
		if (f >= *aF)
		{
			*aF = f;
			memcpy(aX, individual_pos, mN*sizeof(double));
		}
	}
	
	if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
		std::cout << "\tEvolving swarm..." << std::endl;
	
	// main loop
	for (int generation_id(0); generation_id < aMaxNumGenerations; ++ generation_id)
	{
		if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
			std::cout << "\t\tGeneration " << generation_id << std::endl;

		const double dt = 1.0;
		const double inverse_dt = 1.0 / dt;
		
		// evolve the velocities
		for (int individual_id(0); individual_id<mPopSize; ++individual_id)
		{	
			double *individual_pos = mPositions + individual_id*mN;
			double *individual_best_pos = mBestPositions + individual_id*mN;
			double *individual_velocity = mVelocities + individual_id*mN;
			
			const double inertia 		= 0.2 + 0.7/(1.0+static_cast<double>(generation_id));//1.4 * pow(0.9, 1.0+static_cast<double>(generation_id)); //
			const double trust_memory 	= 1.5 * sqrt(randFrom0to1()) * inverse_dt;
			const double trust_best		= 2.5 * sqrt(randFrom0to1()) * inverse_dt;
			
			// inertia term
			dscal_(&mN, &inertia, individual_velocity, &I1);
			
			// confidence in memory
			memcpy(mWorkSpace, individual_best_pos, mN*sizeof(double));
			daxpy_(&mN, &minus_one, individual_pos, &I1, mWorkSpace, &I1);
			daxpy_(&mN, &trust_memory, mWorkSpace, &I1, individual_velocity, &I1);
			
			// confidence in the swarm
			memcpy(mWorkSpace, aX, mN*sizeof(double));
			daxpy_(&mN, &minus_one, individual_pos, &I1, mWorkSpace, &I1);
			daxpy_(&mN, &trust_best, mWorkSpace, &I1, individual_velocity, &I1);
		}
		
		// update positions
		int total_num_variables = mN*mPopSize;
		daxpy_(&total_num_variables, &dt, mVelocities, &I1, mPositions, &I1);
		
		// treat particles outside the bounds
		//#pragma omp parallel for
		for (int id(0); id<mPopSize; ++id)
		{
			int pos_id = mN*id;
			//bool reset_velocity = false;
			for (int bound_id(0); bound_id<mN; ++bound_id)
			{
				double x = mPositions[pos_id];

				double v = mVelocities[pos_id];
				if ( x <= mLowerBound[bound_id] )
				{
					x = mLowerBound[bound_id];
					v = (v > 0.0) ? v : 0.0;
				}
				if ( x >= mUpperBound[bound_id] )
				{
					x = mUpperBound[bound_id];
					v = (v > 0.0) ? 0.0 : v;
				}
				mVelocities[pos_id] = v;
				mPositions[pos_id++] = x;
			}
		}
		
		// update fitness
		for (int individual_id(0); individual_id<mPopSize; ++individual_id)
		{
			double *individual_pos = mPositions + individual_id*mN;
			double f = evaluateLikelihood(individual_pos);
			mFitnesses[individual_id] = f;
			if (f >=  mBestFitnesses[individual_id])
			{
				mBestFitnesses[individual_id] = f;
				memcpy(&mBestPositions[individual_id*mN], individual_pos, mN*sizeof(double));
				if (f >= *aF) 
				{
					*aF = f;
					memcpy(aX, individual_pos, mN*sizeof(double));
					if( mVerbose >= VERBOSE_MORE_INFO_OUTPUT )
						std::cout << "\t\tNew f PSO: " << f << " at step " << generation_id << std::endl;
				}
			}
		}
	}
}
#endif // BOOTSTRAP_PSO




