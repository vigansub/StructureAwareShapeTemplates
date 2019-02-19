#include <limits>
#include <memory.h>
#include <QString>
#include <QtGlobal>
#include "desolver.h"
#include "Scene.h"

//see http://www1.icsi.berkeley.edu/~storn/code.html

#define Element(a,b,c)  a[b*m_dimension+c]
#define RowVector(a,b)  (&a[b*m_dimension])
#define CopyVector(a,b) memcpy((a),(b),m_dimension*sizeof(double))

DESolver::DESolver(Scene *scene)
: m_dimension(0), 
  m_populationSize(0),
  m_generations(0), 
  m_strategy(Strategy_Rand1Exp),
  m_scale(0.7), 
  m_probability(0.5), 
  m_bestEnergy(0.0), 
  m_trialEnergy(0.0),
  m_trialSolution(NULL), 
  m_bestSolution(NULL),
  m_populationEnergy(NULL),
  m_population(NULL),
  m_scene(scene)
{
}

DESolver::~DESolver()
{
  free();
}

void DESolver::free()
{
  if (m_trialSolution) {
    delete m_trialSolution;
  }
  if (m_bestSolution) {
    delete m_bestSolution;
  }
  if (m_populationEnergy) {
    delete m_populationEnergy;
  }
  if (m_population) {
    delete m_population;
  }
  m_trialSolution = m_bestSolution = m_populationEnergy = m_population = NULL;
}

void DESolver::initialize(int dimension, int populationSize, double *min, double *max, Strategy strategy, const double& scale, const double& probability)
{
    free();

    m_dimension       = dimension;
    m_populationSize  = populationSize;
    m_strategy	      = strategy;
    m_scale		        = scale;
    m_probability     = probability;

    m_trialSolution     = new double[m_dimension];
    m_bestSolution      = new double[m_dimension];
    m_populationEnergy	= new double[m_populationSize];
    m_population	    = new double[m_populationSize * m_dimension];
	
	for (int i=0; i < m_populationSize; ++i) 
    {
        for (int j=0; j < m_dimension; ++j) 
        {
          Element(m_population, i, j) = random_uniform(min[j], max[j]);
        }
		m_populationEnergy[i] = std::numeric_limits<double>::max();
	}

    for (int i=0; i < m_dimension; ++i) 
    {
        m_bestSolution[i] = 0.0;
    }

}

bool DESolver::solve(int maxGenerations)
{
    emit bestCurrentSolution(10.0);

    m_bestEnergy = std::numeric_limits<double>::max();
	
    bool foundSolution = false;

    for (m_generations = 0; m_generations < maxGenerations && !foundSolution; ++m_generations) 
    {
		for (int candidate = 0; candidate < m_populationSize; ++candidate) 
        {
            switch (m_strategy) 
            {
		        case Strategy_Best1Exp:
			        Best1Exp(candidate);
			        break;
		        case Strategy_Rand1Exp:
			        Rand1Exp(candidate);
			        break;
		        case Strategy_RandToBest1Exp:
			        RandToBest1Exp(candidate);
			        break;
		        case Strategy_Best2Exp:
			        Best2Exp(candidate);
			        break;
		        case Strategy_Rand2Exp:
			        Rand2Exp(candidate);
			        break;
		        case Strategy_Best1Bin:
			        Best1Bin(candidate);
			        break;
		        case Strategy_Rand1Bin:
			        Rand1Bin(candidate);
			        break;
		        case Strategy_RandToBest1Bin:
			        RandToBest1Bin(candidate);
			        break;
		        case Strategy_Best2Bin:
			        Best2Bin(candidate);
			        break;
		        case Strategy_Rand2Bin:
			        Rand2Bin(candidate);
			        break;
	        }

			m_trialEnergy = energy_function(m_trialSolution, foundSolution);

			if (m_trialEnergy < m_populationEnergy[candidate]) 
            {
				// New low for this candidate
				m_populationEnergy[candidate] = m_trialEnergy;
				CopyVector(RowVector(m_population, candidate), m_trialSolution);

				// Check if all-time low
				if (m_trialEnergy < m_bestEnergy) 
                {
					m_bestEnergy = m_trialEnergy;
					CopyVector(m_bestSolution, m_trialSolution);
				}
			}
		}

    //qDebug("Generation %d: %f", m_generations, m_bestEnergy);
    
    QString solutionString;
    for (int i = 0; i < m_dimension; ++i) 
    {
        solutionString.append(QString::number(m_bestSolution[i], 'f', 2));
        if (i < m_dimension - 1) 
        {
            solutionString.append(", ");
        }
    }

    qDebug("Generation %d: %g [%s]", m_generations, m_bestEnergy, qPrintable(solutionString));

    if(m_scene) m_scene->updateGrammarGeometry();
    //emit bestCurrentSolution(10.0);
  }

	return foundSolution;
}

void DESolver::Best1Exp(int candidate)
{
	int r1, r2;
	select_samples(candidate, &r1, &r2);

	int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i = 0; (random_uniform(0.0,1.0) < m_probability) && (i < m_dimension); ++i) 
    {
		m_trialSolution[n] = m_bestSolution[n] + m_scale * (Element(m_population, r1, n) - Element(m_population, r2, n));
		n = (n + 1) % m_dimension;
	}
}

void DESolver::Rand1Exp(int candidate)
{
	int r1, r2, r3;
	select_samples(candidate,&r1,&r2,&r3);
	
    int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
    
    for (int i=0; (random_uniform(0.0,1.0) < m_probability) && (i < m_dimension); ++i) 
    {
		m_trialSolution[n] = Element(m_population,r1,n)	+ m_scale * (Element(m_population,r2,n)	- Element(m_population,r3,n));
		n = (n + 1) % m_dimension;
	}
}

void DESolver::RandToBest1Exp(int candidate)
{
	int r1, r2;
	select_samples(candidate,&r1,&r2);
	
  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i=0; (random_uniform(0.0,1.0) < m_probability) && (i < m_dimension); ++i) {
		m_trialSolution[n] += m_scale * (m_bestSolution[n] - m_trialSolution[n]) + m_scale * (Element(m_population,r1,n) - Element(m_population,r2,n));
		n = (n + 1) % m_dimension;
	}
}

void DESolver::Best2Exp(int candidate)
{
	int r1, r2, r3, r4;
	select_samples(candidate,&r1,&r2,&r3,&r4);
	
  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i=0; (random_uniform(0.0,1.0) < m_probability) && (i < m_dimension); ++i)	{
		m_trialSolution[n] = m_bestSolution[n] + m_scale * (Element(m_population,r1,n) + Element(m_population,r2,n)	- Element(m_population,r3,n) - Element(m_population,r4,n));
		n = (n + 1) % m_dimension;
	}
}

void DESolver::Rand2Exp(int candidate)
{
	int r1, r2, r3, r4, r5;
	select_samples(candidate,&r1,&r2,&r3,&r4,&r5);
	
  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i=0; (random_uniform(0.0,1.0) < m_probability) && (i < m_dimension); ++i) {
		m_trialSolution[n] = Element(m_population,r1,n)	+ m_scale * (Element(m_population,r2,n) + Element(m_population,r3,n) - Element(m_population,r4,n) - Element(m_population,r5,n));
		n = (n + 1) % m_dimension;
	}
}

void DESolver::Best1Bin(int candidate)
{
	int r1, r2;
	select_samples(candidate,&r1,&r2);
	
  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i=0; i < m_dimension; ++i) {
    if ((random_uniform(0.0,1.0) < m_probability) || (i == (m_dimension - 1))) {
      m_trialSolution[n] = m_bestSolution[n] + m_scale * (Element(m_population,r1,n) - Element(m_population,r2,n));
    }
		n = (n + 1) % m_dimension;
	}
}

void DESolver::Rand1Bin(int candidate)
{
	int r1, r2, r3;
	select_samples(candidate,&r1,&r2,&r3);
	
  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i=0; i < m_dimension; ++i) {
    if ((random_uniform(0.0,1.0) < m_probability) || (i  == (m_dimension - 1))) {
      m_trialSolution[n] = Element(m_population,r1,n)	+ m_scale * (Element(m_population,r2,n)	- Element(m_population,r3,n));
    }
		n = (n + 1) % m_dimension;
	}
}

void DESolver::RandToBest1Bin(int candidate)
{
	int r1, r2;
	select_samples(candidate,&r1,&r2);
	
  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i=0; i < m_dimension; ++i) {
    if ((random_uniform(0.0,1.0) < m_probability) || (i  == (m_dimension - 1))) {
      m_trialSolution[n] += m_scale * (m_bestSolution[n] - m_trialSolution[n]) + m_scale * (Element(m_population,r1,n)- Element(m_population,r2,n));
    }
		n = (n + 1) % m_dimension;
	}
}

void DESolver::Best2Bin(int candidate)
{
	int r1, r2, r3, r4;
	select_samples(candidate,&r1,&r2,&r3,&r4);
	
  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
	for (int i=0; i < m_dimension; i++) {
    if ((random_uniform(0.0,1.0) < m_probability) || (i  == (m_dimension - 1))) {
      m_trialSolution[n] = m_bestSolution[n] + m_scale * (Element(m_population,r1,n) + Element(m_population,r2,n)	- Element(m_population,r3,n) - Element(m_population,r4,n));
    }
		n = (n + 1) % m_dimension;
	}
}

void DESolver::Rand2Bin(int candidate)
{
	int r1, r2, r3, r4, r5;
	select_samples(candidate,&r1,&r2,&r3,&r4,&r5);

  int n = (int)random_uniform(0.0,(double)m_dimension);

	CopyVector(m_trialSolution,RowVector(m_population,candidate));
  for (int i=0; i < m_dimension; i++)  {
    if ((random_uniform(0.0,1.0) < m_probability) || (i  == (m_dimension - 1))) {
      m_trialSolution[n] = Element(m_population,r1,n) + m_scale * (Element(m_population,r2,n)	+ Element(m_population,r3,n) - Element(m_population,r4,n)	- Element(m_population,r5,n));
    }
		n = (n + 1) % m_dimension;
	}
}

void DESolver::select_samples(int candidate, int *r1, int *r2, int *r3, int *r4, int *r5)
{
	if (r1) {
		do {
			*r1 = (int)random_uniform(0.0,(double)m_populationSize);
		} while (*r1 == candidate);
	}

	if (r2)	{
		do {
			*r2 = (int)random_uniform(0.0,(double)m_populationSize);
		}	while ((*r2 == candidate) || (*r2 == *r1));
	}

	if (r3) {
		do {
			*r3 = (int)random_uniform(0.0,(double)m_populationSize);
		}	while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
	}

	if (r4) {
		do {
			*r4 = (int)random_uniform(0.0,(double)m_populationSize);
		} while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
	}

	if (r5)	{
		do {
			*r5 = (int)random_uniform(0.0,(double)m_populationSize);
		}	while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)	|| (*r5 == *r2) || (*r5 == *r1));
	}
}

double DESolver::random_uniform(double minValue, double maxValue)
{
  #define SEED 3
  #define IM1 2147483563
  #define IM2 2147483399
  #define AM (1.0/IM1)
  #define IMM1 (IM1-1)
  #define IA1 40014
  #define IA2 40692
  #define IQ1 53668
  #define IQ2 52774
  #define IR1 12211
  #define IR2 3791
  #define NTAB 32
  #define NDIV (1+IMM1/NTAB)
  #define EPS 1.2e-7
  #define RNMX (1.0-EPS)

	long j;
	long k;
	static long idum;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	double result;

  if (iy == 0) {
		idum = SEED;
  }

	if (idum <= 0) {
    if (-idum < 1) {
      idum = 1;
    } else {
			idum = -idum;
    }

		idum2 = idum;

		for (j=NTAB+7; j>=0; j--) {
			k = idum / IQ1;
			idum = IA1 * (idum - k*IQ1) - k*IR1;
      if (idum < 0) {
        idum += IM1;
      }
      if (j < NTAB) {
        iv[j] = idum;
      }
		}

		iy = iv[0];
	}

	k = idum / IQ1;
	idum = IA1 * (idum - k*IQ1) - k*IR1;

  if (idum < 0) {
    idum += IM1;
  }

	k = idum2 / IQ2;
	idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;

  if (idum2 < 0) {
		idum2 += IM2;
  }

	j = iy / NDIV;
	iy = iv[j] - idum2;
	iv[j] = idum;

  if (iy < 1) {
		iy += IMM1;
  }

	result = AM * iy;

  if (result > RNMX) {
		result = RNMX;
  }

	result = minValue + result * (maxValue - minValue);

	return result ;
}

double DESolverTest::energy_function(double *trial, bool &foundSolution)
{
  // Himmelblau function -- 4 points with minimum value 0

  double result = 0;

  for (int i = 0; i < m_dimension; i += 2) {
    double x = trial[i];
	  double y = trial[i + 1];

    double term1 = x * x + y - 11;
    double term2 = x + y * y - 7;

    result += term1 * term1 + term2 * term2;
  }

  if (m_count++ % m_populationSize == 0) {
    qDebug("%d %lf", m_count / m_populationSize + 1, energy());
  }

  foundSolution = (result < 0.0000000001);
	
	return result;
}

void DESolverTest::run()
{
  #define DIMENSION 8
  #define POPULATIONSIZE 100
  #define MAX_GENERATIONS	1000

  double min[DIMENSION];
	double max[DIMENSION];

	for (int i = 0; i < DIMENSION; ++i) {
		max[i] =  100.0;
		min[i] = -100.0;
	}

  DESolverTest solver;

  solver.initialize(DIMENSION, POPULATIONSIZE, min, max, DESolver::Strategy_Best1Exp, 0.9, 1.0);
	
	qDebug("Calculating...");
	solver.solve(MAX_GENERATIONS);

	const double* solution = solver.solution();

	qDebug("\n\nBest Coefficients:");
  for (int i = 0; i < DIMENSION; ++i) {
    qDebug("[%d]: %lf", i, solution[i]);
  }
}
