#ifndef DESOLVER_H
#define DESOLVER_H

#include "Headers.h"

//see http://www1.icsi.berkeley.edu/~storn/code.html

class Scene;

class DESolver : public QObject
{
    Q_OBJECT

public:

  enum Strategy
  {
    Strategy_Best1Exp, Strategy_Rand1Exp, 
    Strategy_RandToBest1Exp, Strategy_Best2Exp, 
    Strategy_Rand2Exp, Strategy_Best1Bin, 
    Strategy_Rand1Bin, Strategy_RandToBest1Bin, 
    Strategy_Best2Bin, Strategy_Rand2Bin
  };

  DESolver(Scene *scene);
  ~DESolver();
	
  // initialize() must be called before solve to set min, max, strategy etc.
  void initialize(int dimension, int populationSize, double min[], double max[], Strategy strategy,	const double& scale, const double& probability);

  // Solve() returns true if energy() returns true.
  // Otherwise it runs maxGenerations generations and returns false.
  virtual bool solve(int maxGenerations);
	
  inline int dimension()
  { 
    return m_dimension; 
  }
	
  inline int populationsize() 
  { 
    return m_populationSize; 
  }

	// Call these functions after Solve() to get results.
  
  inline const double& energy() const
  { 
    return m_bestEnergy; 
  }

  inline const double* solution() const
  { 
    return m_bestSolution; 
  }

  inline int generations() const
  { 
    return m_generations; 
  }

protected:

  // 'energy' must be overridden for problem to solve
	// testSolution[] is array of size 'm_dimension' for a candidate solution
	// setting 'bAtSolution' = true indicates solution is found
	// and solve() immediately returns true.
	virtual double energy_function(double trial[], bool &foundSolution) = 0;

	int m_dimension;
	int m_populationSize;
	int m_generations;

//private:
public:
    Scene *m_scene;


    void Best1Exp(int candidate);
	void Rand1Exp(int candidate);
	void RandToBest1Exp(int candidate);
	void Best2Exp(int candidate);
	void Rand2Exp(int candidate);
	void Best1Bin(int candidate);
	void Rand1Bin(int candidate);
	void RandToBest1Bin(int candidate);
	void Best2Bin(int candidate);
	void Rand2Bin(int candidate);

    void  free();
    void  select_samples(int candidate, int* r1, int* r2 = 0, int* r3 = 0, int* r4 = 0, int* r5 = 0);

	static double random_uniform(double min, double max);

    Strategy  m_strategy;
	double    m_scale;
    double    m_probability;
	double    m_trialEnergy;
	double    m_bestEnergy;
	double*   m_trialSolution;
	double*   m_bestSolution;
	double*   m_populationEnergy;
	double*   m_population;

signals: 
    void bestCurrentSolution(double solution);

};

class DESolverTest : public DESolver
{
public:

  DESolverTest() : DESolver(nullptr), m_count(0)
  {
  }
	
  static void run();

private:

  double energy_function(double trial[], bool &foundSolution);
	
  int m_count;

};

#endif  // DESOLVER_H
