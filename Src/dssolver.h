#ifndef DSSOLVER_H
#define DSSOLVER_H

class DSSolver 
{
public:

  DSSolver(const int n, const double& tolerance);
  
  ~DSSolver();

  void  solve(const double start[], const double& delta, const int maxEvaluations);

  inline
	int n()
  { 
    return m_n; 
  }

  inline
  const double& energy() const
  { 
    return m_bestEnergy; 
  }

  inline
  const double* solution() const
  {
    return m_bestSolution;
  }

  inline
	int evaluations() const
  { 
    return m_evaluations; 
  }

protected:

  virtual double energy_function(const double trial[], bool& foundSolution) = 0;

  inline
  double energy_function(const double trial[])
  {
    bool foundSolution;
    return energy_function(trial, foundSolution);
  }

  int     m_n;

private:

  void  debug_simplex();

  void    compute_sum();
  double  trial(const int jHighest, const double& factor);

  int     m_p;
  double  m_tolerance;

  double  m_bestEnergy;
  int     m_evaluations;
  
  double* m_bestSolution;
  double* m_simplex;
  double* m_energy;
  double* m_sum;
  double* m_trial;

};

class DSSolverTest : public DSSolver
{
public:

  DSSolverTest(int n, const double& epsilon)
    : DSSolver(n, epsilon)
  {
  }

  static void run_test();

private:

  double energy_function(const double trial[], bool& foundSolution);

};

#endif  // DSSOLVER_H
