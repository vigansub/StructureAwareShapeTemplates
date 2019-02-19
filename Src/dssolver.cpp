#include <limits>
#include <QtGlobal>
#include <QString>
#include <QTime>
#include "dssolver.h"

DSSolver::DSSolver(int n, const double& tolerance)
: m_n(n),
  m_p(n + 1),
  m_tolerance(std::min(tolerance, 1.0e-14))
{
  m_bestSolution = new double[m_n];
  m_simplex = new double[m_p * m_n];
  m_energy = new double[m_p];
  m_sum = new double[m_n];
  m_trial = new double[m_n];
}

DSSolver::~DSSolver()
{
  delete[] m_bestSolution;
  delete[] m_simplex;
  delete[] m_energy;
  delete[] m_sum;
  delete[] m_trial;
}

#define simplex_value(j,i) m_simplex[j*m_n+i]

void DSSolver::solve(const double start[], const double& delta, const int maxEvaluations)
{
  const double TINY = 1.0e-10;

  m_evaluations = 0;

  double startEnergy = energy_function(start);

  for (int j = 0; j < m_p; ++j) {
    for (int i = 0; i < m_n; ++i) {
      simplex_value(j, i) = start[i];
    }
    if (j != 0) {
      if (rand() % 2) {
        simplex_value(j,j - 1) += delta;
      } else {
        simplex_value(j,j - 1) -= delta;
      }
    }
    m_energy[j] = energy_function(&m_simplex[j*m_n]);
  }

  compute_sum();

  do {
    //debug_simplex();

    int jHighest = 1;
    int jNextHighest = 0;
    int jLowest = 0;

    if (m_energy[0] > m_energy[1]) {
      jHighest = 0;
      jNextHighest = 1;
    }

    for (int j = 0; j < m_p; ++j) {
      if (m_energy[j] <= m_energy[jLowest]) {
        jLowest = j;
      }
      if (m_energy[j] > m_energy[jHighest]) {
        jNextHighest = jHighest;
        jHighest = j;
      } else if (m_energy[j] > m_energy[jNextHighest] && j != jHighest) {
        jNextHighest = j;
      }
    }

    for (int i = 0; i < m_n; ++i) {
      m_bestSolution[i] = simplex_value(jLowest,i);
    }
    m_bestEnergy = m_energy[jLowest];

    double tolerance = 2.0 * abs(m_energy[jHighest] - m_energy[jLowest]) / (abs(m_energy[jHighest]) + abs(m_energy[jLowest]) + TINY);
    if (tolerance < m_tolerance) {
      //qDebug("Tolerance limit %f hit after %d evaluations", m_tolerance, m_evaluations);
      break;
    }

    if (m_evaluations > maxEvaluations) {
      //qDebug("Evaluations limit %d hit at tolerance %f", maxEvaluations, tolerance);
      break;
    }

    m_evaluations += 2;

    double energy = trial(jHighest, -1.0);
    if (energy <= m_energy[jLowest]) {
      energy = trial(jHighest, 2.0);
    } else if (energy >= m_energy[jNextHighest]) {
      double energySave = m_energy[jHighest];
      energy = trial(jHighest, 0.5);
      if (energy >= energySave) {
        for (int j = 0; j < m_p; ++j) {
          if (j != jLowest) {
            for (int i = 0; i < m_n; ++i ) {
              simplex_value(j,i) = m_sum[i] = 0.5 * (simplex_value(j,i) + simplex_value(jLowest,i));
            }
            m_energy[j] = energy_function(m_sum);
          }
        }
        m_evaluations += m_n;
        compute_sum();
      }
    } else {
      --m_evaluations;
    }

  } while(true);

  if (startEnergy < m_bestEnergy) {
    for (int i = 0; i < m_n; ++i) {
      m_bestSolution[i] = start[i];
    }
  }
}

double DSSolver::trial(const int jHighest, const double& fac)
{
  double fac1 = (1.0 - fac) / m_n;
  double fac2 = fac1 - fac;

  for (int i = 0; i < m_n; ++i) {
    m_trial[i] = m_sum[i] * fac1 - simplex_value(jHighest,i) * fac2;
  }

  double energy = energy_function(m_trial);

  if (energy < m_energy[jHighest]) {
    m_energy[jHighest] = energy;
    for (int i = 0; i < m_n; ++i) {
      m_sum[i] += m_trial[i] - simplex_value(jHighest,i);
      simplex_value(jHighest,i) = m_trial[i];
    }
  }

  return energy;
}

void DSSolver::compute_sum()
{
  for (int i = 0; i < m_n; ++i) {
    m_sum[i] = 0;
    for (int j = 0; j < m_p; ++j) {
      m_sum[i] += simplex_value(j,i);
    }
  }
}

void DSSolver::debug_simplex()
{
  QString string("Simplex: ");
  for (int j = 0; j < m_p; ++j) {
    QString point("(");
    for (int i = 0; i < m_n; ++i) {
      point.append(QString::number(m_simplex[j * m_n + i]));
      if (i != m_n - 1) 
        point.append(",");
    }
    point.append(")=").append(QString::number(m_energy[j]) + " ");
    string.append(point).append(" ");
  }
  qDebug(qPrintable(string));
}

double DSSolverTest::energy_function(const double trial[], bool& foundSolution)
{
  // Himmelblau function

  double result = 0;

  for (int i = 0; i < m_n; i += 2) {
    const double& x = trial[i];
	  const double& y = trial[i + 1];
    double term1 = x * x + y - 11;
    double term2 = x + y * y - 7;
    result += term1 * term1 + term2 * term2;
  }

  foundSolution = (result < 0.000001);
	
	return result;
}

void DSSolverTest::run_test()
{
  const int DIMENSION = 8;
  const double TOLERANCE = 0.00000001;
  const int MAX_EVALUATIONS =	10000;

  const double START[] = {-1, -1, 
                          -1,  1, 
                           1,  1, 
                           1, -1};
  const double DELTA = 0.5;

  DSSolverTest solver(DIMENSION, TOLERANCE);

  qDebug("Calculating...");

  const int RUNS = 10000;

  QTime time;
  time.start();
  for (int i = 0; i < RUNS; ++i) {
    solver.solve(START, DELTA, MAX_EVALUATIONS);
  }

  qDebug("\nTime for %d test runs: %d msecs\n", RUNS, time.elapsed());

  qDebug("\nEvaluations: %d\n", solver.evaluations());

  qDebug("\nEnergy: %lf\n", solver.energy());

  qDebug("\nBest Coefficients:");
	const double* solution = solver.solution();
  for (int i = 0; i < DIMENSION; ++i) {
    qDebug("[%d]: %lf", i, solution[i]);
  }

  qDebug("\n");
}
