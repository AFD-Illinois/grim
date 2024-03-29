#include "physics.hpp"

/* Pass contravariant vectors */
void fluidElement::tetradConToCoordCon(const array vTetrad[NDIM],
                                       array vCoord[NDIM])
{
  for (int mu = 0; mu < NDIM; mu++)
  {
    vCoord[mu] = zero;
    for (int nu = 0; nu < NDIM; nu++)
    {
      vCoord[mu] += eCon[nu][mu]*vTetrad[nu];
    }
  }
}

/* Pass contravariant vectors */
void fluidElement::coordConToTetradCon(const array vCoord[NDIM],
                                       array vTetrad[NDIM])
{
  for (int mu = 0; mu < NDIM; mu++) 
  {
    vTetrad[mu] = zero;
    for (int nu = 0; nu < NDIM; nu++)
    {
      vTetrad[mu] += eCov[mu][nu]*vCoord[nu];
    }
  }
}

void fluidElement::constructTetrads()
{
  const double SMALL_VECTOR = 1.e-30;
  /* BASIS VECTOR 0 */
  /* Adopt time component parallel to uCon*/
  for (int mu = 0; mu < NDIM; mu++)
  {
    eCon[0][mu] = uCon[mu];
  }
  normalize(eCon[0]);
  
  /* BASIS VECTOR 1 */
  /* Adopt normalize bCon, enforce sanity checks */
  array trial[NDIM];
  for (int mu = 0; mu < NDIM; mu++)
  {
    trial[mu] = bCon[mu];
  }
  
  normalize(trial);
  
  /* Check norm of trial vector */
  array norm = zero;
  for (int mu = 0; mu < NDIM; mu++)
  {
    for (int nu = 0; nu < NDIM; nu++)
    {
      norm += trial[mu]*trial[nu]*geom->gCov[mu][nu];
    }
  }
  
  /* Default to X1 direction for bad trial vectors */
  array nanTrialIndices = where(af::isNaN(norm));
  int numNanIndices = nanTrialIndices.elements();
  if (numNanIndices > 0)
  {
    for (int mu = 0; mu < NDIM; mu++)
    {
      trial[mu](nanTrialIndices) = DELTA(mu, 1);
    }
  }
  
  /* Use trial vector */
  for (int mu = 0; mu < NDIM; mu++)
  {
    eCon[1][mu] = trial[mu];
  }
  
  /* Check to see if trial is parallel to X2 or X3 */
  array X2v[NDIM];
  array X3v[NDIM];
  for (int mu=0; mu < NDIM; mu++)
  {
    X2v[mu] = zero;
    X3v[mu] = zero;
  }
  array trial_X2[NDIM];
  array trial_X3[NDIM];
  for (int mu = 0; mu < NDIM; mu++)
  {
    X2v[mu] = DELTA(mu,2);
    X3v[mu] = DELTA(mu,3);
    trial_X2[mu] = trial[mu];
    trial_X3[mu] = trial[mu];
  }
  projectOut(X2v, trial_X2);
  projectOut(X3v, trial_X3);
  array norm_X2 = zero;
  array norm_X3 = zero;
  for (int mu = 0; mu < NDIM; mu++)
  {
    for (int nu = 0; nu < NDIM; nu++)
    {
      norm_X2 += trial_X2[mu]*trial_X2[nu]*geom->gCov[mu][nu];
      norm_X3 += trial_X3[mu]*trial_X3[nu]*geom->gCov[mu][nu];
    }
  }
  
  array parallelIndices_X2 = where(af::abs(norm_X2) <= SMALL_VECTOR);
  array parallelIndices_X3 = where(af::abs(norm_X3) <= SMALL_VECTOR);
  int numParallelIndices_X2 = parallelIndices_X2.elements();
  int numParallelIndices_X3 = parallelIndices_X3.elements();
  
  /* Project out eCon[0] */
  projectOut(eCon[0], eCon[1]);
  normalize(eCon[1]);
  
  /* BASIS VECTOR 2 */
  /* Adopt X2 unit basis vector */
  for (int mu = 0; mu < NDIM; mu++)
  {
    eCon[2][mu] = DELTA(mu, 2);
  }
  /* If eCon[1] parallel to X2, use X1 */
  if (numParallelIndices_X2 > 0)
  {
    for (int mu = 0; mu < NDIM; mu++)
    {
      eCon[2][mu](parallelIndices_X2) = DELTA(mu, 1);
    }
  }
  
  /* Project out eCon[0] and eCon[1] */
  projectOut(eCon[0], eCon[2]);
  projectOut(eCon[1], eCon[2]);
  normalize(eCon[2]);
  
  /* BASIS VECTOR 3 */
  /* Adopt X3 unit basis vector */
  for (int mu = 0; mu < NDIM; mu++)
  {
    eCon[3][mu] = DELTA(mu, 3);
  }
  /* If eCon[1] parallel to X3, use X1 */
  if (numParallelIndices_X3 > 0)
  {
    for (int mu = 0; mu < NDIM; mu++)
    {
      eCon[3][mu](parallelIndices_X3) = DELTA(mu, 1);
    }
  }
  
  /* Project out eCon[0], eCon[1], and eCon[2] */
  projectOut(eCon[0], eCon[3]);
  projectOut(eCon[1], eCon[3]);
  projectOut(eCon[2], eCon[3]);
  normalize(eCon[3]);
  
  /* Make dual transformation matrix. Reset eCov first. */
  for (int mu = 0; mu < NDIM; mu++)
  {
    for (int nu = 0; nu < NDIM; nu++)
    {
      eCov[mu][nu] = 0.;
    }
  }
  
  /* Lower coordinate basis index
   *   eCov[mu][nu] = eCon[mu][lambda]*gCov[lambda][nu] 
   */
  for (int mu = 0; mu < NDIM; mu++)
  {
    for (int nu = 0; nu < NDIM; nu++)
    {
      for (int lambda = 0; lambda < NDIM; lambda++)
      {
        eCov[mu][nu] += eCon[mu][lambda]*geom->gCov[lambda][nu];
      }
    }
  }
  
  /* Raise tetrad basis index */
  for (int nu = 0; nu < NDIM; nu++)
  {
    eCov[0][nu] *= -1.;
  }
}

void fluidElement::normalize(array vCon[NDIM])
{
  array norm = zero;
  for (int mu = 0; mu < NDIM; mu++)
  {
    for (int nu = 0; nu < NDIM; nu++)
    {
      norm += vCon[mu]*vCon[nu]*geom->gCov[mu][nu];
    }
  }
  
  norm = af::sqrt(af::abs(norm));
  
  for (int mu = 0; mu < NDIM; mu++)
  {
    vCon[mu] /= norm;
  }
}

void fluidElement::projectOut(const array vConB[NDIM], 
                              array vConA[NDIM])
{
  array aDotB = zero;
  array vConBSquare = zero;
  
  for (int mu = 0; mu < NDIM; mu++)
  {
    for (int nu = 0; nu < NDIM; nu++)
    {
      vConBSquare += vConB[mu]*vConB[nu]*geom->gCov[mu][nu];
    }
  }
  
  for (int mu = 0; mu < NDIM; mu++)
  {
    for (int nu = 0; nu < NDIM; nu++)
    {
      aDotB += vConA[mu]*vConB[nu]*geom->gCov[mu][nu];
    }
  }
  
  for (int mu = 0; mu < NDIM; mu++)
  {
    vConA[mu] -= vConB[mu]*aDotB/vConBSquare;
  }
}

/* Only needed for Compton scattering -- maybe should be moved */
void fluidElement::normalizeNull(array vCon[NDIM])
{
  array A = geom->gCov[0][0];
  array B = zero;
  array C = zero;
  
  for (int mu = 1; mu < NDIM; mu++)
  {
    B += 2.*geom->gCov[mu][0]*vCon[mu];
  }
  
  for (int mu = 1; mu < NDIM; mu++)
  {
    for (int nu = 1; nu < NDIM; nu++)
    {
      C += geom->gCov[mu][nu]*vCon[mu]*vCon[nu];
    }
  }
  
  vCon[0] = (-B - af::sqrt(af::abs(B*B - 4.*A*C)))/(2.*A);
}
