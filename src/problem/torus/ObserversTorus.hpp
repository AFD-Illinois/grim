#ifndef GRIM_OBSTORUS_HPP
#define GRIM_OBSTORUS_HPP

#include "../problem.hpp"

void ComputeEnergyIntegrals(fluidElement* elemObs, grid* primObs, geometry* geomObs, const double volElem)
{
  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
  af::seq domainX1 = *primObs->domainX1;
  af::seq domainX2 = *primObs->domainX2;
  af::seq domainX3 = *primObs->domainX3;

  // Integrate baryon mass
  array MassIntegrand = primObs->vars[vars::RHO]*volElem*geomObs->g*elemObs->gammaLorentzFactor/geomObs->alpha;;
  array BaryonMass_af = af::sum(af::flat(MassIntegrand(domainX1, domainX2, domainX3)),0);
  double BaryonMass = BaryonMass_af.host<double>()[0];
  /* Communicate to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  BaryonMass += temp;
	}
    }
  else
    {
      MPI_Send(&BaryonMass, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&BaryonMass,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  // Integrate magnetic energy
  array EMagIntegrand = volElem*elemObs->bSqr*geomObs->g*elemObs->gammaLorentzFactor/geomObs->alpha;
  array EMag_af = af::sum(af::flat(EMagIntegrand(domainX1, domainX2, domainX3)),0);
  double EMag = EMag_af.host<double>()[0];
  /* Communicate to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  EMag += temp;
	}
    }
  else
    {
      MPI_Send(&EMag, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&EMag,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  // Integrate thermal energy
  array EThIntegrand = volElem*primObs->vars[vars::U]*geomObs->g*elemObs->gammaLorentzFactor/geomObs->alpha;
  array ETh_af = af::sum(af::flat(EThIntegrand(domainX1, domainX2, domainX3)),0);
  double ETh = ETh_af.host<double>()[0];
  /* Communicate to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  ETh += temp;
	}
    }
  else
    {
      MPI_Send(&ETh, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&ETh,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD,"Baryon Mass = %e; Magnetic Energy = %e; Thermal Energy = %e;\n",BaryonMass,EMag,ETh);
}

void ComputeMinMaxVariables(fluidElement* elemObs, grid* primObs, geometry* geomObs)
{
  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
  af::seq domainX1 = *primObs->domainX1;
  af::seq domainX2 = *primObs->domainX2;
  af::seq domainX3 = *primObs->domainX3;
  //Find maximum density
  array rhoMax_af = af::max(af::max(af::max(primObs->vars[vars::RHO](domainX1,domainX2,domainX3),2),1),0);
  double rhoMax = rhoMax_af.host<double>()[0];
  /* Communicate rhoMax to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  if(rhoMax < temp)
	    rhoMax = temp;
	}
    }
  else
    {
      MPI_Send(&rhoMax, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&rhoMax,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  
  // Find minimum beta
  const array& bSqr = elemObs->bSqr;
  const array& Pgas = elemObs->pressure;
  array PlasmaBeta = 2.*(Pgas+1.e-13)/(bSqr+1.e-18);
  array BetaMin_af = af::min(af::min(af::min(PlasmaBeta(domainX1,domainX2,domainX3),2),1),0);
  double betaMin = BetaMin_af.host<double>()[0];
  /* Use MPI to find minimum over all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  if( betaMin > temp)
	    betaMin = temp;
	}
    }
  else
    {
      MPI_Send(&betaMin, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&betaMin,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD,"rhoMax = %e; betaMin = %e;\n",rhoMax,betaMin);
}

void ComputeBoundaryFluxes(fluidElement* elemObs, grid* primObs, geometry* geomObs, const double volElem)
{
  int world_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(PETSC_COMM_WORLD, &world_size);
  af::seq domainX1 = *primObs->domainX1;
  af::seq domainX2 = *primObs->domainX2;
  af::seq domainX3 = *primObs->domainX3;

  array MassIntegrand = primObs->vars[vars::RHO]*volElem*geomObs->g*elemObs->uCon[1];
  double MassFlowIn = 0.;
  if(primObs->iLocalStart == 0)
    {
      array MassFlowIn_af = af::sum(af::flat(MassIntegrand(3, domainX2, domainX3)),0);
      MassFlowIn += MassFlowIn_af.host<double>()[0];
    }
  double MassFlowOut = 0.;
  if(primObs->iLocalEnd == primObs->N1)
    {
      array MassFlowOut_af = af::sum(af::flat(MassIntegrand(primObs->N1Local+2, domainX2, domainX3)),0);
      MassFlowOut += MassFlowOut_af.host<double>()[0];
    }
  array Bernoulli = (1.+primObs->vars[vars::U]/primObs->vars[vars::RHO]*params::adiabaticIndex)*elemObs->uCov[0];  
  array UnboundMassIntegrand = MassIntegrand * (Bernoulli < -1.);
  array RelativisticUnboundMassIntegrand = UnboundMassIntegrand * (elemObs->gammaLorentzFactor > 2.);
  double UnboundMassFlowOut = 0.;
  if(primObs->iLocalEnd == primObs->N1)
    {
      array UnboundMassFlowOut_af = af::sum(af::flat(UnboundMassIntegrand(primObs->N1Local+2, domainX2, domainX3)),0);
      UnboundMassFlowOut += UnboundMassFlowOut_af.host<double>()[0];
    }
  double RelativisticUnboundMassFlowOut = 0.;
  if(primObs->iLocalEnd == primObs->N1)
    {
      array RelativisticUnboundMassFlowOut_af = af::sum(af::flat(RelativisticUnboundMassIntegrand(primObs->N1Local+2, domainX2, domainX3)),0);
      RelativisticUnboundMassFlowOut += RelativisticUnboundMassFlowOut_af.host<double>()[0];
    }
  /* Communicate to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  MassFlowIn += temp;
	}
    }
  else
    {
      MPI_Send(&MassFlowIn, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&MassFlowIn,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  /* Communicate to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  MassFlowOut += temp;
	}
    }
  else
    {
      MPI_Send(&MassFlowOut, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&MassFlowOut,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  /* Communicate to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  UnboundMassFlowOut += temp;
	}
    }
  else
    {
      MPI_Send(&UnboundMassFlowOut, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&UnboundMassFlowOut,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);
  /* Communicate to all processors */
  if (world_rank == 0) 
    {
      double temp; 
      for(int i=1;i<world_size;i++)
	{
	  MPI_Recv(&temp, 1, MPI_DOUBLE, i, i, PETSC_COMM_WORLD,MPI_STATUS_IGNORE);
	  RelativisticUnboundMassFlowOut += temp;
	}
    }
  else
    {
      MPI_Send(&RelativisticUnboundMassFlowOut, 1, MPI_DOUBLE, 0, world_rank, PETSC_COMM_WORLD);
    }
  MPI_Barrier(PETSC_COMM_WORLD);
  MPI_Bcast(&RelativisticUnboundMassFlowOut,1,MPI_DOUBLE,0,PETSC_COMM_WORLD);
  MPI_Barrier(PETSC_COMM_WORLD);

  PetscPrintf(PETSC_COMM_WORLD,"MdotIn = %e; MdotOut = %e; ( Unbound = %e; Relativistic = %e)\n",
	      MassFlowIn,MassFlowOut,UnboundMassFlowOut,RelativisticUnboundMassFlowOut);
}

#endif
