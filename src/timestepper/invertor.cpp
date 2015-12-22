#include "timestepper.hpp"

void consToPrim(grid &cons, geometry &geom, grid &primGuess)
  {
    grid primGuessTrial(vars::dof, 0);
    grid consGuess(vars::dof, 0);

    grid residual(vars::dof, 0);

    grid primGuessPlusEps(vars::dof, 0);
    grid consGuessPlusEps(vars::dof, 0);

    double epsilon = 4e-8;
    array jacobianSoA = af::constant(0, 
                                     gridParams::N1Local,
                                     gridParams::N2Local,
                                     gridParams::N3Local,
                                     vars::dof * vars::dof, f64
                                    );
    array bSoA = af::constant(0., 
                              gridParams::N1Local,
                              gridParams::N2Local,
                              gridParams::N3Local,
                              vars::dof, f64
                             );
    array deltaPrimAoS = af::constant(0., 
                                      vars::dof,
                                      gridParams::N1Local,
                                      gridParams::N2Local,
                                      gridParams::N3Local, f64
                                     );

    array stepLengths = af::constant(1.,
                                     gridParams::N1Local,
                                     gridParams::N2Local,
                                     gridParams::N3Local, f64
                                    );

    array residualSoA = af::constant(0.,
                                     gridParams::N1Local,
                                     gridParams::N2Local,
                                     gridParams::N3Local,
                                     vars::dof, f64
                                    );

    array fPrime0 = af::constant(0.,
                                     gridParams::N1Local,
                                     gridParams::N2Local,
                                     gridParams::N3Local,
                                     vars::dof, f64
                                    );

    fluidElement elem(primGuess, geom, locations::CENTER);

    for (int iter=0; iter < 10; iter++)
    {
      elem.set(primGuess, geom, locations::CENTER);
      elem.computeFluxes(geom, 0, consGuess);

      for (int var=0; var < vars::dof; var++)
      {
        residual.vars[var] = consGuess.vars[var] - cons.vars[var];
        residualSoA(span, span, span, var) = residual.vars[var];

        primGuessPlusEps.vars[var]  = primGuess.vars[var];
      }

      double globalL2Norm =  af::norm(af::flat(residualSoA));
      printf("Nonliner iter = %d, error = %.15f\n", iter, globalL2Norm);

      for (int row=0; row < vars::dof; row++)
      {
        primGuessPlusEps.vars[row]  = (1. + epsilon)*primGuess.vars[row]; 

        elem.set(primGuessPlusEps, geom, locations::CENTER);
        elem.computeFluxes(geom, 0, consGuessPlusEps);

        for (int column=0; column < vars::dof; column++)
        {
          array residualPlusEps  = consGuessPlusEps.vars[column];
          array residual         = consGuess.vars[column];

          jacobianSoA(span, span, span, column + vars::dof*row)
            = (residualPlusEps - residual)/(epsilon*primGuess.vars[row]);
        }

        /* reset */
        primGuessPlusEps.vars[row]  = primGuess.vars[row]; 
      }

      /* Solve A x = b */
      for (int var=0; var < vars::dof; var++)
      {
        bSoA(span, span, span, var) = -residual.vars[var];
      }

      array jacobianAoS = af::reorder(jacobianSoA, 3, 0, 1, 2);
      array bAoS        = af::reorder(bSoA, 3, 0, 1, 2);

      for (int k=0; k<gridParams::N3Local; k++)
      {
        for (int j=0; j<gridParams::N2Local; j++)
        {
          for (int i=0; i<gridParams::N1Local; i++)
          {
            array A = af::moddims(jacobianAoS(span, i, j, k), 
                                  vars::dof, vars::dof
                                 );

            deltaPrimAoS(span, i, j, k) = af::solve(A, bAoS(span, i, j, k));
          }
        }
      }

      array deltaPrimSoA = af::reorder(deltaPrimAoS, 1, 2, 3, 0);

//      /* Quartic backtracking :
//       * funcToMinimize \equiv 0.5 * normsL2 * normsL2
//       */
//      array f0      = 0.5 * af::sum(af::pow(residualSoA, 2.), 3);
//      array fPrime0 = -2.*f0;
//      
//      array residualAoS = af::reorder(residualSoA, 3, 0, 1, 2);
//
//      for (int k=0; k<gridParams::N3Local; k++)
//      {
//        for (int j=0; j<gridParams::N2Local; j++)
//        {
//          for (int i=0; i<gridParams::N1Local; i++)
//          {
//            array jac = af::moddims(jacobianAoS(span, i, j, k),
//                                    vars::dof, vars::dof
//                                   );
//            array y = deltaPrimAoS(span, i, j, k);
//            array f = residualAoS(span, i, j, k);
//            af_print(jac);
//            af_print(f);
//            array w = af::matmul(jac, f);
//            
//            //fPrime0(i, j, k, span) = af::sum(f * w, 3);
//          }
//        }
//      }
//
//
//
//      af_print(fPrime0, 10);
//      /* Start with a full step */
//      stepLengths = 1.;
//      for (int lineSearchIter=0; lineSearchIter < 3; lineSearchIter++)
//      {
//        /* 1) First take the full step */
//        for (int var=0; var<vars::dof; var++)
//        {
//          primGuessTrial.vars[var] =  
//            primGuess.vars[var] + stepLengths*deltaPrimSoA(span, span, span, var);
//        } 
//
//        /* ...and then compute the norm */
//        elem.set(primGuessTrial, geom, locations::CENTER);
//        elem.computeFluxes(geom, 0, consGuess);
//        for (int var=0; var<vars::dof; var++)
//        {
//          residual.vars[var] = consGuess.vars[var] - cons.vars[var];
//          residualSoA(span, span, span, var) = residual.vars[var];
//        }
//        array f1 = 0.5 * af::sum(af::pow(residualSoA, 2.), 3);
//
//        /* We have 3 pieces of information:
//         * a) f(0)
//         * b) f'(0) 
//         * c) f(1) 
//         */
//      
//        double alpha    = 1e-4;
//        array condition = f1 > f0 + alpha*fPrime0;
//        array lamda     = -fPrime0/(2.*(f1 - f0 - fPrime0));
//        stepLengths     = 1. - condition + condition*lamda;
//
//        array conditionIndices = where(condition > 0);
////        if (conditionIndices.elements() == 0)
////        {
////          break;
////        }
//      }

      /* stepLengths has now been set */
      stepLengths = 1.;
      for (int var=0; var<vars::dof; var++)
      {
        primGuess.vars[var] = 
          primGuess.vars[var] + stepLengths*deltaPrimSoA(span, span, span, var);
      }
    }

  }
