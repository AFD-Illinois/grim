#include "grim.h"

int main(int argc, char **argv)
{ 
  PetscInitialize(&argc, &argv, NULL, help);
  
  struct timeStepper ts;
  timeStepperInit(&ts);

//  while (ts.t + ts.dt < FINAL_TIME)
//  {
//    timeStep(&ts);
//  }
//  
//  /* One final step */
//  if (ts.t < FINAL_TIME)
//  {
//    ts.dt = FINAL_TIME - ts.t;
//    timeStep(&ts);
//  }
//
//  timeStepperDestroy(&ts);
//
  PetscFinalize();  
  return(0);
}
