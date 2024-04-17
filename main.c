/*
 *  main.c
 *  Projet 2022-2023
 *  Elasticite lineaire plane
 *
 *  Code de calcul
 *
 *  Copyright (C) 2023 UCL-IMMC : Vincent Legat
 *  All rights reserved.
 *
 */

#include "fem.h"

int main(void) {
  femGeo *theGeometry = geoGetGeometry();
  geoMeshRead("../data/mesh.txt");
  femProblem *theProblem = femElasticityRead(theGeometry, "../data/problem.txt");
  femElasticityPrint(theProblem);
   // Rajouté par moi
  //femBandSystemPrintInfos(theProblem->bandSystem);
  double *theSoluce = femElasticitySolve(theProblem);
  int nNodes = theGeometry->theNodes->nNodes;
  printf("Number of nodes = %d\n", nNodes);
  printf("soluce = %lf\n", theSoluce[0]);
  femSolutionWrite(nNodes, 2, theSoluce, "../data/UV.txt");
  femElasticityFree(theProblem);
  geoFree();
  return 0;
}
