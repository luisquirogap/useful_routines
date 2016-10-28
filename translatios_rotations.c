
int centerMass(int *Index, int nindex, int limlower, int limuper,CM *centerm)
{
  int i, count;
  double mtot=0.0, cm[3], vcm[3];

  for(i=0; i<3; i++)
    {
      cm[i]=0.0;
      vcm[i]=0.0;
      (*centerm).cm[i] = 0.0;
      (*centerm).vcm[i] = 0.0;
    }
  
  if(Index==NULL)
    {
      
      printf("Computing CM with interval [%d,%d)\n",limlower,limuper);
      
      mtot = 0.0;
      for(i=limlower; i<limuper; i++)
	{
	  cm[0] = cm[0] + particles[i].pos[0]*particles[i].mass;
	  cm[1] = cm[1] + particles[i].pos[1]*particles[i].mass;
	  cm[2] = cm[2] + particles[i].pos[2]*particles[i].mass;
	  
	  vcm[0] = vcm[0] + particles[i].vel[0]*particles[i].mass;
	  vcm[1] = vcm[1] + particles[i].vel[1]*particles[i].mass;
	  vcm[2] = vcm[2] + particles[i].vel[2]*particles[i].mass;
	  
	  mtot = mtot + particles[i].mass;
	}
      
      (*centerm).cm[0] = cm[0]/mtot;
      (*centerm).cm[1] = cm[1]/mtot;
      (*centerm).cm[2] = cm[2]/mtot;
      
      (*centerm).vcm[0] = vcm[0]/mtot;
      (*centerm).vcm[1] = vcm[1]/mtot;
      (*centerm).vcm[2] = vcm[2]/mtot;
      
    }
  else
    {
      
      printf("Computing CM with index\n");
      mtot = 0.0;
      count = 0;
      for(i=0; i<nindex; i++)
	{
	  cm[0] = cm[0] + particles[Index[i]].pos[0]*particles[Index[i]].mass;
	  cm[1] = cm[1] + particles[Index[i]].pos[1]*particles[Index[i]].mass;
	  cm[2] = cm[2] + particles[Index[i]].pos[2]*particles[Index[i]].mass;
	  
	  vcm[0] = vcm[0] + particles[Index[i]].vel[0]*particles[Index[i]].mass;
	  vcm[1] = vcm[1] + particles[Index[i]].vel[1]*particles[Index[i]].mass;
	  vcm[2] = vcm[2] + particles[Index[i]].vel[2]*particles[Index[i]].mass;
	  
	  mtot = mtot + particles[Index[i]].mass;
	  count++;
	}
      
      (*centerm).cm[0] = cm[0]/mtot;
      (*centerm).cm[1] = cm[1]/mtot;
      (*centerm).cm[2] = cm[2]/mtot;
      
      (*centerm).vcm[0] = vcm[0]/mtot;
      (*centerm).vcm[1] = vcm[1]/mtot;
      (*centerm).vcm[2] = vcm[2]/mtot;
      
      printf("with %d particles\n",count);
      
    }

  (*centerm).M = mtot;  
  
  /*
  printf("RCM : %e %e %e\n",(*centerm).cm[0],(*centerm).cm[1],(*centerm).cm[2]);
  printf("VCM : %e %e %e\n",(*centerm).vcm[0],(*centerm).vcm[1],(*centerm).vcm[2]);
  printf("LCM : %e %e %e\n",(*centerm).lcm[0],(*centerm).lcm[1],(*centerm).lcm[2]);
  printf("M TOTAL : %e\n\n",mtot);
  */
  
  return 0;
}

int compute_angmom(int *Index, int nindex, int limlower, int limuper,CM *centerm)
{
  
  int i,count;
  double lcm[3];
  
  for(i=0; i<3; i++)
    lcm[i]=0.0;
  
  if(Index==NULL)
    {
      
      printf("Computing LCM with interval (%d,%d)\n",limlower,limuper);
      
      
      for(i=limlower; i<limuper; i++)
	{
	  
	  lcm[0] = lcm[0] + (particles[i].pos[1]*particles[i].vel[2] 
			     - particles[i].pos[2]*particles[i].vel[1])*particles[i].mass;
	  lcm[1] = lcm[1] - (particles[i].pos[0]*particles[i].vel[2] 
			     - particles[i].pos[2]*particles[i].vel[0])*particles[i].mass;
	  lcm[2] = lcm[2] + (particles[i].pos[0]*particles[i].vel[1] 
			     - particles[i].pos[1]*particles[i].vel[0])*particles[i].mass;
	  
	}
      
      (*centerm).lcm[0] = lcm[0];
      (*centerm).lcm[1] = lcm[1];
      (*centerm).lcm[2] = lcm[2];
      
    }
  else
    {
      
      printf("Computing LCM with index\n");
      count = 0;
      for(i=0; i<nindex; i++)
	{
	  lcm[0] = lcm[0] + (particles[Index[i]].pos[1]*particles[Index[i]].vel[2] 
			     - particles[Index[i]].pos[2]*particles[Index[i]].vel[1])*particles[Index[i]].mass;
	  lcm[1] = lcm[1] - (particles[Index[i]].pos[0]*particles[Index[i]].vel[2] 
			     - particles[Index[i]].pos[2]*particles[Index[i]].vel[0])*particles[Index[i]].mass;
	  lcm[2] = lcm[2] + (particles[Index[i]].pos[0]*particles[Index[i]].vel[1] 
			     - particles[Index[i]].pos[1]*particles[Index[i]].vel[0])*particles[Index[i]].mass;
	  
	  count++;
	}
      
      (*centerm).lcm[0] = lcm[0];
      (*centerm).lcm[1] = lcm[1];
      (*centerm).lcm[2] = lcm[2];
      
      
      printf("with %d particles\n",count);
            
    }

  
  return 0;
}


int traslate(CM *zero, int *Index, int nindex)
{
  int i;

  for(i=0; i<nindex; i++)
    {
      
      //printf("%d %d\n",i,Index[i]);

      particles[Index[i]].pos[0] = particles[Index[i]].pos[0] - (*zero).cm[0];
      particles[Index[i]].pos[1] = particles[Index[i]].pos[1] - (*zero).cm[1]; 
      particles[Index[i]].pos[2] = particles[Index[i]].pos[2] - (*zero).cm[2]; 
      
      particles[Index[i]].vel[0] = particles[Index[i]].vel[0] - (*zero).vcm[0];  
      particles[Index[i]].vel[1] = particles[Index[i]].vel[1] - (*zero).vcm[1];  
      particles[Index[i]].vel[2] = particles[Index[i]].vel[2] - (*zero).vcm[2];  
      
    }

  return 0;

}

int totalTranslationPlus(CM *zero, int nParticles)
{
  int i;

  for( i=0; i<nParticles; i++ )
    {
      particles[i].pos[X] = particles[i].pos[X] + (*zero).cm[X];
      particles[i].pos[Y] = particles[i].pos[Y] + (*zero).cm[Y]; 
      particles[i].pos[Z] = particles[i].pos[Z] + (*zero).cm[Z]; 
      
      particles[i].vel[X] = particles[i].vel[X] + (*zero).vcm[X];  
      particles[i].vel[Y] = particles[i].vel[Y] + (*zero).vcm[Y];  
      particles[i].vel[Z] = particles[i].vel[Z] + (*zero).vcm[Z];  
    }

  return 0;

}

int totalTranslationMinus(CM *zero, int nParticles)
{
  int i;

  for( i=0; i<nParticles; i++ )
    {
      particles[i].pos[X] = particles[i].pos[X] - (*zero).cm[X];
      particles[i].pos[Y] = particles[i].pos[Y] - (*zero).cm[Y]; 
      particles[i].pos[Z] = particles[i].pos[Z] - (*zero).cm[Z]; 
      
      particles[i].vel[X] = particles[i].vel[X] - (*zero).vcm[X];  
      particles[i].vel[Y] = particles[i].vel[Y] - (*zero).vcm[Y];  
      particles[i].vel[Z] = particles[i].vel[Z] - (*zero).vcm[Z];  
    }

  return 0;

}

//int rotate(CM *vector, int *Index, int nindex)
int rotate(float *lcm, int *Index, int nindex)
{
  
  
  int m;
  double thetar=0.0, betar=0.0, r, theta, beta;
  double M1[3][3], M2[3][3];
  float vec1[3], vec2[3], lcmr[3];
  
  
  /* en el plano y,z  --> rotacion al rededor de x*/
  
  r = sqrt( lcm[1]*lcm[1] + lcm[2]*lcm[2] );
  theta = acos(fabs(lcm[1])/r);
  theta = theta*180.0/M_PI;
    
  if((lcm[1] >= 0) && (lcm[2] >= 0))
    thetar = theta + 270.0;
  
  if((lcm[1] >= 0) && (lcm[2] < 0))
    thetar = 180.0 + (90.0 - theta);
  
  if((lcm[1] < 0) && (lcm[2] < 0))
    thetar = theta + 90.0;
  
  if((lcm[1] < 0) && (lcm[2] >= 0))
    thetar = 90.0 - theta;
  
  printf("theta=%lf thetar=%lf\n",theta,thetar);
  thetar  = thetar*M_PI/180.0;
  
  //rotate about x 
  M1[0][0] = 1.0;         M1[0][1] = 0.0;               M1[0][2] = 0.0;
  M1[1][0] = 0.0;         M1[1][1] = cos(thetar);       M1[1][2] = sin(thetar);
  M1[2][0] = 0.0;         M1[2][1] = -sin(thetar);      M1[2][2] = cos(thetar);
  
  
  for(m=0; m<nindex; m++)
    {
      
      vec1[0] = particles[Index[m]].pos[0];
      vec1[1] = particles[Index[m]].pos[1];
      vec1[2] = particles[Index[m]].pos[2];      
      
      vec2[0] = particles[Index[m]].vel[0];
      vec2[1] = particles[Index[m]].vel[1];
      vec2[2] = particles[Index[m]].vel[2];
      
      particles[Index[m]].pos[0] = M1[0][0]*vec1[0] + M1[0][1]*vec1[1] + M1[0][2]*vec1[2];
      particles[Index[m]].pos[1] = M1[1][0]*vec1[0] + M1[1][1]*vec1[1] + M1[1][2]*vec1[2];
      particles[Index[m]].pos[2] = M1[2][0]*vec1[0] + M1[2][1]*vec1[1] + M1[2][2]*vec1[2];

      particles[Index[m]].vel[0] = M1[0][0]*vec2[0] + M1[0][1]*vec2[1] + M1[0][2]*vec2[2];
      particles[Index[m]].vel[1] = M1[1][0]*vec2[0] + M1[1][1]*vec2[1] + M1[1][2]*vec2[2];
      particles[Index[m]].vel[2] = M1[2][0]*vec2[0] + M1[2][1]*vec2[1] + M1[2][2]*vec2[2]; 
      
    } 
  


  //////////////////////////////////////
  /*          Segunda rotacion        */
  //////////////////////////////////////   
  
  
  lcmr[0] = M1[0][0]*lcm[0] + M1[0][1]*lcm[1] + M1[0][2]*lcm[2];
  lcmr[1] = M1[1][0]*lcm[0] + M1[1][1]*lcm[1] + M1[1][2]*lcm[2];
  lcmr[2] = M1[2][0]*lcm[0] + M1[2][1]*lcm[1] + M1[2][2]*lcm[2];

  /* en el plano x,z --> rota al rededor de y */ 
  r = sqrt( lcmr[0]*lcmr[0] + lcmr[2]*lcmr[2] );
  beta = acos(fabs(lcmr[0])/r); 
  beta = beta*180.0/M_PI;
  
  if((lcmr[0] >= 0.0) && (lcmr[2] >= 0.0))
    betar = 90 - beta;
  
  if((lcmr[0] < 0.0) && (lcmr[2] >= 0.0))
    betar = beta + 270.0;
  
  if((lcmr[0] < 0.0) && (lcmr[2] < 0.0))
    betar = (90.0 - beta) + 180.0;
  
  if((lcmr[0] >= 0.0) && (lcmr[2] < 0.0))
    betar = beta + 90.0;
    
  printf("beta=%lf betar=%lf\n",beta,betar);
  betar = betar*M_PI/180.0;
  
  //rotate about Y
  M2[0][0] = cos(betar);      M2[0][1] = 0.0;            M2[0][2] = -sin(betar);
  M2[1][0] = 0.0;             M2[1][1] = 1.0;            M2[1][2] = 0.0;
  M2[2][0] = sin(betar);      M2[2][1] = 0.0;            M2[2][2] = cos(betar);
  

  for(m=0; m<nindex; m++)
    {
      
      vec1[0] = particles[Index[m]].pos[0];
      vec1[1] = particles[Index[m]].pos[1];
      vec1[2] = particles[Index[m]].pos[2];
      
      vec2[0] = particles[Index[m]].vel[0];
      vec2[1] = particles[Index[m]].vel[1];
      vec2[2] = particles[Index[m]].vel[2]; 
      
      particles[Index[m]].pos[0] = M2[0][0]*vec1[0] + M2[0][1]*vec1[1] + M2[0][2]*vec1[2];
      particles[Index[m]].pos[1] = M2[1][0]*vec1[0] + M2[1][1]*vec1[1] + M2[1][2]*vec1[2];
      particles[Index[m]].pos[2] = M2[2][0]*vec1[0] + M2[2][1]*vec1[1] + M2[2][2]*vec1[2];

      particles[Index[m]].vel[0] = M2[0][0]*vec2[0] + M2[0][1]*vec2[1] + M2[0][2]*vec2[2];
      particles[Index[m]].vel[1] = M2[1][0]*vec2[0] + M2[1][1]*vec2[1] + M2[1][2]*vec2[2];
      particles[Index[m]].vel[2] = M2[2][0]*vec2[0] + M2[2][1]*vec2[1] + M2[2][2]*vec2[2];
      
    } 
  

  return 0;
  
}

int rotationInclinationPosition(double inclination, double position, int nParticles)
{

  int i;

  double auxVector[3], auxVector2[3];

  for( i=0; i<nParticles; i++ )
    {
      
      // Clockwise rotation of host galaxy its angle inclination around X  
      auxVector[X] = particles[i].pos[X];
      auxVector[Y] = particles[i].pos[Y]*cos(inclination) + particles[i].pos[Z]*sin(inclination);
      auxVector[Z] = -particles[i].pos[Y]*sin(inclination) + particles[i].pos[Z]*cos(inclination);
      
      auxVector2[X] = particles[i].vel[X];
      auxVector2[Y] = particles[i].vel[Y]*cos(inclination) + particles[i].vel[Z]*sin(inclination);
      auxVector2[Z] = -particles[i].vel[Y]*sin(inclination) + particles[i].vel[Z]*cos(inclination);

      // Counterclockwise rotation of host galaxy its angle position around Z  
      particles[i].pos[X] = auxVector[X]*cos(position) - auxVector[Y]*sin(position); 
      particles[i].pos[Y] = auxVector[X]*sin(position) + auxVector[Y]*cos(position); 
      particles[i].pos[Z] = auxVector[Z];

      particles[i].vel[X] = auxVector2[X]*cos(position) - auxVector2[Y]*sin(position); 
      particles[i].vel[Y] = auxVector2[X]*sin(position) + auxVector2[Y]*cos(position); 
      particles[i].vel[Z] = auxVector2[Z];
      
    }

  return 0;
}
