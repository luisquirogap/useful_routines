

///////////////////////////////////////////////////////////////////////////////////////
//ROUTINES
///////////////////////////////////////////////////////////////////////////////////////
FILE *fileOpen(char filename[],char mode[]);
int gsl_int_int_sort(int dimension,int *fvector,int *ivector);

int read_gadget1(char filename[])
{
  
  int i, j, type;
  int indexmin, indexmax;
  int nPartWithMass;
  FILE *fdata;
  int dummy;
  
  ////////////////////////////////////////////////////////////////////////
  //READ DATA FILE
  ////////////////////////////////////////////////////////////////////////
  fdata = fileOpen(filename,"r");
  
  ////////////////////////////////////////////////////////////////////////
  //READ HEADER
  ////////////////////////////////////////////////////////////////////////
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&Header,sizeof(io_header),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  ////////////////////////////////////////////////////////////////////////
  //SAVE HEADER
  ////////////////////////////////////////////////////////////////////////

  /*
  sprintf(tmp,"%s_Header.dat",filename);
  fHeader=fopen(tmp,"w");
  
  fprintf(fHeader,"\n");
  fprintf(fHeader,"Read header from %s\n",filename);
  fprintf(fHeader,"\n");
  fprintf(fHeader,"int Npart[6] : \n");
  fprintf(fHeader,"Npart[0] = %d \n",Header.Npart[0]);
  fprintf(fHeader,"Npart[1] = %d \n",Header.Npart[1]);
  fprintf(fHeader,"Npart[2] = %d \n",Header.Npart[2]);
  fprintf(fHeader,"Npart[3] = %d \n",Header.Npart[3]);
  fprintf(fHeader,"Npart[4] = %d \n",Header.Npart[4]);
  fprintf(fHeader,"Npart[5] = %d \n",Header.Npart[5]);
  fprintf(fHeader,"double mass[6] : \n");
  fprintf(fHeader,"mass[0] = %lf \n",Header.mass[0]);
  fprintf(fHeader,"mass[1] = %lf \n",Header.mass[1]);
  fprintf(fHeader,"mass[2] = %lf \n",Header.mass[2]);
  fprintf(fHeader,"mass[3] = %lf \n",Header.mass[3]);
  fprintf(fHeader,"mass[4] = %lf \n",Header.mass[4]);
  fprintf(fHeader,"mass[5] = %lf \n",Header.mass[5]);
  fprintf(fHeader,"double time = %lf \n",Header.time);
  fprintf(fHeader,"double redshift = %lf \n",Header.redshift);
  fprintf(fHeader,"int flag_sfr = %d \n",Header.flag_sfr);
  fprintf(fHeader,"int flag_feedback = %d \n",Header.flag_feedback);
  fprintf(fHeader,"int npartTotal[6] : \n");
  fprintf(fHeader,"npartTotal[0] = %u \n",Header.npartTotal[0]);
  fprintf(fHeader,"npartTotal[1] = %u \n",Header.npartTotal[1]);
  fprintf(fHeader,"npartTotal[2] = %u \n",Header.npartTotal[2]);
  fprintf(fHeader,"npartTotal[3] = %u \n",Header.npartTotal[3]);
  fprintf(fHeader,"npartTotal[4] = %u \n",Header.npartTotal[4]);
  fprintf(fHeader,"npartTotal[5] = %u \n",Header.npartTotal[5]);
  fprintf(fHeader,"int flag_cooling = %d \n",Header.flag_cooling);
  fprintf(fHeader,"int num_files = %d \n",Header.num_files); 
  fprintf(fHeader,"double BoxSize = %lf \n",Header.BoxSize);
  fprintf(fHeader,"double Omega0 = %lf \n",Header.Omega0);
  fprintf(fHeader,"double OmegaLambda = %lf \n",Header.OmegaLambda);
  fprintf(fHeader,"double HubbleParam = %lf \n",Header.HubbleParam);
  fprintf(fHeader,"int flag_stellarage = %d \n",Header.flag_stellarage);
  fprintf(fHeader,"int flag_metals = %d \n",Header.flag_metals);
  fprintf(fHeader,"int npartTotalHighWord[6] : \n");
  fprintf(fHeader,"npartTotalHighWord[0] = %u \n",Header.npartTotalHighWord[0]);
  fprintf(fHeader,"npartTotalHighWord[1] = %u \n",Header.npartTotalHighWord[1]);
  fprintf(fHeader,"npartTotalHighWord[2] = %u \n",Header.npartTotalHighWord[2]);
  fprintf(fHeader,"npartTotalHighWord[3] = %u \n",Header.npartTotalHighWord[3]);
  fprintf(fHeader,"npartTotalHighWord[4] = %u \n",Header.npartTotalHighWord[4]);
  fprintf(fHeader,"npartTotalHighWord[5] = %u \n",Header.npartTotalHighWord[5]);
  fprintf(fHeader,"int flag_entropy_instead_u = %d \n",Header.flag_entropy_instead_u);
  fprintf(fHeader,"Rest size to 256 Bytes = %lu",sizeof(Header.fill));
  */
  
  ////////////////////////////////////////////////////////////////////////
  //TOTAL NUMBER OF PARTICLES
  ////////////////////////////////////////////////////////////////////////
  N_part_total = 0;
  nPartWithMass = 0;
  printf("Reading snapshot %s with:\n",filename);
  for(i=0; i<6; i++)
    {
      N_part[i] = Header.npartTotal[i];
      N_part_total += N_part[i];
      printf("%d particles of type %d\n",N_part[i],i);
      if( Header.mass[i]>0.0 )
	nPartWithMass = nPartWithMass +0;
      else
	nPartWithMass = nPartWithMass + Header.npartTotal[i];
    }
  printf("%d particles in the snapshot\n",N_part_total);
  
  
  ////////////////////////////////////////////////////////////////////////
  //ALLOCATE AND READ
  ////////////////////////////////////////////////////////////////////////
  
  particles = (particulas *)malloc((size_t)N_part_total*sizeof(particulas));
  if(particles == NULL){
    printf("Allocation of particles failed\n");
    exit(0);
  }
  
  if( N_part[0] > 0 )
    {
      gaspro = (gas_properties *)malloc((size_t) N_part[0]*sizeof(gas_properties));
      if(gaspro == NULL){
	printf("Allocation of gaspro failed\n");
	exit(0);
      }
    }
  
  //fprintf(fHeader,"\nRead Blocks:\n");
  
  //*************************************************************************
  //POSITIONS
  //*************************************************************************
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].pos[0],sizeof(float),3,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //*************************************************************************
  //VELOCITIES
  //*************************************************************************
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].vel[0],sizeof(float),3,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //*************************************************************************
  //IDS
  //*************************************************************************
  returnRead = fread(&dummy,sizeof(dummy),1,fdata); 
#ifdef LONGIDS
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].Id,sizeof(unsigned long long),1,fdata);
#else
  for(i=0; i<N_part_total; i++)
    returnRead = fread(&particles[i].id,sizeof(unsigned int),1,fdata);
#endif
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  //*************************************************************************
  //MASSES
  //*************************************************************************

if( nPartWithMass>0  )  
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  
  N_min=N_max=0;
  
  for(j=0;j<=5;j++)
    {
      N_max=N_max+Header.npartTotal[j];
      if( Header.npartTotal[j]>0 )
	{
	  if( Header.mass[j]>0.0 )
	    {
	      for(i=N_min;i<N_max;i++)
		particles[i].mass = Header.mass[j];
	    }
	  else
	    {
	      for(i=N_min;i<N_max;i++)
		returnRead = fread(&particles[i].mass,sizeof(float),1,fdata); 
	    }
	  N_min=N_max;
	}
    }
  
  if( nPartWithMass>0  )  
    returnRead = fread(&dummy,sizeof(dummy),1,fdata);
    
  
  //int N_metal = N_part_total - Header.npartTotal[1];  
 
  //*****************************************************
  //Reading additional properties
  //*****************************************************

  if(Header.npartTotal[0]!=0)
    {
      
      //Read internal energy for particles of gas
      //===============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].U,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //===============================
      
      //Read density for particles of gas
      //===============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].rho,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //===============================
      
#ifdef COOLING
      //Read electron density for particles of gas
      //==============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].Ne,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //==============================
      
      //Read neutral hydrigen density for particles of gas
      //==============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].Nh,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //==============================
#endif
      
      //Read SPH smoothing length for particles of gas
      //===============================
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].h,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //===============================
      
      
      //Read star formation rate for particles of gas
      //=============================
#ifdef SFR
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].sfr,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
      //=============================

  //Read formation time of star for new stars
  //======================================
    }

  if(Header.npartTotal[4]!=0)
    {
      N_min = N_part[0] + N_part[1] + N_part[2] + N_part[3];
      N_max = N_min + N_part[4];
      //   printf("%d %d %d\n",Header.npartTotal[4],N_min,N_max);
#ifdef STELLARAGE
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //printf("Stellar age for %d particles\n",dummy/(int)sizeof(float));
      for(i=N_min;i<N_max;i++)
	returnRead = fread(&particles[i].stellar_age,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      //printf("Stellar age for %d particles\n",dummy/(int)sizeof(float));
#endif
    }
  
  /*
  //Read metallicity for gas and stars
  //===============================
  if((Header.flag_sfr==1) && (metallicity==1))
    {
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<N_metal;i++)
	returnRead = fread(&Z[i],sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
    }
  //===============================
  */

  //Read gravitational potential for all particles
  //===============================
#ifdef OUTPUTPOTENTIAL
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].pot,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
  //===============================
  
  //Read acceleration for all particles
  //===============================
#ifdef OUTPUTACCELERATION
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].acce,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
  //===============================
  
  //Read rate of change of entropic function for gas 
  //===============================
  if( N_part[0] >0 )
    {
#ifdef OUTPUTCHANGEOFENTROPY
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
      for(i=0;i<Header.npartTotal[0];i++)
	returnRead = fread(&gaspro[i].ecr,sizeof(float),1,fdata);
      returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
    }
  //===============================
  
  //Read timestep for all particles
  //=============================== 
#ifdef OUTPUTTIMESTEP
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  for(i=0;i<N_part_total;i++)
    returnRead = fread(&particles[i].timestep,sizeof(float),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
#endif
  //====================================================
  
  //fclose(fHeader);  //OJO, AGREGAR NUMERO Y MASA TOTAL DE PARTICULAS 
  
  //######################################
  //Sorting particles by id
  //######################################
  
  particulas  *aux;
  size_t  *p;
#ifdef LONGIDS
  unsigned long long *ID;
#else
  unsigned int *ID;
#endif
  gas_properties *auxgas;

  
for( type=0; type<6; type++)
  {
    indexmin = 0;
    for(i=0;i<type;i++)
      indexmin = indexmin + N_part[i];
    indexmax = indexmin + N_part[type];
    
    aux = (particulas *)malloc((size_t)N_part[type]*sizeof(particulas));
    if(particles == NULL){
      printf("Allocation of aux failed\n");
      exit(0);
    }
    
    p = (size_t *)malloc((size_t)N_part[type]*sizeof(size_t));
    if(p == NULL){
      printf("Allocation of p failed\n");
      exit(0);
    }
    
#ifdef LONGIDS
    ID = (unsigned long long *)malloc((size_t)N_part[type]*sizeof(unsigned long long));
    if(ID == NULL){
      printf("Allocation of ID failed\n");
      exit(0);
    }
#else
    ID = (unsigned int *)malloc((size_t)N_part[type]*sizeof(unsigned int));
    if(ID == NULL){
      printf("Allocation of ID failed\n");
      exit(0);
    }
#endif

    for( i=indexmin; i<indexmax; i++)
      {  
	ID[i-indexmin] = particles[i].id;  
	aux[i-indexmin] = particles[i];
      }

    gsl_sort_uint_index(p,ID,1,(size_t) N_part[type]); 		

    for( i=indexmin; i<indexmax; i++)
      particles[i] = aux[p[i-indexmin]];   
    
    if( (type == 0) && (N_part[0] > 0) )
      {
	auxgas = (gas_properties *)malloc((size_t) N_part[0]*sizeof(gas_properties));
	if(auxgas == NULL){
	  printf("Allocation of auxgas failed\n");
	  exit(0);
	}
	
	for( i=indexmin; i<indexmax; i++)
	  auxgas[i-indexmin] = gaspro[i];
	
	for( i=indexmin; i<indexmax; i++)
	    gaspro[i] = auxgas[p[i-indexmin]];
	
	free(auxgas);
      }

    
    free(ID);
    free(aux);
    free(p);
   
    
  }
  
  n0 = Header.npartTotal[0];
  n1 = Header.npartTotal[1];
  n2 = Header.npartTotal[2];
  n3 = Header.npartTotal[3];
  n4 = Header.npartTotal[4];
  n5 = Header.npartTotal[5];
  
  for( i=0; i<N_part_total; i++)
    {
      if(i < n0) 
	particles[i].type = 0;
      
      if( (i >= n0) && ( i < (n0+n1)) ) 
	particles[i].type = 1;
      
      if( (i >= (n0+n1)) && (i < (n0+n1+n2)) ) 
	particles[i].type = 2;
      
      if( (i >= (n0+n1+n2) ) && (i < (n0+n1+n2+n3)) ) 
	particles[i].type = 3;
      
      if( (i >= (n0+n1+n2+n3)) && (i < (n0+n1+n2+n3+n4)) ) 
	particles[i].type = 4;
      
      if((i >= (n0+n1+n2+n3+n4)) && (i < (n0+n1+n2+n3+n4+n5)) ) 
	particles[i].type = 5;
    }    
  
  //free(gaspro);  
  
  return 0;
}

int read_header_gadget(char filename[])
{
  FILE *fdata;
  int dummy;
  
  ////////////////////////////////////////////////////////////////////////
  //READ DATA FILE
  ////////////////////////////////////////////////////////////////////////
  fdata = fileOpen(filename,"r");
  
  ////////////////////////////////////////////////////////////////////////
  //READ HEADER
  ////////////////////////////////////////////////////////////////////////
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);
  returnRead = fread(&Header,sizeof(io_header),1,fdata);
  returnRead = fread(&dummy,sizeof(dummy),1,fdata);

  fclose(fdata);

  return 0;
}

int write_gadget1(char *filename, char outfile2[])
{
  
  int i,type;
  FILE *fParam,*fHeader,*fGadget;
  char outfile[500];
  int dummy;
    
  
  ////////////////////////////////////////////////////////////////////////
  // READ DATA FILE
  ////////////////////////////////////////////////////////////////////////
  //filename = argv[1];
  fParam = fileOpen(filename,"r");
  
  ////////////////////////////////////////////////////////////////////////
  // NAME OUTPUT FILE
  ////////////////////////////////////////////////////////////////////////
  //returnRead = fscanf(fParam,"%s",outfile);
  //fData = fileOpen(outfile,"r");
  //printf("\nBuilding Gadget binary with format 1 in:\n%s\n\n",outfile);

  //sprintf(outfile2,"%s.gad",outfile);
  fGadget = fopen(outfile2,"w");

  
  ////////////////////////////////////////////////////////////////////////
  // READ HEADER
  ////////////////////////////////////////////////////////////////////////
  returnRead = fscanf(fParam,"%d %d %d %d %d %d",
	 &Header.Npart[0],&Header.Npart[1],&Header.Npart[2],
	 &Header.Npart[3],&Header.Npart[4],&Header.Npart[5]);
  returnRead = fscanf(fParam,"%lf %lf %lf %lf %lf %lf",
	 &Header.mass[0],&Header.mass[1],&Header.mass[2],
	 &Header.mass[3],&Header.mass[4],&Header.mass[5]);
  returnRead = fscanf(fParam,"%lf",&Header.time);
  returnRead = fscanf(fParam,"%lf",&Header.redshift);
  returnRead = fscanf(fParam,"%d",&Header.flag_sfr);
  returnRead = fscanf(fParam,"%d",&Header.flag_feedback);
  returnRead = fscanf(fParam,"%u %u %u %u %u %u",
	 &Header.npartTotal[0],&Header.npartTotal[1],&Header.npartTotal[2],
	 &Header.npartTotal[3],&Header.npartTotal[4],&Header.npartTotal[5]);
  returnRead = fscanf(fParam,"%d",&Header.flag_cooling);
  returnRead = fscanf(fParam,"%d",&Header.num_files);
  returnRead = fscanf(fParam,"%lf",&Header.BoxSize);
  returnRead = fscanf(fParam,"%lf",&Header.Omega0);
  returnRead = fscanf(fParam,"%lf",&Header.OmegaLambda);
  returnRead = fscanf(fParam,"%lf",&Header.HubbleParam);
  returnRead = fscanf(fParam,"%d",&Header.flag_stellarage);
  returnRead = fscanf(fParam,"%d",&Header.flag_metals);
  returnRead = fscanf(fParam,"%u %u %u %u %u %u",
	 &Header.npartTotalHighWord[0],&Header.npartTotalHighWord[1],
	 &Header.npartTotalHighWord[2],&Header.npartTotalHighWord[3],
	 &Header.npartTotalHighWord[4],&Header.npartTotalHighWord[5]);
 returnRead = fscanf(fParam,"%d",&Header.flag_entropy_instead_u);

 fclose(fParam);
  
 ////////////////////////////////////////////////////////////////////////
 // SAVE HEADER
 ////////////////////////////////////////////////////////////////////////
 sprintf(outfile,"%s_Header_wrote.output",outfile2);
 fHeader=fopen(outfile,"w");
 
 fprintf(fHeader,"\n");
 fprintf(fHeader,"Wrote Header for %s\n",outfile2);
 fprintf(fHeader,"\n");
 fprintf(fHeader,"int Npart[6] : \n");
 fprintf(fHeader,"Npart[0] = %d \n",Header.Npart[0]);
 fprintf(fHeader,"Npart[1] = %d \n",Header.Npart[1]);
 fprintf(fHeader,"Npart[2] = %d \n",Header.Npart[2]);
 fprintf(fHeader,"Npart[3] = %d \n",Header.Npart[3]);
 fprintf(fHeader,"Npart[4] = %d \n",Header.Npart[4]);
 fprintf(fHeader,"Npart[5] = %d \n",Header.Npart[5]);
 fprintf(fHeader,"double mass[6] : \n");
 fprintf(fHeader,"mass[0] = %.10e \n",Header.mass[0]);
 fprintf(fHeader,"mass[1] = %.10e \n",Header.mass[1]);
 fprintf(fHeader,"mass[2] = %.10e \n",Header.mass[2]);
 fprintf(fHeader,"mass[3] = %.10e \n",Header.mass[3]);
 fprintf(fHeader,"mass[4] = %.10e \n",Header.mass[4]);
 fprintf(fHeader,"mass[5] = %.10e \n",Header.mass[5]);
 fprintf(fHeader,"double time = %lf \n",Header.time);
 fprintf(fHeader,"double redshift = %lf \n",Header.redshift);
 fprintf(fHeader,"int flag_sfr = %d \n",Header.flag_sfr);
 fprintf(fHeader,"int flag_feedback = %d \n",Header.flag_feedback);
 fprintf(fHeader,"int npartTotal[6] : \n");
 fprintf(fHeader,"npartTotal[0] = %u \n",Header.npartTotal[0]);
 fprintf(fHeader,"npartTotal[1] = %u \n",Header.npartTotal[1]);
 fprintf(fHeader,"npartTotal[2] = %u \n",Header.npartTotal[2]);
 fprintf(fHeader,"npartTotal[3] = %u \n",Header.npartTotal[3]);
 fprintf(fHeader,"npartTotal[4] = %u \n",Header.npartTotal[4]);
 fprintf(fHeader,"npartTotal[5] = %u \n",Header.npartTotal[5]);
 fprintf(fHeader,"int flag_cooling = %d \n",Header.flag_cooling);
 fprintf(fHeader,"int num_files = %d \n",Header.num_files); 
 fprintf(fHeader,"double BoxSize = %lf \n",Header.BoxSize);
 fprintf(fHeader,"double Omega0 = %lf \n",Header.Omega0);
 fprintf(fHeader,"double OmegaLambda = %lf \n",Header.OmegaLambda);
 fprintf(fHeader,"double HubbleParam = %lf \n",Header.HubbleParam);
 fprintf(fHeader,"int flag_stellarage = %d \n",Header.flag_stellarage);
 fprintf(fHeader,"int flag_metals = %d \n",Header.flag_metals);
 fprintf(fHeader,"int npartTotalHighWord[6] : \n");
 fprintf(fHeader,"npartTotalHighWord[0] = %u \n",Header.npartTotalHighWord[0]);
 fprintf(fHeader,"npartTotalHighWord[1] = %u \n",Header.npartTotalHighWord[1]);
 fprintf(fHeader,"npartTotalHighWord[2] = %u \n",Header.npartTotalHighWord[2]);
 fprintf(fHeader,"npartTotalHighWord[3] = %u \n",Header.npartTotalHighWord[3]);
 fprintf(fHeader,"npartTotalHighWord[4] = %u \n",Header.npartTotalHighWord[4]);
 fprintf(fHeader,"npartTotalHighWord[5] = %u \n",Header.npartTotalHighWord[5]);
 fprintf(fHeader,"int flag_entropy_instead_u = %d \n",Header.flag_entropy_instead_u);
 fprintf(fHeader,"Rest size to 256 Bytes = %lu",sizeof(Header.fill));

 fclose(fHeader);   	       
 
 ////////////////////////////////////////////////////////////////////////
 // TOTAL NUMBER OF PARTICLES
 ////////////////////////////////////////////////////////////////////////
 N_part_total = 0;
 for(i=0; i<6; i++)
   {
     N_part[i] = Header.npartTotal[i];
     N_part_total += N_part[i];
   }
 
 ////////////////////////////////////////////////////////////////////////
 // ALLOCATE AND READ
 ////////////////////////////////////////////////////////////////////////
 /*
 particles = (particulas *)malloc((size_t)N_part_total*sizeof(particulas));
 if(particles == NULL){
   printf("Allocation of particles failed\n");
   exit(0);
 }
 
 U = (float *)malloc((size_t)Header.npartTotal[0]*sizeof(float));
 if(particles == NULL){
   printf("Allocation of U failed\n");
   exit(0);
 }
 */

 /*
 ////////////////////////////////////////////////////////////
 // READING PARTICLES DATA
 ///////////////////////////////////////////////////////////
 
 //////////////////////////////////
 // READING GAS
 //////////////////////////////////
 if( Header.npartTotal[0] != 0)
   {
     for( i=0; i<Header.npartTotal[0]; i++)
       {
#ifdef LONGIDS
	 returnRead = fscanf(fData,"%llu %f %f %f %f %f %f %f %f",
		&particles[i].id,
		&particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
		&particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
		&particles[i].mass,
		&U[i]);
#else
	 returnRead = fscanf(fData,"%u %f %f %f %f %f %f %f %f",
		&particles[i].id,
		&particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
		&particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
		&particles[i].mass,
		&U[i]);
#endif
       }
     
     imin = Header.npartTotal[0];
   }
 
////////////////////////////////////////////////////////
// READING OTHER PARTICLES TYPE
////////////////////////////////////////////////////////
 for( i=imin; i<N_part_total; i++)
   {
#ifdef LONGIDS
     returnRead = fscanf(fData,"%llu %f %f %f %f %f %f %f",
	    &particles[i].id,
	    &particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
	    &particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
	    &particles[i].mass);
#else
     returnRead = fscanf(fData,"%u %f %f %f %f %f %f %f",
	    &particles[i].id,
	    &particles[i].pos[X],&particles[i].pos[Y],&particles[i].pos[Z],
	    &particles[i].vel[X],&particles[i].vel[Y],&particles[i].vel[Z],
	    &particles[i].mass);
#endif
   }
 
 fclose(fData);
*/
 //////////////////////////////////
 // WRITING HEADER
 /////////////////////////////////
 dummy = sizeof(Header);
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 fwrite(&Header,sizeof(Header),1,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING POSITIONS
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].pos,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING VELOCITIES
 ////////////////////////////////////////////////////////////////////////
 dummy = 3*N_part_total*sizeof(float);
 fwrite(&dummy,sizeof(dummy),1,fGadget); 
 for( i=0; i<N_part_total; i++ ) 
   fwrite(&particles[i].vel,sizeof(float),3,fGadget);
 fwrite(&dummy,sizeof(dummy),1,fGadget);

 ////////////////////////////////////////////////////////////////////////
 // WRITING IDs
 ////////////////////////////////////////////////////////////////////////
#ifdef LONGIDS
 dummy = N_part_total*sizeof(unsigned long long);	
#else
 dummy = N_part_total*sizeof(unsigned int);		
#endif

 fwrite(&dummy,sizeof(dummy),1,fGadget);
#ifdef LONGIDS
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].Id,sizeof(unsigned long long),1,fGadget);
#else
 for(i=0; i<N_part_total; i++)
   fwrite(&particles[i].id,sizeof(unsigned int),1,fGadget);
#endif
 fwrite(&dummy,sizeof(dummy),1,fGadget);
 
 ////////////////////////////////////////////////
 // WRITING MASSES
 ///////////////////////////////////////////////
 dummy = 0; 
 for( type=0; type<6; type++)
   if( Header.npartTotal[type]>0 )
     {
       if( Header.mass[type]>0.0 )
	 continue;
       else
	 dummy = dummy + Header.npartTotal[type]*sizeof(float);
     }
 if(dummy>0)
   {
     fwrite(&dummy,sizeof(dummy),1,fGadget);
     
     N_min = N_max = 0;
     for( type=0; type<6; type++)
       {
	 N_max = N_max + Header.npartTotal[type];
	 if( (Header.npartTotal[type]>0) && (Header.mass[type]>0.0) )
	   continue;
	 else
	   {
	     for( i=N_min; i<N_max; i++)
	       {
		 fwrite(&particles[i].mass,sizeof(float),1,fGadget);
	       }
	     N_min = N_max;
	   }
       } 
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 ////////////////////////////////////////////////////////////
 // WRITING INTERNAL ENERGY FOR GAS
 ///////////////////////////////////////////////////////////
 if( Header.npartTotal[0]>0 )
   {
     dummy = Header.npartTotal[0]*sizeof(float);
     fwrite(&dummy,sizeof(dummy),1,fGadget); 
     for( i=0; i<Header.npartTotal[0]; i++ )
      	 fwrite(&U[i],sizeof(float),1,fGadget);
     fwrite(&dummy,sizeof(dummy),1,fGadget);
   }

 printf("Initial conditions in %s\n\n",outfile2);

 fclose(fGadget);   	      
  
 //free(particles);
 //free(U);
 
 return 0;
}


/*
  ======================================================================
  Open a file
  ======================================================================
*/
FILE *fileOpen(char filename[],char mode[])
{
  FILE *f;

  if( !(f=fopen(filename,mode)) ){
    fprintf(stderr,"Error openning '%s' for %s\n",filename,mode);
    exit(1);
  }

  return f;
}

