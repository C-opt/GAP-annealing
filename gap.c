/******************************************************************************
 A template program for developing a GAP solver. Subroutines to read instance
 data and compute the cost of a given solution are included.
 
 This program can also be used to compute the cost and check the feasibility
 of a solution given from a file. The format of a file is:
 for each job j from 1 to n in this order, the index of the agent (the value
 should be given as values from [1, m]) to which j is assigned. For example,
 if n=4 and m=3, and jobs 1, 2, 3 and 4 are assigned to agents 2, 1, 3 and 1,
 respectively, then the data in the file should be as follows:  2 1 3 1.
 
 NOTE: Index i of agents ranges from 0 to m-1, and
 index j of jobs   ranges from 0 to n-1 in the program,
	while in the solution file,
	index i of agents ranges from 1 to m, and
 index j of jobs   ranges from 1 to n in the program.
	Sorry for the confusion.
 
 If you would like to use various parameters, it might be useful to modify
 the definition of struct "Param" and mimic the way the default value of
 "timelim" is given and how its value is input from the command line.
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include "cpu_time.c"

/***** default values of parameters ******************************************/
#define	TIMELIM	300	/* the time limit for the algorithm in seconds */
#define	GIVESOL	0	/* 1: input a solution; 0: do not give a solution */

typedef struct {
  int		timelim;	/* the time limit for the algorithm in secs. */
  int		givesol;	/* give a solution (1) or not (0) */
  /* Never modify the above two lines.  */
  /* You can add more components below. */
} Param;			/* parameters */

typedef struct {
  int	n;	/* number of jobs */
  int	m;	/* number of agents */
  int	**c;	/* cost matrix c_{ij} */
  int	**a;	/* resource requirement matrix a_{ij} */
  int	*b;	/* available amount b_i of resource for each agent i */
} GAPdata;	/* data of the generalized assignment problem */

typedef struct {
  double	timebrid;	/* the time before reading the instance data */
  double	starttime;	/* the time the search started */
  double	endtime;	/* the time the search ended */
  int		*bestsol;	/* the best solution found so far */
  /* Never modify the above four lines. */
  /* You can add more components below. */
  int            cost;          /* the cost for the best solution so far */
  int            penal;         /* the penalty for the best solution so far
                                   (penalty = 0: it satisfies the constraints;
                                   penalty >= 0: it does not satisfy the constraints) */
    
} Vdata;		/* various data often necessary during the search */

/*************************** functions ***************************************/
void copy_parameters(int argc, char *arcv[], Param *param);
void read_instance(GAPdata *gapdata);
void prepare_memory(Vdata *vdata, GAPdata *gapdata);
void free_memory(Vdata *vdata, GAPdata *gapdata);
void read_sol(Vdata *vdata, GAPdata *gapdata);
void recompute_cost(Vdata *vdata, GAPdata *gapdata);
void *malloc_e(size_t size);

/***** check the feasibility and recompute the cost **************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void recompute_cost(Vdata *vdata, GAPdata *gapdata)
{
  int	i, j;		/* indices of agents and jobs */
  int	*rest_b;	/* the amount of resource available at each agent */
  int	cost, penal;	/* the cost; the penalty = the total capacity excess */
  int	temp;		/* temporary variable */
    
    
  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->bestsol[j]] -= gapdata->a[vdata->bestsol[j]][j];
    cost += gapdata->c[vdata->bestsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){penal -= temp;}
  }
  printf("recomputed cost = %d\n", cost);
  if(penal>0){
    printf("INFEASIBLE!!\n");
    printf(" resource left:");
    for(i=0; i<gapdata->m; i++){printf(" %3d", rest_b[i]);}
    printf("\n");
  }
  printf("time for the search:       %7.2f seconds\n",
	 vdata->endtime - vdata->starttime);
  printf("time to read the instance: %7.2f seconds\n",
	 vdata->starttime - vdata->timebrid);
    
  free((void *) rest_b);
}

/* It recomputes the vdata->cost and vdata->penal.
   This function is called whenever a new best solution is found.
*/
void cost_total(Vdata *vdata, GAPdata *gapdata)
{
  int	i,j;		/* indices of agents and jobs */
  int	*rest_b;	/* the amount of resource available at each agent */
  int	cost, penal;	/* the cost; the penalty = the total capacity excess */
  int	temp;		/* temporary variable */
    
    
  rest_b = (int *) malloc_e(gapdata->m * sizeof(int));
  cost = penal = 0;
  for(i=0; i<gapdata->m; i++){rest_b[i] = gapdata->b[i];}
  for(j=0; j<gapdata->n; j++){
    rest_b[vdata->bestsol[j]] -= gapdata->a[vdata->bestsol[j]][j];
    cost += gapdata->c[vdata->bestsol[j]][j];
  }
  for(i=0; i<gapdata->m; i++){
    temp = rest_b[i];
    if(temp<0){
      penal -= temp;
    }
  }
  vdata->penal = penal;
  vdata->cost = cost;
  free((void *) rest_b);
}

/* It performs shifting or swapping of the jobs on a 0.5 porcent of probabilty.
   The third argument (a int random number) must be 0 or 1.
   0 is for swapping and 1 is for shifting jobs. */
void swap_shift(Vdata *vdata, GAPdata *gapdata, int random){ //random = 0 or 1
  
  int temp, a_random_job, another_random_job, a_random_worker;
    
  if (random == 0){ //swap
    
    
    a_random_job = rand () % gapdata->n;
    while((another_random_job = rand() % gapdata->n) == a_random_job) //avoid swapping between the same job
      ;
    //simple swap
    temp = vdata->bestsol[a_random_job];
    vdata->bestsol[a_random_job] = vdata->bestsol[another_random_job];
    vdata->bestsol[another_random_job] = temp;
        
        
  }
  if (random == 1){ //shift
    a_random_job = rand() % gapdata->n;
    a_random_worker = rand() % gapdata->m;
    vdata->bestsol[a_random_job] = a_random_worker;
        
  }
  cost_total(vdata,gapdata);
}

void copy_vdata(Vdata *from, Vdata *to, GAPdata *gapdata){

  int j;
  for(j = 0; j<gapdata->n;j++)
    to->bestsol[j] = from->bestsol[j];
  to->cost = from->cost;
  to->penal = from->penal;

}

/* Annealing algorithm aplied to the GAP problem */
/* Step 1. Initialize the solution (vdata->bestsol[]) by selecting, for each job, the agent that makes the best expenditure,i.e.,the minimal expenditure;
   Step 2. Swap or shift jobs of the best solution so far (using swap_shift) on a 50% probability and recompute the cost and penalty (using cost_total);
   If:
   Case 1. the candidate solution cost < the best solution cost (temp->cost < vdata->cost) and satisfies the constraints (temp->penal == 0) then set the best solution so far as the candidate solution (vdata = temp);
   Case 2. the candidate solution cost >= the best solution cost (temp->cost >= vdata->cost) and does not satisfy the constraints (temp->penal >= 0) then set the best solution so far as the candidate solution (vdata = temp) ON A PROBABILITY OF exp(-(temp->cost - vdata->cost)/t);
   Case 3. the candidate solution cost < the best solution cost (temp->cost < vdata->cost) and does not satisfy the constraints (temp->penal == 0) then set the best solution so far as the candidate solution (vdata = temp) ON A PROBABILITY OF exp(-(temp->penal)/t);
   Case 4. the candidate solution cost >= the best solution cost (temp->cost >= vdata->cost) and does not satisfy the constraints (temp->penal >= 0) then set the best solution so far as the candidate solution (vdata = temp) ON A PROBABILITY OF exp(-(temp->penal + temp->cost - vdata->cost)/t);
   Step 3. If the time limit is not exceeded, return to 2.*/
void annealing(Vdata *vdata, GAPdata *gapdata, Param *param){
    
  int   n_jobs,m_agents;
  int	temp2;		/* temporary variable */
  int iter, i, j;
  int r; //random number
  double t, t_old, t0, t0_old; //heatness
  double  delta, delta1, delta2; //deltas
    
  Vdata *Xold;
  Vdata *optX;
  Vdata *Xnew; /* temporary variable for copying Vdata*/
    
    
    
  /*
    initialization: for each job, take the agent that makes the best expendidure
    i.e., the minimal expenditure.
  */
    
  for(n_jobs = 0; n_jobs<gapdata->n;n_jobs++){
        
    temp2 = gapdata->a[0][n_jobs];
        
    for(m_agents = 0 ; m_agents<gapdata->m ; m_agents++){
      if(gapdata->a[m_agents][n_jobs] <= temp2){
	vdata->bestsol[n_jobs] = m_agents;
	temp2 = gapdata->a[m_agents][n_jobs];
      }
    }
  }
    
  int k = 30;
  t = t_old = 0.1;
  t0 = t0_old = 0.7;
  
  cost_total(vdata,gapdata);
  Xold  = calloc(1,sizeof(Vdata));
  Xold->bestsol = (int *)  malloc_e(gapdata->n * sizeof(int));
  copy_vdata(vdata, Xold, gapdata);
  
  for(i = 1; i <= k; i++){
    while((cpu_time() - vdata->starttime) <= i*param->timelim/k){
      
      //if((cpu_time() - vdata->starttime) > param->timelim){break;}
      
      r = rand()%2;

      Xnew  = calloc(1,sizeof(Vdata));
      Xnew->bestsol = (int *)  malloc_e(gapdata->n * sizeof(int));
      copy_vdata(Xold,Xnew,gapdata);
      
      //perform swap or shift (50%)
      swap_shift(Xnew,gapdata,r);
      
      //case 1
      if((Xnew->cost - Xold->cost < 0) && (Xnew->penal == 0)){
	//cost_total(vdata,gapdata);
	copy_vdata(Xnew, Xold, gapdata);
	printf("searching-cost: %d, opt-cost: %d\n", Xold->cost, vdata->cost);
	if(Xold->cost < vdata->cost)
	  copy_vdata(Xold, vdata, gapdata);
      }
            
      //case 2
      if((Xnew->cost - Xold->cost < 0) && (Xnew->penal >= 0)){
                
	delta = 10*Xnew->penal;
	if(((double) rand()/(RAND_MAX)) <= exp(-delta/t))
	  copy_vdata(Xnew, Xold, gapdata);
	
      }
            
      //case 3
      if((Xnew->cost - Xold->cost >= 0) && (Xnew->penal == 0)){
                
	delta = (Xnew->cost -Xold->cost);
	if(((double) rand()/(RAND_MAX)) <= exp(-delta/t))
	  copy_vdata(Xnew, Xold, gapdata);
      }
            
      //case 4
      if((Xnew->cost - Xold->cost >= 0) && (Xnew->penal >= 0)){
                
	delta1 = Xnew->penal;
	delta2 = (Xnew->cost - Xold->cost);
	delta = delta1 + (100)*delta2;
	if(((double) rand()/(RAND_MAX)) <= exp(-delta/t))
	  copy_vdata(Xnew, Xold, gapdata);	
      }
            
      //remember to free the memory allocated by temp->bestsol and temp IN THIS ORDER.
      free((void *) Xnew->bestsol);
      free((void *) Xnew);
    }
    //copy_vdata(optX, Xold, gapdata);
    if(i == 1)
      copy_vdata(Xold,vdata,gapdata);
    
    if(Xold->penal == 0 && (Xold->cost >= vdata->cost) && i != 1){
      t_old = t;
      t = t + 1;
      printf("(1) t0_old:%lf\tt0:%lf\tt_old:%lf\tt:%lf\n", t0_old, t0, t_old, t);
    }
    if(Xold->penal == 0 && (Xold->cost < vdata->cost) && i != 1){
      copy_vdata(Xold,vdata,gapdata);
      t_old = t;
      t0_old = t0;
      t = (t + t0)/2;
      printf("(2) t0_old:%lf\tt0:%lf\tt_old:%lf\tt:%lf\n", t0_old, t0, t_old, t);
    }
    if(Xold->penal != 0 && (Xold->cost >= vdata->cost)){
      t = t_old;
      t0 = t0_old;
      printf("(3) t0_old:%lf\tt0:%lf\tt_old:%lf\tt:%lf\n", t0_old, t0, t_old, t);
    }
    if(Xold->penal != 0 && (Xold->cost < vdata->cost)){
      t = t_old;
      t0 = t0_old;
      printf("(4) t0_old:%lf\tt0:%lf\tt_old:%lf\tt:%lf\n", t0_old, t0, t_old, t);
    }
  }
  
  free((void *) Xold->bestsol);
  free((void *) Xold);
}


/***** read a solution from STDIN ********************************************/
void read_sol(Vdata *vdata, GAPdata *gapdata)
{
  int	j;		/* index of jobs */
  int	value_read;	/* the value read by fscanf */
  FILE	*fp=stdin;	/* set fp to the standard input */
    
  for(j=0; j<gapdata->n; j++){
    fscanf(fp, "%d", &value_read);
    /* change the range of agents from [1, m] to [0, m-1] */
    vdata->bestsol[j] = value_read - 1;
  }
}

/***** prepare memory space **************************************************/
/***** Feel free to modify this subroutine. **********************************/
void prepare_memory(Vdata *vdata, GAPdata *gapdata)
{
  int j;
  vdata->bestsol = (int *)  malloc_e(gapdata->n * sizeof(int));
  /* the next line is just to avoid confusion */
  for(j=0; j<gapdata->n; j++){vdata->bestsol[j] = 0;}
  vdata->cost = 0;
  vdata->penal = 0;
}

/***** free memory space *****************************************************/
/***** Feel free to modify this subroutine. **********************************/
void free_memory(Vdata *vdata, GAPdata *gapdata)
{
  free((void *) vdata->bestsol);
  free((void *) gapdata->c[0]);
  free((void *) gapdata->c);
  free((void *) gapdata->a[0]);
  free((void *) gapdata->a);
  free((void *) gapdata->b);
}

/***** read the instance data ************************************************/
/***** NEVER MODIFY THIS SUBROUTINE! *****************************************/
void read_instance(GAPdata *gapdata)
{
  int	i, j;		/* indices of agents and jobs */
  int	value_read;	/* the value read by fscanf */
  FILE	*fp=stdin;	/* set fp to the standard input */
    
  /* read the number of agents and jobs */
  fscanf(fp, "%d", &value_read);	/* number of agents */
  gapdata->m = value_read;
  fscanf(fp,"%d",&value_read);		/* number of jobs */
  gapdata->n = value_read;
    
  /* initialize memory */
  gapdata->c    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->c[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->c[i] = gapdata->c[i-1] + gapdata->n;}
  gapdata->a    = (int **) malloc_e(gapdata->m * sizeof(int *));
  gapdata->a[0] = (int *)  malloc_e(gapdata->m * gapdata->n * sizeof(int));
  for(i=1; i<gapdata->m; i++){gapdata->a[i] = gapdata->a[i-1] + gapdata->n;}
  gapdata->b    = (int *)  malloc_e(gapdata->m * sizeof(int));
    
  /* read the cost coefficients */
  for(i=0; i<gapdata->m; i++){
    for(j=0; j<gapdata->n; j++){
      fscanf(fp, "%d", &value_read);
      gapdata->c[i][j] = value_read;
    }
  }
    
  /* read the resource consumption */
  for(i=0; i<gapdata->m; i++){
    for(j=0; j<gapdata->n; j++){
      fscanf(fp, "%d", &value_read);
      gapdata->a[i][j] = value_read;
    }
  }
    
  /* read the resource capacity */
  for(i=0; i<gapdata->m; i++){
    fscanf(fp,"%d", &value_read);
    gapdata->b[i] = value_read;
  }
}

/***** copy and read the parameters ******************************************/
/***** Feel free to modify this subroutine. **********************************/
void copy_parameters(int argc, char *argv[], Param *param)
{
  int i;
    
  /**** copy the parameters ****/
  param->timelim = TIMELIM;
  param->givesol = GIVESOL;
  /**** read the parameters ****/
  if(argc>0 && (argc % 2)==0){
    printf("USAGE: ./gap [param_name, param_value] [name, value]...\n");
    exit(EXIT_FAILURE);}
  else{
    for(i=1; i<argc; i+=2){
      if(strcmp(argv[i],"timelim")==0) param->timelim = atoi(argv[i+1]);
      if(strcmp(argv[i],"givesol")==0) param->givesol = atoi(argv[i+1]);
    }
  }
}

/***** malloc with error check ***********************************************/
void *malloc_e( size_t size ) {
  void *s;
  if ( (s=malloc(size)) == NULL ) {
    fprintf( stderr, "malloc : Not enough memory.\n" );
    exit( EXIT_FAILURE );
  }
  return s;
}


/***** main ******************************************************************/
int main(int argc, char *argv[])
{
  Param		param;		/* parameters */
  GAPdata	gapdata;	/* GAP instance data */
  Vdata		vdata;		/* various data often needed during search */
    
  vdata.timebrid = cpu_time();
  copy_parameters(argc, argv, &param);
  read_instance(&gapdata);
  prepare_memory(&vdata, &gapdata);
  if(param.givesol==1){read_sol(&vdata, &gapdata);}
  vdata.starttime = cpu_time();
    
  /*
    Write your program here. Of course you can add your subroutines
    outside main(). At this point, the instance data is stored in "gapdata".
    gapdata.n	number of jobs n
    gapdata.m	number of agents m
    gapdata.c[i][j]	cost c_{ij}
    gapdata.a[i][j]	resource requirement a_{ij}
    gapdata.b[i]	available amount b_i of resource at agent i
    Note that i ranges from 0 to m-1, and j ranges from 0 to n-1. Note also
    that  you should write, e.g., "gapdata->c[i][j]" in your subroutines.
    Store your best solution in vdata.bestsol, then "recompute_cost" will
    compute its cost and its feasibility. The format of vdata.bestsol is:
    For each job j from 0 to n-1 in this order, the index of the agent
    (the value should be given as values from [0, m-1]) to which j is
    assigned. For example, if n=4 and m=3, and jobs 0, 1, 2 and 3 are
    assigned to agents 1, 0, 2 and 0, respectively, then vdata.bestsol
    should be as follows:
    vdata.bestsol[0] = 1
    vdata.bestsol[1] = 0
    vdata.bestsol[2] = 2
    vdata.bestsol[3] = 0.
    Note that you should write "vdata->bestsol[j]" in your subroutines.
  */
    
  annealing(&vdata, &gapdata, &param);
  vdata.endtime = cpu_time();
  recompute_cost(&vdata, &gapdata);
  free_memory(&vdata, &gapdata);
    
  return EXIT_SUCCESS;
}
