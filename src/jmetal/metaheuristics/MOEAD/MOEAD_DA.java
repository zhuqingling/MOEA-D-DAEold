//  MOEAD.java
//
//  Author:
//       Antonio J. Nebro <antonio@lcc.uma.es>
//       Juan J. Durillo <durillo@lcc.uma.es>
//
//  Copyright (c) 2011 Antonio J. Nebro, Juan J. Durillo
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
// 
//  You should have received a copy of the GNU Lesser General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.

package jmetal.metaheuristics.MOEAD;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

import jmetal.core.*;
import jmetal.util.Distance;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.util.wrapper.XReal;
import jmetal.qualityIndicator.util.MetricsUtil;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Vector;

public class MOEAD_DA extends Algorithm {

  private int populationSize_;
  private double [] utility_;
  /**
   * Stores the population
   */
  private SolutionSet population_;
  private SolutionSet population_around;
  private SolutionSet archive;
  /**
   * Z vector (ideal point)
   */
  double[] z_;
  /**
   * Lambda vectors
   */
  //Vector<Vector<Double>> lambda_ ; 
  double[][] lambda_;
  /**
   * T: neighbour size
   */
  int T_;
  int H_ = 13;
  /**
   * Neighborhood
   */
  int[][] neighborhood_;
  /**
   * delta: probability that parent solutions are selected from neighbourhood
   */
  double delta_;
  double whole_cv;
  double max_cv;
//  int awayfeasible_count=0;
  /**
   * nr: maximal number of solutions replaced by each child solution
   */
  int nr_;
  int TrialFlag=0;
  double f_r;
  double threshold=0.001;
  String functionType_;
  int evaluations_;
  /**
   * Operators
   */
  Operator crossover_;
  Operator mutation_;
  Distance distance = new Distance();
  String dataDirectory_;
  MetricsUtil utils_ = new MetricsUtil();
  /** 
   * Constructor
   * @param problem Problem to solve
   */
  public MOEAD_DA(Problem problem) {
    super (problem) ;
    functionType_ = "_TCHE1";
//    functionType_ = "_NBI";
  } // DMOEA

  public SolutionSet execute() throws JMException, ClassNotFoundException {
    int maxEvaluations;
    
    evaluations_ = 0;
    maxEvaluations = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
    populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
    dataDirectory_ = this.getInputParameter("dataDirectory").toString();
    //System.out.println("POPSIZE: "+ populationSize_) ;
    boolean SBX_flag=((Boolean)this.getInputParameter("SBX")).booleanValue();
    population_ = new SolutionSet(populationSize_);
    population_around=new SolutionSet(populationSize_);
    archive = new SolutionSet(populationSize_);
//    archive = new CrowdingArchive(populationSize_, problem_.getNumberOfObjectives());

    utility_     = new double[populationSize_];
//    T_ = ((Integer) this.getInputParameter("T")).intValue();
//    nr_ = ((Integer) this.getInputParameter("nr")).intValue();
//    delta_ = ((Double) this.getInputParameter("delta")).doubleValue();

    T_ = 20;
    delta_ = 0.9;
    nr_ = 2;
    
/*
    T_ = (int) (0.1 * populationSize_);
    delta_ = 0.9;
    nr_ = (int) (0.01 * populationSize_);
*/
    neighborhood_ = new int[populationSize_][T_];

    z_ = new double[problem_.getNumberOfObjectives()];
    //lambda_ = new Vector(problem_.getNumberOfObjectives()) ;
    lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];

    crossover_ = operators_.get("crossover"); // default: DE crossover
    mutation_ = operators_.get("mutation");  // default: polynomial mutation

    // STEP 1. Initialization
    // STEP 1.1. Compute euclidean distances between weight vectors and find T
    initUniformWeight();
    //for (int i = 0; i < 300; i++)
   // 	System.out.println(lambda_[i][0] + " " + lambda_[i][1]) ;
    
    initNeighborhood();

    // STEP 1.2. Initialize population
    initPopulation();

    // STEP 1.3. Initialize z_
    initIdealPoint(population_);
    int gen = 0;
    double epsilon=max_cv;
    SolutionSet offspring_pop=new SolutionSet(populationSize_);
    // STEP 2. Update
    do {
    	List<Integer> order = tour_selection(10);
        for (int i = 0; i < order.size(); i++) {
    	int n = order.get(i);
        int type;
        double rnd = PseudoRandom.randDouble();

        // STEP 2.1. Mating selection based on probability
        if (rnd < delta_) // if (rnd < realb)    
        {
          type = 1;   // neighborhood
        } else {
          type = 2;   // whole population
        }
        Vector<Integer> p = new Vector<Integer>();
        matingSelection(p, n, 2, type);

        // STEP 2.2. Reproduction
        Solution child=new Solution();
        
        if(SBX_flag) {
        	Solution[] parents = new Solution[2];
        	parents[0] = population_.get(n);//CURRENT
        	parents[1] = population_.get(p.get(1));
            
	        Solution[] offSpring = (Solution[]) crossover_.execute(parents);
	        child=offSpring[0];
        }else {
        	Solution[] parents = new Solution[3];
        	parents[0] = population_.get(p.get(0));
        	parents[1] = population_.get(p.get(1));
        	parents[2] = population_.get(n);//CURRENT
            
            // Apply DE crossover 
	        child = (Solution) crossover_.execute(new Object[]{population_.get(n), parents});
        }
        // Apply mutation
        mutation_.execute(child);

        // Evaluation
        problem_.evaluate(child);
        problem_.evaluateConstraints(child); 
        
        evaluations_++;

        // STEP 2.3. Repair. Not necessary
        // STEP 2.4. Update z_
        updateReference(child);
        // STEP 2.5. Update of solutions
        if(child.getOverallConstraintViolation()==0) {
        	offspring_pop.add(child);
        }
        if(TrialFlag==1) {
        	updateProblem2(child, n, type, epsilon);
        	updateProblem3(child, n, type);
        }else {
        	updateProblem2(child, n, type, epsilon);
        }
      }// for
      gen++;
      if(offspring_pop.size()>10 || (offspring_pop.size()>0 && archive.size()<1)) {
    	  ArchiveUpdate(offspring_pop);
    	  offspring_pop.clear();
      }
      if(TrialFlag==0 || TrialFlag==2) {
    	  int feasible=0;
	    	for (int i = 0; i < populationSize_; i++) {
	    		if(population_.get(i).getOverallConstraintViolation()==0)
	    			feasible=feasible+1;
	    	}
	    	f_r=(double)feasible/(double)populationSize_;
	    	f_r=Math.max(f_r, 0.01);
	    	if(f_r<0.95){
	    		epsilon=(1.0-f_r)*epsilon;//improve diversity for the population
	    	}else if(TrialFlag==0){
	    	  //TrialProcedure
	    	  TrialFlag=1;
	    	  //copy current population
	    	  epsilon=1.0e+30;
	    	}else if(TrialFlag==2) {
	    		epsilon= max_cv;//reset epsilon;*(0.5-0.4/(1.0+Math.exp(-10.0*((double)evaluations_/(double)maxEvaluations-0.5))));
	    	}else {
	    		System.out.println("error and need to stop");
	    		return null;
	    	}
      }
      if(gen%10==0 && TrialFlag!=2)
      {
    	  double whole_old=whole_cv;
          whole_cv=0.0;
          for (int i = 0; i < populationSize_; i++) {
            whole_cv=whole_cv+Math.abs(population_.get(i).getOverallConstraintViolation());
          } // for
          double ROC=(double)Math.abs(whole_cv-whole_old)/(double)whole_cv;
          if(ROC<1.0e-5 && whole_cv>max_cv) {
        	  if(TrialFlag==0) { //a local optima of constraint violation is detected
	        	//TrialProcedure
		    	  TrialFlag=1;
		    	  //copy current population
		    	  epsilon=1.0e+30;
	          }else if(TrialFlag==1){//When the population is stagnated
//	          if(whole_cv-whole_old>1.0e-5) {
//	        	  awayfeasible_count++;
//	          }else {
//	        	  awayfeasible_count=0;
//	          }
//    		  if(change<0.00001 || awayfeasible_count>0.1*(double)maxEvaluations/(double)populationSize_) {//increase f_r by change
		    	TrialFlag=2;
		    	for (int i = 0; i < populationSize_; i++) {
		    		if(Math.abs(population_.get(i).getOverallConstraintViolation())>max_cv)
		    			max_cv=Math.abs(population_.get(i).getOverallConstraintViolation());
		    	}
		    	population_.clear();
		    	for (int i = 0; i < populationSize_; i++) {
		    		population_.add(new Solution(population_around.get(i)));
		    	}
//		    	System.out.println(max_cv);
		    	epsilon=max_cv;
		    	//re-initialize the ideal point
//			    	if(archive.size()>1)
//			    		initIdealPoint(archive);
//			    	else
//			    		initIdealPoint(population_);
		    }
          }
      }
    } while (evaluations_ < maxEvaluations);
    return archive;
  }

  public void ArchiveUpdate(SolutionSet offspring) {
	  SolutionSet union=archive.union(offspring);
	  SolutionSet archive_temp=new SolutionSet(2*populationSize_);
	  archive.clear();
	  
	  //1. find nonDomianted solutions (Corner sort)
	  int j=0;
	  int m=problem_.getNumberOfObjectives();
	  do {
		  //get the best solution according to objective j
		  int best_idx=union.indexBest(new jmetal.util.comparators.ObjectiveComparator(j, false));
		  Solution solution_remove=new Solution(union.get(best_idx));
//		  System.out.println(solution_remove.getObjective(j));//check the solution to remove has the least j-th objective function value
		  archive_temp.add(solution_remove);
		  union.remove(best_idx);
		  j=Math.floorMod(j+1, m);//loop objectives
		  for (int i=0;i<union.size();i++) {
			  Solution u_i=union.get(i);
			  boolean flag=true;
			  for(int k=0;k<m;k++) {
				  if(u_i.getObjective(k)<solution_remove.getObjective(k)) {
					  flag=false;break;
				  }
			  }
			  if(flag) {//remove the solution that dominated by solution_remove
				  union.remove(i);
			  }
		  }
	  }while(union.size()>0);
	  archive_temp.Suppress();
	  //
	  if(archive_temp.size()>=populationSize_) {
		  int l=archive_temp.size();
		  if(m==2) {
			  double[] f_max=new double[m];
			  double[] f_min=new double[m];
			  for(j=0;j<m;j++) {
				  f_max[j]=-1.0e+30;
				  f_min[j]=1.0e+30;
			  }
			  for(int i=0;i<l;i++) {
				  for(j=0;j<m;j++) {
					  if(archive_temp.get(i).getObjective(j)>f_max[j]) {
						  f_max[j]=archive_temp.get(i).getObjective(j);
					  }
					  if(archive_temp.get(i).getObjective(j)<f_min[j]) {
						  f_min[j]=archive_temp.get(i).getObjective(j);
					  }
				  }
			  }
			  while(archive_temp.size()>populationSize_) {
				  //crowding distance assignment
				  l=archive_temp.size();
				  //sort archive according to the first objective value
				  archive_temp.sort(new jmetal.util.comparators.ObjectiveComparator(0, false));//decending order ==false
				  double d_min=1.0e+30;
				  int min_idx=-1;
				  archive_temp.get(0).setFitness(1.0e+30);
				  archive_temp.get(l-1).setFitness(1.0e+30);
				  for(int i=1;i<l-1;i++) {
					  double d=0.0;
					  for(j=0;j<m;j++) {
						  d=d+Math.abs(archive_temp.get(i+1).getObjective(j)-archive_temp.get(i-1).getObjective(j))/(f_max[j]-f_min[j] + Double.MIN_VALUE);
					  }
					  archive_temp.get(i).setFitness(d);
					  if(d<d_min) {
						  d_min=d;
						  min_idx=i;
					  }
				  }
				  //remove the most crowded solution
				  archive_temp.remove(min_idx);
			  }
			  for(int i=0;i<archive_temp.size();i++) {
				  archive.add(archive_temp.get(i));
			  }
		  }else {//calculate the crowding distance with niching method in BCE
			  double [][] normalizedFront ;
			  double [][] PF = archive_temp.writeObjectivesToMatrix();//N*m
			  // STEP 1. Obtain the maximum and minimum values of the Pareto front
			  double [] maximumValue = utils_.getMaximumValues(PF, m);
			  double [] minimumValue = utils_.getMinimumValues(PF, m);
			  
			  // STEP 2. Get the normalized front and true Pareto fronts
			  normalizedFront = utils_.getNormalizedFront(PF,maximumValue,minimumValue);//N*m
			  
			  double [] solutionI, solutionJ;

			  //The matrix of distances
			  double [][] d_m = new double [l][l];
			  double [][] d_t = new double [l][l];
			  //-> Calculate the distances
			  for (int i = 0; i < l; i++){
				  d_m[i][i] = 0.0;
				  d_t[i][i] = 0.0;
			      solutionI = normalizedFront[i];
			      for (j = i + 1; j < normalizedFront.length; j++){
			    	  solutionJ = normalizedFront[j];
			    	  d_m[i][j] = utils_.distance(solutionI,solutionJ);
			    	  d_m[j][i] = d_m[i][j];
			    	  d_t[i][j] = d_m[i][j];
			    	  d_t[j][i] = d_t[i][j];
			      } // for
			  } // for  
			  int k = 2;
			  //int k = (int) Math.sqrt(solutionSet_.size()/2.0);
			  double sum_d3=0.0;
			  for (int i = 0; i < d_t.length; i++) {
			    Arrays.sort(d_t[i]);
			    sum_d3 = sum_d3 + d_t[i][k]; // Calcule de D(i) distance         
			  }// for
			  double r = sum_d3/(double)d_t.length;
			  
			  ArrayList<Integer> [] neighbor = new ArrayList [l];
			  double[] R = new double[l];
			  for(int i=0;i<l;i++) {
				  neighbor[i] = new ArrayList<Integer>();
			  }
			  for(int i=0;i<l;i++){
				  R[i]=1.0;
			  	  for(j=0; j<l;j++){
			  		  if(d_m[i][j] <= r && d_m[i][j]>0.0){
			  			  R[i] = R[i] * d_m[i][j]/r;
			  			  neighbor[i].add(j);//solution j is a neighbor of solution i
			  		  }
			  	  }
			  	  archive_temp.get(i).setFitness(R[i]);
			  }
			  ArrayList<Integer> removed_idx= new ArrayList<Integer>();
			  while(l-removed_idx.size()>populationSize_) {
				  //find the most crowded solution index in R
				  double crowded=1.0e+30;
				  int crowded_idx=-1;
				  for(int i=0;i<l;i++) {
					  if(!removed_idx.contains(i)&&R[i]<crowded) {
						  crowded=R[i];
						  crowded_idx=i;
					  }
				  }
				  //update the neighbors of the removed solution
				  for(int i=0;i<l;i++) {
					  if(neighbor[i].contains(crowded_idx)) {
						  R[i]=R[i]*r/d_m[i][crowded_idx];
						  int idx=neighbor[i].indexOf(crowded_idx);
//						  System.out.println(crowded_idx+"="+neighbor[i].get(idx));
						  neighbor[i].remove(idx);
						  archive_temp.get(i).setFitness(R[i]);
					  }
				  }
				  //add the most crowded index to removed_idx
				  removed_idx.add(crowded_idx);
			  }
			  for(int i=0;i<l;i++) {
				  if(!removed_idx.contains(i)) {
					  archive.add(archive_temp.get(i));
				  }
			  }
		  }
	  }else {//archive_temp.size is less than N
		  for(int i=0;i<archive_temp.size();i++) {
			  archive.add(archive_temp.get(i));
		  }
		  //assign fitness for the archive???
	  }
  }
  public List<Integer> tour_selection(int depth)
  {

	// selection based on utility
	List<Integer> selected = new ArrayList<Integer>();
	List<Integer> candidate = new ArrayList<Integer>();

//	for(int k=0; k<problem_.getNumberOfObjectives(); k++)
//     selected.add(k);   // WARNING! HERE YOU HAVE TO USE THE WEIGHT PROVIDED BY QINGFU (NOT SORTED!!!!)
 
 
	for(int n=0; n<populationSize_; n++)
           candidate.add(n);  // set of unselected weights

	while(selected.size()<(int)(populationSize_/5.0))
	{
	    //int best_idd = (int) (rnd_uni(&rnd_uni_init)*candidate.size()), i2;
           int best_idd = (int) (PseudoRandom.randDouble()*candidate.size());
           //System.out.println(best_idd);
           int i2;
	    int best_sub = candidate.get(best_idd);
           int s2;
           for(int i=1; i<depth; i++)
           {
               i2  = (int) (PseudoRandom.randDouble()*candidate.size());
               s2  = candidate.get(i2);
               //System.out.println("Candidate: "+i2);
               if(utility_[s2]>utility_[best_sub])
               {
                   best_idd = i2;
                   best_sub = s2;
               }
           }
	    selected.add(best_sub);
	    candidate.remove(best_idd);
	}
       return selected;
   }
  public int tour_selection1(int depth)
  {

	// selection based on utility
    int selected     = -1;
	List<Integer> candidate = new ArrayList<Integer>();
 
	for(int n=0; n<archive.size(); n++)
	    candidate.add(n);  // set of unselected weights
	
	int best_idd = (int) (PseudoRandom.randDouble()*candidate.size());
	int i2;
	int best_sub = candidate.get(best_idd);
    int s2;
    for(int i=1; i<depth; i++)
    {
        i2  = (int) (PseudoRandom.randDouble()*candidate.size());
        s2  = candidate.get(i2);
        //System.out.println("Candidate: "+i2);
        if(archive.get(s2).getCrowdingDistance()>archive.get(best_sub).getCrowdingDistance())
        {
            best_idd = i2;
            best_sub = s2;
        }
    }
	selected = best_sub;
    return selected;
   }
  
  public void initUniformWeight() { // init lambda vectors
		int nw = 0;
		if (problem_.getNumberOfObjectives() == 2) {
			lambda_[nw][0] = 0.0;
			lambda_[nw][1] = 1.0;
			nw++;
			lambda_[nw][0] = 1.0;
			lambda_[nw][1] = 0.0;
			nw++;
			for (int n = 1; n < populationSize_-1; n++) {
				double a = 1.0 * n / (populationSize_ - 1);
				lambda_[nw][0] = a;
				lambda_[nw][1] = 1 - a;
				nw++;
			} // for
		} // if
		else {
			int i, j;
			lambda_[nw][0] = 0.0;
			lambda_[nw][1] = 0.0;
			lambda_[nw][2] = 1.0;
			nw++;
			lambda_[nw][0] = 0.0;
			lambda_[nw][1] = 1.0;
			lambda_[nw][2] = 0.0;
			nw++;
			lambda_[nw][0] = 1.0;
			lambda_[nw][1] = 0.0;
			lambda_[nw][2] = 0.0;
			nw++;
			for (i = 0; i <= H_; i++) {
				for (j = 0; j <= H_; j++) {
					if (i + j <= H_) {
						if((i==0 && j==0) || (i==0&&j==H_) || (i==H_&&j==0))
							continue;
						lambda_[nw][0] = (double) (1.0 * i) / H_;
						lambda_[nw][1] = (double) (1.0 * j) / H_;
						lambda_[nw][2] = (double) (1.0 * (H_ - i - j) / H_);
						nw++;
					} // if
				} // for
			} // for
		} // else

		if (nw != populationSize_) {
			System.out.println(nw + "---" + (populationSize_));
			System.out.println("ERROR: population size <> #weights");
			System.exit(0);
		}
		//Apply the WS-transformation on the generated weight vectors
		for (int i=0;i<populationSize_;i++){
			double prod = 1.0, sum = 0.0;
			for (int j=0;j<problem_.getNumberOfObjectives();j++){
				prod = prod * lambda_[i][j];
			}
			if(prod != 0.0){
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					sum = sum + 1.0/lambda_[i][j];
				}
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					lambda_[i][j] = 1.0/lambda_[i][j]/sum;
				}
			}else{
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					sum = sum + 1.0/(lambda_[i][j]+0.0000001);
				}
				for (int j=0;j<problem_.getNumberOfObjectives();j++){
					lambda_[i][j] = 1.0/(lambda_[i][j]+0.0000001)/sum;
				}
			}
		}
	} // initUniformWeight
  
  /**
   * 
   */
  public void initNeighborhood() {
    double[] x = new double[populationSize_];
    int[] idx = new int[populationSize_];

    for (int i = 0; i < populationSize_; i++) {
      // calculate the distances based on weight vectors
      for (int j = 0; j < populationSize_; j++) {
        x[j] = Utils.distVector(lambda_[i], lambda_[j]);
        //x[j] = dist_vector(population[i].namda,population[j].namda);
        idx[j] = j;
      //System.out.println("x["+j+"]: "+x[j]+ ". idx["+j+"]: "+idx[j]) ;
      } // for

      // find 'niche' nearest neighboring subproblems
      Utils.minFastSort(x, idx, populationSize_, T_);
      //minfastsort(x,idx,population.size(),niche);

        System.arraycopy(idx, 0, neighborhood_[i], 0, T_);
    } // for
  } // initNeighborhood

  /**
   * 
   */
  public void initPopulation() throws JMException, ClassNotFoundException {
	max_cv=-1.0e+30;
	whole_cv=0.0;
	int feasible=0;
    for (int i = 0; i < populationSize_; i++) {
      Solution newSolution = new Solution(problem_);
      problem_.evaluate(newSolution);
      problem_.evaluateConstraints(newSolution);
      evaluations_++;
      population_.add(newSolution) ;
      population_around.add(newSolution);
      if(newSolution.getOverallConstraintViolation()==0)
    	  feasible++;
      utility_[i]=1;
      if(Math.abs(newSolution.getOverallConstraintViolation())>max_cv){
    	  max_cv=Math.abs(newSolution.getOverallConstraintViolation());
      }
      whole_cv=whole_cv+Math.abs(newSolution.getOverallConstraintViolation());
    } // for
    f_r=(double)feasible/(double)populationSize_;
  } // initPopulation

  /**
   * 
   */
  void initIdealPoint(SolutionSet pop) throws JMException, ClassNotFoundException {
    for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
      z_[i] = 1.0e+30;
    } // for

    for (int i = 0; i < pop.size(); i++) {
      updateReference(pop.get(i));
    } // for
  } // initIdealPoint

  /**
   * 
   */
  public void matingSelection(Vector<Integer> list, int cid, int size, int type) {
    // list : the set of the indexes of selected mating parents
    // cid  : the id of current subproblem
    // size : the number of selected mating parents
    // type : 1 - neighborhood; otherwise - whole population
    int ss;
    int r;
    int p;

    ss = neighborhood_[cid].length;
    while (list.size() < size) {
      if (type == 1) {
        r = PseudoRandom.randInt(0, ss - 1);
        p = neighborhood_[cid][r];
      //p = population[cid].table[r];
      } else {
        p = PseudoRandom.randInt(0, populationSize_ - 1);
      }
      boolean flag = true;
      for (int i = 0; i < list.size(); i++) {
        if (list.get(i) == p) // p is in the list
        {
          flag = false;
          break;
        }
      }

      //if (flag) list.push_back(p);
      if (flag) {
        list.addElement(p);
      }
    }
  } // matingSelection

  /**
   * 
   * @param individual
   */
  void updateReference(Solution individual) {
    for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
      if (individual.getObjective(n) < z_[n]) {
        z_[n] = individual.getObjective(n);
//        if(z_[n]<0.0) z_[n]=0.0;
      }
    }
  } // updateReference

  /**
   * @param individual
   * @param id
   * @param type
   */
  void updateProblem2(Solution indiv, int id, int type, double epsilon) {
    // indiv: child solution
    // id:   the id of current subproblem
    // type: update solutions in - neighborhood (1) or whole population (otherwise)
    int size;
    int time;

    time = 0;

    if (type == 1) {
      size = neighborhood_[id].length;
    } else {
      size = population_.size();
    }
    int[] perm = new int[size];

    Utils.randomPermutation(perm, size);

    for (int i = 0; i < size; i++) {
      int k;
      if (type == 1) {
        k = neighborhood_[id][perm[i]];
      } else {
        k = perm[i];      // calculate the values of objective function regarding the current subproblem
      }
      double offspring_fit=fitnessFunction(indiv, lambda_[k]);
      double child_constraint=Math.abs(indiv.getOverallConstraintViolation());
      double parent_constraint=Math.abs(population_.get(k).getOverallConstraintViolation());
      if((child_constraint<=epsilon && parent_constraint<=epsilon) || child_constraint==parent_constraint) {
	      double f1, f2;
	      f1 = fitnessFunction(population_.get(k), lambda_[k]);
	      f2 = offspring_fit;
	
	      if (f2 < f1) {
	        population_.replace(k, new Solution(indiv));
	        time++;
	        double delta = (f1 - f2)/(f1+1.0e-6);
			if (delta > threshold) {
				utility_[k] = 1.0;
			}else {
				double uti = 0.95 * (1.0 + delta / threshold) * utility_[k];//0.95 + (0.05 * delta / 0.001) * utility_[k];
				utility_[k] = uti < 1.0 ? uti : 1.0;
			}
	      }
      }else {
    	  double f1=fitnessFunction(population_.get(k), lambda_[k]);
    	  double f2=fitnessFunction(indiv, lambda_[k]);
    	  f1 = f_r*f1+(1.0-f_r)*parent_constraint;
	      f2 = f_r*f2+(1.0-f_r)*child_constraint;
	
	      if (f2 < f1) {
	        population_.replace(k, new Solution(indiv));
	        time++;
	        double delta = (f1 - f2)/f1;
			if (delta > threshold) {
				utility_[k] = 1.0;
			}else {
				double uti = 0.95 * (1.0 + delta / threshold) * utility_[k];//0.95 + (0.05 * delta / 0.001) * utility_[k];
				utility_[k] = uti < 1.0 ? uti : 1.0;
			}
	      }
      }
      // the maximal number of solutions updated is not allowed to exceed 'limit'
      if (time >= nr_) {
        return;
      }
    }
  } // updateProblem
  void updateProblem3(Solution indiv, int id, int type) {
	    int size;
	    int time;

	    time = 0;
	    if (type == 1) {
	      size = neighborhood_[id].length;
	    } else {
	      size = population_around.size();
	    }
	    int[] perm = new int[size];

	    Utils.randomPermutation(perm, size);

	    for (int i = 0; i < size; i++) {
	      int k;
	      if (type == 1) {
	        k = neighborhood_[id][perm[i]];
	      } else {
	        k = perm[i];      // calculate the values of objective function regarding the current subproblem
	      }
	      double f1, f2;
	      double c1, c2;
	      f1 = fitnessFunction0(population_around.get(k), lambda_[k]);
	      f2 = fitnessFunction0(indiv, lambda_[k]);
	      c1 = Math.abs(population_around.get(k).getOverallConstraintViolation());
	      c2 = Math.abs(indiv.getOverallConstraintViolation());
	      double e=0.01;
	      if (e*f2+(1.0-e)*c2<e*f1+(1.0-e)*c1) {
	    	  population_around.replace(k, new Solution(indiv));
	    	  time++;
	      }
	      // the maximal number of solutions updated is not allowed to exceed 'limit'
	      if (time >= 1) {
	        return;
	      }
	    }
	} // updateProblem
//  void updateProblem3(Solution indiv, int id, int type) {
//	  //find the most suitable subproblem for this solution
//	  int min_idx=-1;
//	  double min_f=1.0e+30;
//	  for(int k=0;k<populationSize_;k++) {
//	      double f2 = d2(indiv, lambda_[k]);
//	      if(f2<min_f) {
//	    	  min_idx=k;
//	    	  min_f=f2;
//	      }
//	  }
//	  if(Math.abs(indiv.getOverallConstraintViolation())<Math.abs(population_around.get(min_idx).getOverallConstraintViolation())) {
//		  
//	  }
//	} // updateProblem
  void updateProblem(Solution indiv, int id, int type) {
	    // indiv: child solution
	    // id:   the id of current subproblem
	    // type: update solutions in - neighborhood (1) or whole population (otherwise)
	    int size;
	    int time;

	    time = 0;
//	    int nr=nr_;
	    if (type == 1) {
	      size = neighborhood_[id].length;
	    } else {
	      size = population_.size();
	    }
	    int[] perm = new int[size];

	    Utils.randomPermutation(perm, size);

	    for (int i = 0; i < size; i++) {
	      int k;
	      if (type == 1) {
	        k = neighborhood_[id][perm[i]];
	      } else {
	        k = perm[i];      // calculate the values of objective function regarding the current subproblem
	      }
	      double f1, f2;
	      f1 = fitnessFunction(population_.get(k), lambda_[k]);
	      f2 = fitnessFunction(indiv, lambda_[k]);

	      if (f2 < f1) {
	        population_.replace(k, new Solution(indiv));
	        time++;
	        double delta = (f1 - f2)/(f1+1.0e-6);
			if (delta > threshold) {
				utility_[k] = 1.0;
			}else {
				double uti = 0.95 * (1.0 + delta / threshold) * utility_[k];//0.95 + (0.05 * delta / 0.001) * utility_[k];
				utility_[k] = uti < 1.0 ? uti : 1.0;
			}
	      }
	      // the maximal number of solutions updated is not allowed to exceed 'limit'
	      if (time >= nr_) {
	        return;
	      }
	    }
	} // updateProblem
	double fitnessFunction(Solution individual, double[] lamda) {

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				// double diff = Math.abs(individual.getObjective(n)
				// - this.idealPoint[n]);

				double diff = Math.abs(individual.getObjective(n) - z_[n]);

				double feval;
				if (lamda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lamda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			return maxFun;
		} else if (functionType_.equals("_TCHE2")) {
            double maxFun = -1.0e+30;

            for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
                double diff = Math.abs(individual.getObjective(i) - z_[i]);

                double feval;
                if (lamda[i] == 0) {
                    feval = diff / 0.000001;
                } else {
                    feval = diff / lamda[i];
                }
                if (feval > maxFun) {
                    maxFun = feval;
                }
            } // for
            return maxFun;
        } else if (functionType_.equals("_WSUM")) {

			double sum = 0;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				sum += (lamda[n]) * individual.getObjective(n);
			}
			return sum;

		} // if
		else if (functionType_.equals("_NBI")) {
			int i;
			double d1, d2, nl;
			double theta = 5.0;
			double fin;

			d1 = d2 = nl = 0.0;
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d1 += (individual.getObjective(i) - z_[i]) * lamda[i];
				nl += Math.pow(lamda[i], 2.0);
			}
			d1 = Math.abs(d1) / Math.sqrt(nl);
			if (nl == 0.0) {
				System.out
						.println("ERROR: dived by zero(bad weihgted vector)\n");
				System.exit(0);
			}
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d2 += Math.pow((individual.getObjective(i) - z_[i]) - (d1 * lamda[i]), 2.0);
			}
			d2 = Math.sqrt(d2);
			fin = (d1 + theta * d2);
			return fin;
		}
		else {
			System.out.println("SDMOPSO.fitnessFunction: unknown type "
					+ functionType_);
			return 0;
		}
	} // fitnessEvaluation
	double fitnessFunction0(Solution individual, double[] lamda) {

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				// double diff = Math.abs(individual.getObjective(n)
				// - this.idealPoint[n]);

				double diff = Math.abs(individual.getObjective(n) - 0.5);

				double feval;
				if (lamda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lamda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			return maxFun;
		} else if (functionType_.equals("_TCHE2")) {
            double maxFun = -1.0e+30;

            for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
                double diff = Math.abs(individual.getObjective(i) - z_[i]);

                double feval;
                if (lamda[i] == 0) {
                    feval = diff / 0.000001;
                } else {
                    feval = diff / lamda[i];
                }
                if (feval > maxFun) {
                    maxFun = feval;
                }
            } // for
            return maxFun;
        } else if (functionType_.equals("_WSUM")) {

			double sum = 0;
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				sum += (lamda[n]) * individual.getObjective(n);
			}
			return sum;

		} // if
		else if (functionType_.equals("_NBI")) {
			int i;
			double d1, d2, nl;
			double theta = 5.0;
			double fin;

			d1 = d2 = nl = 0.0;
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d1 += (individual.getObjective(i) - z_[i]) * lamda[i];
				nl += Math.pow(lamda[i], 2.0);
			}
			d1 = Math.abs(d1) / Math.sqrt(nl);
			if (nl == 0.0) {
				System.out
						.println("ERROR: dived by zero(bad weihgted vector)\n");
				System.exit(0);
			}
			for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
				d2 += Math.pow((individual.getObjective(i) - z_[i]) - (d1 * lamda[i]), 2.0);
			}
			d2 = Math.sqrt(d2);
			fin = (d1 + theta * d2);
			return fin;
		}
		else {
			System.out.println("SDMOPSO.fitnessFunction: unknown type "
					+ functionType_);
			return 0;
		}
	} // fitnessEvaluation
	double d2(Solution individual, double[] lamda) {
		int i;
		double d1, d2, nl;
		double theta = 0.9;
		double fin;

		d1 = d2 = nl = 0.0;
		for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
			d1 += (individual.getObjective(i) - z_[i]) * lamda[i];
			nl += Math.pow(lamda[i], 2.0);
		}
		d1 = Math.abs(d1) / Math.sqrt(nl);
		if (nl == 0.0) {
			System.out
					.println("ERROR: dived by zero(bad weihgted vector)\n");
			System.exit(0);
		}
		for (i = 0; i < problem_.getNumberOfObjectives(); i++) {
			d2 += Math.pow((individual.getObjective(i) - z_[i]) - (d1 * lamda[i]), 2.0);
		}
		d2 = Math.sqrt(d2);
		fin = ((1.0-theta)*d1 + theta * d2);
		return fin;
	}
//  double fitnessFunction(Solution individual, double[] lambda) {
//    double fitness;
//    fitness = 0.0;
//
//    if (functionType_.equals("_TCHE1")) {
//      double maxFun = -1.0e+30;
//
//      for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
//        double diff = Math.abs(individual.getObjective(n) - z_[n]);
//
//        double feval;
//        if (lambda[n] == 0) {
//          feval = 0.0001 * diff;
//        } else {
//          feval = diff * lambda[n];
//        }
//        if (feval > maxFun) {
//          maxFun = feval;
//        }
//      } // for
//
//      fitness = maxFun;
//    } // if
//    else {
//      System.out.println("MOEAD.fitnessFunction: unknown type " + functionType_);
//      System.exit(-1);
//    }
//    return fitness;
//  } // fitnessEvaluation
  
 double distance(double [] a, double [] b) {
	    double distance = 0.0;
	    
	    for (int i = 0; i < a.length; i++) {
	      distance += Math.pow(a[i]-b[i],2.0);
	    }
	    return Math.sqrt(distance);
	  } // distance
 double distanceToClosedPoint(double [] point, double [][] front) {
	    double minDistance = distance(point,front[0]);
	    
	    
	    for (int i = 1; i < front.length; i++) {
	      double aux = distance(point,front[i]);
	      if (aux < minDistance) {
	        minDistance = aux;
	      }
	    }
	    
	    return minDistance;
	  } // distanceToClosedPoint
 double cec_IGD(double [][] front, double [][] trueParetoFront, int numberOfObjectives) {
	    double sum = 0.0;
	    int i=0;
	    for (;i<trueParetoFront.length;i++){
	      double[] aNormalizedParetoFront=trueParetoFront[i];
	      sum +=distanceToClosedPoint(aNormalizedParetoFront,front);
	    }
	    double generationalDistance = sum / trueParetoFront.length;
	    
	    return generationalDistance;
	  } // generationalDistance
  
} // MOEAD

