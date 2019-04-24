package jmetal.operators.clone;

import jmetal.util.Configuration;
import jmetal.util.JMException;
import jmetal.util.PseudoRandom;
import jmetal.core.SolutionSet;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;

public class proportionalclone extends Clone {

	public static void printGD(String path,double[] GD){
	    try {
	      /* Open the file */
	      FileOutputStream fos   = new FileOutputStream(path)     ;
	      OutputStreamWriter osw = new OutputStreamWriter(fos)    ;
	      BufferedWriter bw      = new BufferedWriter(osw)        ;               
	      for (int i = 0; i < GD.length; i++) {  
	        bw.write(GD[i]+" ");
	        bw.newLine();        
	      }
	      
	      /* Close the file */
	      bw.close();
	    }catch (IOException e) {
	      Configuration.logger_.severe("Error acceding to the file");
	      e.printStackTrace();
	    }       
	  } // printGD
	/**
	 * @param No
	 */
	private int clonesize;

	public proportionalclone(HashMap<String, Object> parameters) {
		// this.clonesize=clonesize;
		super(parameters);
		if (parameters.get("clonesize") != null)
			//clonesize = (int) parameters.get("clonesize");
			clonesize = Integer.valueOf(parameters.get("clonesize").toString());
	} // proportional clone

	/**
	 * /** Executes the operation
	 * 
	 * @param the
	 *            parent population
	 * @return An object containing the offSprings
	 */
	public Object execute(Object parent) throws JMException {
		SolutionSet parents = (SolutionSet) parent;
		SolutionSet offSpring=new SolutionSet(clonesize);
		double min_distance=0.0;
		double max_distance=1.0;
	     double sum_distance=0.0;
	     int k=0;
	     //find max_distance and min_distance,and set the boundery solution to 2*max_distance
	     //note that parents is sorted by crowding distance, this will easier to understand alg below
		 for (k = 0; k < parents.size(); k++ ) 
		 {
		     if(parents.get(k).getCrowdingDistance()!=Double.POSITIVE_INFINITY)
		     {//k is the first not the boundery solution
				 max_distance=parents.get(k).getCrowdingDistance();
				 min_distance=parents.get(parents.size()-1).getCrowdingDistance();
			     for(int l=0;l<k;l++)
			     {
			    	 parents.get(l).setCrowdingDistance(2*parents.get(k).getCrowdingDistance());
			     }
			     break;
		    }
		  } // for
		 //this situation is that all the parents are boundery solution than set CD to 1.0
		 if(parents.get(0).getCrowdingDistance()==Double.POSITIVE_INFINITY)
		 {
			    for(int l=0;l<parents.size();l++)
			    {
			    	parents.get(l).setCrowdingDistance(1.0);
			    } 
		 }//if all the points are in extreme region.
		 //
		  for (k = 0; k < parents.size(); k++ ) 
		  {
		      sum_distance+=parents.get(k).getCrowdingDistance();
		      parents.get(k).setmaxDistance(max_distance);
		      parents.get(k).setminDistance(min_distance);
		  } // for
		  //begin to clone
		  double[] clones=new double[parents.size()];//clone number of each parent
		  for(k=0;k<parents.size();k++)
		  {
		      clones[k]= Math.ceil(clonesize*parents.get(k).getCrowdingDistance()/sum_distance);
			  if(sum_distance==0)
			  {//all individual are to one point
				  clones[k]=Math.ceil((double)clonesize/parents.size());
				  System.out.print("zeros");
				  System.out.print(clones[k]+" ");
			  }
		  /* for (int l=0;l<clones;l++)
		      {
			     Solution Newsolution=new Solution(parents.get(k));
		       //if(remain>0){
		         offSpring.add(Newsolution);
		         //remain--;
		         //}
		    }*/
		  }
		  printGD("clones",clones);
		  int remain=clonesize;
		  int i=0;
		  for(k=0;k<parents.size();k++)
		  {
			   int age=0;
			   //parents.get(k).setclone_num(k);
			    	for(int l=0;l<clones[k];l++)
			    	{		
			    		if(remain>0)
			    		{
				    	    //Solution Newsolution=new Solution(parents.get(k));
				    	    // Newsolution.age_=Newsolution.age_+age;
					        offSpring.add(parents.get(k));
					        remain--;
					        age++;
			    		}
			    	i++;
			    	}
			    	if(remain==0)
			    		{clones[k] = age;break;}
			    	//parents.get(k).age_+=clones[k];
			    	//parents.get(k).setmutationscale(PseudoRandom.randDouble());
				    if(i>600)
				    {
				    	System.out.print("zeros400");
				    }
			    }

		  /* for(int k=0;k<parents.size();k++){
		    	   Solution Newsolution=new Solution(parents.get(k));
			       offSpring.add(Newsolution);
		    	parents.get(k).age_+=1;
		  }*/
		  printGD("clones",clones);  
	  return offSpring;//*/
	} // execute
}
