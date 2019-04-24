//  Author: Qingling Zhu
//
package jmetal.metaheuristics.MOEAD;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.problems.ProblemFactory;
import jmetal.problems.CTP.*;
import jmetal.problems.LIRCMOP.*;
import jmetal.problems.NCTP.*;
import jmetal.problems.DTLZ.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.Arrays;
import java.util.HashMap;

public class MOEAD_DA_main {
  public static void printGD(String path,int[] GD){
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
   * @param args Command line arguments. The first (optional) argument specifies 
   *      the problem to solve.
   * @throws JMException 
   * @throws IOException 
   * @throws SecurityException 
   * Usage: three options
   *      - jmetal.metaheuristics.moead.MOEAD_main
   *      - jmetal.metaheuristics.moead.MOEAD_main problemName
   *      - jmetal.metaheuristics.moead.MOEAD_main problemName ParetoFrontFile
   * @throws ClassNotFoundException 
 
   */
  public static void main(String [] args) throws JMException, SecurityException, IOException, ClassNotFoundException {

    for(int fun=1;fun<=49;fun++){//15-20,42-49
		int runtimes=30;
		int[] IGDarray=new int[runtimes];
		long[] estimatedTime=new long[runtimes];

		Problem problem=null; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover ; // Crossover operator
		Operator mutation; // Mutation operator

		HashMap parameters; // Operator parameters

		QualityIndicator indicators; // Object to get quality indicators

		// Logger object and file to store log messages

		indicators = null;
		
		for(int i=0;i<runtimes;i++){
		
			if (args.length == 1) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
			} // if
			else if (args.length == 2) {
				Object[] params = { "Real" };
				problem = (new ProblemFactory()).getProblem(args[0], params);
				indicators = new QualityIndicator(problem, args[1]);
			} // if
			else { // Default problem
				if(fun==1){
			  	      problem = new LIRCMOP1("Real");
			  	}
				if(fun==2){
			  	      problem = new LIRCMOP2("Real");
			  	}
				if(fun==3){
			  	      problem = new LIRCMOP3("Real");
			  	}
				if(fun==4){
			  	      problem = new LIRCMOP4("Real");
			  	}
				if(fun==5){
			  	      problem = new LIRCMOP5("Real");
			  	}
				if(fun==6){
			  	      problem = new LIRCMOP6("Real");
			  	}
				if(fun==7){
			  	      problem = new LIRCMOP7("Real");
			  	}
				if(fun==8){
			  	      problem = new LIRCMOP8("Real");
			  	}
				if(fun==9){
			  	      problem = new LIRCMOP9("Real");
			  	}
				if(fun==10){
			  	      problem = new LIRCMOP10("Real");
			  	}
				if(fun==11){
			  	      problem = new LIRCMOP11("Real");
			  	}
				if(fun==12){
			  	      problem = new LIRCMOP12("Real");
			  	}
				if(fun==13){
			  	      problem = new LIRCMOP13("Real");
			  	}
				if(fun==14){
			  	      problem = new LIRCMOP14("Real");
			  	}
				if(fun==15){
			  	      problem = new NCTP1("Real",30);
			  	}
				if(fun==16){
			  	      problem = new NCTP2("Real",30);
			  	}
				if(fun==17){
			  	      problem = new NCTP3("Real",30);
			  	}
				if(fun==18){
			  	      problem = new NCTP4("Real",30);
			  	}
				if(fun==19){
			  	      problem = new NCTP5("Real",30);
			  	}
				if(fun==20){
			  	      problem = new NCTP6("Real",30);
			  	}
				if(fun==21){
			  	      problem = new CF1("Real",10);
			  	}
				if(fun==22){
			  	      problem = new CF2("Real",10);
			  	}
				if(fun==23){
			  	      problem = new CF3("Real",10);
			  	}
				if(fun==24){
			  	      problem = new CF4("Real",10);
			  	}
				if(fun==25){
			  	      problem = new CF5("Real",10);
			  	}
				if(fun==26){
			  	      problem = new CF6("Real",10);
			  	}
				if(fun==27){
			  	      problem = new CF7("Real",10);
			  	}
				if(fun==28){
			  	      problem = new CF8("Real",10);
			  	}
				if(fun==29){
			  	      problem = new CF9("Real",10);
			  	}
				if(fun==30){
			  	      problem = new CF10("Real",10);
			  	}
				if(fun==31){
			  	      problem = new C1_DTLZ1("Real");
			  	}
				if(fun==32){
			  	      problem = new C1_DTLZ3("Real");
			  	}
				if(fun==33){
			  	      problem = new C2_DTLZ2("Real");
			  	}
				if(fun==34){
			  	      problem = new C3_DTLZ1("Real");
			  	}
				if(fun==35){
			  	      problem = new C3_DTLZ4("Real");
			  	}
				if(fun==36){
			  	      problem = new DC1_DTLZ1("Real");
			  	}
				if(fun==37){
			  	      problem = new DC1_DTLZ3("Real");
			  	}
				if(fun==38){
			  	      problem = new DC2_DTLZ1("Real");
			  	}
				if(fun==39){
			  	      problem = new DC2_DTLZ3("Real");
			  	}
				if(fun==40){
			  	      problem = new DC3_DTLZ1("Real");
			  	}
				if(fun==41){
			  	      problem = new DC3_DTLZ3("Real");
			  	}
				if(fun==42){
			  	      problem = new CTP1("Real");
			  	}
				if(fun==43){
			  	      problem = new CTP2("Real");
			  	}
				if(fun==44){
			  	      problem = new CTP3("Real");
			  	}
				if(fun==45){
			  	      problem = new CTP4("Real");
			  	}
				if(fun==46){
			  	      problem = new CTP5("Real");
			  	}
				if(fun==47){
			  	      problem = new CTP6("Real");
			  	}
				if(fun==48){
			  	      problem = new CTP7("Real");
			  	}
				if(fun==49){
			  	      problem = new CTP8("Real");
			  	}
			}

    algorithm = new MOEAD_DA(problem);
    double F=0.5, CR=1.0;
    boolean SBX_flag=true;
    // Algorithm parameters
    if (fun<=12){
	    algorithm.setInputParameter("populationSize",100);
	    algorithm.setInputParameter("maxEvaluations",1500*100);
	    SBX_flag=false;
	    F=0.5;CR=1.0;
    }else if(fun<=14){
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",1500*105);
	    SBX_flag=true;
    }else if(fun<=20){
    	algorithm.setInputParameter("populationSize",300);
	    algorithm.setInputParameter("maxEvaluations",300*500);
	    SBX_flag=true;
    }else if(fun<=27) {
    	algorithm.setInputParameter("populationSize",100);
	    algorithm.setInputParameter("maxEvaluations",100*3000);
	    SBX_flag=false;
	    F=0.5;CR=1.0;
    }else if(fun<=29) {
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",105*3000);
	    SBX_flag=false;
	    F=0.5;CR=1.0;
    }else if(fun==30) {
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",105*3000);
	    SBX_flag=false;
	    F=0.5;CR=1.0;
    }else if(fun==33) {
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",105*250);
	    SBX_flag=true;
    }else if(fun==31||fun==36||fun==38) {
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",105*500);
	    SBX_flag=true;
    }else if(fun==34||fun==35||fun==40) {
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",105*750);
	    SBX_flag=true;
    }else if(fun==32||fun==37||fun==39) {
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",105*1000);
	    SBX_flag=true;
    }else if(fun==41) {
    	algorithm.setInputParameter("populationSize",105);
	    algorithm.setInputParameter("maxEvaluations",105*1500);
	    SBX_flag=true;
    }else if(fun<=49){//CTP
    	algorithm.setInputParameter("populationSize",100);
	    algorithm.setInputParameter("maxEvaluations",100*500);
	    SBX_flag=true;
    }else {
    	System.out.println("error!!");
    }
	
    // Directory with the files containing the weight vectors used in 
    // Q. Zhang,  W. Liu,  and H Li, The Performance of a New Version of MOEA/D 
    // on CEC09 Unconstrained MOP Test Instances Working Report CES-491, School 
    // of CS & EE, University of Essex, 02/2009.
    // http://dces.essex.ac.uk/staff/qzhang/MOEAcompetition/CEC09final/code/ZhangMOEADcode/moead0305.rar
    algorithm.setInputParameter("dataDirectory","C:/Users/qlzhu4/eclipse-workspace");
    
    

    algorithm.setInputParameter("T", 20) ;
    algorithm.setInputParameter("delta", 0.9) ;
    algorithm.setInputParameter("nr", 2) ;
    algorithm.setInputParameter("SBX", SBX_flag);
    
    // Crossover operator 
    parameters = new HashMap() ;
    if(SBX_flag) {
    	parameters.put("probability", 0.9) ;
        parameters.put("distributionIndex", 20.0) ;
        crossover = CrossoverFactory.getCrossoverOperator("SBXCrossover", parameters);
    }else {
    	parameters.put("CR", CR) ;
        parameters.put("F", F) ;
        crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters); 
    }
    
    // Mutation operator
    parameters = new HashMap();
    parameters.put("probability", 1.0/problem.getNumberOfVariables()) ;
    parameters.put("distributionIndex", 20.0) ;
    mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);                    
    
    algorithm.addOperator("crossover",crossover);
    algorithm.addOperator("mutation",mutation);
    
    // Execute the Algorithm
    long initTime = System.currentTimeMillis();
    SolutionSet population = algorithm.execute();
    estimatedTime[i] = System.currentTimeMillis() - initTime;
    
    population.printObjectivesToFile("MOEAD_DA_"+problem.getName()+"_"+problem.getNumberOfObjectives()+"_"+"T"+(i+1));
//    IGDarray[i]=indicators.getCEC_IGD(population);
    System.out.println(problem.getName()+"_"+problem.getNumberOfObjectives()+"_"+"T"+(i+1));
	}
//	printGD("C:\\Users\\qlzhu4\\Desktop\\MOEAD-2P2S\\results\\MOEAD_TrialProcedure\\MOEAD_TrialProcedure_"+problem.getName()+"_"+problem.getNumberOfObjectives()+"_FES",IGDarray);
//	Arrays.sort(IGDarray);
//	TestStatistics sta = null;
//	sta = new TestStatistics(IGDarray);
//	System.out.println(sta.getAverage()+"\t"+sta.getStandardDiviation());
//	System.out.println(IGDarray[runtimes/2]);
	}
  } //main
} // MOEAD_main
