package jmetal.qualityIndicator.fastHypervolume.wfg;

import java.io.IOException;

import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.core.Solution;

public class wfghvIndicator {
	Problem     problem_ ; 
	String filePath_=null;
	SolutionSet PF = null;
	public wfghvIndicator(Problem problem, String paretoFrontFile) {
	    problem_ = problem ;
	    filePath_ = new String(paretoFrontFile);
	  } // Constructor 
	public wfghvIndicator(Problem problem, SolutionSet paretoFront) {
	    problem_ = problem;
	    PF = new SolutionSet(paretoFront.size());
	    for(int i=0;i<paretoFront.size();i++)
	    	PF.add(new Solution(paretoFront.get(i)));
	  } // Constructor 
	public static double hv2point(Point Point1, Point ref){
		  double x=ref.objectives_[0]-Point1.objectives_[0];
			for (int j=1;j<Point1.getNumberOfObjectives();j++){
				 x = x*(ref.objectives_[j]-Point1.objectives_[j]);
			}
			return x;
		}
	public double calculatewfghv() throws IOException{
		Problem p = problem_;
	    Front front = new Front() ;
	    
	    Point referencePoint =null;
	    String problem_name = p.getName();
	    int m=p.getNumberOfObjectives();
    	if (problem_name=="ZDT1" || problem_name=="ZDT2"){ //ZDT1-ZDT4
    		referencePoint = new Point(new double [] {2.0, 2.0});
    	}
    	if (problem_name=="ZDT6"){
    		referencePoint = new Point(new double [] {2.0, 2.0});
    	}
    	if (problem_name=="DTLZ1"){
			double [] tempp = new double [m];
			for (int j=0;j<p.getNumberOfObjectives();j++){
				tempp[j] = 2.0*0.5;
			}
			referencePoint = new Point(tempp);
    	}
    	if (problem_name=="DTLZ2" ||problem_name=="DTLZ3"||problem_name=="DTLZ4"||problem_name=="DTLZ5"||problem_name=="DTLZ6"){
			double [] tempp = new double [m];
			for (int j=0;j<m;j++){
				tempp[j] = 1.1*1.0;
			}
			referencePoint = new Point(tempp);
    	}
    	if (problem_name=="DTLZ7"){
			double [] tempp = new double [m];
			for (int j=0;j<m-1;j++){
				tempp[j] = 1.1*1.0;
			}
			tempp[m-1] = 1.1*2.0*m;
			referencePoint = new Point(tempp);
    	}
    	if (problem_name=="WFG1"||problem_name=="WFG2"||problem_name=="WFG3"||problem_name=="WFG4"||problem_name=="WFG5"||problem_name=="WFG6"||problem_name=="WFG7"||problem_name=="WFG8"||problem_name=="WFG9"){
			double [] tempp = new double [m];
			for (int j=0;j<m;j++){
				tempp[j] = 2.0*(j+1)+1;
			}
			referencePoint = new Point(tempp);
    	}

	    double hv;
	    if(PF.size()>0)
	    	front.loadFront1(PF, referencePoint);
	    else if(filePath_!=null)
	    	front.readFront(filePath_, referencePoint);
	    else
	    	return 0.0;
		
		//NORMALIZATION
//		for (int j=0;j<front.nPoints_;j++)
//			for(int k=0;k<front.getNumberOfObjectives();k++)
//				front.getPoint(j).objectives_[k] = front.getPoint(j).objectives_[k]/referencePoint.objectives_[k];
//		for (int j=0;j<front.getNumberOfObjectives();j++)
//			referencePoint.objectives_[j] = 1.0;
		
		if (front.nPoints_ == 0){
			hv = 0.0;
		}else if(front.nPoints_ == 1){
			hv = hv2point(front.getPoint(0), referencePoint);
		}else if (front.nPoints_ == 2){
			double [] mid = new double [front.getNumberOfObjectives()];
			for (int j=0;j<front.getNumberOfObjectives();j++){
				mid[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(1).objectives_[j]);
			}
			Point midp = new Point(mid);
			hv = hv2point(front.getPoint(0),referencePoint)+hv2point(front.getPoint(1),referencePoint)-hv2point(midp,referencePoint);
			
		}else if (front.nPoints_==3){
			double [] w01 = new double [front.getNumberOfObjectives()];
			double [] w02 = new double [front.getNumberOfObjectives()];
			double [] w12 = new double [front.getNumberOfObjectives()];
			double [] w012 = new double [front.getNumberOfObjectives()];
			for (int j=0;j<front.getNumberOfObjectives();j++){
				w01[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(1).objectives_[j]);
				w02[j] = Math.max(front.getPoint(0).objectives_[j], front.getPoint(2).objectives_[j]);
				w12[j] = Math.max(front.getPoint(1).objectives_[j], front.getPoint(2).objectives_[j]);
			}
			for (int j=0;j<front.getNumberOfObjectives();j++){
				w012[j] = Math.max(w02[j], front.getPoint(1).objectives_[j]);
			}
			Point p01 = new Point(w01);Point p02 = new Point(w02);Point p12 = new Point(w12);Point p012 = new Point(w012);
			hv = hv2point(front.getPoint(0),referencePoint)+hv2point(front.getPoint(1),referencePoint)+hv2point(front.getPoint(2),referencePoint)
					-hv2point(p01,referencePoint)-hv2point(p02,referencePoint)-hv2point(p12,referencePoint)+hv2point(p012,referencePoint);
		}
		else{
			WFGHV wfghv = new WFGHV(referencePoint.getNumberOfObjectives(), front.getNumberOfPoints(), referencePoint) ;
		    hv = wfghv.getHV(front);
		}
		double total_hv=1.0;
		for(int i=0;i<referencePoint.getNumberOfObjectives();i++){
			total_hv=total_hv*referencePoint.getObjectives()[i];
		}
		return hv/total_hv;
  } // CalculateHypervolume
}
