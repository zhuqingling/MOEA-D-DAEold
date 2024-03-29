package jmetal.qualityIndicator.fastHypervolume.wfg;

import java.io.IOException;

import jmetal.core.Problem;
import jmetal.core.*;

public class wfghvCalculator {
	SolutionSet pf_ = null;
	double [][] pfMatrix_ = null;
	Point referencePoint =null;
	public wfghvCalculator(SolutionSet paretoFront,double ref[]) {
	    pf_ = paretoFront;
	    pfMatrix_ = null;
	    referencePoint = new Point(ref);
	  } // Constructor 
	public wfghvCalculator(double [][] pfMatrix) {
		pfMatrix_ = pfMatrix.clone();
		pf_ = null;
	  } // Constructor 
	public wfghvCalculator(double [][] pfMatrix,double [] p) {
		pfMatrix_ = pfMatrix.clone();
		pf_ = null;
		referencePoint = new Point(p);
	  } // Constructor 
	public static double hv2point(Point Point1, Point ref){
		  double x=ref.objectives_[0]-Point1.objectives_[0];
			for (int j=1;j<Point1.getNumberOfObjectives();j++){
				 x = x*(ref.objectives_[j]-Point1.objectives_[j]);
			}
			return x;
		}
	public double calculatewfghv() throws IOException{
		Front front = new Front() ;
//		if(referencePoint==null){
//			referencePoint = front.getReferencePoint();
//		}
	    double hv;
//	    if(pfMatrix_!=null)
//	    	front.loadFront1(pfMatrix_);
	    if(pf_!=null && pf_.size()>0)
	    	front.loadFront(pf_,-1);
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
			WFGHV wfghv = new WFGHV(referencePoint.getNumberOfObjectives(), front.nPoints_, referencePoint) ;
		    hv = wfghv.getHV(front);
		}
		return hv;
  } // CalculateHypervolume
}
