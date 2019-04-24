package jmetal.problems.NCTP;

/**
 * Created by stu on 16/7/8.
 */
public class NCTP13 extends NCTP{

    public NCTP13(String solutionType, Integer decisionNumber) throws ClassCastException{
        super(solutionType,2,decisionNumber,new double[]{-0.2 * Math.PI,0.2,10,1,0.5,1,2,4});
        problemName_ = "NCTP13";
    }

}
