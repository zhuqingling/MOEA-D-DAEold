package jmetal.problems.NCTP_test;

/**
 * Created by stu on 16/7/8.
 */
public class NCTP6 extends NCTP{

    public NCTP6(String solutionType, Integer decisionNumber) throws ClassCastException{
        super(solutionType,2,decisionNumber,new double[]{-0.2 * Math.PI,2,10,1,6,38,-0.5});
        problemName_ = "NCTP6";
    }

}
