package jmetal.problems.NCTP_test;

/**
 * Created by stu on 16/7/8.
 */
public class NCTP12 extends NCTP{

    public NCTP12(String solutionType, Integer decisionNumber) throws ClassCastException{
        super(solutionType,2,decisionNumber,new double[]{-0.2 * Math.PI,2,10,1,6,1,1});
        problemName_ = "NCTP12";
    }

}
