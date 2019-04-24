package jmetal.problems.NCTP_test;

/**
 * Created by stu on 16/7/8.
 */
public class NCTP4 extends NCTP{

    public NCTP4(String solutionType, Integer decisionNumber) throws ClassCastException{
        super(solutionType,2,decisionNumber,new double[]{-0.2 * Math.PI,0.2,10,1,0.5,38,-0.5});
        problemName_ = "NCTP4";
    }

}
