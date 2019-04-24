package jmetal.problems.NCTP;

/**
 * Created by stu on 16/7/8.
 */
public class NCTP3 extends NCTP{

    public NCTP3(String solutionType, Integer decisionNumber) throws ClassCastException{
        super(solutionType,2,decisionNumber,new double[]{-0.2 * Math.PI,2,10,1,6,1,-0.5,6});
        problemName_ = "NCTP3";
    }

}
