//  NCTP.java
/**
 * generate constrained many-objective optimization problems
 * author: wenji li email: wenji_li@126.com
 * data: 8th July, 2016
 */

package jmetal.problems.NCTP;

import jmetal.core.Problem;
import jmetal.core.Solution;
import jmetal.encodings.solutionType.BinaryRealSolutionType;
import jmetal.encodings.solutionType.RealSolutionType;
import jmetal.util.JMException;


public abstract class NCTP extends Problem {

    private double[] constraintParams_;


    public NCTP(String solutionType, int objectiveNumber, int variableNumber, double[] constraintParams) {

        if (constraintParams.length == 8) {
            numberOfConstraints_ = 2;
        } else {
            numberOfConstraints_ = 1;
        }

        numberOfObjectives_ = objectiveNumber;
        numberOfVariables_ = variableNumber;
        constraintParams_ = constraintParams;

        upperLimit_ = new double[numberOfVariables_];
        lowerLimit_ = new double[numberOfVariables_];
        for (int var = 0; var < numberOfVariables_; var++) {
            lowerLimit_[var] = 0.0;
            upperLimit_[var] = 5.0;
        } //for

        if (solutionType.compareTo("BinaryReal") == 0)
            solutionType_ = new BinaryRealSolutionType(this);
        else if (solutionType.compareTo("Real") == 0)
            solutionType_ = new RealSolutionType(this);
        else {
            System.out.println("Error: solution type " + solutionType + " invalid");
            System.exit(-1);
        }
    } // NCTP

    public void evaluate(Solution solution) throws JMException {

        double[] x = new double[numberOfVariables_];
        double[] f = new double[numberOfObjectives_];
        double[] con = new double[numberOfConstraints_];

        for (int i = 0; i < numberOfVariables_; i++) {
            x[i] = solution.getDecisionVariables()[i].getValue();
        }

        double theta = constraintParams_[0];
        double a = constraintParams_[1];
        double b = constraintParams_[2];
        double c = constraintParams_[3];
        double d = constraintParams_[4];
        double e = constraintParams_[5];
        double z_1 = constraintParams_[6];

        double g_x = 0;
//        for (int i = 1; i < numberOfVariables_ - 1; i++) {
//            g_x += 100 * Math.pow((x[i + 1] - x[i] * x[i]), 2) + Math.pow((1 - x[i]), 2);
//        }
        for (int i = 2; i < numberOfVariables_; i++) {
            g_x += 10.0*Math.pow(x[i]-2.5,2.0);//100 * Math.pow((x[i + 1] - x[i] * x[i]), 2) + Math.pow((1 - x[i]), 2);
        }
        f[0] = x[0];
        f[1] = g_x - Math.sqrt(f[0]) + z_1;

        solution.setObjective(0,f[0]);
        solution.setObjective(1,f[1]);

        double temp;
        temp = Math.sin(b * Math.PI * Math.pow((Math.sin(theta) * (f[1] - e) + Math.cos(theta) * f[0]), c));
        con[0] = Math.cos(theta) * (f[1] - e) - Math.sin(theta) * f[0] - a * Math.pow(Math.abs(temp), d);

        if(numberOfConstraints_ == 2){
            double z_2 = constraintParams_[7];
            con[1] = z_2 - 0.73 * f[0] - f[1];
        }

        // set the constraint values
        double total = 0.0;
        int number = 0;
        for (int i = 0; i < numberOfConstraints_; i++) {
            if (con[i] < 0.0) {
                total = total + con[i];
                number++;
//                solution.setConstraint(i,con[i]);
            }
        }
        solution.setOverallConstraintViolation(-total/(double)numberOfConstraints_);
        solution.setNumberOfViolatedConstraint(number);
    }

}
