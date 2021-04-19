package LifeCycleModel;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

public class MarketClearingInterpolation {

    public static void main(String[] args){
        double[] excess_demand = {-3.81E-04, -3.34E-04, 4.00E-05};
        double[] price = {7.23, 7.20, 7.18};

        LinearInterpolator linearInterpolator = new LinearInterpolator();
        PolynomialSplineFunction polynomialSplineFunction = linearInterpolator.interpolate(excess_demand, price);
        System.out.println(polynomialSplineFunction.value(0));
    }
}
