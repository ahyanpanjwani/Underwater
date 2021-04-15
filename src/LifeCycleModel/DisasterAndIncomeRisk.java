package LifeCycleModel;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.javatuples.Quartet;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

public class DisasterAndIncomeRisk {

    MiscFunctions miscFunctions = new MiscFunctions();

    int J = 40;                // lifespan in years

    // preferences
    double gamma = 3;          // CRRA utility term
    double psi = 250;          // bequest utility term
    double xi = 0.57;          // bequest utility term
    double alpha_safe = 0;     // homeownership utility in non-SFHA
    double alpha_risky = 5;    // homeownership utility in SFHA; alpha_risky > alpha_safe
    double d;                  // default penalty (in utile)
    double beta = 0.99;        // discount factor (based on Krishnamurthy)

    // subsidy/transfer due to disaster shock
    double tau;

    // wages and interest rate (risk less)
    double g = 0.02;           // growth rate for wages (annual)
    double y0 = 1;             // base wage
    double r = 0.0132;         // risk-free interest rate (based on Krishnamurthy)
    double[] wage_grid = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/wage_grid.csv", 6);

    // housing related costs
    double maintenance_factor = 0.03;
    double closing_factor = 0.05;
    double moving_factor = 0.01;


    // asset params
    double amin = 1e-3; double amax = 100;
    double[] a_grid = {amin, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5};
    double[] b_grid = {amin, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5};
    //double[] a_grid = {amin, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    //double[] b_grid = {amin, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};

    // disaster risk
    double[] theta_grid = {0.00, 0.22};                     // values for disaster risk
    double[][] P_theta = {{0.99, 0.01}, {0.98, 0.02}};      // transition matrix for disaster Markov chain
    double[][] P_e = {{0.95, 0.05}, {0.5, 0.5}};              // transition matrix for income Markov chain

    int N = 30;                // Term of mortgage, 30 years
    String credit_surface_filepath = "/home/ahyan/Dropbox/Underwater/DataWork/creditsurface/corelogic_credit_surface_2012.csv"; // file path for credit surface
    double[] fico_grid = {500, 600, 650, 700, 750, 800};                   // fico grid
    double[] oltv_grid = {0.50, 0.60, 0.65, 0.70, 0.75, 0.78, 0.80, 0.82, 0.86, 0.88, 0.90, 0.92, 0.94, 0.96, 0.98, 0.99, 1.00}; // oltv grid (used for interpolation as well)
    double[][] credit_surface_data = miscFunctions.credit_surface(credit_surface_filepath, fico_grid, oltv_grid); // credit surface

    // misc functions for calculating mortgage-related terms
    // given (oltv, fico), determine the mortgage rate m given the credit surface data
    double credit_surface(int ioltv, int ifico){
        double[][] cs = credit_surface_data;
        double m = cs[ifico][ioltv];
        return m;
    }

    // calculate mortgage payment x given (oltv, fico) and price of housing (risky vs safe)
    double mortgage_payment(int ioltv, int ifico, double price){
        double m = credit_surface(ioltv, ifico);
        double oltv = oltv_grid[ioltv];
        double x = (m * oltv * price) / (1 - pow(1 + m, -N));
        return x;
    }

    // calculate the contempraneous ltv
    double current_ltv(int ioltv, int ifico, int n, double price) {
        double oltv = oltv_grid[ioltv];
        double initial_balance = oltv * price;
        double m = credit_surface(ioltv, ifico);
        double ltv = (initial_balance * (pow(1 + m, N) - pow(1 + m, n)) / (pow(1 + m, N) - 1)) / price;
        return ltv;
    }

    // calculate the cohort-wide wage (a function of age only)
    double wage(int age, int ie){
        // these figures come from the BLS
        // median weekly wages by age group: https://www.bls.gov/news.release/archives/wkyeng_10182012.htm
        // median weekly wage across the population: https://www.bls.gov/opub/ted/2012/ted_20121019.htm
        // all figures are in weeks (hence I multiply by 52 to get annual numbers)
        // I divide each groups median wage by the population median to normalize and then scale to annual

        double y = 0;
        /*
        age = age + 20;
        if (age >= 20 & age < 25){y = (461d / 758);}
        else if (age >= 25 & age < 35){y = (689d / 758);}
        else if (age >= 35 & age < 45){y = (843d / 758);}
        else if (age >= 45 & age < 55){y = (871d / 758);}
        else if (age >= 55 & age < 65){y = (887d / 758);}
        else if (age >= 65){y = (761d / 758);}
        return y;
         */

        double[] age_grid = {0, 5, 15, 25, 35, 45};
        LinearInterpolator linearInterpolator = new LinearInterpolator();
        PolynomialSplineFunction function = linearInterpolator.interpolate(age_grid, wage_grid);

        if (age <= 45){
            y = function.value(age);
        }else if (age > 45){
            y = wage_grid[5];
        }

        if (ie == 0){
            y = y;
        }else if (ie == 1){
            y = 0.5 * y;
        }

        return y;
    }


    // various utility functions
    // terminal utility function (no utility from housing)
    double utility_terminal(double c, double b){
        return (pow(c, 1- gamma) / (1 - gamma)) + (psi * pow(b + xi, 1 - gamma) / (1 - gamma));
    }

    // utility from living in safe area (q = 1 always)
    double utility_safe(double c, double q){
        return (pow(c, 1- gamma) / (1 - gamma)) + (alpha_safe * q);
    }

    // utility from living in SFHA (q = 1 except for default)
    double utility_risky(double c, double q, double theta){
        return (pow(c, 1- gamma) / (1 - gamma)) + (alpha_risky * (1 - theta) * q);
    }

    // utility from living in a rental (q = 0 always)
    double utility_rental(double c, double q){
        return pow(c, 1- gamma) / (1 - gamma);
    }

    // Safe Homeowner
    // Terminal Safe Homeowner (3.5)
    double[] terminal_safe_homeowner(int age, int ia, int ioltv, int ifico, int n, int ie, double p_safe){

        double y = wage(age, ie);
        double ltv = current_ltv(ioltv, ifico, n, p_safe);
        double a = a_grid[ia];

        double VV = -1e5;
        double CC = 0;
        double BB = 0;

        for (int ib = 0; ib < b_grid.length; ib++){
            double b = b_grid[ib];
            double c = (1 + r) * a + y + (1 - ltv) * p_safe - b;
            double util = utility_terminal(c, b);

            if (c <= 0){util = -1e5;}
            if (util >= VV){
                VV = util;
                BB = b;
                CC = c;
            }
        }

        double OO = Double.NaN;
        double FF = fico_grid[ifico];
        double NN = Double.NaN;
        double action = 3.5;

        double[] output = {age, VV, CC, BB, OO, FF, NN, action};
        return output;
    }

    // Interim Safe Homeowner (3.1)
    double[] interim_safe_homeowner(int age, int ia, int ioltv, int ifico, int n, int ie,
                                        double[][][][] E_rental, double[][][][][][][] E_risky, double[][][][][][] E_safe,
                                        double q, double p_risky, double p_safe){

        double y = wage(age, ie);
        double fico = fico_grid[ifico];
        double current_x = mortgage_payment(ioltv, ifico, p_safe);
        double current_ltv = current_ltv(ioltv, ifico, n, p_safe);
        double a = a_grid[ia];

        // housing related costs
        double maintenance = maintenance_factor * p_safe;
        double closing_cost = closing_factor * p_safe;
        double moving_cost = moving_factor * (p_safe + p_risky) / 2;


        double VV = -1e5; double VV_C = -1e5; double VV_F = -1e5; double VV_M = -1e5; double VV_MR = -1e5; double VV_MS = -1e5; double VV_D = -1e5;
        double CC = 0; double CC_C = 0; double CC_F = 0; double CC_M = 0; double CC_MR = 0; double CC_MS = 0; double CC_D = 0;
        double AA = 0; double AA_C = 0; double AA_F = 0; double AA_M = 0; double AA_MR = 0; double AA_MS = 0; double AA_D = 0;
        double OO = 0; double OO_C = 0; double OO_F = 0; double OO_M = 0; double OO_MR = 0; double OO_MS = 0; double OO_D = 0;
        double FF = fico; double FF_C = fico; double FF_F = fico; double FF_M = fico; double FF_MR = fico; double FF_MS = fico; double FF_D = 500;
        double NN = 0; double NN_C = 0; double NN_F = 0; double NN_M = 0; double NN_MR = 0; double NN_MS = 0; double NN_D = 0;
        double action = 3.0; double action_current = 3.1; double action_refinance = 3.2; double action_move = 3.3; double action_move_to_rental = 3.32; double action_move_to_risky =  3.31; double action_default = 3.4;

        for (int iap = 0; iap < a_grid.length; iap++){
            double a_prime = a_grid[iap];

            // CURRENT: START
            double expected_current = E_safe[age + 1][iap][ioltv][ifico][min(n + 1, 29)][ie];
            double c_current = (1 + r) * a + y - a_prime - current_x - (maintenance_factor * p_safe);
            double util_current = utility_safe(c_current, 1) + beta * expected_current;
            if (c_current <= 0){util_current = pow(-10, 5);}
            if (util_current >= VV_C){
                VV_C = util_current;
                CC_C = c_current;
                AA_C = a_prime;
                OO_C = oltv_grid[ioltv];
                FF_C = fico_grid[ifico];
                NN_C = n + 1;
                action_current = 3.1;
            }
            // CURRENT: END

            // DEFAULT: START
            double expected_default = E_rental[age + 1][iap][0][ie];
            double c_default = (1 + r) * a + y - a_prime - q - moving_cost;
            double util_default = utility_safe(c_default, 1) - d + beta * expected_default;
            if (c_default <= 0){util_default = pow(-10, 5);}
            if (util_default >= VV_D){
                VV_D = util_default;
                CC_D = c_default;
                AA_D = a_prime;
                OO_D = Double.NaN;
                FF_D = 500;
                NN_D = Double.NaN;
                action_default = 3.4;
            }
            // DEFAULT: END

            // MOVE TO RENTAL: START
            double expected_move_to_rent = E_rental[age + 1][iap][ifico][ie];
            double c_move_to_rent = (1 + r) * a + y + (1 - current_ltv) * p_safe - a_prime - q - moving_cost;
            double util_move_to_rent = utility_rental(c_move_to_rent, 0) + beta * expected_move_to_rent;
            if (c_move_to_rent <= 0){util_move_to_rent = pow(-10, 5);}
            if (util_move_to_rent >= VV_MR){
                VV_MR = util_move_to_rent;
                CC_MR = c_move_to_rent;
                AA_MR = a_prime;
                OO_MR = Double.NaN;
                FF_MR = fico;
                NN_MR = Double.NaN;
                action_move_to_rental = 3.32;
            }
            // MOVE TO RENTAL: END



            // REFINANCE: START
            for (int ioltvp = 0; ioltvp < oltv_grid.length; ioltvp++){
                double oltvp = oltv_grid[ioltvp];

                // MOVE TO RISKY: START
                double expected_move_to_risky = E_risky[age + 1][iap][ioltvp][ifico][1][0][ie];
                double new_x_risky = mortgage_payment(ioltvp, ifico, p_risky);
                double c_move_to_risky = (1 + r) * a + y + (1 - current_ltv) * p_safe - a_prime - new_x_risky - (1 - oltvp) * p_risky - moving_cost - closing_cost;
                double util_move_to_risky = utility_risky(c_move_to_risky, 1, 0) + beta * expected_move_to_risky;
                if (c_move_to_risky <= 0){util_move_to_risky = pow(-10, 5);}
                if (util_move_to_risky >= VV_MS){
                    VV_MS = util_move_to_risky;
                    CC_MS = c_move_to_risky;
                    AA_MS = a_prime;
                    OO_MS = oltvp;
                    FF_MS = fico;
                    NN_MS = 2;
                    action_move_to_risky = 3.31;
                }
                // MOVE TO SFHA: END

                double expected_refinance = E_safe[age + 1][iap][ioltvp][ifico][1][ie];
                double new_x_refinance = mortgage_payment(ioltvp, ifico, p_safe);
                double c_refinance = (1 + r) * a + y + (oltvp - current_ltv) * p_safe - a_prime - new_x_refinance
                        - (maintenance_factor * p_safe) - (closing_factor * p_safe);
                double util_refinance = utility_safe(c_refinance, 1) + beta * expected_refinance;
                if (c_refinance <= 0){util_refinance = pow(-10, 5);}
                if (util_refinance >= VV_F){
                    VV_F = util_refinance;
                    CC_F = c_refinance;
                    AA_F = a_prime;
                    OO_F = oltvp;
                    FF_F = fico_grid[ifico];
                    NN_F = 2;
                    action_refinance = 3.2;
                }
            }

        }

        VV_M = max(VV_MS, VV_MR);
        if (VV_M == VV_MS){
            CC_M = CC_MS;
            AA_M = AA_MS;
            OO_M = OO_MS;
            FF_M = FF_MS;
            NN_M = NN_MS;
            action_move = action_move_to_risky;
        }else if (VV_M == VV_MR){
            CC_M = CC_MR;
            AA_M = AA_MR;
            OO_M = OO_MR;
            FF_M = FF_MR;
            NN_M = NN_MR;
            action_move = action_move_to_rental;
        }

        VV = max(max(VV_C, VV_F), max(VV_M, VV_D));
        //VV = max(VV_C, VV_F);
        if (VV == VV_C){
            CC = CC_C;
            AA = AA_C;
            OO = OO_C;
            FF = FF_C;
            NN = NN_C;
            action = action_current;
        }else if (VV == VV_F){
            CC = CC_F;
            AA = AA_F;
            OO = OO_F;
            FF = FF_F;
            NN = NN_F;
            action = action_refinance;
        }else if (VV == VV_M){
            CC = CC_M;
            AA = AA_M;
            OO = OO_M;
            FF = FF_M;
            NN = NN_M;
            action = action_move;
        }else if (VV == VV_D){
            CC = CC_D;
            AA = AA_D;
            OO = OO_D;
            FF = FF_D;
            NN = NN_D;
            action = action_default;
        }


        double[] output = {age, VV, CC, AA, OO, FF, NN, action};
        return output;
    }

    // Rental Market
    // Terminal Renter (2.5)
    double[] terminal_renter(int age, int ia, int ifico, int ie){
        double y = wage(age, ie);
        double a = a_grid[ia];

        double VV = -1e5;
        double CC = 0;
        double BB = 0;

        for (int ib = 0; ib < b_grid.length; ib++){
            double b = b_grid[ib];
            double c = (1 + r) * a + y - b;
            double util = utility_terminal(c, b);

            if (c <= 0){util = -1e5;}
            if (util >= VV){
                VV = util;
                CC = c;
                BB = b;
            }
        }

        double FF = fico_grid[ifico];
        double OO = Double.NaN;
        double NN = Double.NaN;
        double action = 2.5;

        double[] output = {age, VV, CC, BB, OO, FF, NN, action};
        return output;
    }

    // Interim Renter (2.1, 2.31, 2.32)
    double[] interim_renter(int age, int ia, int ifico, int ie,
                            double[][][][] E_rental, double[][][][][][][] E_risky, double[][][][][][] E_safe,
                            double q, double p_risky, double p_safe){

        double y = wage(age, ie);
        double a = a_grid[ia];

        // housing related costs
        double maintenance = maintenance_factor * p_risky;
        double closing_cost = closing_factor * p_risky;
        double moving_cost = moving_factor * p_risky;

        double VV_renter = -1e5; double VV_remain_renter = -1e5; double VV_move_to_risky = -1e5; double VV_move_to_safe = -1e5;
        double CC_renter = 0; double CC_remain_renter = 0; double CC_move_to_risky = 0; double CC_move_to_safe = 0;
        double AA_renter = 0; double AA_remain_renter = 0; double AA_move_to_risky = 0; double AA_move_to_safe = 0;
        double OO_renter = 0; double OO_remain_renter = Double.NaN; double OO_move_to_risky = 0; double OO_move_to_safe = 0;
        double NN_renter = 0; double NN_remain_renter = Double.NaN; double NN_move_to_risky = 2; double NN_move_to_safe = 0;
        double action_renter = 2.0; double action_remain_renter = 2.1; double action_move_to_risky = 2.31; double action_move_to_safe = 2.33;

        for (int iap = 0; iap < a_grid.length; iap++){
            double a_prime = a_grid[iap];

            // STAY RENTER: START
            double expected_renter = E_rental[age + 1][iap][ifico][ie];
            double c_renter = (1 + r) * a + y - q - a_prime;
            double util_renter = utility_rental(c_renter, 0) + beta * expected_renter;

            if (c_renter <= 0){util_renter = -1e5;}
            if (util_renter >= VV_remain_renter){
                VV_remain_renter = util_renter;
                CC_remain_renter = c_renter;
                AA_remain_renter = a_prime;
                OO_remain_renter = Double.NaN;
                action_renter = action_remain_renter;
            }
            // STAY RENTER: END

            for (int ioltv = 0; ioltv < oltv_grid.length; ioltv++){
                double oltv = oltv_grid[ioltv];

                // MOVE TO SAFE: START
                double new_x_safe = mortgage_payment(ioltv, ifico, p_safe);
                double expected_safe = E_safe[age + 1][iap][ioltv][ifico][1][ie];
                double c_safe = (1 + r) * a + y - a_prime - new_x_safe - (1 - oltv) * p_safe - maintenance - closing_cost - moving_cost;
                double util_safe = utility_safe(c_safe, 1) + beta * expected_safe;

                if (c_safe <= 0){util_safe = -1e5;}
                if (util_safe >= VV_move_to_safe){
                    VV_move_to_safe = util_safe;
                    CC_move_to_safe = c_safe;
                    AA_move_to_safe = a_prime;
                    OO_move_to_safe = oltv;
                    action_renter = action_move_to_safe;
                }
                // MOVE TO SAFE: END




                // MOVE TO RISKY: START
                double new_x_risky = mortgage_payment(ioltv, ifico, p_risky);
                double expected_risky = E_risky[age + 1][iap][ioltv][ifico][1][0][ie];
                double c_risky = (1 + r) * a + y - a_prime - new_x_risky - (1 - oltv) * p_risky - maintenance - closing_cost - moving_cost;
                double util_risky = utility_risky(c_risky, 1, 0) + beta * expected_risky;

                if (c_risky <= 0){util_risky = -1e5;}
                if (util_risky >= VV_move_to_risky){
                    VV_move_to_risky = util_risky;
                    CC_move_to_risky = c_risky;
                    AA_move_to_risky = a_prime;
                    OO_move_to_risky = oltv;
                    action_renter = action_move_to_risky;
                }
                // MOVE TO RISKY: END
            }
        }

        VV_renter = max(VV_remain_renter, max(VV_move_to_risky, VV_move_to_safe));
        if (VV_renter == VV_remain_renter){
            CC_renter = CC_remain_renter;
            AA_renter = AA_remain_renter;
            OO_renter = OO_remain_renter;
            NN_renter = NN_remain_renter;
            action_renter = action_remain_renter;
        }else if (VV_renter == VV_move_to_risky){
            CC_renter = CC_move_to_risky;
            AA_renter = AA_move_to_risky;
            OO_renter = OO_move_to_risky;
            NN_renter = NN_move_to_risky;
            action_renter = action_move_to_risky;
        }else if (VV_renter == VV_move_to_safe){
            CC_renter = CC_move_to_safe;
            AA_renter = AA_move_to_safe;
            OO_renter = OO_move_to_safe;
            NN_renter = NN_move_to_safe;
            action_renter = action_move_to_safe;
        }

        double FF_renter = fico_grid[ifico];

        double[] output = {age, VV_renter, CC_renter, AA_renter, OO_renter, FF_renter, NN_renter, action_renter};
        return output;
    }

    // At-Risk Homeowner
    // Terminal At-Risk Homeowner (1.5)
    double[] terminal_risky_homeowner(int age, int ia, int ioltv, int ifico, int n, int itheta, int ie, double p_risky){

        double y = wage(age, ie);
        double a = a_grid[ia];
        double ltv = current_ltv(ioltv, ifico, n, p_risky);
        double theta = theta_grid[itheta];

        double VV = -1e5;
        double CC = 0;
        double BB = 0;
        double OO = Double.NaN;
        double FF = fico_grid[ifico];
        double NN = Double.NaN;
        double action = 1.5;

        for (int ib = 0; ib < b_grid.length; ib++){
            double b = b_grid[ib];
            double c = (1 + r) * a + y + (1 - ltv - theta) * p_risky - b;
            double util = utility_terminal(c, b);

            if (c <= 0){util = -1e5;}
            if (util >= VV){
                VV = util;
                CC = c;
                BB = b;
            }
        }

        double[] output = {age, VV, CC, BB, OO, FF, NN, action};
        return output;
    }

    // Interim SFHA Homeowner (1.1, 1.2, 1.31, 1.32, 1.4)
    double[] interim_risky_homeowner(int age, int ia, int ioltv, int ifico, int n, int itheta, int ie,
                                    double[][][][][][][] E_risky, double[][][][] E_rental, double[][][][][][] E_safe,
                                    double q, double p_risky, double p_safe){
        double y = wage(age, ie);
        double a = a_grid[ia];
        double theta = theta_grid[itheta];
        double oltv = oltv_grid[ioltv];
        double fico = fico_grid[ifico];
        double current_x = mortgage_payment(ioltv, ifico, p_risky);
        double current_ltv = current_ltv(ioltv, ifico, n, p_risky);

        // housing related costs
        double maintenance = maintenance_factor * p_risky;
        double closing_cost = closing_factor * p_risky;
        double moving_cost = moving_factor * (p_safe + p_risky) / 2;

        // transfer/subsidy (contingent on being flooded)
        double transfer = 0;
        if (theta > 0){
            transfer = tau;
        }else if (theta == 0){
            transfer = 0;
        }

        double VV_S = -1e5; double VV_C = -1e5; double VV_F = -1e5; double VV_M = -1e5; double VV_MNS = -1e5; double VV_MR = -1e5; double VV_D = -1e5;
        double CC_S = 0; double CC_C = 0; double CC_F = 0; double CC_M = 0; double CC_MNS = 0; double CC_MR = 0; double CC_D = 0;
        double AA_S = 0; double AA_C = 0; double AA_F = 0; double AA_M = 0; double AA_MNS = 0; double AA_MR = 0; double AA_D = 0;
        double OO_S = 0; double OO_C = 0; double OO_F = 0; double OO_M = 0; double OO_MNS = 0; double OO_MR = 0; double OO_D = 0;
        double NN_S = 0; double NN_C = 0; double NN_F = 0; double NN_M = 0; double NN_MNS = 0; double NN_MR = 0; double NN_D = 0;
        double FF_S = 0; double FF_C = 0; double FF_F = 0; double FF_M = 0; double FF_MNS = 0; double FF_MR = 0; double FF_D = 0;
        double action = 1; double action_C = 1.1; double action_F = 1.2; double action_M = 1.3; double action_MNS = 1.31; double action_MR = 1.32; double action_D = 1.4;

        for (int iap = 0; iap < a_grid.length; iap++){
            double a_prime = a_grid[iap];

            // CURRENT: START
            double expected_current = E_risky[age + 1][iap][ioltv][ifico][min(n + 1, 29)][itheta][ie];
            double c_current = (1 + r) * a + y - a_prime - current_x - theta * p_risky - maintenance + transfer;
            double util_current = utility_risky(c_current, 1, theta) + beta * expected_current;
            if (c_current <= 0){util_current = pow(-10, 5);}
            if (util_current >= VV_C){
                VV_C = util_current;
                CC_C = c_current;
                AA_C = a_prime;
                OO_C = oltv;
                FF_C = fico;
                NN_C = n + 1;
                action_C = 1.1;
            }
            // CURRENT: END

            // DEFAULT: START
            double expected_default = E_rental[age + 1][iap][ifico][ie];
            double c_default = (1 + r) * a + y - a_prime - q - moving_cost;
            double util_default = utility_risky(c_default, 1, theta) - d + beta * expected_default;
            if (c_default <= 0){util_default = pow(-10, 5);}
            if (util_default >= VV_D){
                VV_D = util_default;
                CC_D = c_default;
                AA_D = a_prime;
                OO_D = Double.NaN;
                FF_D = 500;
                NN_D = Double.NaN;
                action_D = 1.4;
            }
            // DEFAULT: END

            // MOVE TO RENTAL: START
            double expected_move_to_rent = E_rental[age + 1][iap][ifico][ie];
            double c_move_to_rent = (1 + r) * a + y + (1 - current_ltv - theta) * p_risky - a_prime - q - moving_cost;
            double util_move_to_rent = utility_rental(c_move_to_rent, 0) + beta * expected_move_to_rent;
            if (c_move_to_rent <= 0){util_move_to_rent = pow(-10, 5);}
            if (util_move_to_rent >= VV_MR){
                VV_MR = util_move_to_rent;
                CC_MR = c_move_to_rent;
                AA_MR = a_prime;
                OO_MR = Double.NaN;
                FF_MR = fico;
                NN_MR = Double.NaN;
                action_MR = 1.32;
            }
            // MOVE TO RENTAL: END

            for (int ioltvp = 0; ioltvp < oltv_grid.length; ioltvp++){
                double oltvp = oltv_grid[ioltvp];

                // MOVE TO NON-SFHA: START
                double expected_move_to_safe = E_safe[age + 1][iap][ioltvp][ifico][1][ie];
                double new_x_safe = mortgage_payment(ioltvp, ifico, p_safe);
                double c_move_to_safe = (1 + r) * a + y + (1 - current_ltv - theta) * p_risky - a_prime - new_x_safe - (1 - oltvp) * p_safe - moving_cost - closing_cost;
                double util_move_to_safe = utility_safe(c_move_to_safe, 1) + beta * expected_move_to_safe;
                if (c_move_to_safe <= 0){util_move_to_safe = pow(-10, 5);}
                if (util_move_to_safe >= VV_MNS){
                    VV_MNS = util_move_to_safe;
                    CC_MNS = c_move_to_safe;
                    AA_MNS = a_prime;
                    OO_MNS = oltvp;
                    FF_MNS = fico;
                    NN_MNS = 2;
                    action_MNS = 1.31;
                }
                // MOVE TO NON-SFHA: END

                // REFINANCE: START
                double expected_refinance = E_risky[age + 1][iap][ioltvp][ifico][1][itheta][ie];
                double new_x_refinance = mortgage_payment(ioltvp, ifico, p_risky);
                double c_refinance = (1 + r) * a + y + (oltvp - current_ltv - theta) * p_risky - a_prime - new_x_refinance - theta * p_risky - maintenance - closing_cost + transfer;
                double util_refinance = utility_risky(c_refinance, 1, theta) + beta * expected_refinance;
                if (c_refinance <= 0){util_refinance = pow(-10, 5);}
                if (util_refinance >= VV_F){
                    VV_F = util_refinance;
                    CC_F = c_refinance;
                    AA_F = a_prime;
                    OO_F = oltvp;
                    FF_F = fico;
                    NN_F = 2;
                    action_F = 1.2;
                }
                // REFINANCE: END
            }
        }

        VV_M = max(VV_MNS, VV_MR);
        if (VV_M == VV_MNS){
            CC_M = CC_MNS;
            AA_M = AA_MNS;
            OO_M = OO_MNS;
            FF_M = FF_MNS;
            NN_M = NN_MNS;
            action_M = action_MNS;
        }else if (VV_M == VV_MR){
            CC_M = CC_MR;
            AA_M = AA_MR;
            OO_M = OO_MR;
            FF_M = FF_MR;
            NN_M = NN_MR;
            action_M = action_MR;
        }

        VV_S = max(max(VV_C, VV_F), max(VV_M, VV_D));
        if (VV_S == VV_C){
            CC_S = CC_C;
            AA_S = AA_C;
            OO_S = OO_C;
            FF_S = FF_C;
            NN_S = NN_C;
            action = action_C;
        }else if (VV_S == VV_F){
            CC_S = CC_F;
            AA_S = AA_F;
            OO_S = OO_F;
            FF_S = FF_F;
            NN_S = NN_F;
            action = action_F;
        }else if (VV_S == VV_M){
            CC_S = CC_M;
            AA_S = AA_M;
            OO_S = OO_M;
            FF_S = FF_M;
            NN_S = NN_M;
            action = action_M;
        }else if (VV_S == VV_D){
            CC_S = CC_D;
            AA_S = AA_D;
            OO_S = OO_D;
            FF_S = FF_D;
            NN_S = NN_D;
            action = action_D;
        }


        double[] output = {age, VV_S, CC_S, AA_S, OO_S, FF_S, NN_S, action};
        return output;
    }

    // Newborn (1.0, 2.0, 3.0)
    double[] newborn(int age, int ia, int ifico, int ie,
                     double[][][][][][][] E_risky, double[][][][] E_rental, double[][][][][][] E_safe,
                     double q, double p_risky, double p_safe){

        double y = wage(age, ie);
        double a = a_grid[ia];
        double fico = fico_grid[ifico];
        double moving_cost = moving_factor * (p_safe + p_risky) / 2;
        double theta = theta_grid[0];


        double VV = -1e5; double VV_S = -1e5; double VV_NS = -1e5; double VV_R = -1e5;
        double CC = 0; double CC_S = 0; double CC_NS = 0; double CC_R = 0;
        double AA = 0; double AA_S = 0; double AA_NS = 0; double AA_R = 0;
        double OO = 0; double OO_S = 0; double OO_NS = 0; double OO_R = 0;
        double FF = 0; double FF_S = 0; double FF_NS = 0; double FF_R = 0;
        double NN = 0; double NN_S = 0; double NN_NS = 0; double NN_R = 0;
        double action = 0; double action_risky = 1.0; double action_safe = 3.0; double action_renter = 2.0;

        for (int iap = 0; iap < a_grid.length; iap++) {
            double a_prime = a_grid[iap];

            // RENTAL: START
            double expected_rent = E_rental[age + 1][iap][ifico][ie];
            double c_rent = a + y - q - a_prime - moving_cost;
            double util_rent = utility_rental(c_rent, 0) + beta * expected_rent;
            if (c_rent <= 0) { util_rent = -1e5; }
            if (util_rent >= VV) {
                VV_R = util_rent;
                CC_R = c_rent;
                AA_R = a_prime;
                OO_R = Double.NaN;
                FF_R = fico;
                NN_R = Double.NaN;
            }
            // RENTAL: END

            for (int ioltv = 0; ioltv < oltv_grid.length; ioltv++){
                double oltv = oltv_grid[ioltv];

                // NEWBORN TO NON-SFHA: START
                double expected_safe = E_safe[age + 1][iap][ioltv][ifico][1][ie];
                double x_safe = mortgage_payment(ioltv, ifico, p_safe);
                double c_safe = a + y - x_safe - (1 - oltv) * p_safe - a_prime - moving_cost - (closing_factor * p_safe);
                double util_safe = utility_safe(c_safe, 1) + beta * expected_safe;
                if (c_safe <= 0){util_safe = -1e5;}
                if (util_safe >= VV_NS){
                    VV_NS = util_safe;
                    CC_NS = c_safe;
                    AA_NS = a_prime;
                    OO_NS = oltv;
                    FF_NS = fico;
                    NN_NS = 2;
                }
                // NEWBORN TO NON-SFHA: END

                // NEWBORN TO SFHA: START
                double expected_risky = E_risky[age + 1][iap][ioltv][ifico][1][0][ie];
                double x_risky = mortgage_payment(ioltv, ifico, p_risky);
                double c_risky = a + y - x_risky - (1 - oltv) * p_risky - a_prime - moving_cost - (closing_factor * p_risky);
                double util_risky = utility_risky(c_risky, 1, theta) + beta * expected_risky;
                if (c_risky <= 0){util_risky = -1e5;}
                if (util_risky >= VV_S){
                    VV_S = util_risky;
                    CC_S = c_risky;
                    AA_S = a_prime;
                    OO_S = oltv;
                    FF_S = fico;
                    NN_S = 2;
                }
                // NEWBORN TO SFHA: END
            }
        }

        VV = max(max(VV_S, VV_NS), VV_R);
        if (VV == VV_S){
            CC = CC_S;
            AA = AA_S;
            OO = OO_S;
            FF = FF_S;
            NN = NN_S;
            action = action_risky;
        }else if (VV == VV_NS){
            CC = CC_NS;
            AA = AA_NS;
            OO = OO_NS;
            FF = FF_NS;
            NN = NN_NS;
            action = action_safe;
        } else if (VV == VV_R){
            CC = CC_R;
            AA = AA_R;
            OO = OO_R;
            FF = FF_R;
            NN = NN_R;
            action = action_renter;
        }

        double[] output = {age, VV, CC, AA, OO, FF, NN, action};
        return output;
    }

    Quartet<List<double[][][][][][][]>, List<double[][][][][][]>, List<double[][][][]>, List<double[][][][]>> solve_household_problems(double p_risky, double p_safe, double q, int verbose){

        if (verbose >= 0){
            System.out.println("Solving for policy functions by backward induction.");
        }

        // order of state variables: (age, a, oltv, fico, n, theta, e)
        // order of output array: (age, VV, CC, AA/BB, OO, FF, NN, action)

        int na = a_grid.length;          // number of asset points
        int noltv = oltv_grid.length;    // number of oltv points
        int nfico = fico_grid.length;    // number of fico points
        int ntheta = theta_grid.length;  // number of disaster risk points
        int ne = 2;                      // number of productivity points (employed or not)

        // Building the state space for safe area (dim: 5)
        Integer[][] matrix_safe = new Integer[5][];
        matrix_safe[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_safe[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_safe[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_safe[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_safe[4] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_safe = new CartesianSet<>(matrix_safe);
        int Z_safe = (int) state_space_safe.getCount();

        // Building the state space for rental (dim: 3)
        Integer[][] matrix_rental = new Integer[3][];
        matrix_rental[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_rental[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_rental[2] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_rental = new CartesianSet<>(matrix_rental);
        int Z_rental = (int) state_space_rental.getCount();

        // Building the state space for at-risk area (dim: 6)
        Integer[][] matrix_risky = new Integer[6][];
        matrix_risky[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_risky[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_risky[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_risky[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_risky[4] = ArrayUtils.toObject(IntStream.range(0, ntheta).toArray());
        matrix_risky[5] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_risky = new CartesianSet<>(matrix_risky);
        int Z_risky = (int) state_space_risky.getCount();

        // Building the state space for newborns (dim: 3)
        Integer[][] matrix_newborn = new Integer[3][];
        matrix_newborn[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_newborn[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_newborn[2] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_newborn = new CartesianSet<>(matrix_newborn);
        int Z_newborn = (int) state_space_newborn.getCount();

        // Initialize grid for value functions, policy functions
        double[][][][][][] V_safe = new double[J][na][noltv][nfico][N][ne];
        double[][][][][][] c_safe = new double[J][na][noltv][nfico][N][ne];
        double[][][][][][] a_safe = new double[J][na][noltv][nfico][N][ne];
        double[][][][][][] oltv_safe = new double[J][na][noltv][nfico][N][ne];
        double[][][][][][] fico_safe = new double[J][na][noltv][nfico][N][ne];
        double[][][][][][] n_safe = new double[J][na][noltv][nfico][N][ne];
        double[][][][][][] action_safe = new double[J][na][noltv][nfico][N][ne];

        double[][][][] V_renter = new double[J][na][nfico][ne];
        double[][][][] c_renter = new double[J][na][nfico][ne];
        double[][][][] a_renter = new double[J][na][nfico][ne];
        double[][][][] oltv_renter = new double[J][na][nfico][ne];
        double[][][][] fico_renter = new double[J][na][nfico][ne];
        double[][][][] n_renter = new double[J][na][nfico][ne];
        double[][][][] action_renter = new double[J][na][nfico][ne];

        double[][][][][][][] V_risky = new double[J][na][noltv][nfico][N][ntheta][ne];
        double[][][][][][][] c_risky = new double[J][na][noltv][nfico][N][ntheta][ne];
        double[][][][][][][] a_risky = new double[J][na][noltv][nfico][N][ntheta][ne];
        double[][][][][][][] oltv_risky = new double[J][na][noltv][nfico][N][ntheta][ne];
        double[][][][][][][] fico_risky = new double[J][na][noltv][nfico][N][ntheta][ne];
        double[][][][][][][] n_risky = new double[J][na][noltv][nfico][N][ntheta][ne];
        double[][][][][][][] action_risky = new double[J][na][noltv][nfico][N][ntheta][ne];

        double[][][][] V_newborn = new double[J][na][nfico][ne];
        double[][][][] c_newborn = new double[J][na][nfico][ne];
        double[][][][] a_newborn = new double[J][na][nfico][ne];
        double[][][][] oltv_newborn = new double[J][na][nfico][ne];
        double[][][][] fico_newborn = new double[J][na][nfico][ne];
        double[][][][] n_newborn = new double[J][na][nfico][ne];
        double[][][][] action_newborn = new double[J][na][nfico][ne];

        double[][][][][][][] Expected_risky;
        double[][][][][][] Expected_safe;
        double[][][][] Expected_rental;

        int age = J - 1; // start at 39 and fo to 0;
        if (verbose == 3) {
            System.out.println("age: " + (J - 1));
        }

        // terminal safe region
        for (int z = 0; z < Z_safe; z++){
            List<Integer> node = state_space_safe.get(z);
            int ia = node.get(0);
            int ioltv = node.get(1);
            int ifico = node.get(2);
            int n = node.get(3);
            int ie = node.get(4);
            double[] terminal_safe_homeowner_output = terminal_safe_homeowner(J - 1, ia, ioltv, ifico, n, ie, p_safe);
            V_safe[J - 1][ia][ioltv][ifico][n][ie] = terminal_safe_homeowner_output[1];
            c_safe[J - 1][ia][ioltv][ifico][n][ie] = terminal_safe_homeowner_output[2];
            a_safe[J - 1][ia][ioltv][ifico][n][ie] = terminal_safe_homeowner_output[3];
            oltv_safe[J - 1][ia][ioltv][ifico][n][ie] = terminal_safe_homeowner_output[4];
            fico_safe[J - 1][ia][ioltv][ifico][n][ie] = terminal_safe_homeowner_output[5];
            n_safe[J - 1][ia][ioltv][ifico][n][ie] = terminal_safe_homeowner_output[6];
            action_safe[J - 1][ia][ioltv][ifico][n][ie] = terminal_safe_homeowner_output[7];
        }

        // terminal renter
        for (int z = 0; z < Z_rental; z++){
            List<Integer> node = state_space_rental.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            int ie = node.get(2);
            double[] terminal_renter_output = terminal_renter(J - 1, ia, ifico, ie);
            V_renter[J - 1][ia][ifico][ie] = terminal_renter_output[1];
            c_renter[J - 1][ia][ifico][ie] = terminal_renter_output[2];
            a_renter[J - 1][ia][ifico][ie] = terminal_renter_output[3];
            oltv_renter[J - 1][ia][ifico][ie] = terminal_renter_output[4];
            fico_renter[J - 1][ia][ifico][ie] = terminal_renter_output[5];
            n_renter[J - 1][ia][ifico][ie] = terminal_renter_output[6];
            action_renter[J - 1][ia][ifico][ie] = terminal_renter_output[7];
        }

        // terminal at-risk region
        for (int z = 0; z < Z_risky; z++){
            List<Integer> node = state_space_risky.get(z);
            int ia = node.get(0);
            int ioltv = node.get(1);
            int ifico = node.get(2);
            int n = node.get(3);
            int itheta = node.get(4);
            int ie = node.get(5);
            double[] terminal_risky_homeowner_output = terminal_risky_homeowner(J - 1, ia, ioltv, ifico, n, itheta, ie, p_risky);
            V_risky[J - 1][ia][ioltv][ifico][n][itheta][ie] = terminal_risky_homeowner_output[1];
            c_risky[J - 1][ia][ioltv][ifico][n][itheta][ie] = terminal_risky_homeowner_output[2];
            a_risky[J - 1][ia][ioltv][ifico][n][itheta][ie] = terminal_risky_homeowner_output[3];
            oltv_risky[J - 1][ia][ioltv][ifico][n][itheta][ie] = terminal_risky_homeowner_output[4];
            fico_risky[J - 1][ia][ioltv][ifico][n][itheta][ie] = terminal_risky_homeowner_output[5];
            n_risky[J - 1][ia][ioltv][ifico][n][itheta][ie] = terminal_risky_homeowner_output[6];
            action_risky[J - 1][ia][ioltv][ifico][n][itheta][ie] = terminal_risky_homeowner_output[7];
        }

        Expected_risky = miscFunctions.Expected_risky_income_risk(J, na, noltv, nfico, N, ntheta, ne, V_risky, P_theta, P_e, J - 1);
        Expected_safe = miscFunctions.Expected_safe_income_risk(J, na, noltv, nfico, N, ne, V_safe, P_e, J - 1);
        Expected_rental = miscFunctions.Expected_rental_income_risk(J, na, nfico, ne, V_renter, P_e, J - 1);

        int j;
        for (j = 2; j < J; j++) {
            int final_j = j;
            if (verbose == 3) {
                System.out.println("age: " + (J - final_j));
            }

            double[][][][][][][] final_Expected_risky = Expected_risky;
            double[][][][][][] final_Expected_safe = Expected_safe;
            double[][][][] final_Expected_rental = Expected_rental;

            IntStream.range(0, Z_safe).parallel().forEach(z -> {
                List<Integer> node = state_space_safe.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                int ie = node.get(4);
                double[] interim_safe_homeowner_output = interim_safe_homeowner(J - final_j, ia, ioltv, ifico, n, ie,
                        final_Expected_rental, final_Expected_risky, final_Expected_safe, q, p_risky, p_safe);
                V_safe[J - final_j][ia][ioltv][ifico][n][ie] = interim_safe_homeowner_output[1];
                c_safe[J - final_j][ia][ioltv][ifico][n][ie] = interim_safe_homeowner_output[2];
                a_safe[J - final_j][ia][ioltv][ifico][n][ie] = interim_safe_homeowner_output[3];
                oltv_safe[J - final_j][ia][ioltv][ifico][n][ie] = interim_safe_homeowner_output[4];
                fico_safe[J - final_j][ia][ioltv][ifico][n][ie] = interim_safe_homeowner_output[5];
                n_safe[J - final_j][ia][ioltv][ifico][n][ie] = interim_safe_homeowner_output[6];
                action_safe[J - final_j][ia][ioltv][ifico][n][ie] = interim_safe_homeowner_output[7];
            });

            IntStream.range(0, Z_rental).parallel().forEach(z -> {
                List<Integer> node = state_space_rental.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                int ie = node.get(2);
                double[] interim_renter_output = interim_renter(J - final_j, ia, ifico, ie,
                        final_Expected_rental, final_Expected_risky, final_Expected_safe, q, p_risky, p_safe);
                V_renter[J - final_j][ia][ifico][ie] = interim_renter_output[1];
                c_renter[J - final_j][ia][ifico][ie] = interim_renter_output[2];
                a_renter[J - final_j][ia][ifico][ie] = interim_renter_output[3];
                oltv_renter[J - final_j][ia][ifico][ie] = interim_renter_output[4];
                fico_renter[J - final_j][ia][ifico][ie] = interim_renter_output[5];
                n_renter[J - final_j][ia][ifico][ie] = interim_renter_output[6];
                action_renter[J - final_j][ia][ifico][ie] = interim_renter_output[7];
            });

            IntStream.range(0, Z_risky).parallel().forEach(z -> {
                List<Integer> node = state_space_risky.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                int itheta = node.get(4);
                int ie = node.get(5);
                double[] interim_risky_homeowner_output = interim_risky_homeowner(J - final_j, ia, ioltv, ifico, n, itheta, ie,
                        final_Expected_risky, final_Expected_rental, final_Expected_safe, q, p_risky, p_safe);
                V_risky[J - final_j][ia][ioltv][ifico][n][itheta][ie] = interim_risky_homeowner_output[1];
                c_risky[J - final_j][ia][ioltv][ifico][n][itheta][ie] = interim_risky_homeowner_output[2];
                a_risky[J - final_j][ia][ioltv][ifico][n][itheta][ie] = interim_risky_homeowner_output[3];
                oltv_risky[J - final_j][ia][ioltv][ifico][n][itheta][ie] = interim_risky_homeowner_output[4];
                fico_risky[J - final_j][ia][ioltv][ifico][n][itheta][ie] = interim_risky_homeowner_output[5];
                n_risky[J - final_j][ia][ioltv][ifico][n][itheta][ie] = interim_risky_homeowner_output[6];
                action_risky[J - final_j][ia][ioltv][ifico][n][itheta][ie] = interim_risky_homeowner_output[7];
            });


            Expected_risky = miscFunctions.Expected_risky_income_risk(J, na, noltv, nfico, N, ntheta, ne, V_risky, P_theta, P_e, J - final_j);
            Expected_safe = miscFunctions.Expected_safe_income_risk(J, na, noltv, nfico, N, ne, V_safe, P_e, J - final_j);
            Expected_rental = miscFunctions.Expected_rental_income_risk(J, na, nfico, ne, V_renter, P_e, J - final_j);

        }

        if (verbose == 3){
            System.out.println("age: " + 0);
        }

        for (int z = 0; z < Z_newborn; z++){
            List<Integer> node = state_space_newborn.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            int ie = node.get(2);
            double[] newborn_output = newborn(0, ia, ifico, ie, Expected_risky, Expected_rental, Expected_safe, q, p_risky, p_safe);
            V_newborn[0][ia][ifico][ie] = newborn_output[1];
            c_newborn[0][ia][ifico][ie] = newborn_output[2];
            a_newborn[0][ia][ifico][ie] = newborn_output[3];
            oltv_newborn[0][ia][ifico][ie] = newborn_output[4];
            fico_newborn[0][ia][ifico][ie] = newborn_output[5];
            n_newborn[0][ia][ifico][ie] = newborn_output[6];
            action_newborn[0][ia][ifico][ie] = newborn_output[7];
        }

        List<double[][][][][][]> safe_policy_functions = new ArrayList<>();
        safe_policy_functions.add(0, V_safe);
        safe_policy_functions.add(1, c_safe);
        safe_policy_functions.add(2, a_safe);
        safe_policy_functions.add(3, oltv_safe);
        safe_policy_functions.add(4, fico_safe);
        safe_policy_functions.add(5, n_safe);
        safe_policy_functions.add(6, action_safe);

        List<double[][][][][][][]> risky_policy_functions = new ArrayList<>();
        risky_policy_functions.add(0, V_risky);
        risky_policy_functions.add(1, c_risky);
        risky_policy_functions.add(2, a_risky);
        risky_policy_functions.add(3, oltv_risky);
        risky_policy_functions.add(4, fico_risky);
        risky_policy_functions.add(5, n_risky);
        risky_policy_functions.add(6, action_risky);

        List<double[][][][]> renter_policy_functions = new ArrayList<>();
        renter_policy_functions.add(0, V_renter);
        renter_policy_functions.add(1, c_renter);
        renter_policy_functions.add(2, a_renter);
        renter_policy_functions.add(3, oltv_renter);
        renter_policy_functions.add(4, fico_renter);
        renter_policy_functions.add(5, n_renter);
        renter_policy_functions.add(6, action_renter);

        List<double[][][][]> newborn_policy_functions = new ArrayList<>();
        newborn_policy_functions.add(0, V_newborn);
        newborn_policy_functions.add(1, c_newborn);
        newborn_policy_functions.add(2, a_newborn);
        newborn_policy_functions.add(3, oltv_newborn);
        newborn_policy_functions.add(4, fico_newborn);
        newborn_policy_functions.add(5, n_newborn);
        newborn_policy_functions.add(6, action_newborn);

        Quartet<List<double[][][][][][][]>, List<double[][][][][][]>, List<double[][][][]>, List<double[][][][]>> output =
                new Quartet<>(risky_policy_functions, safe_policy_functions, renter_policy_functions, newborn_policy_functions);

        return output;
    }

    Quartet<double[][][][][][][], double[][][][][][], double[][][][], double[][][][]> compute_distributions_pre_disaster(Quartet<List<double[][][][][][][]>, List<double[][][][][][]>, List<double[][][][]>, List<double[][][][]>> household_policy_functions, int verbose){

        if (verbose >= 0){
            System.out.println("Computing distributions (pre-disaster).");
        }

        // order of state variables: (age, a, oltv, fico, n, theta, e)
        // order of output array: (age, VV, CC, AA/BB, OO, FF, NN, action)

        int na = a_grid.length;          // number of asset points
        int noltv = oltv_grid.length;    // number of oltv points
        int nfico = fico_grid.length;    // number of fico points
        int ntheta = theta_grid.length;  // number of disaster risk points
        int ne = 2;                      // number of productivity points (employed or not)

        // Building the state space for safe area (dim: 5)
        Integer[][] matrix_safe = new Integer[5][];
        matrix_safe[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_safe[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_safe[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_safe[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_safe[4] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_safe = new CartesianSet<>(matrix_safe);
        int Z_safe = (int) state_space_safe.getCount();

        // Building the state space for rental (dim: 3)
        Integer[][] matrix_rental = new Integer[3][];
        matrix_rental[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_rental[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_rental[2] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_rental = new CartesianSet<>(matrix_rental);
        int Z_rental = (int) state_space_rental.getCount();

        // Building the state space for at-risk area (dim: 6)
        Integer[][] matrix_risky = new Integer[6][];
        matrix_risky[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_risky[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_risky[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_risky[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_risky[4] = ArrayUtils.toObject(IntStream.range(0, ntheta).toArray());
        matrix_risky[5] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_risky = new CartesianSet<>(matrix_risky);
        int Z_risky = (int) state_space_risky.getCount();

        // Building the state space for newborns (dim: 3)
        Integer[][] matrix_newborn = new Integer[3][];
        matrix_newborn[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_newborn[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_newborn[2] = ArrayUtils.toObject(IntStream.range(0, ne).toArray());
        CartesianSet<Integer> state_space_newborn = new CartesianSet<>(matrix_newborn);
        int Z_newborn = (int) state_space_newborn.getCount();


        // unpacking the quartet containing all the policy functions
        List<double[][][][][][][]> risky_homeowner_policy_functions = household_policy_functions.getValue0();
        List<double[][][][][][]> safe_homeowner_policy_functions = household_policy_functions.getValue1();
        List<double[][][][]> renter_policy_functions = household_policy_functions.getValue2();
        List<double[][][][]> newborn_policy_functions = household_policy_functions.getValue3();

        // Initialize grid for value functions, policy functions
        double[][][][][][][] V_risky = risky_homeowner_policy_functions.get(0);
        double[][][][][][][] c_risky = risky_homeowner_policy_functions.get(1);
        double[][][][][][][] a_risky = risky_homeowner_policy_functions.get(2);
        double[][][][][][][] oltv_risky = risky_homeowner_policy_functions.get(3);
        double[][][][][][][] fico_risky = risky_homeowner_policy_functions.get(4);
        double[][][][][][][] n_risky = risky_homeowner_policy_functions.get(5);
        double[][][][][][][] action_risky = risky_homeowner_policy_functions.get(6);

        double[][][][][][] V_safe = safe_homeowner_policy_functions.get(0);
        double[][][][][][] c_safe = safe_homeowner_policy_functions.get(1);
        double[][][][][][] a_safe = safe_homeowner_policy_functions.get(2);
        double[][][][][][] oltv_safe = safe_homeowner_policy_functions.get(3);
        double[][][][][][] fico_safe = safe_homeowner_policy_functions.get(4);
        double[][][][][][] n_safe = safe_homeowner_policy_functions.get(5);
        double[][][][][][] action_safe = safe_homeowner_policy_functions.get(6);

        double[][][][] V_renter = renter_policy_functions.get(0);
        double[][][][] c_renter = renter_policy_functions.get(1);
        double[][][][] a_renter = renter_policy_functions.get(2);
        double[][][][] oltv_renter = renter_policy_functions.get(3);
        double[][][][] fico_renter = renter_policy_functions.get(4);
        double[][][][] n_renter = renter_policy_functions.get(5);
        double[][][][] action_renter = renter_policy_functions.get(6);

        double[][][][] V_newborn = newborn_policy_functions.get(0);
        double[][][][] c_newborn = newborn_policy_functions.get(1);
        double[][][][] a_newborn = newborn_policy_functions.get(2);
        double[][][][] oltv_newborn = newborn_policy_functions.get(3);
        double[][][][] fico_newborn = newborn_policy_functions.get(4);
        double[][][][] n_newborn = newborn_policy_functions.get(5);
        double[][][][] action_newborn = newborn_policy_functions.get(6);

        // tensors to hold the distros
        double[][][][][][][] f_interim_risky_homeowner = new double[J][na][noltv][nfico][N + 1][ntheta][ne];
        double[][][][][][] f_interim_safe_homeowner = new double[J][na][noltv][nfico][N + 1][ne];
        double[][][][] f_interim_renter = new double[J][na][nfico][ne];
        double[][][][] f_newborn = new double[J][na][nfico][ne];

        //====================================//
        // INITIALIZING THE DISTRIBUTION
        // set age at zero in preparation for the forward induction process
        int j = 0;

        // also set ithetap = 0 quasi-globally (for the entirety of this method) since this is pre-Sandy simulation/disaster shocks are off
        int ithetap = 0;

        //fix the initial distribution, f_0, uniform weights (for now)
        double[] asset_distribution = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/asset_distribution.csv", na);
        double[] fico_distribution = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/fico_distribution.csv", nfico);
        double[] wage_distribution = {0.95, 0.05};
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                for (int ie = 0; ie < ne; ie++){
                    f_newborn[j][ia][ifico][ie] = asset_distribution[ia] * fico_distribution[ifico] * wage_distribution[ie] / J;
                }
            }
        }

        // just check that cohort weight for j = 0 adds to 1/J
        double sum = 0;
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                for (int ie = 0; ie < ne; ie++){
                    sum = sum + f_newborn[0][ia][ifico][ie];
                }
            }
        }

        if (verbose == 2){
            System.out.println(j + " , " + sum);
        }

        //====================================//

        // NOW TRANSITIONING NEWBORNS TO HOMEOWNERS AND RENTERS
        // all nested within wage loop since that shock effects everyone regardless of housing status
        // no theta loop as this is pre-Sandy so ithetap = 0 (disaster shocks off)
        for (int iep = 0; iep < ne; iep++){
            int final_iep = iep;
            int final_j = j;
            IntStream.range(0, Z_newborn).parallel().forEach(z -> {
                List<Integer> node = state_space_newborn.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                int ie = node.get(2);
                // newborns -> risky homeowners
                if (action_newborn[final_j][ia][ifico][ie] == 1.0){
                    int iap = ArrayUtils.indexOf(a_grid, a_newborn[final_j][ia][ifico][ie]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_newborn[final_j][ia][ifico][ie]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[final_j][ia][ifico][ie]);
                    int np = (int) n_newborn[final_j][ia][ifico][ie];
                    f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep]
                            = f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep]
                            + (1 * P_e[ie][final_iep] * f_newborn[final_j][ia][ifico][ie]);
                }// newborns -> safe homeowners
                else if (action_newborn[final_j][ia][ifico][ie] == 3.0){
                    int iap = ArrayUtils.indexOf(a_grid, a_newborn[final_j][ia][ifico][ie]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_newborn[final_j][ia][ifico][ie]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[final_j][ia][ifico][ie]);
                    int np = (int) n_newborn[final_j][ia][ifico][ie];
                    f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep]
                            = f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep] + (P_e[ie][final_iep] * f_newborn[final_j][ia][ifico][ie]);
                }// newborns -> renters
                else if (action_newborn[final_j][ia][ifico][ie] == 2.0){
                    int iap = ArrayUtils.indexOf(a_grid, a_newborn[final_j][ia][ifico][ie]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[final_j][ia][ifico][ie]);
                    f_interim_renter[final_j + 1][iap][ificop][final_iep]
                            = f_interim_renter[final_j + 1][iap][ificop][final_iep] + (P_e[ie][final_iep] * f_newborn[final_j][ia][ifico][ie]);
                }
            });
        }

        // check cohort weights for j = 1
        double sum_risky = 0; double sum_safe = 0; double sum_renter = 0;
        for (int ia = 0; ia < na; ia++){
            for (int ioltv = 0; ioltv < noltv; ioltv++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int n = 0; n < N; n++){
                        for (int itheta = 0; itheta < ntheta; itheta++){
                            for (int ie = 0; ie < ne; ie++){
                                sum_risky = sum_risky + f_interim_risky_homeowner[j + 1][ia][ioltv][ifico][n][itheta][ie];
                            }
                        }
                    }
                }
            }
        }
        for (int ia = 0; ia < na; ia++){
            for (int ioltv = 0; ioltv < noltv; ioltv++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int n = 0; n < N; n++){
                        for (int ie = 0; ie < ne; ie++){
                            sum_safe = sum_safe + f_interim_safe_homeowner[j + 1][ia][ioltv][ifico][n][ie];
                        }
                    }
                }
            }
        }
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                for (int ie = 0; ie < ne; ie++){
                    sum_renter = sum_renter + f_interim_renter[j + 1][ia][ifico][ie];
                }
            }
        }
        if (verbose == 2){
            System.out.println((j + 1) + " , " + sum_risky + " , " + sum_safe + " , " + sum_renter + " , " + (sum_risky + sum_safe + sum_renter));
        }
        //====================================//


        // NOW TRANSITION THE j=1 cohort to j=2 and then through to end of life, j = J.
        for (j = 1; j < (J - 1); j++){

            // first, a super loop for productivity. See everyone faces this shock, everything should be nested  within this loop
            for (int iep = 0; iep < ne; iep++){

                // since disaster shock is turned off, ithetap = 0 through out no ithetap loop (ithetap = 0) has been declared at the start of this method
                int final_j = j;
                int final_iep = iep;

                for (int z = 0; z < Z_risky; z++){
                    List<Integer> node = state_space_risky.get(z);
                    int ia = node.get(0);
                    int ioltv = node.get(1);
                    int ifico = node.get(2);
                    int n = node.get(3);
                    int itheta = node.get(4);
                    int ie = node.get(5);
                    int iap = ArrayUtils.indexOf(a_grid, a_risky[final_j][ia][ioltv][ifico][n][itheta][ie]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_risky[final_j][ia][ioltv][ifico][n][itheta][ie]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_risky[final_j][ia][ioltv][ifico][n][itheta][ie]);
                    int np = (int) n_risky[final_j][ia][ioltv][ifico][n][itheta][ie];
                    // giving state of an agent, x \in X, find whether it wants to stay current or refinance and if so then for that
                    // agent find x' and transitioning
                    // loop for homeowners wanting to remain in at-risk region (by staying current(1.1) or refinance(1.2))
                    if (action_risky[final_j][ia][ioltv][ifico][n][itheta][ie] == 1.1 | action_risky[final_j][ia][ioltv][ifico][n][itheta][ie] == 1.2){
                        f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep]
                                = f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep] +
                                (P_e[ie][final_iep] * f_interim_risky_homeowner[final_j][ia][ioltv][ifico][n][itheta][ie]);
                    } // do the same but for risky homeowners wishing to move to the safe region (1.31)
                    else if (action_risky[final_j][ia][ioltv][ifico][n][itheta][ie] == 1.31){
                        f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep]
                                = f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep]
                                + (P_e[ie][final_iep] * f_interim_risky_homeowner[final_j][ia][ioltv][ifico][n][itheta][ie]);
                    } // now do the same but for risky homeowners moving to rental (1.32) or defaulting to rental (1.4)
                    else if (action_risky[final_j][ia][ioltv][ifico][n][itheta][ie] == 1.32 | action_risky[final_j][ia][ioltv][ifico][n][itheta][ie] == 1.4){
                        f_interim_renter[final_j + 1][iap][ificop][final_iep]
                                = f_interim_renter[final_j + 1][iap][ificop][final_iep]
                                + (P_e[ie][final_iep] * f_interim_risky_homeowner[final_j][ia][ioltv][ifico][n][itheta][ie]);
                    }else {
                        System.out.println("undefined action (risky)");
                    }
                }

                // first a parallel for loop to transition agents in at-risk region based on their choices
                //IntStream.range(0, Z_risky).parallel().forEach(z -> { });

                //miscFunctions.pause(1);


                for (int z = 0; z < Z_safe; z++){
                    List<Integer> node = state_space_safe.get(z);
                    int ia = node.get(0);
                    int ioltv = node.get(1);
                    int ifico = node.get(2);
                    int n = node.get(3);
                    int ie = node.get(4);
                    int iap = ArrayUtils.indexOf(a_grid, a_safe[final_j][ia][ioltv][ifico][n][ie]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_safe[final_j][ia][ioltv][ifico][n][ie]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_safe[final_j][ia][ioltv][ifico][n][ie]);
                    int np = (int) n_safe[final_j][ia][ioltv][ifico][n][ie];
                    //System.out.println(j + " , " + ia + " , " + ioltv + " , " + ifico + " , " + n + " , " + ie + " , " + iap + " , " + ioltvp + " , " + ificop + " , " + np + " , " + a_safe[final_j][ia][ioltv][ifico][n][ie] + " , " + oltv_safe[final_j][ia][ioltv][ifico][n][ie] + " , " + V_safe[final_j][ia][ioltv][ifico][n][ie]);
                    if (action_safe[final_j][ia][ioltv][ifico][n][ie] < 3.1){
                        System.out.println(action_safe[final_j][ia][ioltv][ifico][n][ie] );
                    }
                    // first, account for remaining and refinancing in the safe region
                    if (action_safe[final_j][ia][ioltv][ifico][n][ie] == 3.1 | action_safe[final_j][ia][ioltv][ifico][n][ie] == 3.2){
                        f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep]
                                = f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep]
                                + (P_e[ie][final_iep] * f_interim_safe_homeowner[final_j][ia][ioltv][ifico][n][ie]);
                    } // do the same but for safe region residents wishing to move to at-risk region (3.31)
                    else if (action_safe[final_j][ia][ioltv][ifico][n][ie] == 3.31){
                        f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep]
                                = f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep]
                                + (P_e[ie][final_iep] * f_interim_safe_homeowner[final_j][ia][ioltv][ifico][n][ie]);
                    } // now do the same but for safe agent residents wanting to move to rental (3.32) or defaulting to rental (3.4)
                    else if (action_safe[final_j][ia][ioltv][ifico][n][ie] == 3.32 | action_safe[final_j][ia][ioltv][ifico][n][ie] == 3.4){
                        f_interim_renter[final_j + 1][iap][ificop][final_iep]
                                = f_interim_renter[final_j + 1][iap][ificop][final_iep]
                                + (P_e[ie][final_iep] * f_interim_safe_homeowner[final_j][ia][ioltv][ifico][n][ie]);
                    }else {
                        System.out.println("undefined action (safe)");
                    }
                }

                 

                // second parallel for loop to transition agents in safe region to their chosen sectors
                //IntStream.range(0, Z_safe).parallel().forEach(z -> { });

                //miscFunctions.pause(1);

                for (int z = 0; z < Z_rental; z++){
                    List<Integer> node = state_space_rental.get(z);
                    int ia = node.get(0);
                    int ifico = node.get(1);
                    int ie = node.get(2);
                    int iap = ArrayUtils.indexOf(a_grid, a_renter[final_j][ia][ifico][ie]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_renter[final_j][ia][ifico][ie]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_renter[final_j][ia][ifico][ie]);
                    int np = (int) n_renter[final_j][ia][ifico][ie];
                    // first, account for agents wanting to remain in rental
                    if (action_renter[final_j][ia][ifico][ie] == 2.1){
                        f_interim_renter[final_j + 1][iap][ificop][final_iep]
                                = f_interim_renter[final_j + 1][iap][ificop][final_iep]
                                + (P_e[ie][final_iep] * f_interim_renter[final_j][ia][ifico][ie]);
                    } // now renters wishing to become risky homeowners
                    else if (action_renter[final_j][ia][ifico][ie] == 2.31) {
                        f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep]
                                = f_interim_risky_homeowner[final_j + 1][iap][ioltvp][ificop][np][ithetap][final_iep]
                                + (P_e[ie][final_iep] * f_interim_renter[final_j][ia][ifico][ie]);
                    } // lastly renters choosing to become safe homeowners
                    else if (action_renter[final_j][ia][ifico][ie] == 2.33) {
                        f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep]
                                = f_interim_safe_homeowner[final_j + 1][iap][ioltvp][ificop][np][final_iep]
                                + (P_e[ie][final_iep] * f_interim_renter[final_j][ia][ifico][ie]);
                    }else {
                        System.out.println("undefined action (renter)");
                    }
                }
                // third parallel for loop to transition agents in the rental sector to their chosen sectors
                //IntStream.range(0, Z_rental).parallel().forEach(z -> { });

                //miscFunctions.pause(1);


            }

            // check cohort weights for j = 2
            sum_risky = 0; sum_safe = 0; sum_renter = 0;
            for (int ia = 0; ia < na; ia++){
                for (int ioltv = 0; ioltv < noltv; ioltv++){
                    for (int ifico = 0; ifico < nfico; ifico++){
                        for (int n = 0; n < N; n++){
                            for (int itheta = 0; itheta < ntheta; itheta++){
                                for (int ie = 0; ie < ne; ie++){
                                    sum_risky = sum_risky + f_interim_risky_homeowner[j + 1][ia][ioltv][ifico][n][itheta][ie];
                                }
                            }
                        }
                    }
                }
            }
            for (int ia = 0; ia < na; ia++){
                for (int ioltv = 0; ioltv < noltv; ioltv++){
                    for (int ifico = 0; ifico < nfico; ifico++){
                        for (int n = 0; n < N; n++){
                            for (int ie = 0; ie < ne; ie++){
                                sum_safe = sum_safe + f_interim_safe_homeowner[j + 1][ia][ioltv][ifico][n][ie];
                            }
                        }
                    }
                }
            }
            for (int ia = 0; ia < na; ia++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int ie = 0; ie < ne; ie++){
                        sum_renter = sum_renter + f_interim_renter[j + 1][ia][ifico][ie];
                    }
                }
            }
            if (verbose == 2){
                System.out.println((j + 1) + " , " + sum_risky + " , " + sum_safe + " , " + sum_renter + " , " + (sum_risky + sum_safe + sum_renter));
            }

        }

        Quartet<double[][][][][][][], double[][][][][][], double[][][][], double[][][][]> output = new Quartet<>(f_interim_risky_homeowner, f_interim_safe_homeowner, f_interim_renter, f_newborn);
        return output;


    }





    public static void main(String[] args){

        long start_time = System.nanoTime();

        DisasterAndIncomeRisk model = new DisasterAndIncomeRisk();
        Quartet<List<double[][][][][][][]>, List<double[][][][][][]>, List<double[][][][]>, List<double[][][][]>> household_policy_functions = model.solve_household_problems(7.23, 5.78, 0.28, 3);
        Quartet<double[][][][][][][], double[][][][][][], double[][][][], double[][][][]> pre_disaster_distributions = model.compute_distributions_pre_disaster(household_policy_functions, 2);

        long end_time = System.nanoTime();
        System.out.println("run time = " + (end_time - start_time) * 1e-9 / 60 + " mins");
    }

































}
