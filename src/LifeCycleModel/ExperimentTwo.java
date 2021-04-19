package LifeCycleModel;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.javatuples.Quartet;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.IntStream;

import static java.lang.Math.*;

/**
 * Experiment two refers to the change in mass, welfare, and other outcomes as the cost of moving into the SFHA/risky region increases.
 * The idea is to test the role of infrastructure (bridges, roads, broadband) in inducing people to move to certain locations; a higher
 * moving cost (moving to risky regions) is a reduced form way of introducing the idea of poor infrastructure in these regions. How
 * would such a cost affect the decision to move to these areas?
 */

public class ExperimentTwo {

    MiscFunctions miscFunctions = new MiscFunctions();

    int J = 50;                      // lifespan in years

    // preferences
    double gamma = 2;                // CRRA utility term
    double psi = 5;                  // bequest utility term
    double xi = 0.57;                // bequest utility term
    double alpha_non_sfha = 0.2;     // homeownership utility in non-SFHA
    double alpha_sfha = 0.32;        // homeownership utility in SFHA; alpha_sfha > alpha_non_sfha
    double d = 0;                    // default penalty (in utile)
    double beta = 0.97;              // discount factor

    // subsidy/transfer due to disaster shock
    double tau;

    // wages and interest rate (risk less)
    double g = 0.02;           // growth rate for wages (annual)
    double y0 = 1;             // base wage
    double r = 0.03;           // risk-free interest rate
    double[] wage_grid = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/wage_grid.csv", 6);

    // housing related costs
    double price_non_sfha = 0.80 * 7.23;
    double maintenance_factor = 0.03;
    double closing_factor = 0.05;
    double moving_factor = 0.12;

    // this is the main experiment value
    double additional_moving_cost; // the additional cost faced by an agent wishing to move into the SFHA, only faced by newborns and renters

    // construction sector param(s)
    double alpha_construction_sector = 0.6;
    double l_bar = 0.00078;

    // asset params
    double amin = 1e-3; double amax = 100;
    double[] a_grid = {amin, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5};
    double[] b_grid = {amin, 0.1, 0.2, 0.3, 0.4, 0.5, 0.75, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5};
    //double[] a_grid = {amin, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};
    //double[] b_grid = {amin, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};

    // disaster risk
    double[] theta_grid = {0.00, 0.22};         // values for disaster risk
    double[][] Pi = {{0.99, 0.01}, {0.98, 0.02}};   // transition matrix for Markov chain


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

    // calculate mortgage payment x given (oltv, fico) and price of housing (sfha vs non-sfha)
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
    double wage(int age){
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

        return y;
    }

    // construction sector
    double construction_sector(double p_sfha){
        double construction_supply = pow(alpha_construction_sector * p_sfha, alpha_construction_sector / (1 - alpha_construction_sector)) * l_bar;
        return construction_supply;
    }

    // various utility functions
    // terminal utility function (no utility from housing)
    double utility_terminal(double c, double b){
        return (pow(c, 1- gamma) / (1 - gamma)) + (psi * pow(b + xi, 1 - gamma) / (1 - gamma));
    }

    // utility from living in non-SFHA (q = 1 always)
    double utility_non_sfha(double c, double q){
        return (pow(c, 1- gamma) / (1 - gamma)) + (alpha_non_sfha * q);
    }

    // utility from living in SFHA (q = 1 except for default)
    double utility_sfha(double c, double q, double theta){
        return (pow(c, 1- gamma) / (1 - gamma)) + (alpha_sfha * (1 - theta) * q);
    }

    // utility from living in a rental (q = 0 always)
    double utility_rental(double c, double q){
        return pow(c, 1- gamma) / (1 - gamma);
    }

    // Non-SFHA Homeowner
    // Terminal Non-SFHA Homeowner (3.5)
    double[] terminal_non_sfha_homeowner(int age, int ia, int ioltv, int ifico, int n, double p_non_sfha){

        double y = wage(age);
        double ltv = current_ltv(ioltv, ifico, n, p_non_sfha);
        double a = a_grid[ia];

        double VV = -1e5;
        double CC = 0;
        double BB = 0;

        for (int ib = 0; ib < b_grid.length; ib++){
            double b = b_grid[ib];
            double c = (1 + r) * a + y + (1 - ltv) * p_non_sfha - b;
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

    // Interim Non-SFHA Homeowner (3.1)
    double[] interim_non_sfha_homeowner(int age, int ia, int ioltv, int ifico, int n,
                                        double[][][][][] V_non_sfha, double p_non_sfha){

        double y = wage(age);
        double current_x = mortgage_payment(ioltv, ifico, p_non_sfha);
        double current_ltv = current_ltv(ioltv, ifico, n, p_non_sfha);
        double a = a_grid[ia];

        double VV = -1e5; double VV_C = -1e5; double VV_F = -1e5;
        double CC = 0; double CC_C = 0; double CC_F = 0;
        double AA = 0; double AA_C = 0; double AA_F = 0;
        double OO = 0; double OO_C = 0; double OO_F = 0;
        double FF = fico_grid[ifico]; double FF_C = fico_grid[ifico]; double FF_F = fico_grid[ifico];
        double NN = 0; double NN_C = 0; double NN_F = 0;
        double action = 3.0; double action_current = 3.1; double action_refinance = 3.2;

        for (int iap = 0; iap < a_grid.length; iap++){
            double a_prime = a_grid[iap];

            // CURRENT: START
            double expected_current = V_non_sfha[age + 1][iap][ioltv][ifico][min(n + 1, 29)];
            double c_current = (1 + r) * a + y - a_prime - current_x - (maintenance_factor * p_non_sfha);
            double util_current = utility_non_sfha(c_current, 1) + beta * expected_current;
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

            // REFINANCE: START
            for (int ioltvp = 0; ioltvp < oltv_grid.length; ioltvp++){
                double oltvp = oltv_grid[ioltvp];
                double expected_refinance = V_non_sfha[age + 1][iap][ioltvp][ifico][1];
                double new_x_refinance = mortgage_payment(ioltvp, ifico, p_non_sfha);
                double c_refinance = (1 + r) * a + y + (oltvp - current_ltv) * p_non_sfha - a_prime - new_x_refinance
                        - (maintenance_factor * p_non_sfha) - (closing_factor * p_non_sfha);
                double util_refinance = utility_non_sfha(c_refinance, 1) + beta * expected_refinance;
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

        VV = max(VV_C, VV_F);
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
        }




        double[] output = {age, VV, CC, AA, OO, FF, NN, action};
        return output;
    }

    // Rental Market
    // Terminal Renter (2.5)
    double[] terminal_renter(int age, int ia, int ifico){
        double y = wage(age);
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

    // Interim Renter (2.1, 2.3)
    double[] interim_renter(int age, int ia, int ifico, double[][][] V_rental, double[][][][][][] E_sfha, double q, double p_sfha){

        double y = wage(age);
        double a = a_grid[ia];

        // housing related costs
        double maintenance = maintenance_factor * p_sfha;
        double closing_cost = closing_factor * p_sfha;
        double moving_cost = moving_factor * y;


        double VV_R = -1e5; double VV_RR = -1e5; double VV_RMS = -1e5;
        double CC_R = 0; double CC_RR = 0; double CC_RMS = 0;
        double AA_R = 0; double AA_RR = 0; double AA_RMS = 0;
        double OO_R = 0; double OO_RR = Double.NaN; double OO_RMS = 0;
        double NN_R = 0; double NN_RR = Double.NaN; double NN_RMS = 2;
        double action_renter = 2.0; double action_stay_renter = 2.1; double action_move_sfha = 2.3;

        for (int iap = 0; iap < a_grid.length; iap++){
            double a_prime = a_grid[iap];

            // STAY RENTER: START
            double expected_rental = V_rental[age + 1][iap][ifico];
            double c_rental = (1 + r) * a + y - q - a_prime;
            double util_rental = utility_rental(c_rental, 0) + beta * expected_rental;

            if (c_rental <= 0){util_rental = -1e5;}
            if (util_rental >= VV_RR){
                VV_RR = util_rental;
                CC_RR = c_rental;
                AA_RR = a_prime;
                OO_RR = Double.NaN;
                action_renter = action_stay_renter;
            }
            // STAY RENTER: END

            // MOVE TO SFHA: START
            for (int ioltv = 0; ioltv < oltv_grid.length; ioltv++){
                double oltv = oltv_grid[ioltv];
                double new_x = mortgage_payment(ioltv, ifico, p_sfha);
                double expected_sfha = E_sfha[age + 1][iap][ioltv][ifico][1][0];
                double c_sfha = (1 + r) * a + y - a_prime - new_x - (1 - oltv) * p_sfha - maintenance - closing_cost - moving_cost - additional_moving_cost;
                double util_sfha = utility_sfha(c_sfha, 1, 0) + beta * expected_sfha;

                if (c_sfha <= 0){util_sfha = -1e5;}
                if (util_sfha >= VV_RMS){
                    VV_RMS = util_sfha;
                    CC_RMS = c_sfha;
                    AA_RMS = a_prime;
                    OO_RMS = oltv;
                    action_renter = action_move_sfha;
                }
            }
            // MOVE TO SFHA: END
        }

        VV_R = max(VV_RR, VV_RMS);
        if (VV_R == VV_RR){
            CC_R = CC_RR;
            AA_R = AA_RR;
            OO_R = OO_RR;
            NN_R = NN_RR;
            action_renter = action_stay_renter;
        }else if (VV_R == VV_RMS){
            CC_R = CC_RMS;
            AA_R = AA_RMS;
            OO_R = OO_RMS;
            NN_R = NN_RMS;
            action_renter = action_move_sfha;
        }

        double FF_R = fico_grid[ifico];

        double[] output = {age, VV_R, CC_R, AA_R, OO_R, FF_R, NN_R, action_renter};
        return output;
    }

    // SFHA Homeowner
    // Terminal SFHA Homeowner (1.5)
    double[] terminal_sfha_homeowner(int age, int ia, int ioltv, int ifico, int n, int itheta, double p_sfha){

        double y = wage(age);
        double a = a_grid[ia];
        double ltv = current_ltv(ioltv, ifico, n, p_sfha);
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
            double c = (1 + r) * a + y + (1 - ltv - theta) * p_sfha - b;
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
    double[] interim_sfha_homeowner(int age, int ia, int ioltv, int ifico, int n, int itheta,
                                    double[][][][][][] E_sfha, double[][][] V_rental, double[][][][][] V_non_sfha,
                                    double q, double p_sfha, double p_non_sfha){
        double y = wage(age);
        double a = a_grid[ia];
        double theta = theta_grid[itheta];
        double oltv = oltv_grid[ioltv];
        double fico = fico_grid[ifico];
        double current_x = mortgage_payment(ioltv, ifico, p_sfha);
        double current_ltv = current_ltv(ioltv, ifico, n, p_sfha);

        // housing related costs
        double maintenance = maintenance_factor * p_sfha;
        double closing_cost = closing_factor * p_sfha;
        double moving_cost = moving_factor * y;

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
            double expected_current = E_sfha[age + 1][iap][ioltv][ifico][min(n + 1, 29)][itheta];
            double c_current = (1 + r) * a + y - a_prime - current_x - theta * p_sfha - maintenance + transfer;
            double util_current = utility_sfha(c_current, 1, theta) + beta * expected_current;
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
            double expected_default = V_rental[age + 1][iap][ifico];
            double c_default = (1 + r) * a + y - a_prime - q - moving_cost;
            double util_default = utility_sfha(c_default, 1, theta) - d + beta * expected_default;
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
            double expected_move_to_rent = V_rental[age + 1][iap][ifico];
            double c_move_to_rent = (1 + r) * a + y + (1 - current_ltv - theta) * p_sfha - a_prime - q - moving_cost;
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
                double expected_move_to_non_sfha = V_non_sfha[age + 1][iap][ioltvp][ifico][1];
                double new_x_non_sfha = mortgage_payment(ioltvp, ifico, p_non_sfha);
                double c_move_to_non_sfha = (1 + r) * a + y + (1 - current_ltv - theta) * p_sfha - a_prime - new_x_non_sfha - (1 - oltvp) * p_non_sfha - moving_cost - closing_cost;
                double util_move_to_non_sfha = utility_non_sfha(c_move_to_non_sfha, 1) + beta * expected_move_to_non_sfha;
                if (c_move_to_non_sfha <= 0){util_move_to_non_sfha = pow(-10, 5);}
                if (util_move_to_non_sfha >= VV_MNS){
                    VV_MNS = util_move_to_non_sfha;
                    CC_MNS = c_move_to_non_sfha;
                    AA_MNS = a_prime;
                    OO_MNS = oltvp;
                    FF_MNS = fico;
                    NN_MNS = 2;
                    action_MNS = 1.31;
                }
                // MOVE TO NON-SFHA: END

                // REFINANCE: START
                double expected_refinance = E_sfha[age + 1][iap][ioltvp][ifico][1][itheta];
                double new_x_refinance = mortgage_payment(ioltvp, ifico, p_sfha);
                double c_refinance = (1 + r) * a + y + (oltvp - current_ltv - theta) * p_sfha - a_prime - new_x_refinance - theta * p_sfha - maintenance - closing_cost + transfer;
                double util_refinance = utility_sfha(c_refinance, 1, theta) + beta * expected_refinance;
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
    double[] newborn(int age, int ia, int ifico, double[][][][][][] E_sfha, double[][][] V_rental, double[][][][][] V_non_sfha,
                     double q, double p_sfha, double p_non_sfha){
        double y = wage(age);
        double a = a_grid[ia];
        double fico = fico_grid[ifico];
        double moving_cost = moving_factor * y;
        double theta = theta_grid[0];


        double VV = -1e5; double VV_S = -1e5; double VV_NS = -1e5; double VV_R = -1e5;
        double CC = 0; double CC_S = 0; double CC_NS = 0; double CC_R = 0;
        double AA = 0; double AA_S = 0; double AA_NS = 0; double AA_R = 0;
        double OO = 0; double OO_S = 0; double OO_NS = 0; double OO_R = 0;
        double FF = 0; double FF_S = 0; double FF_NS = 0; double FF_R = 0;
        double NN = 0; double NN_S = 0; double NN_NS = 0; double NN_R = 0;
        double action = 0; double action_sfha = 1.0; double action_non_sfha = 3.0; double action_renter = 2.0;

        for (int iap = 0; iap < a_grid.length; iap++) {
            double a_prime = a_grid[iap];

            // RENTAL: START
            double expected_rent = V_rental[age + 1][iap][ifico];
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
                double expected_non_sfha = V_non_sfha[age + 1][iap][ioltv][ifico][1];
                double x_non_sfha = mortgage_payment(ioltv, ifico, p_non_sfha);
                double c_non_sfha = a + y - x_non_sfha - (1 - oltv) * p_non_sfha - a_prime - moving_cost - (closing_factor * p_non_sfha);
                double util_non_sfha = utility_non_sfha(c_non_sfha, 1) + beta * expected_non_sfha;
                if (c_non_sfha <= 0){util_non_sfha = -1e5;}
                if (util_non_sfha >= VV_NS){
                    VV_NS = util_non_sfha;
                    CC_NS = c_non_sfha;
                    AA_NS = a_prime;
                    OO_NS = oltv;
                    FF_NS = fico;
                    NN_NS = 2;
                }
                // NEWBORN TO NON-SFHA: END

                // NEWBORN TO SFHA: START
                double expected_sfha = E_sfha[age + 1][iap][ioltv][ifico][1][0];
                double x_sfha = mortgage_payment(ioltv, ifico, p_sfha);
                double c_sfha = a + y - x_sfha - (1 - oltv) * p_sfha - a_prime - moving_cost - (closing_factor * p_sfha) - additional_moving_cost;
                double util_sfha = utility_sfha(c_sfha, 1, theta) + beta * expected_sfha;
                if (c_sfha <= 0){util_sfha = -1e5;}
                if (util_sfha >= VV_S){
                    VV_S = util_sfha;
                    CC_S = c_sfha;
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
            action = action_sfha;
        }else if (VV == VV_NS){
            CC = CC_NS;
            AA = AA_NS;
            OO = OO_NS;
            FF = FF_NS;
            NN = NN_NS;
            action = action_non_sfha;
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

    // solve household problems
    Quartet<List<double[][][][][][]>, List<double[][][][][]>, List<double[][][]>, List<double[][][]>> solve_household_problems(double p_sfha, double p_non_sfha, String verbose){

        System.out.println("Solving for policy functions by backward induction.");

        double q = 0.28;

        // order of state variables: (age, a, oltv, fico, n, theta)
        // order of output array: (VV, CC, AA/BB, OO, status)

        int na = a_grid.length;          // number of asset points
        int noltv = oltv_grid.length;    // number of oltv points
        int nfico = fico_grid.length;    // number of fico points
        int ntheta = theta_grid.length;  // number of disaster risk points

        // Building the state space for non-SFHA (dim: 4)
        Integer[][] matrix_non_sfha = new Integer[4][];
        matrix_non_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_non_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_non_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_non_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        CartesianSet<Integer> state_space_non_sfha = new CartesianSet<>(matrix_non_sfha);
        int Z_non_sfha= (int) state_space_non_sfha.getCount();

        // Building the state space for rental (dim: 2)
        Integer[][] matrix_rental = new Integer[2][];
        matrix_rental[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_rental[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_rental = new CartesianSet<>(matrix_rental);
        int Z_rental = (int) state_space_rental.getCount();

        // Building the state space for SFHA (dim: 5)
        Integer[][] matrix_sfha = new Integer[5][];
        matrix_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_sfha[4] = ArrayUtils.toObject(IntStream.range(0, ntheta).toArray());
        CartesianSet<Integer> state_space_sfha = new CartesianSet<>(matrix_sfha);
        int Z_sfha= (int) state_space_sfha.getCount();

        // Building the state space for newborns (dim: 2)
        Integer[][] matrix_newborn = new Integer[2][];
        matrix_newborn[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_newborn[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_newborn = new CartesianSet<>(matrix_newborn);
        int Z_newborn = (int) state_space_newborn.getCount();

        // Initialize grid for value functions, policy functions
        double[][][][][] V_non_sfha = new double[J][na][noltv][nfico][N];
        double[][][][][] c_non_sfha = new double[J][na][noltv][nfico][N];
        double[][][][][] a_non_sfha = new double[J][na][noltv][nfico][N];
        double[][][][][] oltv_non_sfha = new double[J][na][noltv][nfico][N];
        double[][][][][] fico_non_sfha = new double[J][na][noltv][nfico][N];
        double[][][][][] n_non_sfha = new double[J][na][noltv][nfico][N];
        double[][][][][] action_non_sfha = new double[J][na][noltv][nfico][N];

        double[][][] V_renter = new double[J][na][nfico];
        double[][][] c_renter = new double[J][na][nfico];
        double[][][] a_renter = new double[J][na][nfico];
        double[][][] oltv_renter = new double[J][na][nfico];
        double[][][] fico_renter = new double[J][na][nfico];
        double[][][] n_renter = new double[J][na][nfico];
        double[][][] action_renter = new double[J][na][nfico];

        double[][][][][][] V_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double[][][][][][] c_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double[][][][][][] a_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double[][][][][][] oltv_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double[][][][][][] fico_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double[][][][][][] n_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double[][][][][][] action_sfha = new double[J][na][noltv][nfico][N][ntheta];

        double[][][] V_newborn = new double[J][na][nfico];
        double[][][] c_newborn = new double[J][na][nfico];
        double[][][] a_newborn = new double[J][na][nfico];
        double[][][] oltv_newborn = new double[J][na][nfico];
        double[][][] fico_newborn = new double[J][na][nfico];
        double[][][] n_newborn = new double[J][na][nfico];
        double[][][] action_newborn = new double[J][na][nfico];

        //double[][][][][][] E_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double[][][][][][] Expected_sfha;

        int age = J - 1; // start at 39 and fo to 0;
        if (verbose.equals("yes")) {
            System.out.println("age: " + (J - 1));
        }

        // terminal non-sfha
        for (int z = 0; z < Z_non_sfha; z++){
            List<Integer> node = state_space_non_sfha.get(z);
            int ia = node.get(0);
            int ioltv = node.get(1);
            int ifico = node.get(2);
            int n = node.get(3);
            double[] terminal_non_sfha_homeowner_output = terminal_non_sfha_homeowner(J - 1, ia, ioltv, ifico, n, p_non_sfha);
            V_non_sfha[J - 1][ia][ioltv][ifico][n] = terminal_non_sfha_homeowner_output[1];
            c_non_sfha[J - 1][ia][ioltv][ifico][n] = terminal_non_sfha_homeowner_output[2];
            a_non_sfha[J - 1][ia][ioltv][ifico][n] = terminal_non_sfha_homeowner_output[3];
            oltv_non_sfha[J - 1][ia][ioltv][ifico][n] = terminal_non_sfha_homeowner_output[4];
            fico_non_sfha[J - 1][ia][ioltv][ifico][n] = terminal_non_sfha_homeowner_output[5];
            n_non_sfha[J - 1][ia][ioltv][ifico][n] = terminal_non_sfha_homeowner_output[6];
            action_non_sfha[J - 1][ia][ioltv][ifico][n] = terminal_non_sfha_homeowner_output[7];
        }

        // terminal renter
        for (int z = 0; z < Z_rental; z++){
            List<Integer> node = state_space_rental.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            double[] terminal_renter_output = terminal_renter(J - 1, ia, ifico);
            V_renter[J - 1][ia][ifico] = terminal_renter_output[1];
            c_renter[J - 1][ia][ifico] = terminal_renter_output[2];
            a_renter[J - 1][ia][ifico] = terminal_renter_output[3];
            oltv_renter[J - 1][ia][ifico] = terminal_renter_output[4];
            fico_renter[J - 1][ia][ifico] = terminal_renter_output[5];
            n_renter[J - 1][ia][ifico] = terminal_renter_output[6];
            action_renter[J - 1][ia][ifico] = terminal_renter_output[7];
        }

        for (int z = 0; z < Z_sfha; z++){
            List<Integer> node = state_space_sfha.get(z);
            int ia = node.get(0);
            int ioltv = node.get(1);
            int ifico = node.get(2);
            int n = node.get(3);
            int itheta = node.get(4);
            double[] terminal_sfha_homeowner_output = terminal_sfha_homeowner(J - 1, ia, ioltv, ifico, n, itheta, p_sfha);
            V_sfha[J - 1][ia][ioltv][ifico][n][itheta] = terminal_sfha_homeowner_output[1];
            c_sfha[J - 1][ia][ioltv][ifico][n][itheta] = terminal_sfha_homeowner_output[2];
            a_sfha[J - 1][ia][ioltv][ifico][n][itheta] = terminal_sfha_homeowner_output[3];
            oltv_sfha[J - 1][ia][ioltv][ifico][n][itheta] = terminal_sfha_homeowner_output[4];
            fico_sfha[J - 1][ia][ioltv][ifico][n][itheta] = terminal_sfha_homeowner_output[5];
            n_sfha[J - 1][ia][ioltv][ifico][n][itheta] = terminal_sfha_homeowner_output[6];
            action_sfha[J - 1][ia][ioltv][ifico][n][itheta] = terminal_sfha_homeowner_output[7];
        }

        Expected_sfha = miscFunctions.Expected_sfha(J, na, noltv, nfico, N, ntheta, V_sfha, Pi, J - 1);

        int j;
        for (j = 2; j < J; j++){
            int final_j = j;
            if (verbose.equals("yes")){
                System.out.println("age: " + (J - final_j));
            }


            double[][][][][][] finalExpected_sfha = Expected_sfha;

            IntStream.range(0, Z_non_sfha).parallel().forEach(z -> {
                List<Integer> node = state_space_non_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                double[] interim_non_sfha_homeowner = interim_non_sfha_homeowner(J - final_j, ia, ioltv, ifico, n, V_non_sfha, p_non_sfha);
                V_non_sfha[J - final_j][ia][ioltv][ifico][n] = interim_non_sfha_homeowner[1];
                c_non_sfha[J - final_j][ia][ioltv][ifico][n] = interim_non_sfha_homeowner[2];
                a_non_sfha[J - final_j][ia][ioltv][ifico][n] = interim_non_sfha_homeowner[3];
                oltv_non_sfha[J - final_j][ia][ioltv][ifico][n] = interim_non_sfha_homeowner[4];
                fico_non_sfha[J - final_j][ia][ioltv][ifico][n] = interim_non_sfha_homeowner[5];
                n_non_sfha[J - final_j][ia][ioltv][ifico][n] = interim_non_sfha_homeowner[6];
                action_non_sfha[J - final_j][ia][ioltv][ifico][n] = interim_non_sfha_homeowner[7];
            });

            IntStream.range(0, Z_rental).parallel().forEach(z -> {
                List<Integer> node = state_space_rental.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                double[] interim_renter = interim_renter(J - final_j, ia, ifico, V_renter, finalExpected_sfha, q, p_sfha);
                V_renter[J - final_j][ia][ifico] = interim_renter[1];
                c_renter[J - final_j][ia][ifico] = interim_renter[2];
                a_renter[J - final_j][ia][ifico] = interim_renter[3];
                oltv_renter[J - final_j][ia][ifico] = interim_renter[4];
                fico_renter[J - final_j][ia][ifico] = interim_renter[5];
                n_renter[J - final_j][ia][ifico] = interim_renter[6];
                action_renter[J - final_j][ia][ifico] = interim_renter[7];
            });


            IntStream.range(0, Z_sfha).parallel().forEach(z -> {
                List<Integer> node = state_space_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                int itheta = node.get(4);
                double[] interim_sfha_homeowner_output = interim_sfha_homeowner(J - final_j, ia, ioltv, ifico, n, itheta, finalExpected_sfha, V_renter, V_non_sfha, q, p_sfha, p_non_sfha);
                V_sfha[J - final_j][ia][ioltv][ifico][n][itheta] = interim_sfha_homeowner_output[1];
                c_sfha[J - final_j][ia][ioltv][ifico][n][itheta] = interim_sfha_homeowner_output[2];
                a_sfha[J - final_j][ia][ioltv][ifico][n][itheta] = interim_sfha_homeowner_output[3];
                oltv_sfha[J - final_j][ia][ioltv][ifico][n][itheta] = interim_sfha_homeowner_output[4];
                fico_sfha[J - final_j][ia][ioltv][ifico][n][itheta] = interim_sfha_homeowner_output[5];
                n_sfha[J - final_j][ia][ioltv][ifico][n][itheta] = interim_sfha_homeowner_output[6];
                action_sfha[J - final_j][ia][ioltv][ifico][n][itheta] = interim_sfha_homeowner_output[7];
            });

            Expected_sfha = miscFunctions.Expected_sfha(J, na, noltv, nfico, N, ntheta, V_sfha, Pi, J - final_j);

        }

        if (verbose.equals("yes")){
            System.out.println("age: " + 0);
        }

        for (int z = 0; z < Z_newborn; z++){
            List<Integer> node = state_space_newborn.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            double[] newborn_output = newborn(0, ia, ifico, Expected_sfha, V_renter, V_non_sfha, q, p_sfha, p_non_sfha);
            V_newborn[0][ia][ifico] = newborn_output[1];
            c_newborn[0][ia][ifico] = newborn_output[2];
            a_newborn[0][ia][ifico] = newborn_output[3];
            oltv_newborn[0][ia][ifico] = newborn_output[4];
            fico_newborn[0][ia][ifico] = newborn_output[5];
            n_newborn[0][ia][ifico] = newborn_output[6];
            action_newborn[0][ia][ifico] = newborn_output[7];
        }

        List<double[][][][][]> non_sfha_policy_functions = new ArrayList<>();
        non_sfha_policy_functions.add(0, c_non_sfha);
        non_sfha_policy_functions.add(1, a_non_sfha);
        non_sfha_policy_functions.add(2, oltv_non_sfha);
        non_sfha_policy_functions.add(3, fico_non_sfha);
        non_sfha_policy_functions.add(4, n_non_sfha);
        non_sfha_policy_functions.add(5, action_non_sfha);
        non_sfha_policy_functions.add(6, V_non_sfha);

        List<double[][][][][][]> sfha_policy_functions = new ArrayList<>();
        sfha_policy_functions.add(0, c_sfha);
        sfha_policy_functions.add(1, a_sfha);
        sfha_policy_functions.add(2, oltv_sfha);
        sfha_policy_functions.add(3, fico_sfha);
        sfha_policy_functions.add(4, n_sfha);
        sfha_policy_functions.add(5, action_sfha);
        sfha_policy_functions.add(6, V_sfha);

        List<double[][][]> renter_policy_functions = new ArrayList<>();
        renter_policy_functions.add(0, c_renter);
        renter_policy_functions.add(1, a_renter);
        renter_policy_functions.add(2, oltv_renter);
        renter_policy_functions.add(3, fico_renter);
        renter_policy_functions.add(4, n_renter);
        renter_policy_functions.add(5, action_renter);
        renter_policy_functions.add(6, V_renter);

        List<double[][][]> newborn_policy_functions = new ArrayList<>();
        newborn_policy_functions.add(0, c_newborn);
        newborn_policy_functions.add(1, a_newborn);
        newborn_policy_functions.add(2, oltv_newborn);
        newborn_policy_functions.add(3, fico_newborn);
        newborn_policy_functions.add(4, n_newborn);
        newborn_policy_functions.add(5, action_newborn);
        newborn_policy_functions.add(6, V_newborn);

        Quartet<List<double[][][][][][]>, List<double[][][][][]>, List<double[][][]>, List<double[][][]>> output =
                new Quartet<>(sfha_policy_functions, non_sfha_policy_functions, renter_policy_functions, newborn_policy_functions);


        return output;
    }

    // compute and forward distributions for the four groups: newborns, renters, non-sfha, and sfha
    Quartet<double[][][][][][], double[][][][][], double[][][], double[][][]> compute_distributions_shocks_off(
            Quartet<List<double[][][][][][]>, List<double[][][][][]>, List<double[][][]>, List<double[][][]>> household_policy_functions, String verbose){

        System.out.println("Computing distributions");

        // re-create state spaces for sfha, non-sfha, renters, and newborns (this helps write cleaner for loops)
        int na = a_grid.length;          // number of asset points
        int noltv = oltv_grid.length;    // number of oltv points
        int nfico = fico_grid.length;    // number of fico points
        int ntheta = theta_grid.length;  // number of disaster risk points

        // Building the state space for SFHA (dim: 5)
        Integer[][] matrix_sfha = new Integer[5][];
        matrix_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_sfha[4] = ArrayUtils.toObject(IntStream.range(0, ntheta).toArray());
        CartesianSet<Integer> state_space_sfha = new CartesianSet<>(matrix_sfha);
        int Z_sfha= (int) state_space_sfha.getCount();

        // Building the state space for non-SFHA (dim: 4)
        Integer[][] matrix_non_sfha = new Integer[4][];
        matrix_non_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_non_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_non_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_non_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        CartesianSet<Integer> state_space_non_sfha = new CartesianSet<>(matrix_non_sfha);
        int Z_non_sfha= (int) state_space_non_sfha.getCount();

        // Building the state space for rental (dim: 2)
        Integer[][] matrix_rental = new Integer[2][];
        matrix_rental[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_rental[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_rental = new CartesianSet<>(matrix_rental);
        int Z_rental = (int) state_space_rental.getCount();

        // Building the state space for newborns (dim: 2)
        Integer[][] matrix_newborn = new Integer[2][];
        matrix_newborn[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_newborn[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_newborn = new CartesianSet<>(matrix_newborn);
        int Z_newborn = (int) state_space_newborn.getCount();

        // unpacking the quartet containing all the policy functions
        List<double[][][][][][]> sfha_homeowner_policy_functions = household_policy_functions.getValue0();
        List<double[][][][][]> non_sfha_homeowner_policy_functions = household_policy_functions.getValue1();
        List<double[][][]> renter_policy_functions = household_policy_functions.getValue2();
        List<double[][][]> newborn_policy_functions = household_policy_functions.getValue3();

        // Initialize grid for value functions, policy functions
        double[][][][][][] c_sfha = sfha_homeowner_policy_functions.get(0);
        double[][][][][][] a_sfha = sfha_homeowner_policy_functions.get(1);
        double[][][][][][] oltv_sfha = sfha_homeowner_policy_functions.get(2);
        double[][][][][][] fico_sfha = sfha_homeowner_policy_functions.get(3);
        double[][][][][][] n_sfha = sfha_homeowner_policy_functions.get(4);
        double[][][][][][] action_sfha = sfha_homeowner_policy_functions.get(5);
        double[][][][][][] V_sfha = sfha_homeowner_policy_functions.get(6);

        double[][][][][] c_non_sfha = non_sfha_homeowner_policy_functions.get(0);
        double[][][][][] a_non_sfha = non_sfha_homeowner_policy_functions.get(1);
        double[][][][][] oltv_non_sfha = non_sfha_homeowner_policy_functions.get(2);
        double[][][][][] fico_non_sfha = non_sfha_homeowner_policy_functions.get(3);
        double[][][][][] n_non_sfha = non_sfha_homeowner_policy_functions.get(4);
        double[][][][][] action_non_sfha = non_sfha_homeowner_policy_functions.get(5);
        double[][][][][] V_non_sfha = non_sfha_homeowner_policy_functions.get(6);

        double[][][] c_renter = renter_policy_functions.get(0);
        double[][][] a_renter = renter_policy_functions.get(1);
        double[][][] oltv_renter = renter_policy_functions.get(2);
        double[][][] fico_renter = renter_policy_functions.get(3);
        double[][][] n_renter = renter_policy_functions.get(4);
        double[][][] action_renter = renter_policy_functions.get(5);
        double[][][] V_renter = renter_policy_functions.get(6);

        double[][][] c_newborn = newborn_policy_functions.get(0);
        double[][][] a_newborn = newborn_policy_functions.get(1);
        double[][][] oltv_newborn = newborn_policy_functions.get(2);
        double[][][] fico_newborn = newborn_policy_functions.get(3);
        double[][][] n_newborn = newborn_policy_functions.get(4);
        double[][][] action_newborn = newborn_policy_functions.get(5);
        double[][][] V_newborn = newborn_policy_functions.get(6);

        // tensors to hold the distros
        double[][][][][][] f_interim_sfha_homeowner = new double[J][na][noltv][nfico][N + 1][ntheta];
        double[][][][][] f_interim_non_sfha_homeowner = new double[J][na][noltv][nfico][N + 1];
        double[][][] f_interim_renter = new double[J][na][nfico];
        double[][][] f_newborn = new double[J][na][nfico];

        //====================================//
        // INITIALIZING THE DISTRIBUTION
        // set age at zero in preparation for the forward induction process
        int j = 0;

        //fix the initial distribution, f_0, uniform weights (for now)
        double[] asset_distribution = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/asset_distribution.csv", na);
        double[] fico_distribution = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/fico_distribution.csv", nfico);
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                f_newborn[j][ia][ifico] = asset_distribution[ia] * fico_distribution[ifico] / J;

            }
        }

        // just check that cohort weight for j = 0 adds to 1/J
        double sum = 0;
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                sum = sum + f_newborn[0][ia][ifico];
            }
        }

        if (verbose.equals("yes")){
            System.out.println(j + " , " + sum);
        }

        //====================================//

        // NOW TRANSITIONING NEWBORNS TO HOMEOWNERS AND RENTERS
        // loop for newborns -> sfha homeowners
        for (int z = 0; z < Z_newborn; z++){
            List<Integer> node = state_space_newborn.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            int ithetap = 0; // fix theta_prime = 0 i.e. no shocks
            if (action_newborn[j][ia][ifico] == 1.0){
                int iap = ArrayUtils.indexOf(a_grid, a_newborn[j][ia][ifico]);
                int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_newborn[j][ia][ifico]);
                int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[j][ia][ifico]);
                int np = (int) n_newborn[j][ia][ifico];
                f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] = f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] +
                        (1 * f_newborn[j][ia][ifico]);
            }
        }

        // loop for newborns -> non-sfha homeowners or renters (in either case there is no role for theta hence this loop is not
        // nested within the theta loop above to avoid double accounting/other such errors).
        for (int z = 0; z < Z_newborn; z++){
            List<Integer> node = state_space_newborn.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            if (action_newborn[j][ia][ifico] == 3.0){
                int iap = ArrayUtils.indexOf(a_grid, a_newborn[j][ia][ifico]);
                int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_newborn[j][ia][ifico]);
                int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[j][ia][ifico]);
                int np = (int) n_newborn[j][ia][ifico];
                f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] = f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] + f_newborn[j][ia][ifico];
            }else if (action_newborn[j][ia][ifico] == 2.0){
                int iap = ArrayUtils.indexOf(a_grid, a_newborn[j][ia][ifico]);
                int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[j][ia][ifico]);
                f_interim_renter[j + 1][iap][ificop] = f_interim_renter[j + 1][iap][ificop] + f_newborn[j][ia][ifico];
            }
        }

        // check cohort weights for j = 1
        double sum_sfha = 0; double sum_non_sfha = 0; double sum_renter = 0;
        for (int ia = 0; ia < na; ia++){
            for (int ioltv = 0; ioltv < noltv; ioltv++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int n = 0; n < N; n++){
                        for (int itheta = 0; itheta < ntheta; itheta++){
                            sum_sfha = sum_sfha + f_interim_sfha_homeowner[j + 1][ia][ioltv][ifico][n][itheta];
                        }
                    }
                }
            }
        }
        for (int ia = 0; ia < na; ia++){
            for (int ioltv = 0; ioltv < noltv; ioltv++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int n = 0; n < N; n++){
                        sum_non_sfha = sum_non_sfha + f_interim_non_sfha_homeowner[j + 1][ia][ioltv][ifico][n];
                    }
                }
            }
        }
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                for (int n = 0; n < N; n++){
                    sum_renter = sum_renter + f_interim_renter[j + 1][ia][ifico];
                }
            }
        }
        if (verbose.equals("yes")){
            System.out.println((j + 1) + " , " + sum_sfha + " , " + sum_non_sfha + " , " + sum_renter + " , " + (sum_sfha + sum_non_sfha + sum_renter));
        }
        //====================================//


        // NOW TRANSITION THE j=1 cohort to j=2 and then through to end of life, j = J.
        for (j = 1; j < (J - 1); j++){

            // loop for homeowners wanting to remain in sfha (by staying current(1.1) or refinance(1.2))
            // this loop is nested within the theta loop to allow for disaster shock (you only face disaster shock in the future if you remain in sfha)
            for (int z = 0; z < Z_sfha; z++){
                List<Integer> node = state_space_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                int itheta = node.get(4);
                int ithetap = 0;
                // giving state of an agent, x \in X, find whether it wants to stay current or refinance and if so then for that
                // agent find x' and transitioning
                if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2){
                    int iap = ArrayUtils.indexOf(a_grid, a_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int np = (int) n_sfha[j][ia][ioltv][ifico][n][itheta];
                    f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] = f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] +
                            (1 * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                }
            }

            for (int z = 0; z < Z_sfha; z++){
                List<Integer> node = state_space_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                int itheta = node.get(4);
                if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31){
                    int iap = ArrayUtils.indexOf(a_grid, a_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int np = (int) n_sfha[j][ia][ioltv][ifico][n][itheta];
                    f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] = f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] +
                            f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                }else if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                    int iap = ArrayUtils.indexOf(a_grid, a_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_sfha[j][ia][ioltv][ifico][n][itheta]);
                    f_interim_renter[j + 1][iap][ificop] = f_interim_renter[j + 1][iap][ificop] + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                }
            }

            for (int z = 0; z < Z_non_sfha; z++){
                List<Integer> node = state_space_non_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                if (action_non_sfha[j][ia][ioltv][ifico][n] == 3.1 | action_non_sfha[j][ia][ioltv][ifico][n] == 3.2){
                    int iap = ArrayUtils.indexOf(a_grid, a_non_sfha[j][ia][ioltv][ifico][n]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_non_sfha[j][ia][ioltv][ifico][n]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_non_sfha[j][ia][ioltv][ifico][n]);
                    int np = (int) n_non_sfha[j][ia][ioltv][ifico][n];
                    f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] = f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] +
                            f_interim_non_sfha_homeowner[j][ia][ioltv][ifico][n];
                }
            }

            for (int z = 0; z < Z_rental; z++){
                List<Integer> node = state_space_rental.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                if (action_renter[j][ia][ifico] == 2.1){
                    int iap = ArrayUtils.indexOf(a_grid, a_renter[j][ia][ifico]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_renter[j][ia][ifico]);
                    f_interim_renter[j + 1][iap][ificop] = f_interim_renter[j + 1][iap][ificop] + f_interim_renter[j][ia][ifico];
                }
            }

            for (int z = 0; z < Z_rental; z++){
                List<Integer> node = state_space_rental.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                int ithetap = 0;
                if (action_renter[j][ia][ifico] == 2.3){
                    int iap = ArrayUtils.indexOf(a_grid, a_renter[j][ia][ifico]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_renter[j][ia][ifico]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_renter[j][ia][ifico]);
                    int np = (int) n_renter[j][ia][ifico];
                    f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] = f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] +
                            (1 * f_interim_renter[j][ia][ifico]);

                }
            }
            // check cohort weights for j = 2
            sum_sfha = 0; sum_non_sfha = 0; sum_renter = 0;
            for (int ia = 0; ia < na; ia++){
                for (int ioltv = 0; ioltv < noltv; ioltv++){
                    for (int ifico = 0; ifico < nfico; ifico++){
                        for (int n = 0; n < N; n++){
                            for (int itheta = 0; itheta < ntheta; itheta++){
                                sum_sfha = sum_sfha + f_interim_sfha_homeowner[j + 1][ia][ioltv][ifico][n][itheta];
                            }
                        }
                    }
                }
            }
            for (int ia = 0; ia < na; ia++){
                for (int ioltv = 0; ioltv < noltv; ioltv++){
                    for (int ifico = 0; ifico < nfico; ifico++){
                        for (int n = 0; n < N; n++){
                            sum_non_sfha = sum_non_sfha + f_interim_non_sfha_homeowner[j + 1][ia][ioltv][ifico][n];
                        }
                    }
                }
            }
            for (int ia = 0; ia < na; ia++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    sum_renter = sum_renter + f_interim_renter[j + 1][ia][ifico];
                }
            }
            if (verbose.equals("yes")){
                System.out.println((j + 1) + " , " + sum_sfha + " , " + sum_non_sfha + " , " + sum_renter + " , " + (sum_sfha + sum_non_sfha + sum_renter));
            }

        }

        Quartet<double[][][][][][], double[][][][][], double[][][], double[][][]> output = new Quartet<>(f_interim_sfha_homeowner, f_interim_non_sfha_homeowner, f_interim_renter, f_newborn);
        return output;
    }

    // compute and forward distributions for the four groups: newborns, renters, non-sfha, and sfha
    Quartet<double[][][][][][], double[][][][][], double[][][], double[][][]> compute_distributions_shocks_on(
            Quartet<List<double[][][][][][]>, List<double[][][][][]>, List<double[][][]>, List<double[][][]>> household_policy_functions, String verbose){

        System.out.println("Computing distributions");

        // re-create state spaces for sfha, non-sfha, renters, and newborns (this helps write cleaner for loops)
        int na = a_grid.length;          // number of asset points
        int noltv = oltv_grid.length;    // number of oltv points
        int nfico = fico_grid.length;    // number of fico points
        int ntheta = theta_grid.length;  // number of disaster risk points

        // Building the state space for SFHA (dim: 5)
        Integer[][] matrix_sfha = new Integer[5][];
        matrix_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_sfha[4] = ArrayUtils.toObject(IntStream.range(0, ntheta).toArray());
        CartesianSet<Integer> state_space_sfha = new CartesianSet<>(matrix_sfha);
        int Z_sfha= (int) state_space_sfha.getCount();

        // Building the state space for non-SFHA (dim: 4)
        Integer[][] matrix_non_sfha = new Integer[4][];
        matrix_non_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_non_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_non_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_non_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        CartesianSet<Integer> state_space_non_sfha = new CartesianSet<>(matrix_non_sfha);
        int Z_non_sfha= (int) state_space_non_sfha.getCount();

        // Building the state space for rental (dim: 2)
        Integer[][] matrix_rental = new Integer[2][];
        matrix_rental[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_rental[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_rental = new CartesianSet<>(matrix_rental);
        int Z_rental = (int) state_space_rental.getCount();

        // Building the state space for newborns (dim: 2)
        Integer[][] matrix_newborn = new Integer[2][];
        matrix_newborn[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_newborn[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_newborn = new CartesianSet<>(matrix_newborn);
        int Z_newborn = (int) state_space_newborn.getCount();

        // unpacking the quartet containing all the policy functions
        List<double[][][][][][]> sfha_homeowner_policy_functions = household_policy_functions.getValue0();
        List<double[][][][][]> non_sfha_homeowner_policy_functions = household_policy_functions.getValue1();
        List<double[][][]> renter_policy_functions = household_policy_functions.getValue2();
        List<double[][][]> newborn_policy_functions = household_policy_functions.getValue3();

        // Initialize grid for value functions, policy functions
        double[][][][][][] c_sfha = sfha_homeowner_policy_functions.get(0);
        double[][][][][][] a_sfha = sfha_homeowner_policy_functions.get(1);
        double[][][][][][] oltv_sfha = sfha_homeowner_policy_functions.get(2);
        double[][][][][][] fico_sfha = sfha_homeowner_policy_functions.get(3);
        double[][][][][][] n_sfha = sfha_homeowner_policy_functions.get(4);
        double[][][][][][] action_sfha = sfha_homeowner_policy_functions.get(5);
        double[][][][][][] V_sfha = sfha_homeowner_policy_functions.get(6);

        double[][][][][] c_non_sfha = non_sfha_homeowner_policy_functions.get(0);
        double[][][][][] a_non_sfha = non_sfha_homeowner_policy_functions.get(1);
        double[][][][][] oltv_non_sfha = non_sfha_homeowner_policy_functions.get(2);
        double[][][][][] fico_non_sfha = non_sfha_homeowner_policy_functions.get(3);
        double[][][][][] n_non_sfha = non_sfha_homeowner_policy_functions.get(4);
        double[][][][][] action_non_sfha = non_sfha_homeowner_policy_functions.get(5);
        double[][][][][] V_non_sfha = non_sfha_homeowner_policy_functions.get(6);

        double[][][] c_renter = renter_policy_functions.get(0);
        double[][][] a_renter = renter_policy_functions.get(1);
        double[][][] oltv_renter = renter_policy_functions.get(2);
        double[][][] fico_renter = renter_policy_functions.get(3);
        double[][][] n_renter = renter_policy_functions.get(4);
        double[][][] action_renter = renter_policy_functions.get(5);
        double[][][] V_renter = renter_policy_functions.get(6);

        double[][][] c_newborn = newborn_policy_functions.get(0);
        double[][][] a_newborn = newborn_policy_functions.get(1);
        double[][][] oltv_newborn = newborn_policy_functions.get(2);
        double[][][] fico_newborn = newborn_policy_functions.get(3);
        double[][][] n_newborn = newborn_policy_functions.get(4);
        double[][][] action_newborn = newborn_policy_functions.get(5);
        double[][][] V_newborn = newborn_policy_functions.get(6);

        // tensors to hold the distros
        double[][][][][][] f_interim_sfha_homeowner = new double[J][na][noltv][nfico][N + 1][ntheta];
        double[][][][][] f_interim_non_sfha_homeowner = new double[J][na][noltv][nfico][N + 1];
        double[][][] f_interim_renter = new double[J][na][nfico];
        double[][][] f_newborn = new double[J][na][nfico];

        //====================================//
        // INITIALIZING THE DISTRIBUTION
        // set age at zero in preparation for the forward induction process
        int j = 0;

        //fix the initial distribution, f_0, uniform weights (for now)
        double[] asset_distribution = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/asset_distribution.csv", na);
        double[] fico_distribution = miscFunctions.read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/fico_distribution.csv", nfico);
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                f_newborn[j][ia][ifico] = asset_distribution[ia] * fico_distribution[ifico] / J;

            }
        }

        // just check that cohort weight for j = 0 adds to 1/J
        double sum = 0;
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                sum = sum + f_newborn[0][ia][ifico];
            }
        }

        if (verbose.equals("yes")){
            System.out.println(j + " , " + sum);
        }

        //====================================//

        // NOW TRANSITIONING NEWBORNS TO HOMEOWNERS AND RENTERS
        // loop for newborns -> sfha homeowners
        for (int ithetap = 0; ithetap < ntheta; ithetap++){
            for (int z = 0; z < Z_newborn; z++){
                List<Integer> node = state_space_newborn.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                if (action_newborn[j][ia][ifico] == 1.0){
                    int iap = ArrayUtils.indexOf(a_grid, a_newborn[j][ia][ifico]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_newborn[j][ia][ifico]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[j][ia][ifico]);
                    int np = (int) n_newborn[j][ia][ifico];
                    f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] = f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] +
                            (Pi[0][ithetap] * f_newborn[j][ia][ifico]);
                }
            }

        }


        // loop for newborns -> non-sfha homeowners or renters (in either case there is no role for theta hence this loop is not
        // nested within the theta loop above to avoid double accounting/other such errors).
        for (int z = 0; z < Z_newborn; z++){
            List<Integer> node = state_space_newborn.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            if (action_newborn[j][ia][ifico] == 3.0){
                int iap = ArrayUtils.indexOf(a_grid, a_newborn[j][ia][ifico]);
                int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_newborn[j][ia][ifico]);
                int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[j][ia][ifico]);
                int np = (int) n_newborn[j][ia][ifico];
                f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] = f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] + f_newborn[j][ia][ifico];
            }else if (action_newborn[j][ia][ifico] == 2.0){
                int iap = ArrayUtils.indexOf(a_grid, a_newborn[j][ia][ifico]);
                int ificop = ArrayUtils.indexOf(fico_grid, fico_newborn[j][ia][ifico]);
                f_interim_renter[j + 1][iap][ificop] = f_interim_renter[j + 1][iap][ificop] + f_newborn[j][ia][ifico];
            }
        }

        // check cohort weights for j = 1
        double sum_sfha = 0; double sum_non_sfha = 0; double sum_renter = 0;
        for (int ia = 0; ia < na; ia++){
            for (int ioltv = 0; ioltv < noltv; ioltv++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int n = 0; n < N; n++){
                        for (int itheta = 0; itheta < ntheta; itheta++){
                            sum_sfha = sum_sfha + f_interim_sfha_homeowner[j + 1][ia][ioltv][ifico][n][itheta];
                        }
                    }
                }
            }
        }
        for (int ia = 0; ia < na; ia++){
            for (int ioltv = 0; ioltv < noltv; ioltv++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int n = 0; n < N; n++){
                        sum_non_sfha = sum_non_sfha + f_interim_non_sfha_homeowner[j + 1][ia][ioltv][ifico][n];
                    }
                }
            }
        }
        for (int ia = 0; ia < na; ia++){
            for (int ifico = 0; ifico < nfico; ifico++){
                for (int n = 0; n < N; n++){
                    sum_renter = sum_renter + f_interim_renter[j + 1][ia][ifico];
                }
            }
        }
        if (verbose.equals("yes")){
            System.out.println((j + 1) + " , " + sum_sfha + " , " + sum_non_sfha + " , " + sum_renter + " , " + (sum_sfha + sum_non_sfha + sum_renter));
        }
        //====================================//


        // NOW TRANSITION THE j=1 cohort to j=2 and then through to end of life, j = J.
        for (j = 1; j < (J - 1); j++){

            // loop for homeowners wanting to remain in sfha (by staying current(1.1) or refinance(1.2))
            // this loop is nested within the theta loop to allow for disaster shock (you only face disaster shock in the future if you remain in sfha)
            for (int ithetap = 0; ithetap < ntheta; ithetap++){
                for (int z = 0; z < Z_sfha; z++){
                    List<Integer> node = state_space_sfha.get(z);
                    int ia = node.get(0);
                    int ioltv = node.get(1);
                    int ifico = node.get(2);
                    int n = node.get(3);
                    int itheta = node.get(4);
                    // giving state of an agent, x \in X, find whether it wants to stay current or refinance and if so then for that
                    // agent find x' and transitioning
                    if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2){
                        int iap = ArrayUtils.indexOf(a_grid, a_sfha[j][ia][ioltv][ifico][n][itheta]);
                        int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_sfha[j][ia][ioltv][ifico][n][itheta]);
                        int ificop = ArrayUtils.indexOf(fico_grid, fico_sfha[j][ia][ioltv][ifico][n][itheta]);
                        int np = (int) n_sfha[j][ia][ioltv][ifico][n][itheta];
                        f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] = f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] +
                                (Pi[itheta][ithetap] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                    }
                }
            }

            for (int z = 0; z < Z_sfha; z++){
                List<Integer> node = state_space_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                int itheta = node.get(4);
                if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31){
                    int iap = ArrayUtils.indexOf(a_grid, a_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int np = (int) n_sfha[j][ia][ioltv][ifico][n][itheta];
                    f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] = f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] +
                            f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                }else if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                    int iap = ArrayUtils.indexOf(a_grid, a_sfha[j][ia][ioltv][ifico][n][itheta]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_sfha[j][ia][ioltv][ifico][n][itheta]);
                    f_interim_renter[j + 1][iap][ificop] = f_interim_renter[j + 1][iap][ificop] + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                }
            }

            for (int z = 0; z < Z_non_sfha; z++){
                List<Integer> node = state_space_non_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                if (action_non_sfha[j][ia][ioltv][ifico][n] == 3.1 | action_non_sfha[j][ia][ioltv][ifico][n] == 3.2){
                    int iap = ArrayUtils.indexOf(a_grid, a_non_sfha[j][ia][ioltv][ifico][n]);
                    int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_non_sfha[j][ia][ioltv][ifico][n]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_non_sfha[j][ia][ioltv][ifico][n]);
                    int np = (int) n_non_sfha[j][ia][ioltv][ifico][n];
                    f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] = f_interim_non_sfha_homeowner[j + 1][iap][ioltvp][ificop][np] +
                            f_interim_non_sfha_homeowner[j][ia][ioltv][ifico][n];
                }
            }

            for (int z = 0; z < Z_rental; z++){
                List<Integer> node = state_space_rental.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                if (action_renter[j][ia][ifico] == 2.1){
                    int iap = ArrayUtils.indexOf(a_grid, a_renter[j][ia][ifico]);
                    int ificop = ArrayUtils.indexOf(fico_grid, fico_renter[j][ia][ifico]);
                    f_interim_renter[j + 1][iap][ificop] = f_interim_renter[j + 1][iap][ificop] + f_interim_renter[j][ia][ifico];
                }
            }

            for (int ithetap = 0; ithetap < ntheta; ithetap++){
                for (int z = 0; z < Z_rental; z++){
                    List<Integer> node = state_space_rental.get(z);
                    int ia = node.get(0);
                    int ifico = node.get(1);
                    if (action_renter[j][ia][ifico] == 2.3){
                        int iap = ArrayUtils.indexOf(a_grid, a_renter[j][ia][ifico]);
                        int ioltvp = ArrayUtils.indexOf(oltv_grid, oltv_renter[j][ia][ifico]);
                        int ificop = ArrayUtils.indexOf(fico_grid, fico_renter[j][ia][ifico]);
                        int np = (int) n_renter[j][ia][ifico];
                        f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] = f_interim_sfha_homeowner[j + 1][iap][ioltvp][ificop][np][ithetap] +
                                (Pi[0][ithetap] * f_interim_renter[j][ia][ifico]);

                    }
                }
            }


            // check cohort weights for j = 2
            sum_sfha = 0; sum_non_sfha = 0; sum_renter = 0;
            for (int ia = 0; ia < na; ia++){
                for (int ioltv = 0; ioltv < noltv; ioltv++){
                    for (int ifico = 0; ifico < nfico; ifico++){
                        for (int n = 0; n < N; n++){
                            for (int itheta = 0; itheta < ntheta; itheta++){
                                sum_sfha = sum_sfha + f_interim_sfha_homeowner[j + 1][ia][ioltv][ifico][n][itheta];
                            }
                        }
                    }
                }
            }
            for (int ia = 0; ia < na; ia++){
                for (int ioltv = 0; ioltv < noltv; ioltv++){
                    for (int ifico = 0; ifico < nfico; ifico++){
                        for (int n = 0; n < N; n++){
                            sum_non_sfha = sum_non_sfha + f_interim_non_sfha_homeowner[j + 1][ia][ioltv][ifico][n];
                        }
                    }
                }
            }
            for (int ia = 0; ia < na; ia++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    sum_renter = sum_renter + f_interim_renter[j + 1][ia][ifico];
                }
            }
            if (verbose.equals("yes")){
                System.out.println((j + 1) + " , " + sum_sfha + " , " + sum_non_sfha + " , " + sum_renter + " , " + (sum_sfha + sum_non_sfha + sum_renter));
            }

        }

        Quartet<double[][][][][][], double[][][][][], double[][][], double[][][]> output = new Quartet<>(f_interim_sfha_homeowner, f_interim_non_sfha_homeowner, f_interim_renter, f_newborn);
        return output;
    }

    double[] compute_market_clearing_and_moments(
            Quartet<List<double[][][][][][]>, List<double[][][][][]>, List<double[][][]>, List<double[][][]>> household_policy_functions,
            Quartet<double[][][][][][], double[][][][][], double[][][], double[][][]> distributions, double price_sfha,
            String compute_moments, String verbose){

        // cohort weight
        double mu_j = 1d / J;

        // re-create state spaces for sfha, non-sfha, renters, and newborns (this helps write cleaner for loops)
        int na = a_grid.length;          // number of asset points
        int noltv = oltv_grid.length;    // number of oltv points
        int nfico = fico_grid.length;    // number of fico points
        int ntheta = theta_grid.length;  // number of disaster risk points

        // Building the state space for SFHA (dim: 5)
        Integer[][] matrix_sfha = new Integer[5][];
        matrix_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        matrix_sfha[4] = ArrayUtils.toObject(IntStream.range(0, ntheta).toArray());
        CartesianSet<Integer> state_space_sfha = new CartesianSet<>(matrix_sfha);
        int Z_sfha= (int) state_space_sfha.getCount();

        // Building the state space for non-SFHA (dim: 4)
        Integer[][] matrix_non_sfha = new Integer[4][];
        matrix_non_sfha[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_non_sfha[1] = ArrayUtils.toObject(IntStream.range(0, noltv).toArray());
        matrix_non_sfha[2] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        matrix_non_sfha[3] = ArrayUtils.toObject(IntStream.range(0, N).toArray());
        CartesianSet<Integer> state_space_non_sfha = new CartesianSet<>(matrix_non_sfha);
        int Z_non_sfha= (int) state_space_non_sfha.getCount();

        // Building the state space for rental (dim: 2)
        Integer[][] matrix_rental = new Integer[2][];
        matrix_rental[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_rental[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_rental = new CartesianSet<>(matrix_rental);
        int Z_rental = (int) state_space_rental.getCount();

        // Building the state space for newborns (dim: 2)
        Integer[][] matrix_newborn = new Integer[2][];
        matrix_newborn[0] = ArrayUtils.toObject(IntStream.range(0, na).toArray());
        matrix_newborn[1] = ArrayUtils.toObject(IntStream.range(0, nfico).toArray());
        CartesianSet<Integer> state_space_newborn = new CartesianSet<>(matrix_newborn);
        int Z_newborn = (int) state_space_newborn.getCount();

        // unpacking the quartet containing all the policy functions
        List<double[][][][][][]> sfha_homeowner_policy_functions = household_policy_functions.getValue0();
        List<double[][][][][]> non_sfha_homeowner_policy_functions = household_policy_functions.getValue1();
        List<double[][][]> renter_policy_functions = household_policy_functions.getValue2();
        List<double[][][]> newborn_policy_functions = household_policy_functions.getValue3();

        // unpacking each list to get policy functions for each state variable
        double[][][][][][] c_sfha = sfha_homeowner_policy_functions.get(0);
        double[][][][][][] a_sfha = sfha_homeowner_policy_functions.get(1);
        double[][][][][][] oltv_sfha = sfha_homeowner_policy_functions.get(2);
        double[][][][][][] fico_sfha = sfha_homeowner_policy_functions.get(3);
        double[][][][][][] n_sfha = sfha_homeowner_policy_functions.get(4);
        double[][][][][][] action_sfha = sfha_homeowner_policy_functions.get(5);
        double[][][][][][] V_sfha = sfha_homeowner_policy_functions.get(6);

        double[][][][][] c_non_sfha = non_sfha_homeowner_policy_functions.get(0);
        double[][][][][] a_non_sfha = non_sfha_homeowner_policy_functions.get(1);
        double[][][][][] oltv_non_sfha = non_sfha_homeowner_policy_functions.get(2);
        double[][][][][] fico_non_sfha = non_sfha_homeowner_policy_functions.get(3);
        double[][][][][] n_non_sfha = non_sfha_homeowner_policy_functions.get(4);
        double[][][][][] action_non_sfha = non_sfha_homeowner_policy_functions.get(5);
        double[][][][][] V_non_sfha = non_sfha_homeowner_policy_functions.get(6);

        double[][][] c_renter = renter_policy_functions.get(0);
        double[][][] a_renter = renter_policy_functions.get(1);
        double[][][] oltv_renter = renter_policy_functions.get(2);
        double[][][] fico_renter = renter_policy_functions.get(3);
        double[][][] n_renter = renter_policy_functions.get(4);
        double[][][] action_renter = renter_policy_functions.get(5);
        double[][][] V_renter = renter_policy_functions.get(6);

        double[][][] c_newborn = newborn_policy_functions.get(0);
        double[][][] a_newborn = newborn_policy_functions.get(1);
        double[][][] oltv_newborn = newborn_policy_functions.get(2);
        double[][][] fico_newborn = newborn_policy_functions.get(3);
        double[][][] n_newborn = newborn_policy_functions.get(4);
        double[][][] action_newborn = newborn_policy_functions.get(5);
        double[][][] V_newborn = newborn_policy_functions.get(6);

        // unpacking the quartet to get the distros
        double[][][][][][] f_interim_sfha_homeowner = distributions.getValue0();
        double[][][][][] f_interim_non_sfha_homeowner = distributions.getValue1();
        double[][][] f_interim_renter = distributions.getValue2();
        double[][][] f_newborn = distributions.getValue3();

        // calculate the aggregate housing demand
        double aggregate_housing_demand_newborns = 0;
        double aggregate_housing_demand_renters = 0;
        double aggregate_housing_supply = 0;
        double aggregate_housing_excess_demand = 0;

        int j = 0;
        for (int z = 0; z < Z_newborn; z++){
            List<Integer> node = state_space_newborn.get(z);
            int ia = node.get(0);
            int ifico = node.get(1);
            // integrate to find the housing demand from newborns
            if (action_newborn[j][ia][ifico] == 1.0){
                aggregate_housing_demand_newborns = aggregate_housing_demand_newborns + (1 * f_newborn[j][ia][ifico]);
            }
        }


        if (verbose.equals("yes")){
            //System.out.println("age , renter , homeowner");
        }

        for (j = 1; j < J; j++){

            for (int z = 0; z < Z_rental; z++){
                List<Integer> node = state_space_rental.get(z);
                int ia = node.get(0);
                int ifico = node.get(1);
                if (action_renter[j][ia][ifico] == 2.3){
                    aggregate_housing_demand_renters = aggregate_housing_demand_renters + (1 * f_interim_renter[j][ia][ifico]);
                    if (f_interim_renter[j][ia][ifico] >= 1){
                        System.out.println("Alert: pdf > 1 (homeonwer/supply)");
                    }
                }
            }

            for (int z = 0; z < Z_sfha; z++){
                List<Integer> node = state_space_sfha.get(z);
                int ia = node.get(0);
                int ioltv = node.get(1);
                int ifico = node.get(2);
                int n = node.get(3);
                int itheta = node.get(4);
                if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31 |
                        action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32 |
                        action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4
                    //        | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.5 // supply from bequeathers (take it out and sub with construction sector
                ){
                    aggregate_housing_supply = aggregate_housing_supply + (1 * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                    if (f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta] >= 1){
                        System.out.println("Alert: pdf > 1 (homeonwer/supply)");
                    }
                }
            }
            if (verbose.equals("yes")){
                //System.out.println(j + " , " + aggregate_housing_demand_renters + " , " + aggregate_housing_supply);
            }

        }

        double housing_supply_construction_sector = construction_sector(price_sfha);
        aggregate_housing_excess_demand = (aggregate_housing_demand_newborns + aggregate_housing_demand_renters)
                - (aggregate_housing_supply + housing_supply_construction_sector);

        System.out.println("==================================");
        System.out.println("housing demand (newborns): " + (aggregate_housing_demand_newborns));
        System.out.println("housing demand (renters): " + (aggregate_housing_demand_renters));
        System.out.println("housing supply (sfha homeowners): " + (aggregate_housing_supply));
        System.out.println("housing supply (construction sector):" + (housing_supply_construction_sector));
        System.out.println("aggregate housing excess demand: " + aggregate_housing_excess_demand);
        System.out.println("==================================");

        if (aggregate_housing_excess_demand > 1e-5){
            System.out.println("Housing market not cleared. Excess demand: " + aggregate_housing_excess_demand);
        }

        double[] output = new double[4];
        if (compute_moments.equals("yes")){

            // CALCULATE MOMENTS OF INTEREST IF INSTRUCTED BY USER
            double defaulters = 0;
            double remainers = 0;
            double refinancers = 0;
            double sfha_homeowners = 0;
            double non_sfha_homeowners = 0;
            double[] consumption_by_age_sfha = new double[J];
            double[] consumption_by_age_non_sfha = new double[J];
            double[] consumption_by_age_rental = new double[J];
            double[] consumption_by_age = new double[J];

            double defaulters_control = 0;
            double defaulters_treatment = 0;
            double remainers_control = 0;
            double remainers_treatment = 0;
            double sfha_homeowners_control = 0;
            double sfha_homeowners_treatment = 0;
            double refinancers_control = 0;
            double refinancers_treatment = 0;
            double movers_rent_control = 0;
            double movers_safe_control = 0;
            double movers_rent_treatment = 0;
            double movers_safe_treatment = 0;

            double welfare_defaulters_control = 0;
            double welfare_defaulters_treatment = 0;
            double welfare_remainers_control = 0;
            double welfare_remainers_treatment = 0;
            double welfare_sfha_homeowners_control = 0;
            double welfare_sfha_homeowners_treatment = 0;
            double welfare_refinancers_control = 0;
            double welfare_refinancers_treatment = 0;
            double welfare_movers_rent_control = 0;
            double welfare_movers_rent_treatment = 0;
            double welfare_movers_safe_control = 0;
            double welfare_movers_safe_treatment = 0;


            for (j = 0; j < J; j++){

                for (int z = 0; z < Z_sfha; z++){
                    List<Integer> node = state_space_sfha.get(z);
                    int ia = node.get(0);
                    int ioltv = node.get(1);
                    int ifico = node.get(2);
                    int n = node.get(3);
                    int itheta = node.get(4);

                    // calculate moments for SFHA overall
                    if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1){
                        remainers = remainers + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                    }
                    if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2){
                        refinancers = refinancers + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                    }
                    if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                        defaulters = defaulters + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                    }
                    if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.0 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1 |
                            action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31 |
                            action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                        sfha_homeowners = sfha_homeowners + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                    }
                    consumption_by_age_sfha[j] = consumption_by_age_sfha[j] + (c_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);

                    // calculate moments for control group within SFHA
                    if (theta_grid[itheta] == theta_grid[0]){
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1){
                            remainers_control = remainers_control +  f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_remainers_control = welfare_remainers_control
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2){
                            refinancers_control = refinancers_control + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_refinancers_control = welfare_refinancers_control
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31){
                            movers_rent_control = movers_rent_control + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_movers_rent_control = welfare_movers_rent_control
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32){
                            movers_safe_control = movers_safe_control + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_movers_safe_control = welfare_movers_safe_control
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                            defaulters_control = defaulters_control + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_defaulters_control = welfare_defaulters_control
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.0 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1 |
                                action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31 |
                                action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                            sfha_homeowners_control = sfha_homeowners_control + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_sfha_homeowners_control = welfare_sfha_homeowners_control
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                    }
                    // calculate moments for treatment group within SFHA
                    else if (theta_grid[itheta] == theta_grid[1]){
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1){
                            remainers_treatment = remainers_treatment +  f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_remainers_treatment = welfare_remainers_treatment
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2){
                            refinancers_treatment = refinancers_treatment + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_refinancers_treatment = welfare_refinancers_treatment
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31){
                            movers_rent_treatment = movers_rent_treatment + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_movers_rent_treatment = welfare_movers_rent_treatment
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32){
                            movers_safe_treatment = movers_safe_treatment + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_movers_safe_treatment = welfare_movers_safe_treatment
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                            defaulters_treatment = defaulters_treatment + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_defaulters_treatment = welfare_defaulters_treatment
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                        if (action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.0 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.1 |
                                action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.2 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.31 |
                                action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.32 | action_sfha[j][ia][ioltv][ifico][n][itheta] == 1.4){
                            sfha_homeowners_treatment = sfha_homeowners_treatment + f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta];
                            welfare_sfha_homeowners_treatment = welfare_sfha_homeowners_treatment
                                    + (V_sfha[j][ia][ioltv][ifico][n][itheta] * f_interim_sfha_homeowner[j][ia][ioltv][ifico][n][itheta]);
                        }
                    }

                }

                for (int z = 0; z < Z_non_sfha; z++){
                    List<Integer> node = state_space_non_sfha.get(z);
                    int ia = node.get(0);
                    int ioltv = node.get(1);
                    int ifico = node.get(2);
                    int n = node.get(3);
                    if (action_non_sfha[j][ia][ioltv][ifico][n] == 3.0 | action_non_sfha[j][ia][ioltv][ifico][n] == 3.1 |
                            action_non_sfha[j][ia][ioltv][ifico][n] == 3.2){
                        non_sfha_homeowners = non_sfha_homeowners + f_interim_non_sfha_homeowner[j][ia][ioltv][ifico][n];
                    }
                    consumption_by_age_non_sfha[j] = consumption_by_age_non_sfha[j] + (c_non_sfha[j][ia][ioltv][ifico][n] * f_interim_non_sfha_homeowner[j][ia][ioltv][ifico][n]);
                }

                for (int z = 0; z < Z_rental; z++){
                    List<Integer> node = state_space_rental.get(z);
                    int ia = node.get(0);
                    int ifico = node.get(1);
                    consumption_by_age_rental[j] = consumption_by_age_rental[j] + (c_renter[j][ia][ifico] * f_interim_renter[j][ia][ifico]);
                }

                consumption_by_age[j] = consumption_by_age_sfha[j] + consumption_by_age_non_sfha[j] + consumption_by_age_rental[j];
            }

            double homeownership_rate = sfha_homeowners + non_sfha_homeowners;
            double homeownership_rate_ratio = non_sfha_homeowners / sfha_homeowners;
            double price_ratio = price_non_sfha / price_sfha;
            double default_rate = (defaulters / sfha_homeowners);
            double current_rate = (remainers / sfha_homeowners);
            double refinance_rate = (refinancers / sfha_homeowners);
            output[0] = homeownership_rate;
            output[1] = homeownership_rate_ratio;
            output[2] = price_ratio;
            output[3] = default_rate;

            if (verbose.equals("yes")){
                System.out.println("defaulters: " + defaulters);
                System.out.println("sfha homeowners: " + sfha_homeowners);
                System.out.println("non-sfha homeowners: " + non_sfha_homeowners);

                System.out.println("homeownership rate: " + homeownership_rate);
                System.out.println("homeownership ratio (non-sfha / sfha): " + homeownership_rate_ratio);
                System.out.println("default rate: " + default_rate);
                System.out.println("current rate: " + current_rate);
                System.out.println("refinance rate: " + refinance_rate);
                System.out.println("consumption by age (sfha): " + Arrays.toString(consumption_by_age_sfha));
                System.out.println("consumption by age (non-sfha): " + Arrays.toString(consumption_by_age_non_sfha));
                System.out.println("consumption by age (renter): " + Arrays.toString(consumption_by_age_rental));
                System.out.println("consumption by age: (overall):" + Arrays.toString(consumption_by_age));
                System.out.println("==================================");

                System.out.println("##################################");
                System.out.println("control mass: " + sfha_homeowners_control);
                System.out.println("control current rate: " + (remainers_control / sfha_homeowners_control));
                System.out.println("control refinance rate: " + (refinancers_control / sfha_homeowners_control));
                System.out.println("control moving rate (rent): " + (movers_rent_control / sfha_homeowners_control));
                System.out.println("control moving rate (safe): " + (movers_safe_control / sfha_homeowners_control));
                System.out.println("control default rate: " + (defaulters_control / sfha_homeowners_control));

                System.out.println("treatment mass: " + sfha_homeowners_treatment);
                System.out.println("treatment current rate: " + (remainers_treatment / sfha_homeowners_treatment));
                System.out.println("treatment refinance rate: " + (refinancers_treatment / sfha_homeowners_treatment));
                System.out.println("treatment moving rate (rent): " + (movers_rent_treatment / sfha_homeowners_treatment));
                System.out.println("treatment moving rate (safe): " + (movers_safe_treatment / sfha_homeowners_treatment));
                System.out.println("treatment default rate: " + (defaulters_treatment / sfha_homeowners_treatment));

                System.out.println("!diff in current rate: " + ((remainers_treatment / sfha_homeowners_treatment) - (remainers_control / sfha_homeowners_control)) * 100);
                System.out.println("!diff in refi rate   : " + ((refinancers_treatment / sfha_homeowners_treatment) - (refinancers_control / sfha_homeowners_treatment)) * 100);
                System.out.println("##################################");

                System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");
                System.out.println("Welfare values");
                System.out.println("sfha control (overall): " + welfare_sfha_homeowners_control);
                System.out.println("remain control: " + welfare_remainers_control);
                System.out.println("refi control: " + welfare_refinancers_control);
                System.out.println("move (rent) control: " + welfare_movers_rent_control);
                System.out.println("move (safe) control: " + welfare_movers_safe_control);
                System.out.println("default control: " + welfare_defaulters_control);
                System.out.println("\n");
                System.out.println("sfha treatment (overall): " + welfare_sfha_homeowners_treatment);
                System.out.println("remain treatment: " + welfare_remainers_treatment);
                System.out.println("refi treatment: " + welfare_refinancers_treatment);
                System.out.println("move (rent) treatment: " + welfare_movers_rent_treatment);
                System.out.println("move (safe) treatment: " + welfare_movers_safe_treatment);
                System.out.println("default treatment: " + welfare_defaulters_treatment);
                System.out.println("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%");

                System.out.println("##################################");
                System.out.println("control default rate: " + (defaulters_control / sfha_homeowners_control));
                System.out.println("treatment default rate: " + (defaulters_treatment / sfha_homeowners_treatment));
                System.out.println("##################################");

            }


        }
        return output;
    }

    double[] calibration(double price_sfha, double[] empirical_moments, String verbose){

        double[] price_ratio_grid = {0.8}; //0.80
        double[] alpha_sfha_grid = {0.42};  //2.2, 0.52
        double[] alpha_non_sfha_grid = {0.2}; //1.9, 0.3
        double[] default_penalty_grid = {0.0}; //0.0707
        double[] transfer_grid = {1.22}; //1.22
        additional_moving_cost = 0.5;

        double SS = 1e5; double AAS = 0; double AANS = 0; double PR = 0; double DD = 0; double TT = 0;

        // Building the state space for calibration (dim: 5)
        Integer[][] matrix_calibration = new Integer[5][];
        matrix_calibration[0] = ArrayUtils.toObject(IntStream.range(0, price_ratio_grid.length).toArray());
        matrix_calibration[1] = ArrayUtils.toObject(IntStream.range(0, alpha_sfha_grid.length).toArray());
        matrix_calibration[2] = ArrayUtils.toObject(IntStream.range(0, alpha_non_sfha_grid.length).toArray());
        matrix_calibration[3] = ArrayUtils.toObject(IntStream.range(0, default_penalty_grid.length).toArray());
        matrix_calibration[4] = ArrayUtils.toObject(IntStream.range(0, transfer_grid.length).toArray());
        CartesianSet<Integer> state_space_calibration = new CartesianSet<>(matrix_calibration);
        int Z_calibration = (int) state_space_calibration.getCount();
        System.out.println("number of simulations: " + Z_calibration);

        Quartet<List<double[][][][][][]>, List<double[][][][][]>, List<double[][][]>, List<double[][][]>> household_policy_functions;
        Quartet<double[][][][][][], double[][][][][], double[][][], double[][][]> distributions_shocks_off;
        double[] model_moments_shocks_off = new double[empirical_moments.length];

        Quartet<double[][][][][][], double[][][][][], double[][][], double[][][]> distributions_shocks_on;
        double[] model_moments_shocks_on = new double[empirical_moments.length];

        for (int z = 0; z < Z_calibration; z++){
            List<Integer> node = state_space_calibration.get(z);
            price_non_sfha = price_ratio_grid[node.get(0)] * price_sfha;
            alpha_sfha = alpha_sfha_grid[node.get(1)];
            alpha_non_sfha = alpha_non_sfha_grid[node.get(2)];
            d = default_penalty_grid[node.get(3)];
            tau = transfer_grid[node.get(4)];
            household_policy_functions = solve_household_problems(price_sfha, price_non_sfha, "no");
            distributions_shocks_off = compute_distributions_shocks_off(household_policy_functions, "no");
            model_moments_shocks_off = compute_market_clearing_and_moments(household_policy_functions, distributions_shocks_off, price_sfha,"yes", "yes");

            double sse = 0;
            for (int i = 0; i < model_moments_shocks_off.length; i++){
                sse = sse + pow(model_moments_shocks_off[i] - empirical_moments[i], 2);
            }
            System.out.println("==================================");
            System.out.println("moments of interest (model / empirical)");
            System.out.println("homeownership rate: " + model_moments_shocks_off[0] + " / "  + empirical_moments[0]);
            System.out.println("homeownership ratio: " + model_moments_shocks_off[1] + " / "  +  + empirical_moments[1]);
            System.out.println("price ratio: " + model_moments_shocks_off[2] + " / "  + empirical_moments[2]);
            System.out.println("foreclosure rate: : " + model_moments_shocks_off[3] + " / "  + empirical_moments[3]);
            System.out.println("result: [" + sse + " , " + alpha_sfha + " , " + alpha_non_sfha + " , " + (price_non_sfha / price_sfha) + " , " + d + "]");
            System.out.println("==================================");

            System.out.println("## Turning Shocks On ##");
            distributions_shocks_on = compute_distributions_shocks_on(household_policy_functions, "no");
            model_moments_shocks_on = compute_market_clearing_and_moments(household_policy_functions, distributions_shocks_on, price_sfha,"yes", "yes");

            System.out.println("==================================");
            System.out.println("moments of interest (model / empirical)");
            System.out.println("homeownership rate: " + model_moments_shocks_on[0] + " / "  + empirical_moments[0]);
            System.out.println("homeownership ratio: " + model_moments_shocks_on[1] + " / "  +  + empirical_moments[1]);
            System.out.println("price ratio: " + model_moments_shocks_on[2] + " / "  + empirical_moments[2]);
            System.out.println("foreclosure rate: : " + model_moments_shocks_on[3] + " / "  + empirical_moments[3]);
            System.out.println("result: [" + sse + " , " + alpha_sfha + " , " + alpha_non_sfha + " , " + (price_non_sfha / price_sfha) + " , " + d + " , " + tau +  "]");
            System.out.println("==================================");
            System.out.println("++++++++++++++++++++++++++++++++++");
            System.out.println("++++++++++++++++++++++++++++++++++");
            System.out.println("++++++++++++++++++++++++++++++++++");
            System.out.println("\n");
            System.out.println("\n");
            System.out.println("\n");
            System.out.println("\n");


            if (sse <= SS){
                SS = sse;
                AAS = alpha_sfha;
                AANS = alpha_non_sfha;
                PR = price_non_sfha / price_sfha;
                DD = d;
                TT = tau;
            }


        }



        double[] output = {SS, AAS, AANS, PR, DD, TT};
        return output;
    }

    public static void main (String[] args){

        long start_time = System.nanoTime();

        ExperimentTwo model = new ExperimentTwo();

        double[] empirical_moments = new MiscFunctions().read_vector("/home/ahyan/Dropbox/Underwater/DataWork/calibration/moments_of_interest.csv", 4);
        System.out.println(Arrays.toString(model.calibration(7.23, empirical_moments, "no")));


        long end_time = System.nanoTime();
        System.out.println("run time = " + (end_time - start_time) * 1e-9 / 60 + " mins");

    }

}
