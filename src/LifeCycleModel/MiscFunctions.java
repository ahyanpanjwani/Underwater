package LifeCycleModel;


import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.interpolation.UnivariateInterpolator;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class MiscFunctions {

    double[] linspace(int nx, double xmin, double xmax){

        double[] agrid = new double[nx];
        agrid[0] = xmin;
        double astep = (xmax - xmin) / (nx - 1);

        for (int i = 1; i < agrid.length; i++){
            agrid[i] = agrid[i - 1] + astep;
        }

        return agrid;
    }

    // matrix for the full credit surface matrix
    double[][] credit_surface(String filePathForCS, double[] fico_grid, double[] oltv_grid_interpolate){

        // this is the oltv grid used in R to estimate the empirical credit surface; it is different from the oltv grid which
        // is part of the agents' state space. oltv_grid_data is used for interpolation over a finer grid, oltv_grid_interpolate
        // and the latter is used in the main model (IdiosyncraticRiskWithSpace.java).
        double[] oltv_grid_data = {0.50, 0.60, 0.70, 0.80, 0.90, 1.00};

        // first read the empirical credit surface developed by R. The R scrip that produces the csv that is being read here is stored in
        // "/home/ahyan/Dropbox/Underwater/DataWork/creditsurface/empirical_credit_surface.R"
        List<String[]> rowList2 = new ArrayList<String[]>();
        try(BufferedReader br = new BufferedReader(new FileReader(filePathForCS))){
            String line;
            while((line = br.readLine()) != null){
                String[] lineItems = line.split(",");
                rowList2.add(lineItems);
            }
            br.close();
        }catch (Exception e){e.printStackTrace();}

        // create a matrix to hold the empirical credit surface, this matrix is used for interpolation
        double[][] creditSurface = new double[fico_grid.length][oltv_grid_data.length];

        // core loops to map empirical credit surface from list form in the csv to a 2x2 matrix
        for (int ifico = 0; ifico < fico_grid.length; ifico++){
            for (int ioltv = 0; ioltv < oltv_grid_data.length; ioltv++){
                creditSurface[ifico][ioltv] = Double.parseDouble(rowList2.get(ioltv + (oltv_grid_data.length * ifico))[2]);
            }
        }

        // moving on to the interpolation step
        int noltv_interpolate = oltv_grid_interpolate.length;
        double[][] credit_surface_interpolate = new double[fico_grid.length][noltv_interpolate];

        UnivariateInterpolator univariateInterpolator = new LinearInterpolator();

        for (int ifico = 0; ifico < fico_grid.length; ifico++){

            UnivariateFunction univariateFunction = univariateInterpolator.interpolate(oltv_grid_data, creditSurface[ifico]);

            for (int ioltv = 0; ioltv < noltv_interpolate; ioltv++){
                credit_surface_interpolate[ifico][ioltv] = univariateFunction.value(oltv_grid_interpolate[ioltv]);
            }
        }


        return credit_surface_interpolate;
    }

    // grid for oltv (it calls the credit surface csv and makes a grid for the oltvs for fico=500)
    double[] oltv_grid(String filePathForCS, int noltv){
        double[] lgrid = new double[noltv];
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(filePathForCS));
            String row;
            int it = 0;
            while ((row = bufferedReader.readLine()) != null & it < noltv) {
                String[] data = row.split(",");
                lgrid[it] = Double.parseDouble(data[1]);
                it++;
            }
        }catch (IOException e){e.printStackTrace();}
        return lgrid;
    }

    double[][][][][][] Expected_sfha(int J, int na, int noltv, int nfico, int N, int ntheta,
                                     double[][][][][][] V_sfha, double[][] transition_matrix, int age){
        double[][][][][][] Expected_sfha = new double[J][na][noltv][nfico][N][ntheta];
        double expected = 0d;
        for (int ia = 0; ia < na; ia++){
            for (int ioltv = 0; ioltv < noltv; ioltv++){
                for (int ifico = 0; ifico < nfico; ifico++){
                    for (int n = 0; n < N; n++){
                        for (int itheta = 0; itheta < ntheta; itheta++){
                            for (int ithetap = 0; ithetap < ntheta; ithetap++){
                                expected = expected + transition_matrix[itheta][ithetap] * V_sfha[age][ia][ioltv][ifico][n][ithetap];
                            }
                            Expected_sfha[age][ia][ioltv][ifico][n][itheta] = expected;
                            expected = 0;
                        }
                    }
                }
            }
        }
        return Expected_sfha;
    }

    double[] read_vector(String file_path, int number_of_grid_points){
        double[] x_grid = new double[number_of_grid_points];
        try {
            BufferedReader bufferedReader = new BufferedReader(new FileReader(file_path));
            String row;
            int it = 0;
            while ((row = bufferedReader.readLine()) != null & it < number_of_grid_points) {
                String[] data = row.split(",");
                x_grid[it] = Double.parseDouble(data[1]);
                it++;
            }
        }catch (IOException e){e.printStackTrace();}
        return x_grid;
    }
}
