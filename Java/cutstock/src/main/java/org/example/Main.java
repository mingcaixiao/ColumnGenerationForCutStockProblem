package org.example;



import com.google.ortools.linearsolver.MPConstraint;

import com.google.ortools.linearsolver.MPObjective;

import com.google.ortools.linearsolver.MPSolver;

import com.google.ortools.linearsolver.MPVariable;

import com.google.ortools.Loader;



import java.io.*;

import java.util.ArrayList;

import java.util.HashMap;





class Data {

    Double materialLen = 0.0;

    ArrayList<Double> partLen = new ArrayList<>();

    ArrayList<Integer> partDemand = new ArrayList<>();

    //pattern list ,key:part index , value:part frequency

    ArrayList<HashMap<Integer, Integer>> patterns = new ArrayList<>();



    @Override

    public String toString() {

        return "Data{" + "\n" +

                "   materialLen=" + materialLen + "\n" +

                "   PartLen=" + partLen + "\n" +

                "   partDemand=" + partDemand + "\n" +

                '}';

    }



    public void parseFromFile(String fileName) {

        try {

            BufferedReader in = new BufferedReader(new FileReader(fileName));

            String materialLenStr = in.readLine();

            this.materialLen = Double.valueOf(materialLenStr);

            String partLenStr = in.readLine();

            partLenStr = partLenStr.replace(']', ' ');

            partLenStr = partLenStr.replace('[', ' ');

            partLenStr = partLenStr.trim();

            String[] partLenStrList = partLenStr.split(",");

            for (String len : partLenStrList) {

                this.partLen.add(Double.valueOf(len.trim()));

            }

            String partDemandStr = in.readLine();

            partDemandStr = partDemandStr.replace(']', ' ');

            partDemandStr = partDemandStr.replace('[', ' ');

            partDemandStr = partDemandStr.trim();

            String[] partDemandStrList = partDemandStr.split(",");

            for (String num : partDemandStrList) {

                this.partDemand.add(Integer.valueOf(num.trim()));

            }

        } catch (IOException e) {

            System.out.println(e);

        }

    }





}



class StockCutter {

    private Data data;



    public StockCutter(Data data) {

        this.data = data;

    }



    public void run() {

        long start = System.currentTimeMillis();

        //create initial pattern

        for (int i = 0; i < data.partLen.size(); ++i) {

            HashMap<Integer, Integer> pattern = new HashMap<>();

            pattern.put(i, Double.valueOf(data.materialLen / data.partLen.get(i)).intValue());

            data.patterns.add(pattern);

        }

        while (true) {

            ArrayList<Double> dualValues = solveMasterProblem();

            double reduceCost = solveSubProblem(dualValues);

            System.out.printf("reduce cost is %f%n", reduceCost);

            if (reduceCost > -1e-3) break;

        }

        solveFinalProblem();

        System.out.printf("time is %d ms\n", System.currentTimeMillis() - start);

    }



    private double solveSubProblem(ArrayList<Double> dualValues) {

        MPSolver solver = MPSolver.createSolver("SCIP");

        ArrayList<MPVariable> x = new ArrayList<>();

        for (int j = 0; j < data.partLen.size(); ++j) {

            MPVariable var = solver.makeIntVar(0.0, data.materialLen / data.partLen.get(j), "x_" + j);

            x.add(var);

        }

        MPConstraint cons = solver.makeConstraint(0.0, data.materialLen, "cons");

        for (int j = 0; j < data.partLen.size(); ++j) {

            cons.setCoefficient(x.get(j), data.partLen.get(j));

        }

        MPObjective objective = solver.objective();

        for (int j = 0; j < data.partLen.size(); ++j) {

            objective.setCoefficient(x.get(j), dualValues.get(j));

        }

        objective.setMaximization();

        final MPSolver.ResultStatus resultStatus = solver.solve();

        if (resultStatus == MPSolver.ResultStatus.OPTIMAL) {

            HashMap<Integer, Integer> pattern = new HashMap<>();

            for (int i = 0; i < x.size(); ++i) {

                int value = (int) Math.round(x.get(i).solutionValue());

                if (value != 0) pattern.put(i, value);

            }

            data.patterns.add(pattern);

        } else {

            System.err.println("no optimal solution!");

        }

        return 1 - objective.value();

    }



    private ArrayList<Double> solveMasterProblem() {

        MPSolver solver = MPSolver.createSolver("GLOP");

        double infinity = Double.POSITIVE_INFINITY;

        ArrayList<MPVariable> x = new ArrayList<>();

        for (int j = 0; j < data.patterns.size(); ++j) {

            MPVariable var = solver.makeNumVar(0.0, infinity, "x_" + j);

            x.add(var);

        }

        ArrayList<MPConstraint> cons = new ArrayList<>();

        for (int i = 0; i < data.partLen.size(); ++i) {

            MPConstraint constraint = solver.makeConstraint(data.partDemand.get(i), infinity, "c_" + i);

            for (int j = 0; j < data.patterns.size(); ++j) {

                if (data.patterns.get(j).containsKey(i))

                    constraint.setCoefficient(x.get(j), data.patterns.get(j).get(i));

            }

            cons.add(constraint);

        }

        MPObjective objective = solver.objective();

        for (MPVariable mpVariable : x) {

            objective.setCoefficient(mpVariable, 1);

        }

        objective.setMinimization();

//        String lp = solver.exportModelAsLpFormat();

//        try {

//            BufferedWriter out = new BufferedWriter(new FileWriter("master.lp"));

//            out.write(lp);

//            out.close();

//        } catch (IOException e) {

//            System.out.println(e);

//        }



        final MPSolver.ResultStatus resultStatus = solver.solve();

        ArrayList<Double> dualValues = new ArrayList<>();

        if (resultStatus == MPSolver.ResultStatus.OPTIMAL) {

            for (MPConstraint c : cons) {

                dualValues.add(c.dualValue());

            }

            System.out.printf("relax value is %f%n", objective.value());

        } else {

            System.err.println("no optimal solution!");

        }

        return dualValues;

    }



    private void solveFinalProblem() {

        MPSolver solver = MPSolver.createSolver("SCIP");

        double infinity = Double.POSITIVE_INFINITY;

        ArrayList<MPVariable> x = new ArrayList<>();

        for (int j = 0; j < data.patterns.size(); ++j) {

            MPVariable var = solver.makeIntVar(0.0, infinity, "x_" + j);

            x.add(var);

        }

        for (int i = 0; i < data.partLen.size(); ++i) {

            MPConstraint constraint = solver.makeConstraint(data.partDemand.get(i), infinity, "c_" + i);

            for (int j = 0; j < data.patterns.size(); ++j) {

                if (data.patterns.get(j).containsKey(i))

                    constraint.setCoefficient(x.get(j), data.patterns.get(j).get(i));

            }

        }

        MPObjective objective = solver.objective();

        for (MPVariable mpVariable : x) {

            objective.setCoefficient(mpVariable, 1);

        }

        objective.setMinimization();

        final MPSolver.ResultStatus resultStatus = solver.solve();

        if (resultStatus == MPSolver.ResultStatus.OPTIMAL) {

            ArrayList<Integer> solution = new ArrayList<>();

            for (MPVariable mpVariable : x) {

                solution.add((int) Math.round(mpVariable.solutionValue()));

            }

            System.out.println("solution is");

            double totalPartLen = 0.0;

            for (int i = 0; i < data.partLen.size(); ++i) {

                totalPartLen += data.partLen.get(i) * data.partDemand.get(i);

            }

            for (int j = 0, index = 0; j < solution.size(); ++j) {

                if (solution.get(j) != 0) {

                    System.out.printf("pattern_%d:", index);

                    data.patterns.get(j).forEach((lenIndex, num) -> System.out.printf("%.2f*%d ", data.partLen.get(lenIndex), num));

                    System.out.printf(" %d\n", solution.get(j));

                    ++index;

                }

            }

            System.out.printf("used roll is %d\n", (int) objective.value());

            System.out.printf("utilization is %.2f %%\n", totalPartLen / (data.materialLen * objective.value()) * 100);

        } else {

            System.err.println("no optimal solution!");

        }

    }

}



public class Main {

    public static void main(String[] args) {

        if (args.length != 1) {

            System.out.println("usage:java -jar ./cutstock-1.0-SNAPSHOT-jar-with-dependencies.jar filename");

            System.exit(1);

        } else {

            Loader.loadNativeLibraries();

            Data data = new Data();

            data.parseFromFile(args[0]);

            System.out.println(data);

            StockCutter stockCutter = new StockCutter(data);

            stockCutter.run();

        }



    }

}