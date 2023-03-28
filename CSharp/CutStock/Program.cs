using System;
using System.IO;
using System.Security.AccessControl;
using System.Text;
using Google.OrTools.LinearSolver;
class Data
{
    public double mateiralLen;
    public List<double> partLen = new List<double>();
    public List<int> partDemand = new List<int>();
    //patterns,key:part index in partLen and partDemand,value:num of part
    public List<Dictionary<int, int>> patterns = new List<Dictionary<int, int>>();
    public void ReadData(String filename)
    {
        String[] text = File.ReadAllLines(filename, Encoding.UTF8);
        this.mateiralLen = Convert.ToDouble(text[0]);
        String partLenTxt = text[1].Replace('[', ' ').Replace(']', ' ').Trim();
        string[] partLenList = partLenTxt.Split(",");
        foreach (string str in partLenList)
        {
            partLen.Add(Convert.ToDouble(str));
        }
        String partDemandTxt = text[2].Replace('[', ' ').Replace(']', ' ').Trim();
        string[] partDemandList = partDemandTxt.Split(",");
        foreach (string str in partDemandList)
        {
            partDemand.Add(Convert.ToInt32(str));
        }
    }
    public void PrintData()
    {
        Console.WriteLine($"materialLen:{mateiralLen}");
        Console.WriteLine($"part len:");
        partLen.ForEach(Console.WriteLine);
        Console.WriteLine($"part demand:");
        partDemand.ForEach(Console.WriteLine);
    }
}
class CutStock
{
    private Data data = new Data();
    public CutStock(Data data)
    {
        this.data = data;
    }
    private double SolveSubProblme(List<double> dualValues)
    {
        Solver solver = Solver.CreateSolver("SCIP");
        List<Variable> x = new List<Variable>();
        for (int j = 0; j < data.partLen.Count; ++j)
        {
            Variable var = solver.MakeIntVar(0.0, double.PositiveInfinity, $"x_{j}");
            x.Add(var);
        }
        Constraint constraint = solver.MakeConstraint(0, data.mateiralLen, $"cons");
        for (int j = 0; j < x.Count; ++j)
        {
            constraint.SetCoefficient(x[j], data.partLen[j]);
        }
        Objective objective = solver.Objective();
        for (int j = 0; j < x.Count; ++j)
        {
            objective.SetCoefficient(x[j], dualValues[j]);
        }
        objective.SetMaximization();
        Solver.ResultStatus resultStatus = solver.Solve();
        if (resultStatus != Solver.ResultStatus.OPTIMAL)
        {
            Console.WriteLine("master problem no optimal solution!");
            return 1;
        }
        else
        {
            Dictionary<int, int> newPattern = new Dictionary<int, int>();
            for (int j = 0; j < x.Count; j++)
            {
                if (x[j].SolutionValue() != 0)
                {
                    newPattern.Add(j, Convert.ToInt32(x[j].SolutionValue()));
                }
            }
            data.patterns.Add(newPattern);
            return 1 - objective.Value();
        }
    }
    private List<double> SolveMasterProblem()
    {
        Solver solver = Solver.CreateSolver("CLP");
        if (solver == null)
        {
            return new List<double>();
        }
        List<Variable> x = new List<Variable>();
        for (int j = 0; j < data.patterns.Count; ++j)
        {
            Variable var = solver.MakeNumVar(0.0, double.PositiveInfinity, $"x_{j}");
            x.Add(var);
        }
        List<Constraint> cons = new List<Constraint>();
        for (int i = 0; i < data.partLen.Count; ++i)
        {
            Constraint constraint = solver.MakeConstraint(data.partDemand[i], double.PositiveInfinity, $"cons_{i}");
            for (int j = 0; j < x.Count; ++j)
            {
                if (data.patterns[j].ContainsKey(i))
                    constraint.SetCoefficient(x[j], data.patterns[j][i]);
            }
            cons.Add(constraint);
        }
        Objective objective = solver.Objective();
        for (int j = 0; j < x.Count; ++j)
        {
            objective.SetCoefficient(x[j], 1);
        }
        objective.SetMinimization();
        //string lp = solver.ExportModelAsLpFormat(false);
        //File.WriteAllText("master.lp", lp);
        Solver.ResultStatus resultStatus = solver.Solve();
        List<double> dualValues = new List<double>();
        if (resultStatus != Solver.ResultStatus.OPTIMAL)
        {
            Console.WriteLine("master problem no optimal solution!");
            return dualValues;
        }
        else
        {
            Console.WriteLine($"relax value is {objective.Value()}");
            for (int i = 0; i < cons.Count; ++i)
            {
                dualValues.Add(cons[i].DualValue());
            }
            return dualValues;
        }
    }
    private void SolveFinalProblem()
    {
        Solver solver = Solver.CreateSolver("SCIP");
        List<Variable> x = new List<Variable>();
        for (int j = 0; j < data.patterns.Count; ++j)
        {
            Variable var = solver.MakeIntVar(0.0, double.PositiveInfinity, $"x_{j}");
            x.Add(var);
        }
        List<Constraint> cons = new List<Constraint>();
        for (int i = 0; i < data.partLen.Count; ++i)
        {
            Constraint constraint = solver.MakeConstraint(data.partDemand[i], double.PositiveInfinity, $"cons_{i}");
            for (int j = 0; j < x.Count; ++j)
            {
                if (data.patterns[j].ContainsKey(i))
                    constraint.SetCoefficient(x[j], data.patterns[j][i]);
            }
        }
        Objective objective = solver.Objective();
        for (int j = 0; j < x.Count; ++j)
        {
            objective.SetCoefficient(x[j], 1);
        }
        objective.SetMinimization();
        Solver.ResultStatus resultStatus = solver.Solve();
        if (resultStatus != Solver.ResultStatus.OPTIMAL)
        {
            Console.WriteLine("master problem no optimal solution!");
        }
        else
        {
            List<int> solution = new List<int>();
            for (int i = 0; i < x.Count; ++i)
            {
                solution.Add(Convert.ToInt32(Math.Round(x[i].SolutionValue())));
            }
            double totalPartLen = 0.0;
            for (int i = 0; i < data.partLen.Count; ++i)
            {
                totalPartLen += data.partLen[i] * data.partDemand[i];
            }
            Console.WriteLine("solution:");
            for (int i = 0, index = 0; i < solution.Count; ++i)
            {
                if (solution[i] != 0)
                {
                    Console.Write($"pattern_{index}:");
                    foreach (KeyValuePair<int, int> kvp in data.patterns[i])
                    {
                        Console.Write($"{data.partLen[kvp.Key]}*{kvp.Value} ");
                    }
                    Console.Write(solution[i]);
                    Console.WriteLine();
                    ++index;
                }
            }
            Console.WriteLine($"used roll is {Convert.ToInt64(objective.Value())}");
            Console.WriteLine("utilization:{0:P}", totalPartLen / (Convert.ToInt64(objective.Value()) * data.mateiralLen));
        }
    }
    public void Run()
    {
        //create initial patterns
        for (int i = 0; i < data.partLen.Count; ++i)
        {
            Dictionary<int, int> pattern = new Dictionary<int, int>();
            pattern.Add(i, Convert.ToInt32(Math.Floor(data.mateiralLen / data.partLen[i])));
            data.patterns.Add(pattern);
        }
        while (true)
        {
            List<double> dualValues = SolveMasterProblem();
            double reduceCost = SolveSubProblme(dualValues);
            System.Console.WriteLine($"reduce cost is {reduceCost}");
            if (reduceCost > -1e-3) break;
        }
        SolveFinalProblem();
    }
}
class Program
{
    static void Main(String[] args)
    {
        if (args.Length != 1)
        {
            Console.WriteLine("usage:CutStock filename");
        }
        else
        {
            var stopWatch = new System.Diagnostics.Stopwatch();
            stopWatch.Start();
            Data d = new Data();
            d.ReadData(args[0]);
            d.PrintData();
            CutStock cutStock = new CutStock(d);
            cutStock.Run();
            stopWatch.Stop();
            Console.WriteLine("used time is {0} ms", stopWatch.ElapsedMilliseconds);
        }
    }
}