import sys

from ortools.linear_solver import pywraplp

import time





class Data:

    def __init__(self):

        self.partLen = []

        self.partDemand = []

        self.materialLen = 0

        self.patterns: list[dict[int, int]] = []  # key:len index,value:num





def read_array_from_str(string: str) -> list[float]:

    string = string.removesuffix('\n')

    string = string.removesuffix(']')

    string = string.removeprefix('[')

    string_list = string.split(',')

    d: list[float] = [float(string_list[i]) for i in range(len(string_list))]

    return d





def read_data(filename: str) -> Data:

    d: Data = Data()

    with open(filename, "r") as f:

        lines = f.readlines()

        d.materialLen = float(lines[0])

        d.partLen = read_array_from_str(lines[1])

        d.partDemand = read_array_from_str(lines[2])

    return d





def print_data(d: Data) -> None:

    print(f'material len is {d.materialLen}')

    print(f'part len is  {d.partLen}')

    print(f'part demand is {d.partDemand}')





def solve_restricted_master_problem(d: Data) -> list[float]:

    solver = pywraplp.Solver.CreateSolver('CLP')

    if not solver:

        return []

    x_list = []

    for (pattern, j) in zip(d.patterns, range(len(d.patterns))):

        x = solver.NumVar(0, solver.infinity(), f'x_{j}')

        x_list.append(x)

    cons_list = []

    for i in range(len(d.partLen)):

        constraint = solver.RowConstraint(d.partDemand[i], solver.infinity(), f'c_{i}')

        for j in range(len(d.patterns)):

            if i in d.patterns[j]:

                constraint.SetCoefficient(x_list[j], d.patterns[j][i])

        cons_list.append(constraint)



    solver.Minimize(solver.Sum([x_list[j] for j in range(len(d.patterns))]))

    # with open("master.lp", "w") as out_f:

    #     lp_text = solver.ExportModelAsLpFormat(obfuscated=False)

    #     out_f.write(lp_text)

    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:

        print(f'relax value is {solver.Objective().Value():.5f}')

        return [cons_list[i].dual_value() for i in range(len(d.partLen))]

    else:

        return []





def solve_sub_problem(d: Data, dual_values: list[float]) -> float:

    """

    solve sub problem to find new pattern

    :param d: input data

    :param dual_values: master problem dual values

    :return: reduce cost

    """

    solver = pywraplp.Solver.CreateSolver('SCIP')

    if not solver:

        return 1

    x = []

    # Create decision variables

    for i in range(len(d.partLen)):

        var = solver.IntVar(0, d.materialLen / d.partLen[i], 'x_%i' % i)

        x.append(var)



    # Create constraint

    solver.Add(sum([x[i] * d.partLen[i] for i in range(len(d.partLen))]) <= d.materialLen)

    solver.Maximize(sum([x[i] * dual_values[i] for i in range(len(d.partLen))]))

    # with open("sub.lp", "w") as out_f:

    #     lp_text = solver.ExportModelAsLpFormat(obfuscated=False)

    #     out_f.write(lp_text)

    status = solver.Solve()

    if status == pywraplp.Solver.OPTIMAL:

        item = {}

        for i in range(len(d.partLen)):

            if x[i].solution_value() != 0:

                item[i] = x[i].solution_value()

        d.patterns.append(item)

        # Calculates reduce cost

        return 1 - solver.Objective().Value()

    else:

        return 1





def solve_final_master_problem(d: Data):

    solver = pywraplp.Solver.CreateSolver('SCIP')

    if not solver:

        return

    x_list = []

    for (pattern, j) in zip(d.patterns, range(len(d.patterns))):

        x = solver.IntVar(0, solver.infinity(), f'x_{j}')

        x_list.append(x)

    cons_list = []

    for i in range(len(d.partLen)):

        constraint = solver.RowConstraint(d.partDemand[i], solver.infinity(), f'c_{i}')

        for j in range(len(d.patterns)):

            if i in d.patterns[j]:

                constraint.SetCoefficient(x_list[j], d.patterns[j][i])

        cons_list.append(constraint)



    solver.Minimize(solver.Sum([x_list[j] for j in range(len(d.patterns))]))

    status = solver.Solve()

    solution = {}

    if status == pywraplp.Solver.OPTIMAL:

        for j in range(len(d.patterns)):

            if x_list[j].solution_value() != 0:

                solution[j] = x_list[j].solution_value()

        print(f'used roll num is {solver.Objective().Value()}')

        print('plan is:')

        for (pattern_index, num) in solution.items():

            for len_index, item_num in d.patterns[pattern_index].items():

                print(f'{d.partLen[len_index]}*{int(item_num)}', end=' ')

            print(int(num))

        total_part_len = sum([d.partLen[i] * d.partDemand[i] for i in range(len(d.partLen))])

        print(f'utilization is {total_part_len / (solver.Objective().Value() * d.materialLen) * 100:.2f}%')

    else:

        return





def cut_stock(d: Data):

    # init some patterns

    for i in range(len(d.partLen)):

        pattern = {i: int(d.materialLen / d.partLen[i])}

        d.patterns.append(pattern)

    while True:

        dual_values = solve_restricted_master_problem(d)

        reduce_cost = solve_sub_problem(d, dual_values)

        print(f'reduce cost is {reduce_cost}')

        if reduce_cost > -1e-3:

            break

    solve_final_master_problem(d)





if __name__ == '__main__':

    if len(sys.argv) != 2:

        print('usage:python cutstock.py filename')

        exit(1)

    else:

        data = read_data(sys.argv[1])

        print_data(data)

        start = time.time()

        cut_stock(data)

        print(f'time used is {time.time() - start:.5f} s')

