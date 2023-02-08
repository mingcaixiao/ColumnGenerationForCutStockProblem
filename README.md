# cut stock problem (use column generation method to solve)
## cut stock problem 
[cut stock problem](https://en.wikipedia.org/wiki/Cutting_stock_problem)
In operations research, the cutting-stock problem is the problem of cutting standard-sized pieces of stock material, such as paper rolls or sheet metal, into pieces of specified sizes while minimizing material wasted. It is an optimization problem in mathematics that arises from applications in industry. In terms of computational complexity, the problem is an NP-hard problem reducible to the knapsack problem. The problem can be formulated as an integer linear programming problem.

## column generation
[column generation](https://en.wikipedia.org/wiki/Column_generation)
Column generation or delayed column generation is an efficient algorithm for solving large linear programs.

In cutting stock problem,Num cutting pattern will increases exponentially as  num of different part increases.

So exploiting all cutting patterns in advance is not practical.

Column generation is a very efficient method to solve this problem.
## Depend on
[Google OR-Tools](https://github.com/google/or-tools):Google's Operations Research tools.

I use ortools to solve master linear programing,sub and final interger programing.

To **build** this project,must install google or-tools first

## Directory
- root
    - src:contains cpp file
    - data:data file

Data files contains three lines.The first line is roll width.
the second line is part length list.the last line is part demand list.

## Build
This project use cmake as build system.
You can build like this:
```
mkdir build
cd build
cmake ..
cmake --build .
```

After build,
run `cutstock filename` to execute.

filename is file in directory "root/data"
