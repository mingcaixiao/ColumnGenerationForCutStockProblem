# cut stock problem (use column generation method to solve)
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