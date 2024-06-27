# qdags
 
To compile:
```
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib program.cpp -o program -lsdsl -ldivsufsort -ldivsufsort64 
```
-std=c++11 -O3 -DNDEBUG -I /opt/homebrew/include -L /opt/homebrew/lib -lsdsl -ldivsufsort -ldivsufsort64

```
g++ -std=c++11 -O3 -DNDEBUG -I ~/include -L ~/lib -lsdsl -ldivsufsort -ldivsufsort64 program.cpp -o program  
```

To run main 
```
./main ./data/prop-direct-P197 ./data/prop-direct-P800 ./data/prop-direct-P1366 
```

# Tests
In the queries folder we have all the queries to be tested grouped by the type of test: 
- partial: 
  - dfuds:
    - fixedQueue
    - nonFixedQueue
  - louds:
    - fixedQueue
    - nonFixedQueue
- ranked:
    - dfuds:
        - fixedQueue
        - nonFixedQueue
    - louds:
        - fixedQueue
        - nonFixedQueue
- lqdags
- traditional algorithm (Arroyuelo et al. 2022)

The script files to run the tests are in the runqueries folder. On each subfolder we have:
- Compile script: it compiles all the queries.
- Run script: it runs all the queries and put the outputs in the respective output folder.
- All the run queries to be tested with different data. The files are named as runqueries-[type of query]-bfs-sorted.sh

## To run a particular test
1. First go to the runqueries folder and choose the subfolder you want to test. For example
` cd runqueries/partial/dfuds/fixedQueue`
2. Run the compile script (the firs)

    For GNU CC: `./compilePartialDfudsFixedQueue.sh`

    For Clang: `./compilePartialDfudsFixedQueueM1.sh`

3. Run the run script
`./runPartialDfudsFixed.sh`

## To run all the tests
1. Go to the runqueries folder
`cd runqueries`
2. Run the compile script

    For GNU CC `./compileAllAlgorithms.sh`

    For Clang `./compileAllAlgorithmsM1.sh`

3. Run all scripts
`./runAllAlgorithms.sh`