import matplotlib.pyplot as plt
import pathlib
import os.path


colors = ['b', 'g', 'r', 'c', 'm', 'y', '#C0C0C0', '#800000', '#FFA500', '#FFFFF0', '#F0E68C', '#FF69B4' , '#7FFFD4', '#FFDAB9' , '#8FBC8F' , '#DB7093' , '#FFC0CB', '#87CEFA' , '#FFE4B5', '#40E0D0', '#D8BFD8', '#228B22' , '#FF7F50', '#FFA07A' , '#B0E0E6']

benchmark_dir = pathlib.Path(__file__).parent.resolve()

### RETRIEVE VALUES
if os.path.isfile(str(benchmark_dir)+"/statsFiles/tmpQueryTime"):
    with open(str(benchmark_dir)+"/statsFiles/tmpQueryTime", "r") as fh:
        qsizes = fh.readline().rstrip().split("q")[1:] 
        vals = [0] * len(qsizes) * 2 * 95

        for c, line in enumerate(fh):
            vals[c] = float(line.rstrip())


    ### QUERY TIMES
    queryTimes = [[list(), list()] for _ in qsizes]
    def ms_to_millionsPerSec(n):
        n /= 1000 #ms to s
        return 1/n / 10 #how many 100k queries in 1s then how many millions
    for i in range(len(qsizes)):
        queryTimes[i][0] = list(map(ms_to_millionsPerSec, vals[i*190 : i*190+190-1 : 2] ))
        queryTimes[i][1] = list(map(ms_to_millionsPerSec, vals[i*190+1 : i*190+190 : 2] ))


    loads = range(1, 96)
    for i in range(len(qsizes)):
        plt.plot(loads, queryTimes[i][0], label = "q="+qsizes[i]+"-positiveQuery", color=colors[i])
        plt.plot(loads, queryTimes[i][1], '--', label = "q="+qsizes[i]+"-negativeQuery", color=colors[i])

    plt.xlabel("Load factor (%)")
    plt.ylabel("Millions insertions / s")
    plt.title(f"Query perfs for q={', '.join(qsizes)}")
    plt.legend()
    plt.show()
else:
    print("no query times found")


if os.path.isfile(str(benchmark_dir)+"/statsFiles/tmpInsertTime"):
    ### RETRIEVE VALUES
    with open(str(benchmark_dir)+"/statsFiles/tmpInsertTime", "r") as fh:
        infos = fh.readline().rstrip().split(";")
        insertTimes = [0] * 95

        for c, line in enumerate(fh):
            insertTimes[c] = float(line.rstrip())


    ### INSERT TIMES
    def ms_to_millionsPerSec(n):
        unit = int(infos[1][4:])
        n /= 1000 #ms to s
        return 1 / n  * unit / 100000

    insertTimes = list(map(ms_to_millionsPerSec, insertTimes))



    plt.plot(range(1, 96), insertTimes)

    plt.xlabel("Load factor (%)")
    plt.ylabel("Millions insertions / s")
    plt.title(f"Insertion perfs for q={infos[0][1:]}")
    plt.show()
else:
    print("no query times found")