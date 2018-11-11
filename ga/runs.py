import os
import time

runs = 2
xover_op = ["edge-recombination"]
xover_rates = [0.4]
mut_op = ["insertion"]
mut_rates = [0.25]
uniform_xover_rates = [0.5]

start = time.time()

for x_o in xover_op:
    for x in xover_rates:
        for m_o in mut_op:
            for m in mut_rates:
                for u in uniform_xover_rates:
                    for i in range(runs):
                        os.system("./ga params opfiles fx.tsp " + str(x) + " " + str(m) + " " + str(u) + " " + str(
                            x_o) + " " + str(m_o))

end = time.time()
duration = end - start
avg_dur = duration / runs
print("total duration: " + str(duration))
print("avg duration per run: " + str(avg_dur))