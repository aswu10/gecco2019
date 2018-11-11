import os
import time

runs = 1
parsimony = [1]
xover_rates = [0.8]
mut_rates = [0.01]

start = time.time()

for p in parsimony:
    for x in xover_rates:
        for m in mut_rates:
            for i in range(runs):
                os.system("./ga params_tsp opfiles fx.tsp " + str(x) + " " + str(m) + " " + str(p))
                                
end = time.time()
duration = end - start
avg_dur = duration / 2
print("total duration: " + str(duration))
print("avg duration per run: " + str(avg_dur))