import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd


if __name__ == "__main__":
    file = os.path.join('', 'index-comparison.csv')
    output_file = 'index-comparison.pdf'
    df = pd.read_csv(file, delimiter=',')
    df.avg_query_time = df.avg_query_time.round(2)

    fig, axs = plt.subplots(nrows=2, ncols=2, figsize=(8,8))

    dfs = df.groupby(['benchmark_name'])

    benchmarks = df['benchmark_name'].unique()
    j = 0
    i =0
    for benchmark in benchmarks:
        df_b = df[df['benchmark_name'] == benchmark]
        x = df_b['index_name']
        y = df_b['avg_query_time']
        b = axs[i,j].bar(x,y)
        axs[i,j].bar_label(b,padding=12,fontsize=8)
        axs[i,j].set_title(benchmark)
        axs[i,j].set_ylabel("Average Query Time (ms)")
        axs[i,j].margins(y=0.2)
        print(df_b)
        j+=1
        if j == 2:
            j = 0
            i += 1
    fig.suptitle('Index Comparison for Range Query')
    fig.savefig(os.path.join('', output_file), bbox_inches='tight')
    #plt.show()