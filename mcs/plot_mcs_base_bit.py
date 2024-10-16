import os,sys
import argparse

import random

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
except (ImportError, RuntimeError):
    print("COULD not import matplotlib")

import seaborn as sns
import pandas as pd
from matplotlib import pyplot

def plot_fraction(input_csv, outfile, linewidth = 2.5): #, palette, tools, linewidth = 2.5):
    matplotlib.rcParams.update({'font.size': 18})
    # sns.set(font_scale=1.5)
    sns.set(font_scale=1.6) 
    sns.set_style("whitegrid")
    indata = pd.read_csv(input_csv)
    g = sns.relplot(data=indata, x="B", y="Exp_B", hue="k", style="w", linewidth = linewidth, kind="line") #, #dashes = dashes,
        #col="dataset") # hue_order = tools, # hue="datastructure", style="datastructure",
        # col_wrap=3, col_order=["SIM1", "SIM2", "SIM4"], palette=palette)
        # col_order=["SIM3"],  palette=palette)
    # ax = sns.lineplot(data=indata, x="k", y="unique", hue="datastructure", style="chr", palette = sns.color_palette()[:7])
    # axes = g.axes
    g.set_titles("Genome size 2^16 (repeats)")
    g.set_axis_labels("Bit space B (2^x)", "Fraction of B for base")
    # g.set(ylim=(94, 99), xticks=[50,75,100,150,200,250,300,500])
    # g.set_xticks([12, 16, 20, 24, 28])
    g.set_xticklabels(labels=[12, 16, 20, 24, 28]) #rotation=60, 
    g.tight_layout()
    # g.set(ylim=(95, 100))
    plt.savefig(outfile)
    plt.close()

def main(args):
    sns.set_style("whitegrid")

    # palette = {
    # 'minimap2': 'tab:blue',
    # 'bwa_mem': 'tab:orange',
    # 'strobealign_v0140': 'tab:green',
    # # "strobealign_v0120_opt" : 'pink',
    # "strobealign_mcs_022a721" : 'black'
    # }
    # # tools = ["minimap2", "bwa_mem", "strobealign_v071", "strobealign_v0120_opt", "strobealign_multicontext"]
    # tools = ["minimap2", "bwa_mem", "strobealign_v0140", "strobealign_mcs_022a721"]

    plot_fraction(args.csv, args.outfile) #, palette, tools, args.type)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calc identity", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('csv', type=str, default= "", help='results file')
    parser.add_argument('outfile', type=str,  help='outfile.')
    args = parser.parse_args()



    if len(sys.argv)==1:
        parser.print_help()
        sys.exit()

    main(args)