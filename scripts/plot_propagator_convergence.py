#!/usr/bin/env python

import sys
import numpy as np
import csv
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

#######################################
# Convergence Analysis
#######################################

prop_list = sys.argv[1:]

if not len(prop_list):
    print ('TOO FEW ARGUMENTS: Please pass the error files as command line parameters');

fig = plt.figure()
ax = fig.gca()


for p in prop_list:

    with open(p, 'rb') as f:
        data = csv.reader(row for row in f if (not row.startswith('#') and row.strip()));
        meta = data.next();
        name, coefs, T = meta;
        print; print(meta);
        conv = []
        for row in data:
            conv.append(row);

        Dt = np.array(conv)[:,0]
        err = np.array(conv)[:,1]

        ax.loglog(Dt, err, '-o', label=name+' (' + coefs + ' splitting)')

lgd = ax.legend(bbox_to_anchor=(1.02, 1), loc=2, borderaxespad=0.)
# ax.set_xlim([5e1,1e5])
# ax.set_ylim(view[2:])
# ax.ticklabel_format(style="sci", scilimits=(0, 0), axis="y")
ax.set_title(r"$L_2$ error vs. number of step size (T=" + T + ")");
ax.grid(True);
ax.set_xlabel(r"step size $\Delta t$");
ax.set_ylabel(r"$L_2$ error $\frac{\| u (x) - u_* (x) \|}{\| u_* (x) \|}$")
fig.savefig("error_analysis.png",bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.show()
plt.close(fig)
