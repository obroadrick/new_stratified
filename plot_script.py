"""
This script is used to plot the first round sizes in round_size.json.

Oliver Broadrick 2020
"""

import json
import pprint
import matplotlib.pyplot as plt

# open json file
with open('round_sizes.json') as json_file:
    data = json.load(json_file)
    pprint.pprint(data)

# parse file for desired data
new_method_round_sizes = []
suite_r2bravo_round_sizes = []
fractional_margins = []
comparisons = data['audits'][0]['comparisons']
for audit in data['audits']:
    fractional_margin = audit['fractional_margin']
    new_method_round_size = audit['round_size_new_method']
    suite_r2bravo_round_size = audit['round_size_suite_r2bravo']

    fractional_margins.append(fractional_margin)
    new_method_round_sizes.append(new_method_round_size)
    suite_r2bravo_round_sizes.append(suite_r2bravo_round_size)

# make a pretty plot of the data
fig = plt.figure(figsize=(20,10))
fig.suptitle('First Round Sizes for Varying Margins (with '+str(comparisons)+' comparisons)', fontsize=20)
ax = fig.add_subplot(111)#numrows, numcols, num of subplot being referenced
ax.plot(fractional_margins, new_method_round_sizes, color='b', marker='o', label='New Method (Sick Math)')
ax.plot(fractional_margins, suite_r2bravo_round_sizes, color='r', marker='x', label='SUITE w/ R2 Bravo')
ax.set_xlabel('Fractional Margin', fontsize=20)
ax.set_ylabel('First Round Size (achieves 90% stopping probability)', fontsize=20)
plt.legend(loc='upper left', fontsize=20)
plt.setp(ax.get_xticklabels(), fontsize=18)
plt.setp(ax.get_yticklabels(), fontsize=18)
plt.show()


