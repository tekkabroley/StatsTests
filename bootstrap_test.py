import pandas as pd

from matplotlib import pyplot as plt

from numpy import sqrt
from scipy import stats

# option to pass the filepath as a commandline param
from sys import argv

from sys import exit

######################## Dependencies
# path to directory containing file
work_path = ""

# name of file
filename = ""

# path to file
filepath = work_path + filename

# value of variant for each group
control_group = "control"
test_group = "test"

# the position of each column in the file
id_col = 0
variant_col = 1
metric_col = 2

# the number of resamples to take for bootstrapping
resamples = 1000

######################## Data Import
# csv should be formatted to: id, variant, metric
# load csv into dataframe

if len(argv) > 1:
    filepath = argv[1]
else:
    filepath = work_path + filename

try:
    dataframe = pd.read_csv(filepath)
except Exception as ex:
    print("caught execption: {}".format(ex))
    print('exiting script...')
    exit()


columns = list(dataframe.columns)

########################## Descriptive Stats

id = dataframe[columns[id_col]]
rows = id.count()
print("population size: %s" % rows)

variant_name = columns[variant_col]
variant = dataframe[variant_name]

# mask for each variant
control_mask = variant.map(lambda v: v == control_group)
test_mask = variant.map(lambda v: v == test_group)

# create metric data sets for each group
metric_name = columns[metric_col]
control_metric = dataframe[control_mask][metric_name]
test_metric = dataframe[test_mask][metric_name]

# count number of rows in each group
control_pop = control_metric.count()
test_pop = test_metric.count()

print("number of rows of data in population for each group")
print("control: %s test: %s" % (control_pop, test_pop))

# calculate proportion of test cohort in dataframe
test_proportion = float(test_pop)/rows
test_percent = test_proportion*100
print("percentage of population in test group %s" % test_percent)

# calculate descriptive stats
control_mean = control_metric.mean()
control_var = control_metric.var()
print("conrol mean: %s std: %s" % (control_mean, sqrt(control_var)))

test_mean = test_metric.mean()
test_var = test_metric.var()
print("test mean: %s std: %s" % (test_mean, sqrt(test_var)))

# observed difference between test mean and control mean
observed_difference = test_mean - control_mean
print("observed difference between test and control: %s" % observed_difference)

# percent uplift in test over control
control_metric_sum = control_metric.sum()
test_metric_sum = test_metric.sum()

percent_uplift = ((test_metric_sum - control_metric_sum)/(control_metric_sum))*100
#print("percent uplift in test over control: %s" % percent_uplift)

# percent uplifts in test mean over control mean
percent_uplift_mean = (observed_difference/control_mean)*100
print("percent uplift in test mean over control mean: %s" % percent_uplift_mean)


############################### Hypothesis Testing
samples = pd.Series()
binomial_dist = stats.binom(1, test_proportion)

for i in range(resamples):
	cohort_assignment = []
	sample_data = dataframe.copy()
	random_binom = binomial_dist.rvs(rows)
	for j in range(rows):
		if random_binom[j] == 0:
			cohort_assignment.append(control_group)
		else:
			cohort_assignment.append(test_group)

	sample_data[variant_name] = cohort_assignment

	sample_test_mask = sample_data[variant_name].map(lambda v: v == test_group)
	sample_test_mean = sample_data[sample_test_mask][metric_name].mean()

	sample_control_mask = sample_data[variant_name].map(lambda v: v == control_group)
	sample_control_mean = sample_data[sample_control_mask][metric_name].mean()

	sample_difference = sample_test_mean - sample_control_mean

	samples = samples.append(pd.Series(sample_difference))


samples_data = pd.DataFrame({"values": samples})

# output histogram of difference of sample means
samples_data.hist(column='values', bins=100)

# p value calculation
samples_len = samples.count()
print(float(samples.map(lambda s: 1 if abs(observed_difference) >= abs(s) else 0).sum())/samples_len)
p_value = 1 - float(samples.map(lambda s: 1 if abs(observed_difference) >= abs(s) else 0).sum())/samples_len
print("p value: %s" % p_value)

plt.show()
