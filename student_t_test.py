import pandas as pd
from numpy import sqrt
from scipy import stats

# option to pass the filepath as a commandline param
from sys import argv

from sys import exit

##################### Dependencies
# path to directory containing file
work_path = ""

# name of file
filename = ""

# value of variant for each group
control_group = 20181001
test_group = 20181201

# the position of each column in the file
id_col = 0
variant_col = 1
metric_col = 2

# set the alpha level for all tests used in this script
alpha = 0.05

##################### Data Import
# csv should be formatted to: id, variant, metric
# load csv into dataframe

# path to file
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

# lead column should be an id
ids = dataframe[columns[id_col]]
rows = ids.count()
print("population size: %s" % rows)

##################### Descriptive Stats
# mask for each variant
control_mask = dataframe[columns[variant_col]] == control_group
test_mask = dataframe[columns[variant_col]] == test_group

# create data sets for each group
metric = columns[metric_col]
control_metric = dataframe[control_mask][metric]
test_metric = dataframe[test_mask][metric]

# count number of rows in each group
control_pop = control_metric.count()
test_pop = test_metric.count()
print("number of rows of data in population for each group")
print("control: %s test: %s" % (control_pop, test_pop))

# calculate proportion of test cohort in data_frame
percent_test = (test_pop/rows)*100
print("percentage of population in test group %s" % percent_test)

# calculate descriptive stats
control_mean = control_metric.mean()
control_var = control_metric.var()
print("conrol mean: %s std: %s" % (control_mean, sqrt(control_var)))

test_mean = test_metric.mean()
test_var = test_metric.var()
print("test mean: %s std: %s" % (test_mean, sqrt(test_var)))

# percent uplifts in test mean over control mean
percent_uplift_mean = ((test_mean - control_mean)/control_mean)*100
print("percent uplift in %s mean over %s mean: %s" % (test_group, control_group, percent_uplift_mean))

#########################  Hypothesis Testing
# compute pearsonr test for h_A:r_test != r_control
pearsonr_obj = stats.pearsonr(control_metric, test_metric)
#print("correlation coef: %s p-value: %s" % (pearsonr_obj[0], pearsonr_obj[1]))
p_rtest = pearsonr_obj[1]
is_correlated = p_rtest <= alpha
print("are groups correlated? %s" % is_correlated)

# compute flinger's test for h_A: sig^2_test != sig^2_control
fligner_obj = stats.fligner(control_metric, test_metric, center="mean")
p_fligner = fligner_obj[1]
is_var_equal = p_fligner > alpha
print("is variance of groups equal? %s" % is_var_equal)

# compute student t test for h_A: mu_test != mu_control
if is_correlated:
    p_ttest = stats.ttest_rel(control_metric, test_metric)[1]
else:
    p_ttest = stats.ttest_ind(control_metric, test_metric, equal_var=is_var_equal)[1]

print("t test p value: %s" % p_ttest)

# output test results
if p_ttest <= alpha:
  print("reject null hypothesis, means are not equal")
else:
  print("fail to reject null hypothesis")
