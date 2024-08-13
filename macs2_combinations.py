import pandas as pd
from itertools import product
import sys
# Load the data from the file (replace 'input_file.csv' with your actual file path)
#input_file = '/lrlhps/genomics/prod/lgm/dna_editing/BN24-11550_APOC3_Top8_SiteSeq/user_data/project_data/pipeline_info/atacseq.samplesheet.csv'
input_file = sys.argv[1]
df = pd.read_csv(input_file, sep=',')

# Separate experimental and control samples
experiments = df[df['treatmenttype'] == 'trt']['sample']
controls = df[df['treatmenttype'] == 'control']['sample']

# Generate combinations and corresponding labels
combinations = []
for exper, control in product(experiments, controls):
    #comparison = f"{exper}_vs_{control}"
    comparison = "{}_vs_{}".format(exper, control)
    combinations.append([exper, control, comparison])

# Convert the combinations to a DataFrame
#result_df = pd.DataFrame(combinations, columns=['Experiment', 'Control', 'Comparison'])
result_df = pd.DataFrame(combinations)

# Display the result
#print(result_df)

# Optionally, save to a CSV file
output_file = 'macs2_callpeak_combo.csv'
result_df.to_csv(output_file, index=False, sep=',', header = False)
