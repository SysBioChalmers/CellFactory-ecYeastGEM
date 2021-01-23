# -*- coding: utf-8 -*-

import pandas as pd
# read the product information
chemical_info = pd.read_excel("ComplementaryData/chemicals_info.xlsx")

# input readme file
with open('README.md') as f:
    lines = f.readlines()
part1 =lines[0:20]
index0 = lines.index("## Installation\n")
part2 = lines[index0:]


newline1 = '* Product list:\n'
product_template ='| Name | Formula | KEGG ID | CHEBI ID | Group | class | Gene target |\n'
newline2 = '|:-------:|:--------------:|:---------:|:----------:|:-----:|:-----:|:-----:|\n'
products_info = []
for i, x in chemical_info.iterrows():
    print(i, x)
    newline0 = product_template.replace("Name", x["Name"]).replace("Formula", str(x["Formula"])).replace("KEGG ID", str(x["KEGG ID"])).\
        replace("CHEBI ID", str(x["CHEBI ID"])).replace("Group", x["Group"]).replace("class", x["class"])
    # query the gene target prediction result
    newline0 = newline0.replace("Gene target", "Link")
    products_info.append(newline0)

# combine the product information as a readme file
new_lines = part1 + [newline1] + ['\n'] + [product_template] + [newline2] + products_info + part2

out = open('README.md', 'w')
for line in new_lines:
    print(line)
    out.write(line)
out.close()


















