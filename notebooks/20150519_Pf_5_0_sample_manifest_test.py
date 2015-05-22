# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

%run standard_imports.ipynb

# <codecell>

rob_metadata_fn = '/Users/rpearson/src/github/malariagen/methods-dev/meta/Pf/PGV4_mk5.xlsx'
sampling_date_fn = '/Users/rpearson/src/github/malariagen/methods-dev/meta/Pf/20150519/all_sampling_date.txt'
output_fn = '/Users/rpearson/src/github/malariagen/methods-dev/meta/Pf/20150519/PGV4_mk5_sampling_date.xlsx'

# <codecell>

rob_metadata = petl.io.xlsx.fromxlsx(rob_metadata_fn, 'PGV4.0')
sample_date = petl.io.csv.fromtsv(sampling_date_fn)

# <codecell>

rob_metadata.display(2)

# <codecell>

sample_date.tail(2)

# <codecell>

output_table = (rob_metadata
    .leftjoin(sample_date, lkey='Sample', rkey='ox_code')
)

# <codecell>

(output_table
    .selectne('sampling_date', 'NULL')
    .selectne('sampling_date', None)
    .selectne('Year', '')
    .selectne('Year', None)
    )

# <codecell>

(output_table
    .selectne('sampling_date', 'NULL')
    .selectne('sampling_date', None)
    .selectne('Year', '')
    .selectne('Year', None)
    .addfield('dates_concat', lambda rec: str(rec['Year']) + '_' + rec['sampling_date'][0:4])
    .valuecounts('dates_concat')
    .displayall()
)
print(len(output_table
    .selectne('sampling_date', 'NULL')
    .selectne('sampling_date', None)
    .selectne('Year', '')
    .selectne('Year', None)
)-1)

# <codecell>

(output_table
    .selectne('Year', '')
    .selectne('Year', None)
    .valuecounts('Year')
    .displayall()
)
print(len(output_table
    .selectne('Year', '')
    .selectne('Year', None)
    )-1
)

# <codecell>

(output_table
    .selectne('sampling_date', 'NULL')
    .selectne('sampling_date', None)
    .addfield('dates_concat', lambda rec: rec['sampling_date'][0:4])
    .valuecounts('dates_concat')
    .displayall()
)
print(len(output_table
    .selectne('sampling_date', 'NULL')
    .selectne('sampling_date', None)
)-1)

# <codecell>

petl.io.xlsx.toxlsx(output_table, output_fn, 'PGV4.0')

# <codecell>


