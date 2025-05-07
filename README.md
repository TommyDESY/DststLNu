# DststLNu
You first need to generate a sample. The procedure described here assumes that you can use basf2.
1) Create the initial root data frame: `python3 gen.py [name of dec file]`
To use gen.py, you need a full basf2 release
2) Process the data frame: `python3 inclusivesl_generator.py -input [input root file] -output [output root file]`
3) Transform the root file in a pandas data frame by running create_df_pandas.py

To do the fit, you can run the `fit_BLR.ipynb` notebook. This is still work in progress. Start with smaller samples as Hammer tends to be quite slow especially on larger samples.
