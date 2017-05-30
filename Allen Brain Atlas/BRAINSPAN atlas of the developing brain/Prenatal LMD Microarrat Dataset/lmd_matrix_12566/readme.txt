==============================================
BrainSpan: Atlas of the Developing Human Brain
==============================================

This data set contains expression levels measured by microarray profiling (see the whitepaper at
www.brainspan.org).

expression_matrix.csv

    Contains normalized expression values arranged by (row, column) =
    (probe, samples).  The first column is the probe ID.

pa_call.csv 

    Contains a present/absent flag which indicates whether the probe's
    expression is well above background.  It is set to 1 when both of the
    following conditions are met.
    
        1) The 2-sided t-test p-value is lower than 0.01, (indicating the mean
           signal of the probe's expression is significantly different from the
           corresponding background).
        2) The difference between the background subtracted signal and the
           background is significant (> 2.6 * background standard deviation).

rows_metadata.csv

    Metadata for the probes listed in the same order as the rows in expression_matrix.csv.

columns_metadata.csv

    Metadata for the samples listed in the same order as the columns in
    expression_matrix.csv.
