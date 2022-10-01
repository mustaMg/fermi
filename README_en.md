# Search-of-Extended-Emissions-on-GRB-s
My work about swift's GRB data analysis at TUBITAK STAR project.

[Türkçe](https://github.com/mustafagumustas/Search-of-Extended-Emissions-on-GRB-s/blob/main/README.md)

## What is Gamma Ray Burst (GRB)?
High energy event that happens out of the Milky Way. They are separated by their time that is called T90. T90, %90 of the time that all the photons comes to the detector. Long GRB are the one that has T90>2 and has soft spectrum, the short ones are T90<2 and has hard spectrum.

Some of the GRB's not fits in these two categories exactly. They have an after spike later on the main spike. The reason that we are doing this is once we understand better of these hybrid GRB's our understanding about long GRB's going to get better.


## EE Search
I chose 10 GRB's for this project. These GRB's detected by Swift satellite in the type of "light curve" fits files. From that database, I get these types of events:
- Single energy channel (15-150 eV), 1 sec resolution, background extracted
- Single energy channel, 64 msec resolution, background extracted
- Four energy channel (15-25, 25-50, 50-100, 100-150 eV), 1 sec resolution, background extracted
- Four energy channel, 64 msec resolution, background extracted LC file
- Single energy channel, 64 msec resolution, unweighted LC file


## Filtering Data Based on Morphological Criteria
- The maximum count rate has to be in the first 5 second after the trigger. 
- The count rates remain below 30 per cent (or 40 per cent for some BAT events with low peak rate of <11k count/s) of the peak count rate for at least 50 per cent of the rest of the duration after the peak time until +5s.

## Coding
The criteria on upper page applied to the code and the result visualized with marking EE's on the graph if that event has it.
