READ ME:

The project aims to prepare case study for synergyfinder plus manuscript and test the updated Rpackage synergyfinder (http://bioconductor.org/packages/release/bioc/html/synergyfinder.html).


Issues noticed so far:

1) Warning message issue: It is quite annoying from user point of view to see tons of warning message when calling CalculateSensitivity and CalculateSynergy function: NaN Produced, Convergence failed, The model was not fitted. Speed issue: 
2) *CalculateSensitivity* and *CalculateSynergy* can be reaaaally slow when the number of block is large. 
3) Plotting issue:
    Barometer plot tick can be improve. Currently when the hand is overlapping with the number areas it does not look good.

4) Help page issues:
    (1)  The arugment *data* required in the *CalculateSynergy* function can either be object returned by *ReshapeDate* or*CalculateSensitivity*, depending on what is the need for subsequent analysis. For example the subsequent analysis is to PlotMultiDrugSurface onto the 
    (2) *PlotSensitiveSynergy* have to use the object returned first by CalculateSensitivity then by CalculateSynergy. The current argument defination is wrong. 
    
