## Revision from second submission

Victoria Wimmer asked for us to make one change regarding writing to userspace.

Commented out line that calls fwrite in: 

inst/raw/create_variance_from_additive_lookup_table_step1.R

since this violates CRAN policies.

No other changes.

## Revision from first submission

We got some comments from Beni Altmann, they have been addressed in the following manner:

1) Added single quotes in Description field and removed linebreaks.
2) Added reference to one article in the description field.
3) Made authors as contributors with ctb and added deCODE/Amgen as cph and fnd.
4) Replaced all occurrences T and F with TRUE and FALSE.

Hope this covers all issues, thanks for reviewing the package!

## Test environments

* rhub with devtools::check_rhub
* local OS X install, R 4.1.2
* win builder with devtools::check_win_devel

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release and first release for Audunn S. Snaebjarnarson. 
Gudmundur Einarsson has published two packages to CRAN.
