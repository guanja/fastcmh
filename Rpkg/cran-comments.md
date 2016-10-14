##Resubmission 4 - 13 September 2016

Resubmission/bug fix

**Last fastcmh submission for a long time - promise!**

Two changes:
1) In DESCRIPTION, the line
    Depends: R (>= 3.3.0), bindata
is included to ensure R version 3.3.0 (at least) is used.

In previous submission, the "Depends" was split over two lines, and then
the dependence on the R version was somehow lost.


2) A small typo in an output message of demofastcmh() has been corrected.

Checks:
------------------------------------
Successful checks on:
OSX 10.10.5, R 3.3.1
winbuilder
Linux Mint 18 (Xfce, Ubuntu 16.04), R 3.3.1
Solaris 11.3, R 3.3.0



========================================================================



##Resubmission 3 - 07 September 2016
Resubmission to fix installation failures on Solaris and r-oldrel-windows.

Three changes:
1) 
Put back:
#include<sstream>

and went back to converting using ostringstream. Just #include<fstream> was causing issues.

2)
In fdr.cpp, in function ‘double computeApproxHarmonicLower(long long int)’:
now first cast long long to double before computing log

3)
Now require R version at least 3.3.0 - just to make sure everything works.


Checks:
------------------------------------
Successful checks on:
Solaris 11.3, R 3.3.0
Linux Mint 18 (Ubuntu 16.04), R 3.3.1
winbuilder
OSX 10.10.5, R 3.2.2 and R 3.3.1



========================================================================




##Resubmission 2 - 05 September 2016
This is a resubmission in order to deal with Rf_errors which were generated in the last submission.

Error was something of the form: 
"error: 'RF_error' is not a member of 'std::codecvt_base'   if (__result == codecvt_base::error) ..." etc

Error was fixed by commenting out lines:  
#include<fstream>
#include<sstream>
in fastcmh_cpp.cpp

It will not install on R 3.1.2 (some Rcpp errors), so minimum version has been set to R 3.2.2 (I am certain it works on this version, and R 3.3.0, 3.3.1).



## Resubmission
This is a resubmission. In this version I have:

* Created src/Makevars and src/Makevars.win files which contain  
>> CXX_STD = CXX11
, which I hope will deal with the warnings about the type long long.


* The examples in the function 
>> makefastcmhdata
(previously makeSampleData) do not write to file now - before the examples wrote to the home directory, which was an error.

