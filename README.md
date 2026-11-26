Matlab code for detection of change in polarimetric SAR data
in the covariance matrix formulation with auxiliary code.
See my homepage https://people.compute.dtu.dk/alan under Publications.

Here some of the functions are mentioned.


**freadenvisar**
reads covariance matrix polarimetric
SAR data (with an ENVI-like header file) from disk
and stores them in the order needed for the change
detection programs.
SAR images are stored in 3D arrays: rows by columns by matrix elements.
You may need to change the file names in
the code. Also, some files are written in row major order
and some in column major order, so there may be an issue
with transposing the images read.
There is also an endian issue.\
%\
SAR input images are images of the (row-wise) upper-right elements of
the covariance (or coherency) matrix, e.g. for full polarimetry:
a vector with nine real (i.e., non-complex) elements, namely\
[ShhShh* Re(ShhShv*) Im(ShhShv*) Re(ShhSvv*) Im(ShhSvv*)
 ShvShv* Re(ShvSvv*) Im(ShvSvv*) SvvSvv*]\
%\
NOTE : freadenvisar needs the corresponding image header file generated automatically by ENVI. The ENVI header file must have the same name as the ENVI image file + the '.hdr' extension.
Example hdr-file:\
%\
samples = 1024\
lines   = 1024\
bands   = 6\
data type = 4\
interleave = bsq\
%\
Although there are 9 non-complex elements, there are
6 matrix elements only read from 3 files with real numbers
and 3 files with complex numbers, hence "bands = 6" in the hdr-file.\
%\
To be used with for example **wc_Q** and **wcRjl** below, the 3D SAR image arrays, say,
SAR1, SAR2 and SAR3 from 3 time points, must be concatenated along a fourth dimension resulting in a 4D array,
rows by columns by matrix elements by time point:\
SAR = cat(4,SAR1,SAR2,SAR3);

**rho_hhvv**
calculates the complex correlation
between HH and VV polarizations (including the phase
difference) for full polarimetry, reflection symmetry, and
dual polarimetry with HH and VV.

**wishart_det**
calculates determinants of complex
variance-covariance matrices and flags singularity for full
polarimetry (including reflection symmetry and diagonal
only),and dual polarimetry (including diagonal only).

**wishart_change_Q**
calculates the logarithm of the test statistic Q and the
p-value to test for equality of k ≥ 2 covariance matrices
in the full and dual polarimetry cases and of single band
power data (any polarization).

**wishart_change_Rj**
calculates the logarithm of the test statistics Q and R<sub>j</sub>, and the
associated p-values to test for equality of covariance
matrices in full and dual polarimetry data as well as in
single band power data, k ≥ 2. **wishart_change_Rj**
calculates ln(Q) as the sum of the ln(R<sub>j</sub>). This code also exists in
a version which reads, processes,and writes data from/to
disk line-by-line making it useful for large datasets that do
not fit into computer memory.

**wcRjl**
is a wrapper which calls **wishart_change_Rj**
several times, first with data from all
time points (which generates output for the second column
of Table 1), then with data from all time points
omitting the first (which generates output for the third
column of Table 1), then with data from all time points
omitting the 2 first (which generates output for the fourth
column of Table 1), etc. The output is a structure wc
which contains:
- wc(ℓ).lnR (lnQ and lnRj; lnQ is stored in
  lnR(:,:,1),lnR2 inlnR(:,:,2),etc.),
- wc(ℓ).f (f),
- wc(ℓ).rho (ρ and ρj),
- wc(ℓ).omega2 (ω<sub>2</sub> and ω<sub>2j</sub>),and
- wc(ℓ).P (p-values, stored in the same fashion as
  lnQ and lnRj).
For ℓ= 1, all time points are considered
for ℓ= 2, all time points omitting the first are considered;
for ℓ= 3, all time points omitting the two first are considered,
and so on.\
The structure may optionally contain:
- wc(ℓ).sumMat (the sum of X<sub>i</sub> for i=1,...,j)
and
- wc(ℓ).logDet (the logarithm of the determinant
  of this sum mentioned in the section on updating of
  change).
For large data sets the structure wc takes up a lot of
computer memory. You may want to output wc.P only.
Also, you may want to store this in 4- as opposed to 8-byte
precision.

**wc_Q**
like **wishart_change_Q** calculates the logarithm of the test statistic Q and the
p-value to test for equality of k ≥ 2 covariance matrices
in the full and dual polarimetry cases and of single band
power data (any polarization).
Like **wcRjl**, **wc_Q** outputs results to a struct wc, see above.

**wishart_change_Rj_updt**
given data from k time points **wishart_change_Rj_updt** based on
new data in X<sub>k+1</sub>, calculates ln(R<sub>k+1</sub>), p-values, and
auxiliary variables (f, etc.) for full and dual polarimetry
data as well as for single band power data.

**wcRjl_updt**
like **wcRjl** is a wrapper. **wcRjl_updt** calls **wishart_change_Rj_updt**
several times, first with either data from all time points
or the sum of data from all time points (which updates
output for the second column of Table 1 with an extra
row followed by an update of lnQ based on Equation 3),
then with either data from all time points or the sum
of data from all time points, omitting the first (which
updates output for the third column of Table 1 with
an extra row followed by an update of lnQ based on
Equation 3), then with either data from all time points
or the sum of data from all time points, omitting the
two first (which updates output for the fourth column
of Table 1 with an extra row followed by an update of
lnQ based on Equation 3), and so on. A version called
wcRjl_data_updt that handles the situation with no
storage of the sum of X<sub>i</sub> for i=1,...,j
(and the logarithm of the
determinant of this sum) mentioned in the section on
updating is also available.

**wcWhen**
examines the change probabilities (P) in the
structure wc output by **wcRjl** and outputs all
occurrences of change (given a change probability
significance level).

**wc_Ptable**
calculates the no-change probabilities
(1−P) from the structure wc output by **wcRjl** for a
region of interest (ROI) and outputs a table such as Table 1
(or Table 2). It also plots histograms for the no-change
probabilities in the ROI. These histograms are flat for a
ROI with no change (the distribution of 1−P is uniform).
You may input a ROI image or a color image on which to
draw a ROI polygon.

**wc_lnRtable**
calculates −2ρ ln(Q) and −2ρ<sub>j</sub> ln(R<sub>j</sub>)
from the structure wc output by **wcRjl** for a ROI and
outputs a table such as Table 1. It also plots sample
histograms and theoretical distributions for the these
quantities in the ROI.

\begin{itemize}
\item
{\tt \bfseries wishart\_rand} generates one or more SAR scenes in the covariance matrix formulation with specified ENL and \bm{\Sigma}.
To obtain Wishart distributed data we must multiply the generated data by ENL.
\item
{\tt \bfseries wishart\_rand\_patch} generates a  SAR scene in the covariance matrix formulation consisting of a number of patches all with the same specified size and ENL.
As for {\tt \bfseries wishart\_rand} above
to obtain Wishart distributed data we must multiply the generated data by ENL.
Each patch has an individually specified \bm{\Sigma} which is read from a pixel in a SAR image either by ``clicking the mouse'' or by specifying row and column number
(remember that each pixel in a  SAR image in the covariance matrix formulation consists of a Hermitian covariance matrix).
An auxiliary Matlab function {\tt \bfseries polsar\_find\_coord} can be used to find row and column numbers for a number of potentially interesting pixels.
An interesting pixel is for example one with
maximum or minimum $S_{hh}S_{hh}^*$ (divided by the trace),
%maximum or minimum $S_{hh}S_{vv}^*$ (divided by the trace),
maximum or minimum trace,
or
maximum or minimum determinant etc.
\item
{\tt \bfseries wishart\_rand\_scene} reads a SAR scene in the covariance matrix formulation and generates one or more new realisations of the original scene, all with the same size as the original scene based on \bm{\Sigma} read from each pixel in the original scene.
As for {\tt \bfseries wishart\_rand} above
to obtain Wishart distributed data we must multiply the generated data by ENL.
\item
{\tt \bfseries wishart\_rand\_segment} reads one or more SAR scene(s) and an image with segments
of that (those) scene(s),
one segment only also allowed.
It then generates as many scenes as read from the original SAR input, where covariance matrices in all pixels in all segments are replaced by generated realisations of the mean covariance matrix for the original segments.
\end{itemize}

