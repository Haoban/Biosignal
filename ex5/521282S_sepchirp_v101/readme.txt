02-04-2016 V.1.01 SEPchirp

This package contains the software (m-function source code) and sample test data for polynomial phase chirp signal based SEP parameterization method presented in

Väyrynen E, Noponen K, Vipin A, Thow X, Al-Nashash H, Kortelainen J, All A. Automatic parametrization of somatosensory evoked potentials with chirp modeling. IEEE Transactions on Neural Systems and Rehabilitation Engineering. Accepted.

Please, refer to this paper when using the package for e.g. scientific purposes.

The methods requires a Particle Swarm Optimization (PSO) toolbox for parameter optimization (available for download here: http://www.mathworks.com/matlabcentral/fileexchange/7506-particle-swarm-optimization-toolbox).

A license to the software is granted under the BSD License terms, see LICENSE.

Files included in the package (for version changelog see updated files):

readme.txt - This text.

LICENSE - License file.

chirp_features.m - A function using PSO to fit the chosen chirp model to SEP recordings (v. 1.01)

Chirpmodel.m - A chirp generator to produce single chirps using the chosen model

objFun.m - An objective function used with the PSO to produce linearly weighed mean squared error measures using the chirpmodel.m (v. 1.01)

datasample.mat - A sample data file containing a synthetic SEP recording sample template with 100 sequential repeated measures (AWGN with 5dB SNR) in the used data format.


