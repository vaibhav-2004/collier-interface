# Collier interface for one-loop scalar integrals

Fortran driver and Mathematica wrapper for A0, B0, C0, D0 via COLLIER,
with an independent validation harness.

## Validation

fortran/b0ref.f90 implements B0 analytically, independently of Collier.
fortran/reference.txt holds 15 test points covering unequal and equal
masses, below and above threshold, plus spacelike, on-shell,
zero-momentum and fully massless cases.

Collier agrees with the analytic result at all 15 points to ~15 digits.

## Build

Validate the plumbing first, with no Collier involved:

    gfortran -O2 -o b0ref fortran/b0ref.f90
    ./fortran/check.sh ./b0ref

Then Collier, from https://collier.hepforge.org/ :

    tar xzf collier-1.2.8.tar.gz
    cd COLLIER-1.2.8/build && cmake .. && make -j4 && cd ../..

Then the driver:

    export COL=$PWD/COLLIER-1.2.8
    gfortran -O2 -I$COL/modules fortran/loopdrv.f90 -o loopdrv \
             -L$COL -lcollier -Wl,-rpath,$COL
    ./fortran/check.sh ./loopdrv

## Usage

    ./loopdrv B0 p10 m02 m12 mu2
    ./loopdrv C0 p10 p21 p20 m02 m12 m22 mu2
    ./loopdrv D0 p10 p21 p32 p30 p20 p31 m02 m12 m22 m32 mu2

pij means (p_i - p_j)^2. For D0 the first four momenta are the external
legs and the last two are the diagonals.

## Conventions

SetDeltaUV_cll(0), SetDeltaIR_cll(0,0), mu^2 = SetMuUV2 = SetMuIR2.

Always state mu. B0 is UV divergent and the massless-gluon D0 is IR
divergent, so both finite parts depend on it. C0 with a massive internal
line is finite in both senses and genuinely mu-independent.
