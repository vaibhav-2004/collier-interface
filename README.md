# collier-interface

Command-line driver and Mathematica wrapper for the one-loop scalar
integrals `A0`, `B0`, `C0`, `D0`, computed with
[COLLIER](https://collier.hepforge.org/), plus an independent validation
harness.

Built for a $gg \to t\bar{t}$ spin-density-matrix calculation, but nothing
here is specific to that process.

The driver returns the **full Laurent structure** of each integral:

$$I \;=\; \frac{c_2}{\epsilon^2} \;+\; \frac{c_1}{\epsilon} \;+\; \mathrm{finite}$$

in one call — six numbers, not two. The poles are needed because the IR
divergences of the virtual correction have to cancel against real
emission.

---

## What is in the repo

```
collier-interface/
├── README.md                   you are here
├── fortran/
│   ├── b0ref.f90               B0 computed analytically, no Collier needed
│   ├── loopdrv.f90             the driver: A0, B0, C0, D0 via Collier
│   ├── reference.txt           15 validated B0 values
│   ├── reference_c0d0.txt      10 validated C0 / D0 values
│   ├── check.sh                compares a 2-column executable (b0ref)
│   ├── check_b0.sh             compares the 6-column driver on B0
│   ├── check_c0d0.sh           compares the 6-column driver on C0 / D0
│   └── poles.sh                pole structure of the gg -> ttbar catalogue
└── mathematica/
    └── collier_mathematica.m   Mathematica wrapper, with caching
```

**Start with `b0ref`.** It computes the same `B0` analytically and takes
the same arguments, so you can validate the whole chain — compiler,
argument passing, output parsing, reference values — *before* Collier is
involved. Then if `loopdrv` misbehaves you know the problem is Collier and
not the plumbing.

---

## How the validation works

Three references, and they are not equally strong. Worth knowing which is
which.

| integral | reference | independent? | covers |
|---|---|---|---|
| `B0` finite | analytic derivation | **yes** | everything, 15 points |
| `B0` pole | $c_1 = 1$ by construction | **yes** | all 15 points |
| `C0`, `D0` finite | Feynman-parameter quadrature | **yes** | below threshold only, 10 points |
| `C0`, `D0` above threshold | COLI vs DD | no — internal to Collier | the rest |

The `B0` values come from a hand derivation, so Collier matching them is a
genuine cross-check rather than a tautology. The `C0` / `D0` values come
from direct quadrature of

$$C_0 = -\!\int_\triangle \frac{dx\,dy}{\Delta}, \qquad
  D_0 = \int_\triangle \frac{dx\,dy\,dz}{\Delta^2}$$

which is still independent of Collier, but only converges where the
integrand has no pole inside the region — i.e. below threshold. Every
`C0` / `D0` reference value is therefore real.

**The `B0` grid covers all four mass/threshold cells:**

| | below threshold | above threshold |
|---|---|---|
| equal masses | ✓ | ✓ |
| unequal masses | ✓ | ✓ |

plus spacelike, on-shell, zero-momentum, fully massless, just-above-
threshold, and the exact mass-exchange symmetry.

---

## Requirements

- **gfortran** ≥ 4.8 (Collier is Fortran 2003)
- **cmake**
- **Mathematica**, only for the wrapper

---

## Install

### Linux

```bash
sudo apt install gfortran cmake        # Debian / Ubuntu
# or: sudo dnf install gcc-gfortran cmake
```

### macOS

Apple's `gcc` is clang and has no Fortran compiler, so use Homebrew:

```bash
brew install gcc cmake
```

Homebrew names it `gfortran-13` (or similar), **not** plain `gfortran` —
substitute that name below, or make an alias.

### Windows

Native Windows is not supported; Collier expects a Unix toolchain. Use
**WSL2**:

```powershell
wsl --install -d Ubuntu-22.04
```

then follow the Linux instructions inside it. Note which user you are
(`whoami`) — you need it for the Mathematica section if Mathematica itself
runs on Windows.

---

## Build

### Step 1 — validate the plumbing, no Collier (2 minutes)

```bash
git clone https://github.com/vaibhav-2004/collier-interface.git
cd collier-interface

gfortran -O2 -o b0ref fortran/b0ref.f90
cp fortran/reference.txt fortran/check.sh .
./check.sh ./b0ref
```

Expect **15 OK**. If you get them, your compiler works, arguments pass
correctly, the output format parses, and the reference values are right.
Everything that remains is Collier.

If this fails, fix it before going further.

### Step 2 — build Collier

```bash
wget https://collier.hepforge.org/downloads/collier-1.2.8.tar.gz
tar xzf collier-1.2.8.tar.gz
cd COLLIER-1.2.8
mkdir -p build && cd build
cmake ..
make -j4
cd ../..
```

Then find what it produced — the paths vary between versions:

```bash
find . -name 'libcollier*'
find . -name 'collier.mod'
```

### Step 3 — build the driver

```bash
export COL=$PWD/COLLIER-1.2.8

gfortran -O2 -I$COL/modules fortran/loopdrv.f90 -o loopdrv \
         -L$COL -lcollier -Wl,-rpath,$COL

cp fortran/reference_c0d0.txt fortran/check_b0.sh fortran/check_c0d0.sh .
./check_b0.sh   ./loopdrv      # 15 points, and c1 = 1 on every one
./check_c0d0.sh ./loopdrv      # 10 points
```

On **macOS** shared libraries work differently. If it builds but will not
run, drop `-Wl,-rpath` and set the path at runtime:

```bash
gfortran-13 -O2 -I$COL/modules fortran/loopdrv.f90 -o loopdrv -L$COL -lcollier
export DYLD_LIBRARY_PATH=$COL:$DYLD_LIBRARY_PATH
```

---

## Using the driver

```
./loopdrv A0 m02 mu2
./loopdrv B0 p10 m02 m12 mu2
./loopdrv C0 p10 p21 p20 m02 m12 m22 mu2
./loopdrv D0 p10 p21 p32 p30 p20 p31 m02 m12 m22 m32 mu2
```

Prints six numbers:

```
Re(fin)  Im(fin)   Re(c1)  Im(c1)   Re(c2)  Im(c2)
```

All arguments are **squared** quantities, and `pij` means
$(p_i - p_j)^2$.

> **`D0` argument order is the easy mistake.** The first **four** momenta
> are the external legs; the **last two** are the diagonals (the
> Mandelstams). They are not interleaved, and a wrong order gives a
> plausible-looking number with no error.

Build the invariants from real four-momenta, never by hand — the six `D0`
arguments are not independent, and inconsistent values put Collier into a
degenerate configuration where its two internal algorithms disagree.

### Sanity rules for the poles

| | $1/\epsilon$ | $1/\epsilon^2$ |
|---|---|---|
| `B0` | exactly 1 | 0 |
| `A0(m^2)` | $m^2$ | 0 |
| `C0`, all internal lines massive | 0 | 0 |
| `D0`, two massless internal lines | present | present |
| `D0`, one massless internal line | 0 | 0 |

and every $c_2$ should come out **purely real** — the leading
soft-collinear pole carries no phase.

`./poles.sh ./loopdrv` prints this table for the whole $gg \to t\bar t$
catalogue at a physical phase-space point.

---

## Using it from Mathematica

### Step 1 — point it at the driver

Edit the top of `mathematica/collier_mathematica.m`:

```mathematica
colExe = "/your/path/to/loopdrv";
colWSL = False;
```

`colWSL` is only about where **Mathematica** runs:

| Mathematica | driver | `colWSL` |
|---|---|---|
| Linux or macOS | same machine | `False` |
| inside WSL | inside WSL | `False` |
| **Windows** | **in WSL** | **`True`** |

Unsure? Run `$OperatingSystem` — `"Unix"` means `False`.

**If you set `True`**, also edit the `wsl` call inside `colRun`:

```mathematica
If[colWSL, {"wsl", "-d", "Ubuntu-22.04", "-u", "yourname", "-e", colExe},
           {colExe}],
```

Without this, Windows enters WSL as the *default* user, who may not be
able to read your home directory.

### Step 2 — load it

```mathematica
Get["/your/path/to/collier-interface/mathematica/collier_mathematica.m"]
```

### Step 3 — check it works

```mathematica
colValidate[]
```

Expect **ALL 15 AGREE**. If everything comes back `$Failed`, the path or
the `colWSL` setting is wrong; the raw driver output is printed, which
usually says why.

### Step 4 — use it

```mathematica
mu2 = 172.5^2;

B0c[62500., 125.25^2, 172.5^2, mu2]            (* finite part *)
B0cF[62500., 125.25^2, 172.5^2, mu2]           (* <|"fin","c1","c2"|> *)

C0c[29929., 1000000., 29929., 29929., 0., 29929., mu2]

D0c[0., 0., 29929., 29929., 1000000., -391997.844,
    0., 0., 29929., 29929., mu2]

colPoleTable[]                                  (* the whole catalogue *)
```

> **Use the cached `B0c` / `B0cF` and friends**, not `B0col` / `B0Full`.
> Each scalar integral depends only on the phase-space point, not on which
> helicity amplitude you are building, so the same number is requested many
> times — and without caching each request spawns a new process.

### Building arguments from momenta

```mathematica
kin = colLoopArgs[{500, 0, 0, 500}, {500, 0, 0, -500},
                  {500, 461.6046, 29.5343, -78.2301}];
kin["s"]     (* 1000000     *)
kin["t"]     (* -391997.844 *)
kin["mt2"]   (* 29929       *)
```

---

## Conventions

The driver evaluates each point three times, with

```fortran
(DeltaUV, DeltaIR1, DeltaIR2) = (0,0,0), (1,1,0), (0,0,1)
```

and subtracts, which isolates the Laurent coefficients. UV and IR single
poles are lumped together in $c_1$; set `SetDeltaUV_cll` and
`SetDeltaIR_cll` independently if you need them apart.

The pole convention is
$X = \Delta_\epsilon + X^{\rm fin}$ with
$\Delta_\epsilon = 1/\epsilon - \gamma_E + \ln 4\pi$.

> **Always state $\mu$ when quoting a number.** `B0` is UV divergent, and a
> box with massless internal lines is IR divergent, so both finite parts
> depend on it — the shift equals minus the $1/\epsilon$ coefficient. A
> `C0` with a massive internal line is finite in both senses and genuinely
> $\mu$-independent.

Different libraries split the $\Delta_\epsilon$ pole differently. When
comparing against another implementation, compare **differences** of two
values — the convention cancels.

---

## Known limitations

- `C0` / `D0` **above threshold** rest only on the COLI-vs-DD comparison,
  which is internal to Collier. Quadrature cannot reach there (the
  integrand has a pole inside the region). LoopTools would close this.
- Collier's DD branch aborts on some degenerate configurations — all-zero
  external momenta, exactly on-shell legs — with `fatal error in DD: -10`.
  COLI handles them, which is why mode 1 is the default.
- `RunProcess` spawns a process per call. Fine for validating and for a
  few hundred points; a phase-space scan wants LibraryLink.

---

## Credit

COLLIER is by A. Denner, S. Dittmaier and L. Hofer —
<https://collier.hepforge.org/>. It is **not** included in this
repository; the build instructions download it. Cite their papers if you
use it.
