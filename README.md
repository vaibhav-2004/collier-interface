# collier-interface

Command-line driver and Mathematica wrapper for the one-loop scalar
integrals `A0`, `B0`, `C0`, `D0`, computed with
[COLLIER](https://collier.hepforge.org/), plus an independent validation
harness.

Built for a $gg \to t\bar{t}$ spin-density-matrix calculation, but nothing
here is specific to that process.

---

## What is in the repo

```
collier-interface/
├── README.md                      you are here
├── fortran/
│   ├── b0ref.f90                  B0 computed analytically, no Collier needed
│   ├── loopdrv.f90                the real driver: A0, B0, C0, D0 via Collier
│   ├── reference.txt              15 validated B0 values
│   └── check.sh                   compares an executable against reference.txt
├── mathematica/
│   └── collier_mathematica.m      Mathematica wrapper, with caching
└── examples/
    └── scan_B0.sh                 walks B0 across a threshold
```

**Start with `b0ref`.** It computes the same `B0` analytically and takes
the same arguments as the real driver, so you can validate the whole
chain — compiler, argument passing, output format, reference values —
*before* Collier is involved. Then if `loopdrv` misbehaves you know the
problem is Collier and not the plumbing.

---

## The validation

`fortran/reference.txt` holds 15 `B0` values covering every case that can
go wrong:

| case | why it is in the list |
|---|---|
| unequal masses, below threshold | the common case; result is real |
| unequal masses, above threshold | imaginary part appears |
| equal masses, below threshold | separate closed form, different branch |
| equal masses, above threshold | that branch on the other side |
| spacelike, $p^2 < 0$ | prefactor sits on the branch cut |
| on shell, $p^2 = m^2$ | degenerate case |
| zero momentum | the root formula divides by zero here |
| fully massless | checkable by hand: $2 - \ln(p^2/\mu^2) + i\pi$ |
| just above threshold | tiny imaginary part, tests the edge |
| mass exchange symmetry | $B_0(p^2,m_0^2,m_1^2) = B_0(p^2,m_1^2,m_0^2)$ |

The reference values come from an independent analytic derivation, not
from Collier — so Collier agreeing with them is a genuine cross-check
rather than a tautology. All 15 agree to ~15 digits.

---

## Requirements

- **gfortran** ≥ 4.8 (Collier is Fortran 2003)
- **cmake**
- **Mathematica**, only if you want the wrapper

---

## Install

### Linux

```bash
sudo apt install gfortran cmake        # Debian/Ubuntu
# or: sudo dnf install gcc-gfortran cmake
```

### macOS

Apple's `gcc` is clang and has no Fortran compiler, so use Homebrew:

```bash
brew install gcc cmake
```

Homebrew names it `gfortran-13` (or similar), **not** plain `gfortran` —
substitute that name in every command below, or make an alias.

### Windows

Native Windows is not supported; Collier expects a Unix toolchain. Use
**WSL2**:

```powershell
wsl --install -d Ubuntu-22.04
```

then open the Ubuntu terminal and follow the Linux instructions. Note
which user you are (`whoami`) — you will need it for the Mathematica
section if Mathematica itself runs on Windows.

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

Expect **15 OK**. If you get them, your compiler works, arguments are
being passed correctly, the output format parses, and the reference
values are right. Everything that remains is Collier.

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

`make` takes a couple of minutes. Then find what it produced:

```bash
find . -name 'libcollier*'
find . -name 'collier.mod'
```

Note both paths — the next step needs them.

### Step 3 — build the driver

```bash
export COL=$PWD/COLLIER-1.2.8

gfortran -O2 -I$COL/modules fortran/loopdrv.f90 -o loopdrv \
         -L$COL -lcollier -Wl,-rpath,$COL

./check.sh ./loopdrv
```

Expect **15 OK** again — now from Collier rather than the analytic
formula.

On **macOS**, shared libraries work differently. If it builds but will
not run, drop the `-Wl,-rpath` and set the path at runtime instead:

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

Prints `Re` and `Im` of the finite part. All arguments are **squared**
quantities, and `pij` means $(p_i - p_j)^2$.

> **`D0` argument order is the easy mistake.** The first **four** momenta
> are the external legs; the **last two** are the diagonals (the
> Mandelstams). They are not interleaved, and a wrong order gives a
> plausible-looking number with no error.

Example — a $t$-channel box at $\sqrt{s} = 1$ TeV, $m_t = 173$,
$\cos\theta = -0.167$:

```bash
./loopdrv D0 0 0 29929 29929 1000000 -391840.887 0 0 29929 29929 29756.25
```

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

**If you set `True`**, also edit `colRun` (around line 58):

```mathematica
If[colWSL, {"wsl", "-e", colExe}, {colExe}],
```

becomes

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

B0c[62500., 125.25^2, 172.5^2, mu2]

C0c[29929., 1000000., 29929., 29929., 0., 29929., mu2]

D0c[0., 0., 29929., 29929., 1000000., -391840.887,
    0., 0., 29929., 29929., mu2]
```

> **Use `B0c`, `C0c`, `D0c`, not `B0col`, `C0col`, `D0col`.** The `c`
> versions cache. Each scalar integral depends only on the phase-space
> point, not on which helicity amplitude you are building — so the same
> number is requested many times, and without caching each request spawns
> a new process.

### Building arguments from momenta

```mathematica
kin = colLoopArgs[{500, 0, 0, 500}, {500, 0, 0, -500},
                  {500, 461.6046, 29.5343, -78.2301}];
kin["s"]     (* 1000000     *)
kin["t"]     (* -391840.887 *)
kin["mt2"]   (* 29929       *)
```

Use this rather than typing invariants by hand. The six `D0` arguments
are **not independent** — they must come from four real four-momenta.
Inconsistent values put Collier into a degenerate configuration where its
two internal algorithms (COLI and DD) disagree with each other.

### Other functions

```mathematica
colSelfTest[]    (* C0 mu-independence, D0 IR behaviour *)
colTiming[]      (* per-call cost, cached vs uncached   *)
colClearCache[]  (* forget memoised values              *)
```

---

## Conventions

The driver sets

```fortran
SetDeltaUV_cll(0d0)        ! UV pole dropped
SetDeltaIR_cll(0d0, 0d0)   ! IR poles dropped
SetMuUV2_cll(mu2)
SetMuIR2_cll(mu2)          ! both, not just UV
```

so what you get is the **finite part** in the convention
$X = \Delta_\epsilon + X^{\rm fin}$ with
$\Delta_\epsilon = 1/\epsilon - \gamma_E + \ln 4\pi$.

> **Always state $\mu$ when quoting a number.** `B0` is UV divergent, and
> a box with massless internal lines is IR divergent, so both finite
> parts depend on it. A `C0` with a massive internal line is finite in
> both senses and genuinely $\mu$-independent — `colSelfTest[]`
> demonstrates the difference.

Different libraries split the $\Delta_\epsilon$ pole differently. When
comparing against another implementation, compare **differences** of two
values — the convention cancels.

---

## Troubleshooting

| symptom | cause | fix |
|---|---|---|
| `Can't open module file 'collier.mod'` | wrong `-I` path | `find $COL -name collier.mod` |
| `undefined reference to b0_cll_` | wrong `-L`/`-l` | `nm libcollier.so \| grep -i b0` |
| `error while loading shared libraries` | runtime path | add `-Wl,-rpath,$COL`, or `DYLD_LIBRARY_PATH` on macOS |
| Fortran syntax errors | gfortran too old | needs ≥ 4.8 |
| `check.sh` shows all zeros | banner parsing | the script already filters it; check the executable runs alone |
| every value off by the same real constant | $\Delta_\epsilon$ convention | harmless — compare differences |
| imaginary parts have the wrong sign | $i\varepsilon$ convention | check `Im B0 > 0` above threshold |
| only spacelike points wrong, by $2\pi i$ | branch cut on the prefactor | compare against `b0ref` first |
| COLI and DD disagree | inconsistent kinematics | build invariants from real momenta |

---

## Platforms

Tested on Linux (WSL2, Ubuntu 22.04, gfortran 11.4). macOS should work
but is untested, as is Apple Silicon. Native Windows is not supported.

---

## Credit

COLLIER is by A. Denner, S. Dittmaier and L. Hofer —
<https://collier.hepforge.org/>. It is **not** included here; the build
instructions download it. Cite their papers if you use it.
