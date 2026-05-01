# Installing Claude Code and Gemini CLI on macOS

This guide covers installing two command-line AI coding assistants on a Mac:
1. **Claude Code** (Anthropic)
2. **Gemini CLI** (Google)

---

## Prerequisites

Open the **Terminal** app (Applications → Utilities → Terminal, or press `Cmd + Space` and type "Terminal").

### 1. Install Homebrew (if not already installed)

Homebrew is the standard package manager for macOS.

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

Verify:

```bash
brew --version
```

### 2. Install Node.js (required for both CLIs)

Both tools need Node.js 18 or newer.

```bash
brew install node
```

Verify:

```bash
node --version
npm --version
```

---

## Part 1: Install Claude Code

### Option A — Native installer (recommended)

```bash
curl -fsSL https://claude.ai/install.sh | bash
```

### Option B — npm

```bash
npm install -g @anthropic-ai/claude-code
```

> Do **not** use `sudo npm install -g`. If you hit permission errors, fix your npm prefix or use the native installer above.

### Verify and launch

```bash
claude --version
```

To start Claude Code in a project folder:

```bash
cd ~/path/to/your/project
claude
```

On first launch you'll be prompted to log in with your Anthropic account (Claude Pro/Max subscription) or an API key from https://console.anthropic.com/.

### Useful commands inside Claude Code

- `/help` — show help
- `/login` — switch accounts
- `/config` — open settings
- `/exit` — quit

---

## Part 2: Install Gemini CLI

### Option A — Homebrew (recommended)

```bash
brew install gemini-cli
```

### Option B — npm

```bash
npm install -g @google/gemini-cli
```

### Verify and launch

```bash
gemini --version
```

To start Gemini CLI:

```bash
cd ~/path/to/your/project
gemini
```

On first launch, choose a sign-in method:
- **Sign in with Google** (free tier with personal Google account)
- **Gemini API key** from https://aistudio.google.com/apikey — set with:
  ```bash
  export GEMINI_API_KEY="your_key_here"
  ```
  Add that line to `~/.zshrc` to make it permanent, then run `source ~/.zshrc`.

### Useful commands inside Gemini CLI

- `/help` — show help
- `/auth` — change auth method
- `/quit` — exit

---

## Updating later

```bash
# Claude Code (if installed via the native installer)
claude update

# Claude Code (if installed via npm)
npm update -g @anthropic-ai/claude-code

# Gemini CLI (Homebrew)
brew upgrade gemini-cli

# Gemini CLI (npm)
npm update -g @google/gemini-cli
```

---

## Uninstalling

```bash
# Claude Code
npm uninstall -g @anthropic-ai/claude-code
# or, if installed via native installer:
rm -rf ~/.claude ~/.local/bin/claude

# Gemini CLI
brew uninstall gemini-cli
# or
npm uninstall -g @google/gemini-cli
```

---

## Troubleshooting

- **`command not found` after install** — close and reopen Terminal, or run `source ~/.zshrc`.
- **npm EACCES permission errors** — don't use `sudo`. Run `npm config set prefix ~/.npm-global` and add `export PATH=~/.npm-global/bin:$PATH` to `~/.zshrc`.
- **Apple Silicon (M1/M2/M3/M4)** — both tools support arm64 natively; no Rosetta needed.

---

# Part 3: Reproduce the `xmr-evolve` `python_fortran` Test Suite on macOS

This section reproduces **only** the `python_fortran/` path of
`saisuryadv/xmr-evolve` — the Python + Fortran MR³-GK bidiagonal SVD
implementation. It does **not** build the C/CLAPACK side of the repo.

The repo ships a prebuilt `libxmr.so` for **Linux x86_64 only**, so on a Mac
you must rebuild it from the Fortran sources before any tests will run.

## 1. Install all macOS system dependencies

Open Terminal and install everything in order.

### 1a. Xcode Command Line Tools (clang, make, git)

```bash
xcode-select --install
```

### 1b. Homebrew (if not already installed)

```bash
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

### 1c. Compilers and numerical libraries

The upstream README lists Linux/apt packages
(`gfortran liblapack-dev libblas-dev`). The macOS equivalents:

```bash
brew install gcc lapack openblas python@3.11 git
```

- `gcc` — Homebrew's gcc package ships **gfortran** (the Apple-provided
  toolchain has clang but no Fortran compiler).
- `lapack` and `openblas` — provide `-llapack` and `-lblas` that `build.sh`
  links against.
- `python@3.11` — the runtime for the test driver.

### 1d. Verify everything is on PATH

```bash
gfortran --version
gcc --version
python3 --version
brew --prefix lapack
brew --prefix openblas
```

## 2. Clone the repo

```bash
cd ~
git clone https://github.com/saisuryadv/xmr-evolve.git
cd xmr-evolve/python_fortran
```

(No hardcoded path issues in this directory — `build.sh` uses
`cd "$(dirname "$0")"`, so anywhere you clone is fine.)

## 3. Set up the Python environment

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install numpy scipy
```

(`numpy` and `scipy` are the only Python runtime deps used by `evaluate.py`,
`mr3_gk.py`, and the test scripts.)

## 4. Make `build.sh` macOS-friendly

The script as written assumes Linux conventions. Two small edits before you
run it:

### 4a. Point the linker at Homebrew's LAPACK/BLAS

Open `build.sh` and find the final `gcc -shared` link line. Change:

```
  -lgfortran -llapack -lblas -lm
```

to:

```
  -L"$(brew --prefix lapack)/lib" -L"$(brew --prefix openblas)/lib" \
  -lgfortran -llapack -lblas -lm
```

### 4b. (Apple Silicon only) tell gcc where libgfortran lives

If `-lgfortran` fails to resolve at link time, also add:

```
  -L"$(brew --prefix gcc)/lib/gcc/current"
```

before `-lgfortran`. The exact path is whatever
`brew --prefix gcc`/lib/gcc/<version>/ contains `libgfortran.dylib`.

## 5. Build `libxmr.so` from source

```bash
bash build.sh
```

You should see lines like `Compiled NN unmodified XMR files`, then
`Linking libxmr.so`, then `Built libxmr.so (NNNNNN bytes)`. The shipped
Linux `libxmr.so` is overwritten with your fresh macOS build (a Mach-O
shared library — ctypes loads it fine despite the `.so` extension).

## 6. Run the test suite

From inside `python_fortran/` with the venv active:

```bash
# Full 379-test suite (definitive evaluation)
python3 evaluate.py

# Faster 270-test medium suite (skips n=400 and STCollection)
python3 evaluate.py --medium
```

Per-test output looks like:

```
  pattern_name    n=NNNN  res=X.XXX  ortU=X.XXX  ortV=X.XXX  t=X.XXXXs  PASS/FAIL
```

Pass thresholds: `res ≤ 7.0 nε`, `ortU ≤ 5.0 nε`, `ortV ≤ 5.0 nε`.

**Expected reproduction target: 379/379 (100%) with score 92.58.** This is
the current state with the default build (`USEBLOCKS=.TRUE.` in
`xmr_src/dlaxrs_stat.f`), per `ablations.md`. The 330/379 figure quoted at
the top of `python_fortran/README.md` is from an older phase of the project
and is out of date — trust `ablations.md`. If you build with
`USEBLOCKS=.FALSE.`, you'll get 373/379 instead.

### Optional: extended 224-test suite

```bash
python3 test_dense_to_bidiag.py            # full
python3 test_dense_to_bidiag.py --quick    # sizes {10, 100} only
python3 test_dense_to_bidiag.py --part1    # dense-to-bidiag only
python3 test_dense_to_bidiag.py --part2    # paper tests only
```

Expected: **219/224 passing (≈97.8%)**.

### Optional: baseline comparison

```bash
python3 run_baselines.py
```

Compares MR³-GK against DBDSQR (and TGK+DSTEMR where it doesn't crash).

## 7. Troubleshooting

- **`OSError: ... libxmr.so: invalid ELF header`** — you're trying to load
  the shipped Linux binary on macOS. Rerun `bash build.sh` to rebuild it.
- **`ld: library 'lapack' not found`** — add the `-L"$(brew --prefix lapack)/lib"`
  flag from step 4a; on Apple Silicon, Homebrew lives under `/opt/homebrew`
  and isn't in the default linker search path.
- **`ld: library 'gfortran' not found`** — apply the fix in step 4b; locate
  `libgfortran.dylib` with `find $(brew --prefix gcc) -name 'libgfortran*'`.
- **`gfortran: command not found`** — `brew install gcc` (the gcc bottle
  bundles gfortran; Apple's toolchain doesn't).
- **`ModuleNotFoundError: numpy`** — activate the venv:
  `source .venv/bin/activate`, then `pip install numpy scipy`.
- **All tests FAIL with huge residuals** — likely an ABI mismatch from
  dropping the `-fno-second-underscore` flag in step 4. Keep the original
  `FFLAGS` in `build.sh` unchanged; only edit the link line.
- **Pass rate is 373/379 instead of 379/379** — your build has block
  factorization disabled. Check `xmr_src/dlaxrs_stat.f` for
  `USEBLOCKS = .TRUE.` (the shipped default) and rebuild with
  `bash build.sh`.
