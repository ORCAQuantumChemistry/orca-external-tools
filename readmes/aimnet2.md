# AIMNet2 calculator — options and usage

This document covers AIMNet2-specific options for OET's `oet_aimnet2`
wrapper. For installation, generic ORCA wiring, and server / client
basics see the main [README](../README.md).

## Default model

Default model is `aimnet2` (= `aimnet2-wb97m-d3_0`, ωB97M-D3 trained).
Other models are selectable via `-m`. For non-covalent or screening
work, `aimnet2-2025` (B97-3c trained, faster) is a good alternative;
for general thermochemistry and reaction barriers, `aimnet2` remains
the recommended default per upstream guidance. For open-shell or
charged systems, switch to `aimnet2-nse` (see warning below).

> **Open-shell users — read this.** The default `aimnet2` (and
> `aimnet2-2025`, `aimnet2-b973c-d3*`, `aimnet2-rxn*`) are
> closed-shell-trained. Passing `mult ≠ 1` (e.g. an ORCA input with
> `* xyzfile 0 2 OH.xyz`) is **silently accepted** and produces a
> spin-restricted energy, not a true UKS-equivalent open-shell value.
> For genuine open-shell or charged-species energetics, use
> `-m aimnet2-nse`. OET emits a one-time `UserWarning` at calculator
> setup when this combination is detected.
>
> No AIMNet2 family does broken-symmetry, so a closed-shell-singlet
> biradical (`mult = 1` with two unpaired electrons coupled
> antiferromagnetically — carbenes, nitrenes, stretched singlet σ
> bonds) is not modelable correctly by **any** model in the suite,
> NSE included. Treat such systems with a multireference method.

## Choosing a model

**WARNING — energies from different model families are NOT comparable.**
Different families were trained on different reference data, with
different functionals, with or without solvation. Mixing models within
a single workflow (e.g. optimize with `aimnet2-rxn`, compute SP energies
with `aimnet2-2025`) produces wrong reaction energies of order tens of
kcal/mol — comparable to or larger than typical reaction barriers, with
no failure indication. The aimnet runtime emits a one-time
`UserWarning` when two families are constructed in the same Python
process; this is informational, not a bug.

Recommended pattern: TS search in `aimnet2-nse` → SP energies in
`aimnet2-nse`. Never optimize in one family and compute SP in another.
The exception is `aimnet2-pd`, which can be used end-to-end for
Pd-catalyzed THF-solvated workflows but must not be mixed with any
other family at all.

| Workflow | Recommended model | Notes |
|----------|-------------------|-------|
| Closed-shell, neutral, equilibrium geometry, gas-phase | `aimnet2` | ωB97M-D3 trained |
| Closed-shell, neutral, non-covalent / screening / large systems | `aimnet2-2025` | B97-3c trained; faster but barriers typically 3-5 kcal/mol off vs the ωB97M-D3 default. Fine for relative ranking and screening, not absolute kinetics. |
| Open-shell or charged | `aimnet2-nse` | NSE: handles arbitrary charge / multiplicity. 14-element coverage (H, B, C, N, O, F, Si, P, S, Cl, As, Se, Br, I). Covers single-reference DFT regimes. NOT reliable for biradicals (two unpaired electrons on different atoms with weak coupling), near-degenerate spin states, or stretched singlet bonds (e.g. homolytic dissociation on a singlet surface). |
| Reactive (TS, NEB, IRC), broad element set | `aimnet2-nse` | Same caveats as above; covers reactive open-shell trajectories. |
| Reactive, neutral H/C/N/O only, faster | `aimnet2-rxn` | NEUTRAL only (charge ≠ 0 raises `ValueError`). H/C/N/O only (other elements raise `ValueError`). Coulomb cutoff locked at 4.6 Å (other values fire an upstream `UserWarning`; OET also warns at setup time). |
| Pd-containing (catalysis, organometallic) | `aimnet2-pd` | B97-3c with **implicit CPCM/THF solvation BAKED IN** — energies are NOT gas-phase. Replaces As with Pd vs default; AsR3 ligands (arsine ligands on Pd, e.g. Pd(AsPh3)4-style complexes) are not supported. |

## Reproducibility

Aliases (`aimnet2`, `aimnet2-2025`, …) may be repointed in future
aimnet releases. For bit-stable reproducibility across versions, pin
to canonical keys: `-m aimnet2-wb97m-d3_0` rather than `-m aimnet2`.

## Performance flags

| Flag | Default | Effect |
|------|---------|--------|
| `--compile` | False | **Server mode only.** Compile the model with `torch.compile` for faster subsequent runs after the first call (introduces an additional latency of 10–60 s on the first call). Recompiles on shape changes (catastrophic in NEB) and is incompatible with Hessians, thus **do not use with NEB / OptTS / IRC**. For server mode, set `TORCHINDUCTOR_CACHE_DIR=/persistent/path` to keep the warm cache across worker restarts. |
| `--nb-threshold N` | 120 | Adaptive neighbor-list batch size. |
| `--ensemble-member {0,1,2,3}` | 0 | Define the member of the method ensemble (default 0). OET runs ONE model per call. For production-quality accuracy, run each member (`0`,`1`,`2`,`3`) separately and average the energies / gradients outside ORCA; upstream recommends ensemble averaging for production calculations. |

## Numerical precision

AIMNet2 internally evaluates in float32. For large systems (energy magnitudes
above ~1000 eV), absolute energies are precise to ~1e-4 eV ≈ 4e-6 Eh. ORCA's
default `TolE=5e-6 Eh` is at this noise floor; if optimization stalls on
`Energy convergence not reached`, loosen with `! TightOpt` → `%scf TolE 1e-5 end`
or accept slightly looser convergence.

## GPU server deployment

Each `oet_server aimnet2 -d cuda` worker holds the model resident on GPU
until evicted. `torch.cuda.empty_cache()` returns memory to PyTorch's allocator,
not to the OS — `nvidia-smi` shows steady-state usage equal to the sum of
resident workers' models. Plan worker count by VRAM/model-size, not per-call
peak. For multi-worker GPU deployments, set `CUDA_VISIBLE_DEVICES` per worker
or use `--workers 1`. (Server-side worker-cache eviction with `release()` is
in a sibling PR.)

> **`aimnet2-rxn` family CAUTION**: the Coulomb cutoff is locked at 4.6 Å
> in training. Always pass `--coulomb-cutoff 4.6` when using `-m aimnet2-rxn*`.
> OET emits a `UserWarning` at setup time if a different value is used; the
> upstream `UserWarning` from `aimnet` also fires, and your electrostatics
> are physically suspect outside the trained cutoff.

## Long-range Coulomb flags

| Flag | Default | Effect |
|------|---------|--------|
| `--coulomb {auto,on,off}` | `auto` | Force on/off the long-range Coulomb module. `auto` defers to model. |
| `--coulomb-method {simple,dsf,ewald}` | unset | Override the model's method. |
| `--coulomb-cutoff R` | 15.0 | Cutoff in Å (used by dsf/ewald). Rejected without `--coulomb-method`. **For `aimnet2-rxn` family, pass `--coulomb-cutoff 4.6`** (training-frozen value). |

## Dispersion (DFT-D3) flags

| Flag | Default | Effect |
|------|---------|--------|
| `--dispersion {auto,on,off}` | `auto` | Force on/off the D3 module. |
| `--dftd3-cutoff R` | unset | Override D3 cutoff in Å. |
| `--dftd3-smoothing-fraction f` | unset | Override D3 smoothing fraction. |

## Special flags

| Flag | Default | Effect |
|------|---------|--------|
| `--download-only` | unset | Only download model files and exit. |

## Examples

CLI smoke (default model, water gradient):

```bash
oet_aimnet2 input.extinp.tmp
```

Condensed phase, DSF Coulomb at 12 Å:

```bash
oet_aimnet2 input.extinp.tmp --coulomb-method dsf --coulomb-cutoff 12.0
```

ORCA OptTS with the NSE model on a charged radical (no `--compile`):

```text
! OptTS
%method
  ProgExt "/path/to/oet_aimnet2"
  Ext_Params "-m aimnet2-nse -d cuda"
end
* xyzfile -1 2 ts_radical.xyz
```

ORCA NEB-CI with the RXN model (neutral H/C/N/O only; no `--compile`,
`--coulomb-cutoff 4.6`):

```text
! NEB-CI
%neb Product "product.xyz" NImages 8 end
%method
  ProgExt "/path/to/oet_aimnet2"
  Ext_Params "-m aimnet2-rxn --coulomb-method dsf --coulomb-cutoff 4.6"
end
* xyzfile 0 1 reactant.xyz
```
